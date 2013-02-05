!==================================================================================================
!
!
!
!rdefect = defect tolerance bond length
!
!IMPORTANT
!fort.136 = Tagged H30 (includes the proton transfer)
!fort.138 = More than one H3O listing
!fort.145 = Tagging process
!fort.198 = Absolute index of proton transfer
!fort.199 = Absolute index of Singular H3O
!
!==================================================================================================
 PROGRAM Hbond_main
   !
   USE main_variables
   USE defect_module,      only : init_defect, find_defect
   USE defect_module,      only : ndefect, noH3O, dblH3O, multiH3O, OHcount
   USE hbonds_module,      only : init_hbonds, find_hbonds, find_distance_delta
   USE hbonds_module,      only : donate_single, donate_dbl, donate_total
   USE hbonds_module,      only : accept_single, accept_dbl, accept_total
   USE hbonds_module,      only : ndelta_below, ndelta_above, delta
   USE rdf_module,         only : count_pairs, final_rdf
   USE gen_lib
   USE omp_lib
   !
   implicit none
   !
   integer                 :: na,            &  !Atomic index
                              ns,            &  !Speices index
                              ncount=0,      &  !number of steps index
                              readstep,      &  !read the nfi from *.pos 
                              i,j
   
   real(DP)                :: time              !time step

   integer                 :: start_time, end_time, total_time
  
   !
   !Start Time
   !$ start_time = OMP_get_wtime()   
   !
   !Read-in and Initialize main variables
   CALL INIT()
   !
   !Initialize find_defect variables
   call init_defect()
   !
   !Initialize find_Hbond variables
   call init_hbonds()
   !
   !******************************************************************
   Main_loop: do 
      !
      !-----read *.pos file-----
      read(1,*, iostat=ierror) readstep, time
      !
      !-----Check for end of *.pos file-----
      if (ierror < 0) then 
         write(*,*) ''
         write(*,*) ' End of File Reached' 
         exit
      endif
      !
      !Read in current position values
      do ns=1,ntsp,1
         do na=1,nsp(ns),1
            read(1,*) (tau(j,na,ns),j=1,3)   
         enddo
      enddo
      !
      !-----read *cel file------
      read(2, *) dummy, dummy
      do i=1,3,1
         read(2, *) (aprim(i,j),j=1,3)
      enddo
      !
      !----Check for start condition------
      if (readstep < stepstart) then
         cycle Main_loop
      endif
      !
      !calculate the inverse
      call invert(aprim, apinv, omega)
      !
      !region dimensions for THIS loop
      box  = mic*omega**(1./3.)
      Hbox = box/2.
      !
      !Convert to scaled positions 
      do ns=1,ntsp,1
         do na=1,nsp(ns),1
            call r_to_s( tau(1:3,na,ns), stau(1:3,na,ns), apinv )
         enddo
      enddo
      !
      !Reset the H to O reference array
      Href(:,:) = 0
      !
      !Find Defects
      CALL find_defect(ncount, readstep)
      !
      if (.not. skip) then
        !Find Hbonds
       CALL Find_Hbonds(ncount, readstep)
        !
        !Calculate the rdf for all cases 
        if (calc_rdf)  then
           !
           Call count_pairs(1, ngdr)                                     !Full
           !
           !delta RDFs only calculated for single defects
           if (ndefect == 1) then
              if (delta < delta_below) Call count_pairs(2, ngdr_below)   !delta < 0.1
              if (delta > delta_above) Call count_pairs(3, ngdr_above)   !delta > 0.5
            endif
        endif
         !
         !Update Number of Counts (Only if Skip is not .true.)
         ncount = ncount + 1
      endif
      ! 
      !-----Check for stop condition in *.pos file-------
      if(readstep == stepstop) then
         write(*,*) ''
         write(*,'(2X, "Final step reached : ",2X, I9)') stepstop
         write(*,*) ''
         exit
      endif !
      !
   enddo Main_loop
   !******************************************************************
   !
   !Finalize the rdf for all cases
   if (calc_delta) then
      Call final_rdf(1, ngdr,   ncount-noH3O-dblH3O-multiH3O)     !Full
      Call final_rdf(2, ngdr_below, ndelta_below)                 !delta < delta_below
      Call final_rdf(3, ngdr_above, ndelta_above)                 !delta > delta_above
   endif
   !
   close(1)
   close(2)
   close(13)
   close(150)
   if (calc_rdf) then
      close(161)
      close(162)
      close(160)
   endif
  

   write(*,'(2X,"---Configuration Totals---",/)')
   write(*,'(2X, "Total number of steps: ", I9,/)') (ncount)
   write(*,*) ' Total number of configurations with NO H3Oplus:       ', noH3O
   Write(*,*) ' Total number of configurations with TWO H3Oplus:      ', dblH3O
   write(*,*) ' Total number of configurations with MULTI H3Oplus:    ', multiH3O
   write(*,*) ' Total number of configurations with OH Ions (Skipped):', OHcount
   write(*,*) ''
   if (calc_delta) then
      write(*,'(2X, "Total number of deltas < ",f4.2," :", 29X, I6)') delta_below, ndelta_below 
      write(*,'(2X, "Total number of deltas > ",f4.2," :", 29X, I6)') delta_above, ndelta_above 
      write(*,*)
   endif


   write(*,'(2X,"---Hbond Averages---",/)')
   !
   !Single H3O 
   write(*,'(3X,"Single Defect")')
   write(*,'(3X,"Accepted Bonds: ", f5.3)') dble(accept_single)/(dble(ncount-noH3O-dblH3O-multiH3O))
   write(*,'(3X,"Donated Bonds:  ", f5.3)') dble(donate_single)/(dble(ncount-noH3O-dblH3O-multiH3O))
   write(*,*) ''
   !
   !Proton Transfer Cases (for H3O: noH3O and dblH3O each counts for TWO!!)
   write(*,'(3X,"Proton Transfer")')
   write(*,'(3X,"Accepted Bonds: ", f5.3)') dble(accept_dbl)/(dble(2*dblH3O+2*noH3O))
   write(*,'(3X,"Donated Bonds:  ", f5.3)') dble(donate_dbl)/(dble(2*dblH3O+2*noH3O))
   write(*,*)
   !
   !Total
   write(*,'(3X,"Total Defect Accepted Bonds:  ", f5.3)') dble(accept_total)/dble(ncount+dblH3O+noH3O)
   write(*,'(3X,"Total Defect Donated Bonds:   ", f5.3)') dble(donate_total)/dble(ncount+dblH3O+noH3O)
   write(*,*)


   !$ end_time = OMP_get_wtime()
   !$ total_time  = end_time - start_time
   !$ write(*,'(1X,"---------------------------------")')
   !$ write(*,'(1X, "Total Computing Time: " I6, 1X, "s")') total_time
   !$ write(*,'(1X,"---------------------------------")')
   !$ write(*,*)


 END PROGRAM Hbond_main
