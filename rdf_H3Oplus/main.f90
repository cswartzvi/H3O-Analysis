!==================================================================================================
!
!
!atom1 can only be O*
!
!
!Output:
!fort.36 = Tagged H3O+ (includes singular proton trsanfer)
!fort.38 = More than one H3O+ list
!fort.45 = Tagging process
!fort.98 = Absolute index of proton transfers
!fort.99 = Absolute index of Singular H3O+
!
!
!==================================================================================================
 PROGRAM new_rdf
   !
   USE rdf_lib
   USE omp_lib
   !
   implicit none
   !
   integer, parameter   :: namax = 20                 !max number of species   
   real(DP), parameter  :: pi    = 4.d0*atan(1.d0),& !Pi
                           ao    = 0.52917720859_DP
   !
   integer              :: binNum,     &  !number of bins for g(r)
                           ntsp,       &  !number of total atomic species, limit namax
                           nsp(namax), &  !number of each atomic species     
                           nat=0,      &  !total number of atoms
                           mic,        &  !Multiple Supercells 
                                             !  mic = 1, Only the main supercell
                                             !  mic = 2, translated in the pos direction
                                             !  mic = 3, translated in both the pos and neg directions
                                             !  Also, scaled length of region
                           mic3,          &  !number of extended cells, mic**3 (line )
                           stepstart,     &  !nfi that we will start reading *.pos 
                           stepstop,      &  !nfi that we will stop reading *.pos AFTER (optional)
                           interval,      &  !interval between nfi in *pos file
                           nsteps,        &  !total number of steps
                           ncount=0,      &  !number of steps index
                           i, j, k,       &  !general indexes
                           offsetI,       &  !initial offset
                           offsetF,       &  !final offset
                           atom1,         &  !first atom type (IMPORTANT: atom1 with always be the O* or the H'
                                             !  See Tuckerman Acc. Chem. Res.
                                             !  2006, 39, pg 151-158 for more
                                             !  information
                           atom2,         &  !first atom type
                           n,             &  !general index (extended supercell)
                           m,             &  !secondary general index
                           ns,            &  !number of species index
                           na,            &  !numbnfier of atoms index
                           readstep,      &  !read the nfi from *.pos 
                           ierror,        &  !error index
                           nbin,          &  !bin index
                           pcount_end,    &  !first limit used in the pair counting process
                           pcount_start,  &  !second limit used in the pair counting process
                           Ostar=0,       &  !O* of H3O+ 
                           Ostar2=0,       &  !O* of H3O+ 
                           Hprime1=0,     &  !H'1 of H3O+ 
                           Hprime2=0,     &  !H'2 of H3O+ 
                           Hprime3=0,     &  !H'3 of H3O+ 
                           Hprime21=0,     &  !H'1 of H3O+ (Shared Proton) 
                           Hprime22=0,     &  !H'2 of H3O+ (Shared Proton)
                           Hprime23=0,     &  !H'3 of H3O+ (Shared Proton)
                           Ocount,        &  !O counter
                           Hcheck,        &  !double H check, this is to make
                                             !  sure that we are not counting more then one
                                             !  H3o+ ions
                           dblH3O=0,      &  !counter for configurations with TWO H3o+ 
                           noH3O=0,       &  !counter for configurations with NO H3o+
                           multiH3O=0,    &  !coutner for configurations with 2> H3o+
                           H1,            &  !H1 counter
                           H2,            &  !H2 counter
                           H3,            &  !H3 counter
                           H4                !H4 counter (used to make sure there are not more then 3 asscoaited)
                           
   !
   real(DP), allocatable    :: tau(:,:,:),    &  !atomic positions (dim, index, atomic-species)
                               stau(:,:,:),   &  !scaled atomic positions (dim,index, atomic-species)
                               staux(:,:,:),  &  !multiple scaled atomic positions (dim, index, atomic-species)
                               ngdr(:),       &  !main pair correlation (unnormalized g(r)) array
                               ngdrt(:)          !Temp pair correlation (unnormalized g(r)) array
   !
   real(DP)             :: gofr,          &  !FINAL g(r) print-out variable (not an array, print-out values)  
                           r,             &  !total distance (not an array, print-out value)
                           aprim(3,3),    &  !lattice prim vectors
                           apinv(3,3),    &  !inverse of prim vectors, for scaled positions
                           box,           &  !length of the region (includes multiple supercells)
                           Hbox,          &  !length of half the region (including multiple supercells)
                                             ! Recall: we can only compute up to
                                             ! L/2 of the total region, which is
                                             ! why we need multiple super cells
                           res,           &  !resolution (size) of bins  binNum/Hbox (number/length)
                           omega,         &  !volume of cell
                           omegaT=0.0,    &  !total volume of the cell
                           dens,          &  !raw particle density
                           dx,            &  !diff in x-value
                           dy,            &  !diff in y-value
                           dz,            &  !diff in z-value
                           sdist(3),      &  !square of the components distance in scaled coordinates
                           rdist(3),      &  !square of the components distance in real coordinates
                           r2,            &  !total squared distance
                           norm=1.0,      &  !normalizing factor (used to normalize final g(r), if there are multiple atoms  
                           time,          &  !timestep
                           dmic,          &  !convert mic to real
                           vshell,        &  !volume of the infinitesimal radial shell
                           rcut,          &  !Cut-off radius (in Angstrom!!)
                           dummy             !a dummy variable

   character(len=30)       :: filePos, fileCel, fileOut
   logical                 :: H3Oskip, H3Ofound 
   integer                 :: start_time, end_time, total_time
   !
   namelist /rdfinput/ filePos, fileCel, fileOut, binNum, ntsp, nsp, mic, &
                 & stepstart, stepstop, atom1, atom2, norm, rcut
   !
   !
   !Start Time
   !$ start_time = OMP_get_wtime()   
   !
   !Read in parameters from standard input
   !
 
   read(*, nml=rdfinput)
   !
   !
   !******************************************************************
   !Housekeeping (to be moved to a subroutine)
   !******************************************************************
   !Total number of atoms
   do ns=1,ntsp
      nat = nat + nsp(ns)
   enddo
   !
   !Check for multiple Supercells
   select case (mic)
      case(1)
         offsetI = 0
         offsetF = 0
      case default
         write(*,*) ' ERROR: incorrect mic!! Only mic=1 supported!'
         stop
   end select
   mic3 = mic**3
   !
   !Allocate variables 
   allocate( tau(3,(mic**3*nat),ntsp)    )  
   allocate( stau(3,(mic**3*nat),ntsp)   )
   allocate( staux(3,(mic**3*nat),ntsp)  )
   allocate( ngdr(binNum), ngdrt(binNum) )
   !---------------------------------------------
   !
   !
   !---------------------------------------------
   !I/O
   !---------------------------------------------
   !Open Files
   open(unit=1, file=(trim(filePos)), status='old')
   open(unit=2, file=(trim(fileCel)), status='old')
   open(unit=3, file=(trim(fileOut)), status='unknown')
   !******************************************************************
   !
   !
   !******************************************************************
   Main_loop: do 
      !
      !-----read *.pos file-----
      read(1,*, iostat=ierror) readstep, time
      !
      !-----Check for end of *.pos file-----
      if (ierror < 0) then
         write(*,*) ' End of File Reached'
         write(*,*) ' Total number of Steps: ', ncount
         exit
      endif
      !
      !----Check for start condition------
      if (readstep < stepstart) then
         cycle
      endif
      !
      !-----read *cel file------
      read(2, *) dummy, dummy
      do i=1,3,1
         read(2, *) (aprim(i,j),j=1,3)
      enddo
      !
      !calculate the inverse
      call invert(aprim, apinv, omega)
      omegaT = omegaT + omega
      !
      !region dimensions for THIS loop
      box  = mic*omega**(1./3.)
      Hbox = box/2.
      res  = binNum/Hbox
      !
      !Read in current position values
      do ns=1,ntsp,1
         do na=1,nsp(ns),1
            read(1,*) (tau(j,na,ns),j=1,3)   
         enddo
      enddo
      !
      !Convert to scaled positions 
      do ns=1,ntsp,1
         do na=1,nsp(ns),1
            call r_to_s( tau(1:3,na,ns), stau(1:3,na,ns), apinv )
         enddo
      enddo
      !
      !Zero all temp g(r)
      do n=1,binNum,1
         ngdrt(n) = 0
      enddo
      !
      !ADD OPTIONAL PERIODIC
      !
      !Extended Multiple supercell !!!!! NOT IMPLEMENTED
      do ns=1,ntsp,1
         n=0
         do i=offsetI,offsetF
            do j=offsetI,offsetF
               do k=offsetI,offsetF
                  do na=1,nsp(ns),1
                  n = n +1
                  staux(1,n,ns) = stau(1,na,ns) + i
                  staux(2,n,ns) = stau(2,na,ns) + j
                  staux(3,n,ns) = stau(3,na,ns) + k
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      !
      !---------Find H3O+--------------
      !  fort.45 is the H3O+ tagging processing
      !  fort.36 tagged H3O+ list
      !  fort.38 multiple H3O+ list
      write(45,*) 'Configuration :' , (ncount+1), readstep
      Hcheck = 0 
      Ostar = 0
      Ostar2 = 0
      Hprime1= 0
      Hprime2= 0
      Hprime3= 0
      Hprime21= 0
      Hprime22= 0
      Hprime23= 0
      H3Ofound = .false.
      H3Oskip = .false.
      Oloop: do i=1,nsp(1),1
         !Hydrogen counters
         H1 = 0
         H2 = 0
         H3 = 0
         H4 = 0
         !
         Hloop: do j=1,nsp(2),1
            !
            dx = staux(1,i,1) - staux(1,j,2)
            dy = staux(2,i,1) - staux(2,j,2)
            dz = staux(3,i,1) - staux(3,j,2)      
            !
            !Periodic adjusted distance
            sdist(1) = dx - NINT( dx/(dble(mic)) ) * dble(mic)
            sdist(2) = dy - NINT( dy/(dble(mic)) ) * dble(mic)
            sdist(3) = dz - NINT( dz/(dble(mic)) ) * dble(mic)
            !
            !convert from scaled to real
            call s_to_r( sdist, rdist, aprim )
            r2 = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
            !
            !Check for a H3O+
            if (r2 < (rcut/ao)**2 ) then
               !How many Hydrogen bonds are already assocated with the O*
               if ( (H1==0) .and. (H2==0) .and. (H3==0) .and. (H4==0) ) then
                  H1 = j
                  write(45,*) 'Hydrogen1 ', i, j, ' ', (sqrt(r2)*ao)
               elseif ( (H2==0) .and. (H3==0) .and. (H4==0) ) then
                  H2 = j
                  write(45,*) 'Hydrogen2 ', i, j, ' ', (sqrt(r2)*ao)
               elseif ( (H3==0) .and. (H4==0) ) then
                  H3 = j
                  write(45,*) 'Hydrogen3 ', i, j, ' ', (sqrt(r2)*ao)
               elseif ( (H4==0) )then
                  !A warning message will be generated below
                  H4 = j
                  write(45,*) 'Hydrogen4', i, j, ' ', (sqrt(r2)*ao)
                  write(*,*) ' WARNING: Four Hydrogens within Rcut :', readstep
               endif
            endif
            !
         enddo Hloop
         !
         !Check for a singular H3O+
         if ( (H1 /= 0) .and. (H2 /= 0) .and. (H3 /= 0) .and. (H4 == 0) ) then
            !
            write (45, *) 'H3O+ found!!!'
            Ostar = i
            Hprime1 = H1
            Hprime2 = H2
            Hprime3 = H3
            H3Ofound = .true.
            !
            Hcheck = Hcheck + 1
            !
            !Multiple H3O+ check
            if (Hcheck == 1) then
               Ostar2   = i
               Hprime21 = H1 
               Hprime22 = H2
               Hprime23 = H3
            endif
            !
         endif
         !
      enddo Oloop
      !
      !--------------------------------------------------------------
      !Check the number of the H3O+ present
      !--------------------------------------------------------------
      if ( Hcheck == 0 ) then
         !
         write(36,'(I9, " H3O+ not found !!!")')readstep 
         ! 
         noH3O = noH3O + 1
         !
         write(38,*) '0 H3O+ in configuration ', readstep 
         !
      !Only one H3O+ 
      elseif (Hcheck == 1) then
         write (36, '(I9 " Ostar: " I7, " H1: ", I5, " H2: ", I5, " H3: ", I5)') &
               readstep, Ostar, Hprime1, Hprime2, Hprime3
      !
      !Two H3O+
      elseif ( Hcheck == 2 )  then
         !If there are two H3O+ in the configuration, Proton Transfer 
         write (36, '(I9, " Ostar: " I7, " H1: ", I5, " H2: ", I5, " H3: ", I5, " (PT in progress)")')&
                   readstep, Ostar, Hprime1, Hprime2, Hprime3
         write (36, '(I9, " Ostar2: " I6, " H1: ", I5, " H2: ", I5, " H3: ", I5, " (PT in progress)")')& 
                   readstep, Ostar2, Hprime21, Hprime22, Hprime23
         write(38, *) '2 H3O+ in configuration', readstep
         dblH3O = dblH3O + 1
         H3Oskip = .true.
         !
         !create fort.98
         write(98,*) readstep, Ostar, Ostar2
         !
         !
         !More then two H3O+, NOT GOOD!
      elseif ( Hcheck > 2 ) then
         !If three H3O+ are found, 
         write(36,'(I9, " 3+ H3O+ in configuration !!!")') readstep
         write(38,*) '3+ H3O+ in configuration', readstep, '!!!'
         multiH3O = multiH3O+1
         H3Oskip = .true.  
      endif
      !--------------------------------------------------------------
      !
      !
      !---------Count PAIRS------------------------------------------
      !Index atoms to be used in the pair correlation
      !NOTE: atom1 is always either O* or H', atom2 are the other atoms Ow, Hw    
      if (atom1 == 1 ) then                                    
         n = Ostar
      elseif (atom1 == 2 ) then
         n = Hprime1
         !write(*,*) ' ERROR: atom1 cannot equal 2 (There are three Hprimes in H3O+)'
         !stop
      else
         Write(*,*) ' ERROR: atom1 is defined incorrectly!!!'
         stop
      endif
      !
      !
      if (H3Ofound .and. .not. H3Oskip)  then
      
            !Print out absolute index of O* and H'
            write(99,*) readstep, Ostar
            !
            do m=1,nsp(atom2),1
               !
               !Skip this cylce if atom2 is Ow and the current atom is the O*
               if ( (atom2 ==1) .and. (m == Ostar) ) cycle
               !coordinate distance
               dx = staux(1,n,atom1) - staux(1,m,atom2)
               dy = staux(2,n,atom1) - staux(2,m,atom2)
               dz = staux(3,n,atom1) - staux(3,m,atom2)
               !
               !Periodic adjusted distance
               sdist(1) = dx - NINT( dx/(dble(mic)) ) * dble(mic)
               sdist(2) = dy - NINT( dy/(dble(mic)) ) * dble(mic)
               sdist(3) = dz - NINT( dz/(dble(mic)) ) * dble(mic)
               !
               !convert from scaled to real
               call s_to_r( sdist, rdist, aprim )
               r2 = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               !establish which bin this distance is in
               nbin = nint(sqrt(r2)*res)
               !
               !if we are at the end of binNum
               if (nbin > binNum) nbin = binNum
               !
               !update the temp g(r)
               ngdrt(nbin) = ngdrt(nbin) + 1 
            enddo

      endif
      !
      !update true unnormalized, raw pair correlation
      do n=1,binNum,1
         ngdr(n) = ngdr(n) + ngdrt(n)
      enddo
      !
      !-----Check for stop condition in *.pos file-------
      if(readstep == stepstop) then
         write(*,*) ''
         write(*,*) ' Final step reached : ', stepstop
         write(*,*) ' Total number of steps: ', ncount+1
         write(*,*) ''
         exit
      endif
      !
      ncount = ncount + 1
      !
   enddo Main_loop
   !
   close(1)
   close(2)
   !******************************************************************
  
   write(*,*) ''
   write(*,*) ' Total number of configurations with NO H3O+:    ', noH3O
   Write(*,*) ' Total number of configurations with TWO H3O+:   ', dblH3O
   write(*,*) ' Total number of configurations with MULTI H3O+: ', multiH3O
   write(*,*) ''


   !Total region dimensions, averaged over all loops (it should be the same as
   ! individual loop region dimensions if the *.cel file is constant)
   omega = omegaT/dble(ncount)
   box = mic*omega**(1./3.)
   Hbox = box/2.0
   res = binNum/ Hbox
   !
   !
   do n=1,(binNum-1),1
      !
      !construct some details of normalization and distance  
      r = dble(n)/res
      r2 = r**2.0
      vshell = 4.0 * pi * r2 / res
      !if(atom1 /= atom2) then
      !   dens = dble(nsp(atom2))/omega
      !else
      !   dens = dble(nsp(atom2) - 1)/omega
      !endif
      dens = dble(nsp(atom2) +1)/omega
      !
      !Calculate the final, normalized, value of g(r)! Please note that the
      !normalization constant (ncount*mic3*norm*nsp(atom1)*dens*vshell)...
      !  ncount*mic3 = number of steps and number of extended shells
      !  norm = external normalization constant (default is 1.0)
      !  nsp(atom1)*dens = number of pairs
      !  vshell = volume of the infinitesimal shell
      gofr = ngdr(n)/ &     
      &( dble(ncount*mic3*norm*dens*vshell) )

      write(3,*) (r*ao), gofr
   enddo


   !$ end_time = OMP_get_wtime()
   !$ total_time  = end_time - start_time
   !$ write(*,'(A, I5, A)') 'Total Computing Time: ', total_time, ' s'




 END PROGRAM new_rdf
