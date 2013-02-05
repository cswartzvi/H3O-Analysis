Module Hbonds_module
   !
   use main_variables,        only : DP, shell_max, numHref, Hdefect_num, calc_delta, delta_donate
   implicit none        
   !
   integer                 :: Oshell_count,  &  !Total Number of Oxygen mol inside the Hbond_cut radius
                                                !Warning: This establishes a NEW index set for the 
                                                !Oxygens inside Oshell
                              donate_single, &  !Donated Hbonds, Single H3O
                              donate_dbl,    &  !Donated Hbonds, Double H3O
                              donate_total,  &  !total donated Hbonds
                              accept_single, &  !Accepted Hbonds, Single H3O
                              accept_dbl,    &  !Accepted Hbonds, Double H3O
                              accept_total,  &  !total accepted Hbonds
                              ndelta_below,  &  !Number of configurations with delta < delta_below
                              ndelta_above      !Number of configurations with delta > delta_above
   !
   real(DP),allocatable    :: sOshell(:,:),  &  !Oxygens inside the Hbond_cut length (dim, Oshell-index) (assume index < 20)
                              sHshell(:,:,:)    !Hydrogens of Oxygens inside the HBond_cut length (dim, Oshell-index, H#)
                                                !Warning: The indexes used above are for the Oshell index within the Oshell radius
   !
   real(DP)                :: delta             !distance between atoms !  Contains
   !
   Contains
      !
      !Initialize the hbond variables
      Subroutine init_hbonds
         !
         use main_variables,        only : box
         !
         implicit none
         !
         !Hbonds counters
         donate_single = 0
         donate_dbl = 0
         donate_total = 0
         accept_single = 0
         accept_dbl = 0
         accept_total = 0 
         !
         !Number of deltas (delta < 0.1 and delta > 0.5)
         ndelta_below = 0
         ndelta_above = 0
         delta = box
         !
         allocate( sOshell(3, shell_max)               )
         allocate( sHshell(3, shell_max,numHref)             )
         !
         sOshell(:,:) = 0.0
         sHshell(:,:,:) = 0.0
         !
         return
      End Subroutine init_hbonds
      !
      Subroutine Find_Hbonds(ncount, readstep)
      !
      use main_variables
      use defect_module,        only : Ostar, Ostar2, ndefect, freeH
      use gen_lib
      !
      !
      implicit none
      !
      integer,intent(in)   :: readstep          !current production step
      integer,intent(in)   :: ncount            !current number of steps, incremented when completed
      !
      real(DP)             :: rdist(3),      &  !square of the components distance in real coordinates
                              r2                !squared distanced
      !
      real(DP)             :: r2AB, r2AC, r2BC,   &    !for the law of cosines variables
                              cos_angle
      !
      integer              :: na,ns,no,            &  !index: atoms, species, oxygen-shell
                              Oshell_count_Ostar,  &  !For 2 H3O, in the term of the number of
                                                      !Oxygens in the Hbond_cut radius, which is the Ostar
                              temp_donate,         &  !local donated Hbonds, for each configuration
                              temp_accept,         &  !local accepted Hbonds, for each configuration
                              Ostar_x,             &  !executable Oxygen index
                              i,j,k    
      !
      !
      !Initialize the delta distance for this configuration
      Call Find_distance_delta(1, ncount, readstep)
      !
      !
      !----------------------------------------------------------------------------- 
      !Single defect
      !----------------------------------------------------------------------------- 
      if (ndefect == 1)  then
         !
         temp_donate = 0
         temp_accept = 0
         !
         !Reset the Oshell_count
         Oshell_count = 0
         sOshell(:,:) = 0.0_DP
         sHshell(:,:,:) = 0.0_DP
         !
         !Loop over all Oxygen atoms (cycle when Ostar) find
         !All atoms within Hbond length
         Oloop: do na=1,nsp(1),1
            !
            if (na == Ostar) cycle
            !
            !Oxygen-Oxygen Distance
            CALL get_rdist(stau(1:3,Ostar,1),stau(1:3,na,1),rdist,mic,aprim)
            r2 = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
            !
            if (r2 < (Hbond_cut/ao)**2 ) then 
               ! 
               !Current, running, number of of atoms in first Oshell 
               Oshell_count = Oshell_count + 1
               !
               !Create sOshell
               sOshell(:,Oshell_count) = stau(1:3,na,1)
               !
               !Create sHshell, using Href(O-index, H#)
               !because only ONE defect, found no possible zeros in Href(index, 1:numHref)
               Hloop: do j=1,numHref,1
                  !
                  !cycle if Href is 0 (no Hydrogen association) 
                  if (Href(na,j) == 0) cycle Hloop 
                  !
                  !Assign Hydrogen to the sHshell array
                  sHshell(:,Oshell_count,j) = stau(1:3,Href(na,j),2)
                  !
               enddo Hloop
               ! 
            endif
            !
         enddo Oloop
         !
         write(13,*) readstep, Oshell_count
         !
         !------------------------------------------
         !Loop over all donated H bonds (for each Hprime), check angle
         !write to coord-num.dat-> readstep -- accepted-Hbonds donated-Hbonds
         !r2AB = OstarHprime
         !r2BC = HprimeOw
         !r2AC = OstarOw
         Hprime: do i=1,Hdefect_num,1
            !
            CALL  get_rdist(stau(1:3,Ostar,1),stau(1:3,Href(Ostar,i),2),rdist,mic,aprim)
            r2AB = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
            do j=1,Oshell_count,1
               !
               CALL  get_rdist(stau(1:3,Href(Ostar,i),2),sOshell(1:3,j),rdist,mic,aprim)
               r2BC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               CALL  get_rdist(stau(1:3,Ostar,1),sOshell(1:3,j),rdist,mic,aprim)
               r2AC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               cos_angle = ( (r2AB + r2BC) - r2AC )/(2.0*sqrt(r2AB*r2BC))
               ! 
               !confirm Hbond
               if (cos_angle <= ref_cos) then
                  temp_donate = temp_donate + 1
                  !
                  !Calc deltas 
                  if (calc_delta .and. delta_donate) then
                     Call Find_distance_delta(2, ncount, readstep, Oat1=stau(1:3,Ostar,1), &
                        & Hat=stau(1:3,Href(Ostar,i),2), Oat2=sOshell(1:3,j) )
                  endif
                  !
               endif  
               !
            enddo
         enddo Hprime
         !------------------------------------------
         ! 
         ! 
         !------------------------------------------
         !Loop over all  accepted H bonds, check angle
         !write to coord-num.dat -> readstep -- accepted-Hbonds donated-Hbonds
         !r2AB = OstarHw
         !r2BC = HwOw
         !r2AC = OstarOw
         do i=1,Oshell_count,1
            do j=1,2,1  !loop over each Hydrogen of normal water molecules
               ! 
               CALL  get_rdist(stau(1:3,Ostar,1),sHshell(1:3,i,j),rdist,mic,aprim)
               r2AB = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               ! 
               CALL  get_rdist(sHshell(1:3,i,j),sOshell(1:3,i),rdist,mic,aprim)
               r2BC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               CALL  get_rdist(stau(1:3,Ostar,1),sOshell(1:3,i),rdist,mic,aprim)
               r2AC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               cos_angle = ( (r2AB + r2BC) - r2AC )/(2.0*sqrt(r2AB*r2BC))
               ! 
               if (cos_angle <= ref_cos) then
                  temp_accept = temp_accept + 1
                  !
                  !Calc deltas 
                  if (calc_delta .and. .not. delta_donate) then
                     Call Find_distance_delta(2, ncount, readstep, Oat1=stau(1:3,Ostar,1), &
                        & Hat=sHshell(1:3,i,j), Oat2=sOshell(1:3,i) )
                  endif
                  !
               endif  
               !
            enddo
         enddo
         !------------------------------------------
         !         
         !Finalize the delta distance for this configuration
         if (calc_delta) Call Find_distance_delta(3, ncount, readstep)
         ! 
         write(177,*) readstep, '---', temp_donate, temp_accept, '|', Ostar
         ! 
         donate_single   =  donate_single + temp_donate
         accept_single   =  accept_single + temp_accept
         donate_total    =  donate_total + temp_donate
         accept_total    =  accept_total + temp_accept
         !
      endif
!----------------------------------------------------------------------------- 
!
!
!----------------------------------------------------------------------------- 
! Proton Transfer (Similar to the above but with Ostar and Ostar2)
!----------------------------------------------------------------------------- 
      if (ndefect == 0 .or. ndefect == 2) then 
         !Two cases of Proton transfer
         ! 1- No H3O, Hydrogen in the middle of a transfer
         ! 2- Two H3O, Hydrogen on the edge of BOTH covalent radii
         !         
         do no=1,2,1
            !
            temp_donate = 0
            temp_accept = 0
            !
            select case (no)
               case(1)
                  Ostar_x = Ostar
               case(2) 
                  Ostar_x = Ostar2
               case default
            end select
            !
            !Reset the Oshell_count
            Oshell_count = 0
            sOshell(:,:) = 0.0_DP
            sHshell(:,:,:) = 0.0_DP
            !
            !Loop over all Oxygen atoms (cycle when Ostar) find
            !All atoms within Hbond length
            Oloop2: do na=1,nsp(1),1
               !
               !cycle if we are looking at the current ]3O 
               if (na == Ostar_x) cycle Oloop2
               !
               !Oxygen-Oxygen Distance
               CALL get_rdist(stau(1:3,Ostar_x,1),stau(1:3,na,1),rdist,mic,aprim)
               r2 = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               if (r2 < (Hbond_cut/ao)**2 ) then 
                  !    
                  Oshell_count = Oshell_count + 1
                  !
                  !Create sOshell
                  sOshell(:,Oshell_count) = stau(1:3,na,1)
                  !
                  !Create sHshell, using Href(Oshell-index, H#)
                  Hloop2: do j=1,numHref,1
                     !
                     if (Href(na,j) == 0 ) cycle Hloop2
                     !
                     sHshell(:,Oshell_count,j) = stau(1:3,Href(na,j),2)
                     !
                  enddo Hloop2
                  !
                  !In this outer do loop we are looping over the TWO defects.
                  !While we are considering one of the defects the other only has
                  !has to be tagged in this new Oshell-index
                  if (Ostar_x == Ostar .and. na == Ostar2) then
                     Oshell_count_Ostar = Oshell_count
                  elseif (Ostar_x == Ostar2 .and. na == Ostar) then 
                     Oshell_count_Ostar = Oshell_count
                  endif
                  !
               endif
               !
            enddo Oloop2
            ! 
            write(13,*) readstep, Oshell_count
            !
            !------------------------------------------
            !Loop over all donated H bonds, check angle
            !write to coord-num.dat -> readstep -- accepted-Hbonds donated-Hbonds
            !r2AB = OstarHprime
            !r2BC = HprimeOw
            !r2AC = OstarOw
            Hprime2: do i=1,Hdefect_num,1
               !
               CALL  get_rdist(stau(1:3,Ostar_x,1),stau(1:3,Href(Ostar_x,i),2),rdist,mic,aprim)
               r2AB = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               do j=1,Oshell_count,1
                  !
                  CALL  get_rdist(stau(1:3,Href(Ostar_x,i),2),sOshell(1:3,j),rdist,mic,aprim)
                  r2BC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
                  !
                  CALL  get_rdist(stau(1:3,Ostar_x,1),sOshell(1:3,j),rdist,mic,aprim)
                  r2AC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
                  !
                  cos_angle = ( (r2AB + r2BC) - r2AC )/(2.0*sqrt(r2AB*r2BC))
                  ! 
                  if (cos_angle <= ref_cos) then
                     temp_donate = temp_donate + 1
                  endif  
                  !
               enddo
            enddo Hprime2
            !------------------------------------------
            ! 
            !
            !------------------------------------------
            !Loop over all accepted H bonds, check angle
            !write to coord-num.dat -> readstep -- accepted-Hbonds donated-Hbonds
            !r2AB = OstarHw
            !r2BC = HwOw
            !r2AC = OstarOw
            do i=1,Oshell_count,1
               do j=1,2,1  !loop over each Hydrogen 
                  !
                  !Because we are dealing with multiple defects (in the Oshell
                  !index) we need to be sure to skip 
                  if (j == Hdefect_num .and. i == Oshell_count_Ostar) cycle
                  ! 
                  CALL  get_rdist(stau(1:3,Ostar_x,1),sHshell(1:3,i,j),rdist,mic,aprim)
                  r2AB = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
                  !
                  CALL  get_rdist(sHshell(1:3,i,j),sOshell(1:3,i),rdist,mic,aprim)
                  r2BC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
                  !
                  CALL  get_rdist(stau(1:3,Ostar_x,1),sOshell(1:3,i),rdist,mic,aprim)
                  r2AC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
                  !
                  cos_angle = ( (r2AB + r2BC) - r2AC )/(2.0*sqrt(r2AB*r2BC))
                  ! 
                  if (cos_angle <= ref_cos) then
                     temp_accept = temp_accept + 1
                  endif  
                  !
               enddo
            enddo
            !------------------------------------------
            !
            !
#if 0
            !------------------------------------------
            !for the free Hydrogen Atom and the other Ostar on the second defect
            !use the Oshell_count_Ostar to locate this other Ostar within the 
            !Oshell-index--> Considered a DONATED Bond for EACH Ostar
            !r2AB = Ostar_x-freeH
            !r2BC = freeH-Oshell_count_Ostar
            !r2AC = Ostar_x-Oshell_count_Ostar
            CALL  get_rdist(stau(1:3,Ostar_x,1),stau(1:3,freeH,2),rdist,mic,aprim)
            r2AB = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
            !
            CALL  get_rdist(stau(1:3,freeH,2),sOshell(1:3,Oshell_count_Ostar),rdist,mic,aprim)
            r2BC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
            !
            CALL  get_rdist(stau(1:3,Ostar_x,1),sOshell(1:3,Oshell_count_Ostar),rdist,mic,aprim)
            r2AC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
            !
            cos_angle = ( (r2AB + r2BC) - r2AC )/(2.0*sqrt(r2AB*r2BC))
            ! 
            if (cos_angle <= ref_cos) then
               temp_donate = temp_donate + 1
            endif  
            !------------------------------------------
#ndif
            !
            !
            select case (no)
               case(1)
                  write(177,*) readstep, 'PT1', temp_donate, temp_accept, '|', Ostar, Ostar2
               case(2) 
                  write(177,*) readstep, 'PT2', temp_donate, temp_accept, '|', Ostar, Ostar2
               case default
            end select
            ! 
            donate_dbl   =  donate_dbl + temp_donate
            accept_dbl   =  accept_dbl + temp_accept
            donate_total =  donate_total + temp_donate
            accept_total =  accept_total + temp_accept
            !
         enddo
         !        
      endif
      !----------------------------------------------------------------------------- 
      !
      return

   End Subroutine Find_Hbonds
   !
   !This Subroutine is run after Hbonds have been found but ONLY for ndefect = 1
   Subroutine Find_distance_delta(cal, ncount, readstep, Oat1, Oat2, Hat)
      !
      use main_variables,        only : stau, box, aprim, mic, ao
      use main_variables,        only : delta_above, delta_below
      use defect_module,         only : Ostar, ndefect
      use gen_lib
      !
      implicit none
      !
      integer,intent(in)   :: cal               !Calculation type:
                                                ! 1-Initialize the deltas ! (Beginning of find Hbonds) 
                                                ! 2-Delta calculation
                                                ! 3-Finalize the calculation
      integer,intent(in)   :: readstep          !current production step
      integer,intent(in)   :: ncount            !current number of steps, incremented when completed
      real(DP),intent(inout), optional    :: Oat1(3),Oat2(3),Hat(3)  !Coordinates of the Oxygens and 
                                                                  !Hydrogens involved in PT transfer
      ! 
      integer              :: i,j
      !
      !save these values for the next loop through
      real(DP),save        :: roh1, roh2        !FINAL O to H distances for those molecules hydrogen bonded to the H3O
      !
      real(DP)             :: troh1, troh2,  &  !TEMP O to H distances for those molecules hydrogen bonded to the H3O 
                              rdist(3)          !component distance
      !
      if (cal == 1) then
         !reset the delta counter, and distances
         delta = box
         roh1 = box
         roh2 = box
      endif
      !
      if (cal ==2 ) then
         !
         CALL get_rdist(Oat1(1:3),Hat(1:3),rdist,mic,aprim)              
         troh1 = sqrt(rdist(1)**2 + rdist(2)**2 + rdist(3)**2)*ao                                     
         !
         CALL get_rdist(Oat2(1:3),Hat(1:3),rdist,mic,aprim)                   
         troh2 = sqrt(rdist(1)**2 + rdist(2)**2 + rdist(3)**2)*ao
         !
         !
         if ( abs(troh1 - troh2) < delta ) then 
            delta = abs(troh1 - troh2)
            roh1 = troh1
            roh2 = troh2
         endif
         !
         !
      endif
      !
      if (cal == 3) then
         write(150,'(1X,I9,5X,f7.4, " | ", 2(f7.4, 5X))') readstep, delta, roh1, roh2 
         !
         !Counters for special cases of delta
         if (delta < delta_below) then
            ndelta_below = ndelta_below + 1
         elseif (delta > delta_above) then
            ndelta_above = ndelta_above + 1
         endif
      endif
      !
      return
      !
   End Subroutine find_distance_delta
   !
END Module Hbonds_module
