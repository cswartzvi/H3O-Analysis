Module defect_module
   ! 
   use main_variables,        only : DP, numHref
   implicit none
   !save
   !
   !
   real(DP)             :: roh1,roh2,     &  !For 2 H3O, the TRUE distance between the O's and shared H 
                           PTdelta,         &  !For 2 H3O, the difference of the Oxygen, shared Proton distances
                           OO_rdist          !For 2 H3O, the distance between the Oxygens

   integer              :: Ostar,         &  !O* of H3O See Tuckerman
                           prevOstar,     &  !O* of H3O from last Single H3O
                           Ostar2,        &  !O* of H3O 
                           freeH,         &  !Unreferenced Hydrogen atoms
                           Ocount,        &  !O counter
                           ndefect,        &  !double H check, this is to make
                                             !  sure that we are not counting 2
                                             !  H3O atoms 
                           dblH3O,         &  !counter for configurations with TWO H3O
                           noH3O,          &  !counter for configurations with NO H3O
                           multiH3O,       &  !coutner for configurations with 3>= H3O
                           OHcount            !Number of configurations with OH molecules
   !
   !
   contains
      !
      Subroutine init_defect
         !
         use main_variables
         implicit none
         !
         !defect counters
         dblH3O = 0
         noH3O = 0
         multiH3O = 0
         OHcount = 0
         !
         !initialize the prevOstar
         prevOstar = 0 
         !
         return
      End Subroutine init_defect
      !
      Subroutine find_defect(ncount, readstep)
         !
         ! 
         use main_variables,           only : stau, Href, nsp, rdefect, ao, aprim, mic, box, &
                                              Hbond_cut, skip
         use gen_lib,                  only : get_rdist
         ! 
         implicit none
         ! 
         integer,intent(in)   :: readstep          !current production step
         integer,intent(in)   :: ncount            !current number of steps, incremented when completed
         ! 
         real(DP)             :: rdist(3),      &  !square of the components distance in real coordinates
                                 r2                !total squared distance
         !
         integer              :: H1,H2,H3,H4,   &  !Counters for the Hydrogen around O
                                 i,j,k             !general index
         !
         !
         !
         !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
         ! Associated Hydrogens with Oxygens 
         !======================================================================
         ! fort.145 is the H3O tagging processing
         ! fort.136 tagged H3O list
         write(145,*) 'Configuration :' , (ncount+1), readstep
         !
         !Reset ALL counters
         ndefect = 0 
         Ostar = 0
         Ostar2 = 0
         !
         !reset skip flag
         skip = .FALSE.
         !
         !Loop over all Oxygens
         Oloop: do i=1,nsp(1),1
            !Hydrogen counters, set to zero for each Oxygen
            H1 = 0
            H2 = 0
            H3 = 0
            H4 = 0
            !
            Hloop: do j=1,nsp(2),1
               !
               !Get the real distance in each coordinate direction
               CALL get_rdist(stau(1:3,i,1),stau(1:3,j,2),rdist,mic,aprim)
               r2 = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               if (r2 < (rdefect/ao)**2 ) then
                  !How many Hydrogen bonds are already associated with the O*
                  if ( (H1==0) .and. (H2==0) .and. (H3==0) .and. (H4==0) ) then
                     H1 = j 
                     write(145,*) 'Hydrogen1 ', i, j, ' ', (sqrt(r2)*ao)
                  elseif ( (H2==0) .and. (H3==0) .and. (H4==0) ) then
                     H2 = j 
                     write(145,*) 'Hydrogen2 ', i, j, ' ', (sqrt(r2)*ao)
                  elseif ( (H3==0) .and. (H4==0) ) then
                     H3 = j 
                     write(145,*) 'Hydrogen3 ', i, j, ' ', (sqrt(r2)*ao)
                  elseif ( (H4==0) )then !Error
                     !A warning message will be generated below
                     H4 = j 
                     write(145,*) 'Hydrogen4', i, j, ' ', (sqrt(r2)*ao)
                  endif
               endif
               !
            enddo Hloop
            !
            !---------------------------------------------------
            ! Check to see how many H atoms are associated 
            ! with this Oxygen
            !
            !  1 H (Error)    : H1  = 0 H2  =0 H3  =0 H4  =0 
            !  OH             : H1 /= 0 H2  =0 H3  =0 H4  =0
            !  Normal Water   : H1 /= 0 H2 /=0 H3  =0 H4  =0 
            !  H3O            : H1 /= 0 H2 /=0 H3 /=0 H4  =0
            !  4+H (Error)    : H1 /= 0 H2 /=0 H3 /=0 H4 /=0
            !
            !---------------------------------------------------
            !1-H (Error): Check for a Single O Atom (Error!!)
            if ( (H1 == 0) .and. (H2 == 0) .and. (H3 == 0) .and. (H4 == 0) ) then
               !
               write(*,*) ' ERROR: NO Hydrogens within Rcut:', readstep
               skip = .TRUE.
               !
            !OH: Check for a Single OH molecule (Double PT)
            elseif ( (H1 /= 0) .and. (H2 == 0) .and. (H3 == 0) .and. (H4 == 0) ) then
               !
               OHcount = OHcount + 1 
               write(*,*) ' WARNING: OH molecule (skipping configuration for now) :', OHcount, readstep
               skip = .TRUE.
               !
            !H2O: Check for a Single H2O molecule 
            elseif ( (H1 /= 0) .and. (H2 /= 0) .and. (H3 == 0) .and. (H4 == 0) ) then
               !
               !Set reference array
               Href(i,1:2) = (/ H1, H2 /)
               !
            !H3Check for a singular H3O
            elseif ((H1 /= 0) .and. (H2 /= 0) .and. (H3 /= 0) .and. (H4 == 0)  ) then
               !
               write (145, *) 'H3O found!!!'
               !
               !Set Ostar
               Ostar = i
               !
               !Set reference array
               Href(i,1:3) = (/ H1, H2, H3 /)
               !
               ndefect = ndefect + 1
               !
               !Multiple H3O check, provides storage for second H3O
               if (ndefect == 1) then
                  Ostar2 = i
               endif
               !
            elseif ((H1 /= 0) .and. (H2 /= 0) .and. (H3 /= 0) .and. (H4 /= 0)) then 
               !
               write(*,*) ' ERROR: Four Hydrogens (or more) within Rcut :', readstep
               skip = .TRUE.
               !
            endif
            !
         enddo Oloop
         !======================================================================
         !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         !
         !
         !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
         ! Check the number of the H3O present (using ndefect)
         !======================================================================
         if( .not. skip) then
            !
            !--------------------------------------
            ! No H3O Proton Transfer (no Ostar's defined yet)
            if ( ndefect == 0 ) then
               !Important with H3O, the "normal" transfer is associated with the NO H3O!!!
               !
               write(136,'(I9, " H3O not found (PT in Progress)")') readstep
               noH3O = noH3O + 1
               !
               !
               !Find the Free Hydrogen
               Call find_freeH(1, readstep, freeH)
               !
               !Find the two atoms participating in the proton transfer
               !Determine the Difference in Distance between free (shared) proton 
               Call calc_PTdelta(readstep, Ostar, Ostar2, PTdelta, roh1, roh2, freeH) 
               !
               !Print out absolute index of O* during a PT
               write(198,'(I10,I5, I5, 5X,f6.3, 5X, 2(f6.3))') readstep, Ostar,  Ostar2, PTdelta, roh1, roh2
               !               
            !--------------------------------------
            !Only one H3O (One Ostar defined)
            elseif (ndefect == 1) then
               !
               write (136, '(I9 " Ostar: " I7, " H1: ", I5, " H2: ", I5, " H3: ", I5)') &
                     & readstep, Ostar, Href(Ostar,1), Href(Ostar,2), Href(Ostar,3)
               !
               !Print out absolute index of O* 
               write(199,'(I10, I5)') readstep, Ostar
               !
            !--------------------------------------
            !Two H3O (One Ostar and Ostar2)
            elseif ( ndefect == 2 )  then
               !
               !Alternative form of a proton transfer: Hydrogen on the edge of
               !two covalent radii
               !
               !If there are two H3O in the configuration
               write (136, '(I9, " Ostar: " I7, " H1: ", I5, " H2: ", I5, " H3: ", I5, " (PT in progress)")')&
                     readstep, Ostar, Href(Ostar,1), Href(Ostar,2), Href(Ostar,3)
               write (136, '(I9, " Ostar2: " I6, " H1: ", I5, " H2: ", I5, " H3: ", I5, " (PT in progress)")')&
                     readstep, Ostar2, Href(Ostar2,1), Href(Ostar2,2), Href(Ostar2,3)
               !
               dblH3O = dblH3O + 1
               !
               !Find the Shared Hydrogen
               Call find_freeH(2, readstep, freeH)
               !
               !Find the two atoms participating in the proton transfer
               !Determine the Difference in Distance between shared proton 
               Call calc_PTdelta(readstep, Ostar, Ostar2, PTdelta, roh1, roh2, freeH) 
               !
               !Print out absolute index of O* during a PT
               write(198,'(I10,I5, I5, 5X,f6.3, 5X, 2(f6.3))') readstep, Ostar,  Ostar2, PTdelta, roh1, roh2
               !
               !More then two H3O, NOT GOOD!
            elseif ( ndefect > 2 ) then
               !
               !If three H3O are found
               write(136,*) readstep, '3+ H3O in configuration !!!'
               !
               multiH3O = multiH3O+1
               !
            endif
            !--------------------------------------
            !
         else
            if (H1 == 0) then 
               !
               write(136,*) readstep, 'SKIPPED due to errors: Bare Oxygen'
               !
            elseif (H1 /= 0 .and. H2 == 0) then
               !
               write(136,*) readstep, 'SKIPPED due to errors: OH Ion'
               !
            elseif( H4 /= 0) then
               !
               write(136,*) readstep, 'SKIPPED due to errors: H4O Ion'
               !
            endif
         endif
         !======================================================================
         !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         !
         return
         !
   End Subroutine find_defect
   !
   !
   !
   Subroutine find_freeH(val, readstep, Hat)
      !Finds the Hydrogen that is NOT associated with an Oxygen molecule OR,
      !alternatively, finds the Hydrogen that is shared (in the middle of the
      !covalent radii) of two O atoms along with those two atoms 
      !
      use main_variables,           only : DP
      use main_variables,           only : Href, nsp ,nat
      !
      implicit none
      !
      integer, intent(in)              :: val      !val=1 Finds free Hydrogen
                                                   !  between H2Os
                                                   !val=2 Finds Shared Hydrogen
                                                   !  between H3O
      integer, intent(in)              :: readstep !current step number
      integer, intent(inout)           :: Hat      !Unassociated Hydrogen index (local)
      !
      integer                       :: i,j,m,n     !indexes
      !
      select case (val)
         case (1)
            !No defect found, only one unassociated Hydrogen atom (Normal Proton Transfer)
            !
            do i=1,nsp(2),1
               if (All(Href /= i)) then
                  Hat = i
                  return
               endif
            enddo
            !
            write (*,*) ' Error: no unassociated Hydrogen Found', readstep
            !
            return
         case(2)
            !Two defects found, Hydrogen is in the middle of the of covalent radii 
            !
            !Loop over the number of Oxygens
            do i=1,(nsp(1)-1)
               do j=(i+1),nsp(1)
                  !
                  !Loop Over possible H
                  do m=1,numHref,1
                     do n=1,numHref,1
                        !
                        if ( (Href(i,m) == Href(j,n)) .and. (Href(i,m) /= 0) ) then
                           Hat = Href(i,m)
                        endif
                        !
                     enddo
                  enddo
                  !
               enddo
            enddo

         case default
      end select 
      !
   End Subroutine find_freeH

   Subroutine calc_PTdelta(readstep, O, O2, del, r1, r2, Hat)
      !
      Use main_variables,           only : DP
      Use main_variables,           only : stau, nsp, Href
      Use main_variables,           only : box, mic, aprim, Hbond_cut, ao
      use gen_lib,                  only : get_rdist
      !
      !
      implicit none
      !
      !
      integer, intent(in)           :: readstep             !Current Step Number
      integer,intent(out)           :: O, O2                !Ostar and Ostar2 index number
      real(DP), intent(inout)       :: del,              &  !PTdelta between to r1 and r2
                                       r1,r2                !distances
      integer, intent(in)           :: Hat                  !Unassociated Hydrogen index (local)
      !
      integer                       :: i,j,m,n
      !
      real(DP)                      :: r1t, r2t,         &  !temp distances
                                       rr,               &  !distance sqaured
                                       rdist(3)             !component
      !
      !Reset the temp variables
      r1   = box
      r2   = box
      !
      do i=1,nsp(1),1
         !
         !Check to see what two O atoms are closet to the free (shared) Hydrogen
         CALL get_rdist(stau(1:3,Hat,2),stau(1:3,i,1),rdist,mic,aprim)
         rr = sqrt(rdist(1)**2 + rdist(2)**2 + rdist(3)**2)*ao
         !
         !If inside the H2O-H2O hydrogen bond sphere (elminate some extra work)
         if ( rr <= Hbond_cut ) then
            !
            !Find the 2 lowest distance ( r1 lower then r2)
            if (rr < r2) then
               if (rr < r1) then 
                  !
                  !Store (shift) old r1's
                  O2 = O
                  r2 = r1
                  !
                  !Update new
                  O = i
                  r1 = rr
               else
                  O2 = i
                  r2 = rr
               endif
            endif

         endif
      enddo
      !
      !Set delta for the smallest distances 
      del = abs(r1-r2) 
      !
      return
      !
   End Subroutine calc_PTdelta

End Module defect_module
