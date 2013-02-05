
MODULE main_variables
   !
   implicit none
   ! 
   integer, parameter   :: DP  = SELECTED_REAL_KIND(15,99) 
   !
   integer, parameter   :: namax = 2            !!!!Has only been tested with water and Ions!!!!
   integer, parameter   :: numHref = 3          !Number of Hydrogens in Href
   integer, parameter   :: Hdefect_num = 3      !number of Hdyreogens on the defect
   integer, parameter   :: shell_max = 20       !max number of Oxygen assumed around the defect
   !
   real(DP), parameter  :: pi    = 4.d0*atan(1.d0),& !Pi
                           ao    = 0.52917720859_DP
   !
   !
   !
   integer              :: ntsp,          &  !number of total atomic species, limit namax
                           nsp(namax),    &  !number of each atomic species     
                           mic=1,         &  !Scaled length of the box region, (for future use)
                           stepstart,     &  !nfi that we will start reading *.pos 
                           stepstop,      &  !nfi that we will stop reading *.pos AFTER (optional)
                           interval,      &  !interval between nfi in *pos file (iprint)
                           binNum,        &  !Total Number of Bins for the g(r)
                           atom1,         &  !Atom 1 in the g(r)
                           atom2,         &  !Atom 2 in the g(r)
                           norm,          &  !external normaization for g(r)
                           nat=0,         &  !total number of atoms
                           ierror            !error index
   !
   integer,allocatable      ::Href(:,:),      &  !Hydrogens assoication with Oxygens (O-index, ref1-ref2-ref3 )
                              ngdr(:),        &  !the final g(r) array ngdr(binNum) for full rdf
                              ngdr_below(:),  &  !the final g(r) array ngdr(binNum) for delta < delta_below
                              ngdr_above(:)      !the final g(r) array ngdr(binNum) for delta > delta_above
   !                        
   ! 
   real(DP), allocatable    :: tau(:,:,:),    &  !atomic positions (dim, index, atomic-species)
                               stau(:,:,:),   &  !scaled atomic positions (dim,index, atomic-species)
                               staux(:,:,:)      !multiple scaled atomic positions (dim, index, atomic-species)
                               
   ! 
   real(DP)             :: aprim(3,3),    &  !lattice prim vectors
                           apinv(3,3),    &  !inverse of prim vectors, for scaled positions
                           box,           &  !length of the region 
                           Hbox,          &  !length of half the region 
                           omega=0.0,     &  !volume of cell
                           rdefect,       &  !Defect cut-off distance (in Angstrom!!)
                           Hbond_cut,     &  !Hbond cut-off distance
                           Hbond_angle,   &  !Hbond Angle cut-off
                           ref_cos,       &  !Cosine of the above Hbond_angle
                           delta_above,   &  !Interval for the deltas above
                           delta_below,   &  !Interval for deltas below
                           dummy             !a dummy variable
   !
   character(len=30)       :: filePos, fileCel, filerdf_root, filerdf_full, &
                              filerdf_above, filerdf_below, dum
   ! 
   logical                 :: calc_delta, delta_donate,calc_rdf, full_rdf
   !
   logical                 :: skip                                   !if .true. Skip this Configuration
                              
   !   
   !
   namelist /system/ filePos, fileCel, ntsp, nsp, stepstart, stepstop, &
                     rdefect, calc_delta, calc_rdf

   namelist /hbonds/  Hbond_cut, Hbond_angle

   namelist /ddelta/ delta_donate, delta_above, delta_below
                  
   Namelist /rdf/ atom1, atom2, binNum, norm, filerdf_root
   ! 
   !
   CONTAINS
      !
      !Initialize the Main Variables, Open fiels and read them
      subroutine init
         !
         implicit none
         !
         integer           :: ns       !speices index
         !
         !
         !defaults
         filePos    = 'cp.pos'
         fileCel    = 'cp.cel'
         ntsp       = 0
         nsp(:)     = 0
         stepstart  = -10
         rdefect    = 1.1655
         calc_delta = .true.
         delta_donate = .true. !For H3O only
         delta_above = 0.5
         delta_below = 0.1
         calc_rdf   = .true.
         !hbonds
         Hbond_cut  = 3.5
         Hbond_angle= 150.0
         !rdf
         filerdf_root   = 'rdf'
         atom1          = 1
         atom2          = 1
         binNum         = 100
         norm           = 1
         !
         !
         !Read in parameters from standard input
         read(*, nml=system)
         !
         !deltas always used in the rdf
         if (calc_rdf) calc_delta = .true.
         !
         read(*, nml=hbonds)
         !
         if (calc_delta) read(*,nml=ddelta)
         !
         if (calc_rdf) read(*,nml=rdf)
         !
         !Check values
         !ADD ERROR CHECK!!! 
         !
         !Open Files
         open(unit=1,   file=(trim(filePos)), status='old')
         open(unit=2,   file=(trim(fileCel)), status='old')
         open(unit=13,  file='nearest_O.dat', status='unknown')
         open(unit=177, file='coord-num.dat', status='unknown')
         if (calc_delta) open(unit=150, file='H3O-deltas.dat', status='unknown')
         if (calc_rdf) then
            write(dum,'(f3.1)') delta_below
            filerdf_below = TRIM(filerdf_root)//'-'//TRIM(dum)//'.dat'
            open(unit=161, file=filerdf_below, status='unknown')
            !
            write(dum,'(f3.1)') delta_above
            filerdf_above = TRIM(filerdf_root)//'-'//TRIM(dum)//'.dat'
            open(unit=162, file=filerdf_above, status='unknown')
            !
            filerdf_full = TRIM(filerdf_root)//'-full.dat'
            open(unit=160, file=filerdf_full, status='unknown')
            !
         endif
         !open(unit=198, file='O-index-PT.txt', status='unknown')
         !open(unit=199, file='O-index-noPT.txt', status='unknown')
         !
         !Total number of atoms
         do ns=1,ntsp
            nat = nat + nsp(ns)
         enddo
         !
         !Set the cos_angle problem
         ref_cos = cos(Hbond_angle*(pi/180.0))
         !
         !Allocate variables 
         allocate( tau(3,(mic**3*nat),ntsp)    )  
         allocate( stau(3,(mic**3*nat),ntsp)   )
         allocate( staux(3,(mic**3*nat),ntsp)  )
         allocate( Href(nsp(1),numHref)        )  
         if (calc_rdf) then 
            allocate( ngdr(binNum)) 
            allocate( ngdr_below(binNum)) 
            allocate( ngdr_above(binNum)) 
            ngdr(:) = 0.0_DP
            ngdr_below(:) = 0.0_DP
            ngdr_above(:) = 0.0_DP
         endif
         !
         Href(:,:) = 0
         !
         call hello_print()
         ! 
         return
      end subroutine init
      !
      !Print out welcome message, echo the variables 
      Subroutine hello_print
         !
         implicit none
         !
         character(len=30)    :: dummy
         integer     :: i
         !
         write(*,*)
         write(*,'(1X,"--------------------------------------------------")')
         write(*,'(1X,"|      Solvation Structures of H3Oplus           |")')
         write(*,'(1X,"|                                                |")')
         write(*,'(1X,"|           Charles W. Swartz VI                 |")')
         write(*,'(1X,"--------------------------------------------------")')
         write(*,*)
         write(*,'(1X,"---Read-in values---")')
         !
         !Set-up
         write(*,'(/,2X,"Set-up")')
         write(*,'(2X,"Starting Step: ",15X,I10)') stepstart
         write(*,'(2X,"Ending Step: ",17X,I10)') stepstop
         write(*,'(2X,"Defect Radius: ",19X,f6.4)') rdefect
         write(*,'(2X, "Calculate Deltas:",22X, L)') calc_delta
         write(*,'(2X, "Calculate RDF:", 25X, L)') calc_rdf
         !
         !System
         write(*,'(/,2X,"System")')
         write(*,'(2X,"Position File: ", 13X, A12)') Trim(filePos)
         write(*,'(2X,"Cell File: ", 14X, a15)') trim(fileCel)
         write(*,'(2X,"Number of Species: ",18X, I3)') ntsp
         do i=1,ntsp,1
            write(*,'(2X,"Species", I3,": ",23X, I5)') i, nsp(i)
         enddo
         !
         !Hbonds
         write(*,'(/,2X,"Hydrogen Bonds")')
         write(*,'(2X,"H-Bond Radius: ",19X,f6.4)') Hbond_cut
         write(*,'(2X,"H-Bond Angle: ",20X,f6.2)') Hbond_angle
         write(*,*)
         !
         !Deltas
         if(calc_delta) then
            write(*,'(2X,"Deltas")')
            if (delta_donate) then
               write(*,'(2X,"Calculated for Donated Bonds")')
            else
               write(*,'(2X,"Calculated for Accepted Bonds")')
            endif
            write(*,'(2X,"Delta Above:  ",21X,f5.2)') delta_above
            write(*,'(2X,"Delta Below:  ",21X,f5.2)') delta_below
            write(*,*)
         endif
         !rdf
         if (calc_rdf) then
            write(*,'(2X,"Radial Distribution Functions")')
            write(*, '(2X,"Atom1: ", 31X,I2 )' )atom1
            write(*, '(2X,"Atom2: ", 31X,I2 )' )atom2
            write(*, '(2X,"BinNum: ", 27X,I5 )') binNum
            write(*,*)
         endif
         !Output
         write(*,'(1X,"---Output Files---",/)')
         write(*,'(2X,"Coordinate (Hbond) File:           coord-num.dat")')
         write(*,'(2X,"Non Proton Transfer Defect Index:  fort.199")')
         write(*,'(2X,"Proton Transfer Defect Index:      fort.198")')
         if (calc_delta) then
            write(*,'(2X,"H3O delta File:                    H3O-deltas.dat")')
         endif
         if(calc_rdf) then
            write(*,'(2X,"RDF Full:                          ", a)') filerdf_full
            write(*,'(2X,"RDF delta < ", f3.1,": ", 18X, a)') delta_below, filerdf_below
            write(*,'(2X,"RDF delta > ", f3.1,": ", 18X, a)') delta_above, filerdf_above
         endif
         write(*,'(2X,"Nearest Oxygen File:               nearest_O.dat")')
         write(*,'(2X,"Tagged Defect List:                fort.136")')
         write(*,'(2X,"Defect Tagging Process:            fort.145")')
         write(*,*)
         write(*,'(1X,"-------------------")')
         write(*,'(1X,"---Program START---")')
         write(*,'(1X,"-------------------")')
         write(*,*)
      End Subroutine hello_print

END MODULE
