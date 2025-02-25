!-------------------------------------------------------------------------------
!Program To calculate the local density of states around an OH^- defect
!
!Written By Charles Swartz Fri July 1 2012, converted to parallel on 
!Fri Aug 3 2012
!
!
!-------------------------------------------------------------------------------
!NOTE: This Program in sdesigned to read the .xsf file produed by the 
!cppp.x program. It will NOT read the .xsf files from the PP code because
!the PP generates the atomic positions with letters NOT numbers ~ Needs to be fixed
!
!-------------------------------------------------------------------------------
!INPUT(s):
!1)input_ldos.dat    !Namelist: For Input 
!2)eig.dat           !Single columun eig-value
!
!OUTPUT(s):
!1)intKS.dat 
!-------------------------------------------------------------------------------


 Program main
   !
   USE ks_lib
   USE omp_lib
   !
   implicit none
   !
   include 'mpif.h'
   !
   !Constants
   real(DP), parameter        :: ao=0.52917721_DP
   !Atom speices derived data types
   type atom_positions
      integer                 :: Num
      real(DP), allocatable   :: r(:,:)
   end type atom_positions
   !
   type(atom_positions) :: tau(2) !!! tau(1) => O, tau(2) => H !!!
   !
   !
   integer                 :: i, j, k, n, m, &
                              Nat,     &  !Total Number of Atoms
                              Ostar1,  &  !Index of the Ostar1    
                              Ostar2,  &  !Index of the Ostar2
                              Hshared, &  !Index of the shared Proton
                              O_count, &  !Counting index of the Ostars
                              Ostar_x, &  !Current Index of the Ostars
                              atype,   &  !Atomic Nummber of the Atom
                              grid(3), &  !FFT grid Dimensions
                              ntsp,    &  !Total Number of Atoms
                              ns,      &  !State index
                              na,      &  !Atomic index
                              nbnd,    &  !Number of states
                              nE,      &  !Number of Energy Divisions
                              ierr,    &  !Error int
                              shared_count,     &  !Index of nearest H lists
                              DeltaR_index,     &  !Index of nearest H lists
                              near_list1(10),   &  !Nearest H list Ostar1
                              near_list2(10)      !Nearest H list Ostar2

   !
   real(DP), allocatable   :: rho(:,:,:),       &  !ix, iy, iz, value at point
                              space(:,:,:,:),   &  !ix, iy, iz, point with it's coordinates
                              val(:),           &  !local DOS array
                              eig(:)               !Eigenvalues READ in from eig.dat
   !
   real(DP)                :: valt,       &  !Temp value
                              norm,       &  !Normalization temp value
                              aprim(3,3), &  !Prim matrix
                              alat,       &  !Lattice Constant
                              temp(3),    &  !General temp vector
                              dr(3),      &  !discretization units
                              origin(3),  &  !Origin of the FFT grid
                              omega,      &  !Volume of the Unit Cell
                              disp(3),    &  !Displacment vectors
                              radius,     &  !Radius of the integration shell
                              radius2,    &  !Radius sqaures
                              shared_lim, &  !Min Distance
                              R1,         &  !R1 Distance (From Ostar1, Shared Proton)
                              R2,         &  !R2 Distance (From Ostar2, Shared Proton)
                              DeltaR,     &  ! = |R2 - R1|, Current Value
                              DeltaR_MIN, &  ! = |R2 - R1|, MIN value = Hshared 
                              dist2,      &  !Distance Squared
                              midpt(3),   &  !Midpoint
                              ldos,       &  !FINAL DOS Values
                              degauss,    &  !Parameter for the guassian
                              wk=2.0_DP,  &  !Weight for the integration
                              E, dE,      &  
                              E_start,    &
                              E_stop,     &  !Engery val, increment, start and stop  
                              ans,        &  !tranfer value amoung processors
                              start_time, end_time, &
                              near_dist1(10),       & !Nearest H Dist (Ostar1)
                              near_dist2(10)          !Nearest H Dist (Ostar2)

   !
   !----------------------------------------------
   !For parallel execution
   !----------------------------------------------
   integer                 :: myid,       &  !Absoulte Index of the processor 0 .... (nproc -1)
                              me,         &  !relative index 
                              nproc,      &  !total number of processors
                              per_proc,   &  !Operations per processor (ideal)
                              remain,     &  !diffence in per_proc and nproc
                              np,         &  !processr index
                              status(MPI_STATUS_SIZE)
   !----------------------------------------------

   !
   !----------------------------------------------
   !Varaiables fo reading in rho over the FFT grid
   !----------------------------------------------
   integer                 :: ix, iy, iz, i1, i2, i3,       &
                              ind_x(6), ind_y(6), ind_z(6), &
                              dummy_int, ncount
   !----------------------------------------------

   character(len=256)      :: dummy_array, buffer, outfile
   character(len=3)        :: x1
   character(len=100)      :: fileroot    !Root of the input file
   character(len=256)      :: output, eigfile

   logical                 :: shared_proton

   NAMELIST /inputpp/ nbnd, radius, dE, degauss, E_start, E_stop, shared_proton, &
                        output, eigfile

   Call MPI_INIT(ierr)
   Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
   Call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
   if (myid == 0) start_time = MPI_Wtime()
   !
   !Number of Atoms
   tau(1)%Num = 63 
   tau(2)%Num = 127
   allocate ( tau(1)%r(tau(1)%Num,3), tau(2)%r(tau(2)%Num,3) )
   !
   !nbnd=256
   !radius = 1.0_DP 
   !dE=0.01_DP
   !degauss=0.09524_DP
   !E_start = -30.0_DP
   !E_stop=0.0_DP
   !
   !
   !-------------------------------------------
   !Read in Values
   !-------------------------------------------
   if (myid == 0 ) then
      !
      !Read in command line
      !Input: COMMAND LINE ARUGMENTS! 1) outfile 2) Ostar 
      CALL getarg(1, fileroot)
      CALL getarg(2, buffer)
      read(buffer,*) Ostar1
      CALL getarg(3, buffer)
      read(buffer,*) Ostar2
      !
      !Read in the input file    
      !open(unit=21, file='input_ldos.dat', IOSTAT=ierr)
      !if (ierr /= 0) then
      !   write(6,*) '  ERROR: Opening input_ldos.dat!'
      !   stop
      !endif
      read(5,nml=inputpp,IOSTAT=ierr)
      if (ierr /= 0) then
         write(6,*) '  ERROR: Incorrect format for the namelist!!'
         stop 
      endif
      close(21)
   endif
   !
   !Broadcast all values to all nodes
   Call MPI_BCAST(fileroot, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(output, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(eigfile, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(Ostar1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(Ostar2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(nbnd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(radius, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(dE, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(degauss, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(E_start, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(E_stop, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(shared_proton, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   !
   if (ierr /= 0 ) then
      print *, '  MPI_BCAST error ', myid
   endif
   !-------------------------------------------
   !
   !-------------------------------------------
   !Initialize
   !-------------------------------------------
   allocate ( val(nbnd), eig(nbnd) )
   radius2 = radius**2
   nE=abs(E_stop - E_start)/dE + 1 !add one for zero
   !
   !Number jobs per processor
   per_proc = Ceiling(nbnd/real(nproc))
   remain = per_proc * nproc - nbnd
   !
   if(myid == 0) then
      !read in eig.dat MUST BE prsent!
      open(unit=2,file=eigfile, IOSTAT=ierr)
      if(ierr /= 0) then
         write(*,*) '  ERROR: The eig.dat file is missing or corrupted'
         stop
      endif
      do i=1,nbnd,1
         read(2,*) eig(i)
      enddo
      close(2)

      write(6,*) 'Eigfile: ', eigfile
      write(6,*) 'output:  ', output
   endif
   Call MPI_BCAST(eig, nbnd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   !-------------------------------------------


   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
   !Start Main Processor loop
   !-------------------------------------------
   proc_loop: do np=1,per_proc,1

      me = myid + (np-1) * nproc + 1
      if (me > nbnd) exit

      write(6,*) ' KS State : ', me

      valt = 0.0_DP   
      ldos = 0.0_DP
      norm = 0.0_DP

      write(x1, '(I0)'), me 
      outfile = TRIM(fileroot)//TRIM(x1)//'.xsf'
      open(unit=1, file=TRIM(outfile),iostat=ierr)
      if (ierr /= 0) then
         print *, '  Error: Opening ', outfile
         stop
      endif

      !---------------------------------
      !Read in the position file
      !---------------------------------
      read(1,*) dummy_array
      read(1,*) dummy_array
      !
      !Read in the aprim matrix 
      do i=1,3,1
         read(1,*) ( aprim(i,j),j=1,3) 
         aprim(i,1:3)= aprim(i,1:3) 
      enddo
      !
      read(1,*) dummy_array
      read(1,*) Nat, dummy_int
      !
      !Read in the atomic positions
      n=0
      m=0
      do k=1,Nat,1
         !
         read(1,*) atype, (temp(i),i=1,3) 
         !
         if(atype == 1) then
            n = n +1
            tau(2)%r(n,1:3) = temp(1:3)
         endif
         !
         if(atype == 8) then
            m = m + 1
            tau(1)%r(m,1:3) = temp(1:3)
         endif  
         !
      enddo
      !
      !Atom number Check
      na = n + m
      if (na /= Nat) then
         print *, '  ERROR:Problem reading Atom numbers'
         stop
      endif
      !
      !Read in the FFT Grid section
      read(1,*) dummy_array
      read(1,*) dummy_array
      read(1,*) dummy_array
      read(1,*) (grid(i),i=1,3)
   
      !Read origin and second aprim matrix
      read(1,*) (origin(i),i=1,3) 
      do i=1,3
         read(1,*) (temp(j), j=1,3)
      enddo

      omega = grid(1)*grid(2)*grid(3)

      if (np == 1) then
         allocate( rho(grid(1), grid(2), grid(3)), space(grid(1), grid(2), grid(3), 3) )
      end if

      !Set up space increments 
      do i=1,3,1
         dr(i) = aprim(i,i)/(grid(i) -1)
      enddo


      !Read/create the space and rho files
      ncount=0
      do iz=1, (grid(3))
         do iy=1,(grid(2))
            do ix=1,(grid(1))
               !
               !SPACE ARRAY
               space(ix,iy,iz,1:3) = (/(ix-1)*dr(1), (iy-1)*dr(2), (iz-1)*dr(3) /)
               !
               !RHO ARRAY
               if (ncount.lt.6) then
                  ncount = ncount + 1
               else
                  !read(1,'(6(1pe15.7))') (rho(ind_x(i),ind_y(i),ind_z(i)),i=1,6)
                  read(1,*) (rho(ind_x(i),ind_y(i),ind_z(i)),i=1,6)
                  ncount=1
               endif
               ind_x(ncount) = ix
               ind_y(ncount) = iy
               ind_z(ncount) = iz
            enddo
         enddo
      enddo
      !read(1,'(6(1pe15.7))') (rho(ind_x(i),ind_y(i),ind_z(i)),i=1,ncount)
      read(1,*) (rho(ind_x(i),ind_y(i),ind_z(i)),i=1,ncount)


      !For Shared protons
      !Loop over all H atoms and find those within a radius of 2.0 Angstroms
      if (shared_proton) then
         !Shared Limit in Angstrom
         shared_lim = 2.0

         !For shared proton's
         do i =1,10
            near_list1(i) = 100
            near_list2(i) = 100
         enddo


         do O_count=1,2,1

            if (O_count ==1) Ostar_x = Ostar1           
            if (O_count ==2) Ostar_x = Ostar2           
            shared_count = 0

            do i=1,tau(2)%Num

               do j=1,3,1
                  disp(j) = tau(2)%r(i,j) - tau(1)%r(Ostar_x,j)
                  disp(j) = disp(j) - NINT(disp(j)/aprim(j,j))*aprim(j,j)
               enddo

               dist2 = disp(1)**2 + disp(2)**2 + disp(3)**2
               !
               !Integration, yay!
               if ( dist2 < shared_lim**2 ) then
                  shared_count = shared_count + 1

                  !Check to make sure the share list isn't greater then 10
                  if(shared_count > 10 ) then
                     write(*,*) '  ERROR: Share List greater then 10!!!'
                     stop
                  endif
   
                  !Write The share index list
                  if(O_count == 1) then
                     near_list1(shared_count) = i
                     near_dist1(shared_count) = sqrt(dist2)
                  elseif(O_count == 2) then
                     near_list2(shared_count) = i
                     near_dist2(shared_count) = sqrt(dist2)
                  endif
               endif
            enddo

            if(myid == 0) then
            
               if(O_count == 1) then

                  !Nearest H list 1
                  open(unit=10, file='nearest_H_list1', status='unknown')
                  do i=1,10,1
                     if (near_list1(i) /= 100) &
                        write(10,*) near_list1(i), near_dist1(i) 
                  enddo
                  close(10)

               elseif(O_count == 2) then

                  !Nearest H list 2
                 open(unit=11, file='nearest_H_list2', status='unknown')
                  do i=1,10,1
                     if (near_list2(i) /= 100) &
                        write(11,*) near_list2(i), near_dist2(i)
                  enddo
                  close(11)

               endif
            endif
         enddo

         
         !find Hshared!!!
         do O_count=1,2,1
            shared_count = 0 
            DeltaR_MIN = 100

            if (O_count == 1) then
            
               !From the first Nearest list
      
               do i=1,10,1

                  if( near_list1(i) /= 100) then
                     

                     !Find R1
                     do j=1,3,1
                        disp(j) = tau(2)%r(near_list1(i),j) - tau(1)%r(Ostar1,j)
                        disp(j) = disp(j) - NINT(disp(j)/aprim(j,j))*aprim(j,j)
                     enddo
                     R1 = sqrt(disp(1)**2 + disp(2)**2 + disp(3)**2)

                     !Find R2
                     do j=1,3,1
                        disp(j) = tau(2)%r(near_list1(i),j) - tau(1)%r(Ostar2,j)
                        disp(j) = disp(j) - NINT(disp(j)/aprim(j,j))*aprim(j,j)
                     enddo
                     R2 = sqrt(disp(1)**2 + disp(2)**2 + disp(3)**2)


                     DeltaR = abs(R1 - R2)

                     if (DeltaR < DeltaR_MIN) then
                        DeltaR_MIN = DeltaR
                        DeltaR_index = near_list1(i)
                     endif
                  endif
               enddo
 
            else
            
               !From the Second Nearest List
               
               do i=1,10,1

                  if( near_list2(i) /= 100) then


                     !Find R1
                     do j=1,3,1
                        disp(j) = tau(2)%r(near_list2(i),j) - tau(1)%r(Ostar1,j)
                        disp(j) = disp(j) - NINT(disp(j)/aprim(j,j))*aprim(j,j)
                     enddo
                     R1 = sqrt(disp(1)**2 + disp(2)**2 + disp(3)**2)

                     !Find R2
                     do j=1,3,1
                        disp(j) = tau(2)%r(near_list2(i),j) - tau(1)%r(Ostar2,j)
                        disp(j) = disp(j) - NINT(disp(j)/aprim(j,j))*aprim(j,j)
                     enddo
                     R2 = sqrt(disp(1)**2 + disp(2)**2 + disp(3)**2)


                     DeltaR = abs(R1 - R2)

                     if (DeltaR < DeltaR_MIN) then
                        DeltaR_MIN = DeltaR
                        DeltaR_index = near_list2(i)
                     endif
                  endif
               enddo


            endif
         enddo

         Hshared = DeltaR_index
      
         if (myid == 0) then
            write(*,*) '  Shared Proton:', Hshared
            write(*,*) (tau(2)%r(Hshared,j),j=1,3)
         endif
         !Non mid-point
         midpt(:) = (/tau(2)%r(Hshared,1),tau(2)%r(Hshared,2),tau(2)%r(Hshared,3) /)

      else

         !Set mid point
         !midpt(:) = (/   0.5*(tau(1)%r(Ostar1,1) + tau(2)%r(Hprime,1)), &
         !                0.5*(tau(1)%r(Ostar1,2) + tau(2)%r(Hprime,2)), &
         !                0.5*(tau(1)%r(Ostar1,3) + tau(2)%r(Hprime,3))  /)

         !Non mid-point
         midpt(:) = (/tau(1)%r(Ostar1,1),tau(1)%r(Ostar1,2),tau(1)%r(Ostar1,3) /)

      endif



   
      !---------------------------------
      !Loop over all point and see 
      !which ones are within the radius
      !---------------------------------
      do iz=1, (grid(3))
         do iy=1,(grid(2))
            do ix=1,(grid(1))
               !
               !Normalization
               norm = norm + rho(ix,iy,iz)
               !
               do i=1,3,1
                  disp(i) = space(ix,iy,iz,i) - midpt(i) 
                  disp(i) = disp(i) - NINT(disp(i)/aprim(i,i))*aprim(i,i)
               enddo
               !
               dist2 = disp(1)**2 + disp(2)**2 + disp(3)**2
               !
               !Integration, yay!
               if ( dist2 < radius2 ) then
                  !
                  valt = valt + rho(ix,iy,iz)
                  !
               endif
               !
            enddo
         enddo
      enddo
      !
      norm = norm/(grid(1)*grid(2)*grid(3))
      !norm = norm
      !valt = valt*dr(1)*dr(2)*dr(3)
      valt = valt*dr(1)*((4.0/3.0)*3.14*radius**3)
      !valt = valt/(grid(1)*grid(2)*grid(3))
      !valt = valt/(omega*( ((4.0/3.0)*3.14*radius**3)/(aprim(1,1)*aprim(2,2)*aprim(3,3))) )
      !
      !write(6,*) 'myid: ', myid, 'Norm: ', norm, 'val:', valt
      !
      !
      if (np == per_proc) deallocate(tau(1)%r, tau(2)%r, rho, space )
      close(1)

      !Send to Root process
      if (myid /= 0) then
         Call MPI_SEND(valt, 1, MPI_DOUBLE_PRECISION, 0, me, MPI_COMM_WORLD, ierr)
      endif

      !Receive by Root Process
      if(myid == 0 ) then
         !
         val(me) = valt
         !
         if (np == per_proc) then !For the last loop, may have unfilled procs
            !
            do i =1,(nproc-1-remain),1 !from the other processors
               Call MPI_RECV(ans, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                              MPI_COMM_WORLD, status, ierr)
               val(status(MPI_TAG)) = ans
            enddo
            !
         else
            !
            do i =1,(nproc-1),1 !from the other processors
               Call MPI_RECV(ans, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                              MPI_COMM_WORLD, status, ierr)
               val(status(MPI_TAG)) = ans
            enddo
            !
         endif
      endif
      !
   enddo proc_loop
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   !End Main Processor loop
   !----------------------------------------------


   if(myid == 0)then
      !Open Main output file
      open(unit=3, file=output, IOSTAT=ierr, STATUS='UNKNOWN')
      if (ierr /= 0 ) then
         write (*,*) '  ERROR: There was some issue opeing the intKS.dat file'
         stop
      endif
      !
      !Open raw output file (no broadening)
      open(unit=4, file=output//raw, IOSTAT=ierr, STATUS='UNKNOWN')
      if (ierr /= 0 ) then
         write (*,*) '  ERROR: There was some issue opeing the intKS.dat file'
         stop
      endif
      !
      !Write the raw file (Easier to find peak states) before broadening
      do ns=1, nbnd,1
         write(4,*) ns, val(ns)
      enddo
      close(4)
      !
      !Final packagin of the ldos and E values
      do i=1,nE,1
         !
         ldos = 0.0_DP
         E = E_start + (i-1)*dE
         !
         do ns=1,nbnd,1
            ldos = ldos +  wk * w0gauss ( (E - eig(ns))/ degauss) * val(ns)
         enddo
         !
         ldos = ldos/degauss
         !
         write(3,*) E , ldos/omega
         !
      enddo
      !
      end_time = MPI_Wtime()
      close(3)
      Write(*,*) '    Total time: ', end_time - start_time
   endif
   !
   deallocate ( val, eig )
   
   Call MPI_FINALIZE(ierr)
   
 CONTAINS
   
   real(DP) function w0gauss (x)

         implicit none

         REAL(DP), PARAMETER  :: pi                = 3.14159265358979323846_DP
         REAL(DP), PARAMETER  :: sqrtpi            = 1.77245385090551602729_DP
         REAL(DP), PARAMETER  :: sqrtpm1           = 1.0_DP / sqrtpi        

         real(DP) :: x
         ! output: the value of the function
         ! input: the point where to compute the function

         real(DP) :: arg

         arg = min (200.d0, x**2)
         w0gauss = exp ( - arg) * sqrtpm1

   end function w0gauss 

End Program main
