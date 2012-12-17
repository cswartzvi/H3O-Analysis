 Program main

   USE ks_lib
   USE omp_lib

   implicit none
   
   include 'mpif.h'

   type atom_positions
      integer                 :: Num
      real(DP), allocatable   :: r(:,:)

   end type atom_positions

   type(atom_positions) :: tau(2) !1 is O 2 is H


   integer                 :: i, j, k, n, m, Nat, Ostar, Hprime, ncount, &
                              grid(3), ntsp

   !For parallel
   integer                 :: myid, me, nproc, per_proc, remain, &
                                 status(MPI_STATUS_SIZE), np

   real(DP), allocatable   :: rho(:,:,:), space(:,:,:,:), & !ix, iy, iz, point with it's coordinates
                              val(:), eig(:), rho_1d(:)

   real(DP)                :: valt, aprim(3,3), temp(3), dr(3), origin(3), &
                              tempphi(6), omega, disp(3), radius, dist2, &
                              midpt(3), radius2, ldos, degauss, wk=2.0_DP,&
                              E, dE, alat, dummy_real, ans, start_time, &
                              end_time

   !For reading rho 
   integer                 :: ix, iy, iz, i1, i2, i3, ns, E_start, E_stop, &
                              ind_x(6), ind_y(6), ind_z(6), nbnd, nE, &
                              ierr, dummy_int, atype, na, grid_tot, &
                              loopx, loopy, loopz


   character(len=256)       :: dummy_array
   character(len=3)        :: x1
   character(len=256)      :: buffer, outfile
   character(len=100)      :: fileroot

   Call MPI_INIT(ierr)
   Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
   Call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
   if (myid == 0) start_time = MPI_Wtime()
   !-------------------------------------------
   !Input: COMMAND LINE ARUGMENTS! 1) outfile 2) Ostar 3)Hprim
   tau(1)%Num = 64 
   tau(2)%Num = 127
   nbnd=256
   radius = 1.0_DP 
   dE=0.01_DP
   degauss=0.09524_DP
   E_start = -30
   E_stop= 20
   !-------------------------------------------
   !-------------------------------------------

   !-------------------------------------------
   !Initialize
   allocate ( val(nbnd), eig(nbnd) )
   radius2 = radius**2
   nE=abs(E_stop - E_start)/dE + 1 !add one for zero

   !Number jobs per processor
   per_proc = Ceiling(nbnd/real(nproc))
   remain = per_proc * nproc - nbnd

   !-------------------------------------------
   !-------------------------------------------


   !-------------------------------------------
   !Read in Values
   !-------------------------------------------
   !read in eig.dat MUST BE prsent!
   if (myid == 0 ) then
      open(unit=2,file='eig.dat', IOSTAT=ierr)
      if(ierr /= 0) then
         write(*,*) 'The eig.dat file is missing or corrupted'
         stop
      endif
      do i=1,nbnd,1
         read(2,*) eig(i)
      enddo
      close(2)

      CALL getarg(1, fileroot)
      CALL getarg(2, buffer)
      read(buffer,*) Ostar
      CALL getarg(3, buffer)
      read(buffer,*) Hprime
   
   endif

   Call MPI_BCAST(fileroot, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(eig, nbnd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(Ostar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(Hprime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   if (ierr /= 0 ) then
      print *, '   MPIT_BCAST error ', myid
   endif

   !-------------------------------------------

   proc_loop: do np=1,per_proc,1

      me = myid + (np-1) * nproc + 1
      if (me > nbnd) exit

      write(6,*) '    KS State : ', me

      valt = 0.0_DP   
      ldos = 0.0_DP

      write(x1, '(I0)'), me 
      outfile = TRIM(fileroot)//TRIM(x1)
      open(unit=1, file=TRIM(outfile),iostat=ierr)
      if (ierr /= 0) then
         print *, 'Error: Opening ', outfile
         stop
      endif

      !Read in the position file
      !read(1,*) dummy_array
      read(1,*) (grid(i), i=1,3), dummy_int, dummy_int, dummy_int, &
                  Nat, ntsp
      read(1,*) dummy_int, alat, dummy_real, dummy_real, dummy_real, &
                  dummy_real, dummy_real
      read(1,*) dummy_array
      read(1,*) dummy_array
      read(1,*) dummy_array


      !Only the first time through, ON EACH PROCESSOR
      if (np == 1) then
         !Find Total Mesh and Allocate the temp array for reading
         grid_tot = grid(1)*grid(2)*grid(3)
         allocate(rho_1d(grid_tot))
         !Atoms dervied data type
         allocate(tau(1)%r(tau(1)%Num, 3), tau(2)%r(tau(2)%Num, 3) )
         !rho and space grids
         allocate( rho(grid(1), grid(2), grid(3)), space(grid(1), grid(2), grid(3), 3) )

         !Set up space with matching indexes as the 
         do i=1,3,1
            aprim(i,i) = alat * 0.52917721_DP
            dr(i) = aprim(i,i)/(grid(i) -1)
         enddo
      end if


      !Read in the atomic positions
      n=1
      m=1
      do k=1,Nat,1
         
         read(1,*) na, (temp(j), j=1,3), atype

         if(atype == 2) then
            tau(2)%r(n,:) = temp * alat * 0.52917721_DP 
            n = n +1
         endif

         if(atype == 1) then
            tau(1)%r(m,:) = temp * alat * 0.52917721_DP
            m = m + 1
         endif  

      enddo

      if (na /= Nat) then
         print *, '  Problem reading Atom numbers'
         stop
      endif


      !Read/create the space and rho files AND 
      read(1,*) (rho_1d(i), i=1,grid_tot)

      if (myid == 0 .and. np == 1 ) then
         do i=1, grid_tot
            write(189,*) rho_1d(i)
         enddo
      endif
     
      !Convert From temp (1D) vector to rho (3D) Array 
      !Quantum Espresso outputs the DATA file as 1D vectors, converts to a
      !3D array in chegens.f90
      do i=1,grid_tot
         loopz = INT( ABS(i - 1)/(grid(1)*grid(2)) ) + 1
         loopy = INT( ABS( (i - 1) - (loopz - 1)*grid(1)*grid(2) ) / grid(1) ) + 1
         loopx = i - (loopz-1)*grid(1)*grid(2) - (loopy-1) * grid(1)
         rho(loopx, loopy, loopz) =  rho_1d(i)
      enddo


      !Set mid point
      midpt(:) = (/   0.5*(tau(1)%r(Ostar,1) + tau(2)%r(Hprime,1)), &
                      0.5*(tau(1)%r(Ostar,2) + tau(2)%r(Hprime,2)), &
                      0.5*(tau(1)%r(Ostar,3) + tau(2)%r(Hprime,3))  /)

      !Non mid-point
      !midpt(:) = (/tau(2)%r(Ostar,1),tau(2)%r(Ostar,2),tau(2)%r(Ostar,3) /)


      
      !Loop over all point and see which ones are within the radius
      do iz=1, (grid(3))
         do iy=1,(grid(2))
            do ix=1,(grid(1))

               !SPACE ARRAY
               space(ix,iy,iz,1:3) = (/(ix-1)*dr(1), (iy-1)*dr(2), (iz-1)*dr(3) /)


               do i=1,3,1
                  disp(i) = space(ix,iy,iz,i) - midpt(i) 
                  disp(i) = disp(i) - NINT(disp(i)/aprim(i,i))*aprim(i,i)
               enddo

               dist2 = disp(1)**2 + disp(2)**2 + disp(3)**2

               !Integration, yay!
               if ( dist2 < radius2 ) then


                  valt = valt + rho(ix,iy,iz)

               endif

            enddo
         enddo
      enddo

      valt = valt*dr(1)*dr(2)*dr(3)


      if (np == per_proc) deallocate(tau(1)%r, tau(2)%r, rho, rho_1d, space )
      close(1)

      !Send to Root process
      if (myid /= 0) then
         Call MPI_SEND(valt, 1, MPI_DOUBLE_PRECISION, 0, me, MPI_COMM_WORLD, ierr)
      endif

      !Receive by Root Process
      if(myid == 0 ) then

         val(me) = valt

         if (np == per_proc) then !For the last loop, may have unfilled procs

            do i =1,(nproc-1-remain),1 !from the other processors
               Call MPI_RECV(ans, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                              MPI_COMM_WORLD, status, ierr)
               val(status(MPI_TAG)) = ans
            enddo

         else

            do i =1,(nproc-1),1 !from the other processors
               Call MPI_RECV(ans, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                              MPI_COMM_WORLD, status, ierr)
               val(status(MPI_TAG)) = ans
            enddo

         endif
      endif

 enddo proc_loop


 if(myid == 0)then
   open(unit=3, file='intKS.dat', IOSTAT=ierr, STATUS='OLD')
   if (ierr /= 0 ) then
      write (*,*) 'There was some issue opeing the intKS.dat file'
      stop
   endif

   do i=1,nE,1

      ldos = 0.0_DP
      E = E_start + (i-1)*dE
   
      do ns=1,nbnd,1
         ldos = ldos +  wk * w0gauss ( (E - eig(ns))/ degauss) * val(ns)
      enddo

      ldos = ldos/degauss

      write(3,*) E , ldos

   enddo
   
   end_time = MPI_Wtime()
   close(3)
   Write(*,*) '     Total time: ', end_time - start_time
 endif
   
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
