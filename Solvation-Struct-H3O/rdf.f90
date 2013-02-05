MODULE rdf_module
   !
   use main_variables,        only : DP
   !
   implicit none
   !
   Contains
      !
      !---------Count PAIRS-----------
      Subroutine count_pairs(cal, ngdr_x)
         !
         use Hbonds_module,   only : delta
         use main_variables,  only : binNum, atom1, atom2, stau, nsp, box, aprim, mic
         use defect_module,   only : Ostar, ndefect
         use gen_lib,         only : get_rdist
         !
         implicit none
         !
         integer,intent(in)                  :: cal            !Calculation type 1-full 2-delta < 0.1, 3-delta > 0.5
         integer,allocatable, intent(inout)  :: ngdr_x(:)      !Current Total g(r) array, 
                                                               !Recall there are 3 types of RDF (full, < 0.1, > 0.5)
         !
         integer                 :: n,m,        &  !indexes
                                    nbin           !Bin Number Index
         !
         real(DP), allocatable   :: ngdrt(:)       !temp g(r)
         !
         real(DP)                :: Hbox,       &  !length of half the region (including multiple supercells)
                                                   !Recall: we can only  compute up to L/2 of the  total region
                                    res,        &  !resolution (size) of bins binNum/Hbox (number/length) 
                                    r2,         &  !Distance calculation
                                    rdist(3)       !component distance
         !
         !
         Hbox = box/2.
         res  = binNum/Hbox
         allocate(ngdrt(binNum))
         ngdrt(:) = 0.0_DP
         !
         !Index atoms to be used in the pair correlation
         !NOTE: atom1 is always either O* or H', atom2 are the other
         !atoms Ow, Hw    
         if (atom1 == 1 ) then
            n = Ostar
         elseif (atom1 == 2 ) then
!            n = Hprime
         else
            Write(*,*) ' ERROR: atom1 is defined incorrectly!!!'
            stop
         endif
         !
         !
         if ( ndefect == 1 ) then 
         
            do m=1,nsp(atom2),1
               !Skip this cylce if atom2 is Ow and the c
                  if ( (atom2 ==1) .and. (m == Ostar) ) cycle
                  !  
                  !calculate distance 
                  CALL get_rdist(stau(1:3,n,atom1),stau(1:3,m,atom2),rdist,mic,aprim) 
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
         !
         !update true unnormalized, raw pair correlation
         do n=1,binNum,1
            ngdr_x(n) = ngdr_x(n) + ngdrt(n)
         enddo

         deallocate(ngdrt)
   End Subroutine count_pairs
   !
   Subroutine final_rdf(cal, ngdr_x, total_count)
      !
      use main_variables,        only : DP, ngdr, nsp, atom1, atom2, mic, omega, box
      use main_variables,        only : ngdr, atom1, atom2, binNum, norm
      use main_variables,        only : ao, pi
      use hbonds_module,         only : delta
      !
      implicit none
      !
      integer,intent(in)                  :: cal            !Calculation type 1-full 2-delta < 0.1, 3-delta > 0.5
      integer,allocatable, intent(inout)  :: ngdr_x(:)      !Current Total g(r) array, 
      integer,intent(in)                  :: total_count    !total number of data points in the current g(r)
      !
      integer              :: iunit,      &  !File number for the particular rdf, determined by cal
                              n              !general index
      !
      real(DP)             :: gofr,       &  !Final value of the g(r) 
                              res,        &  !resolution (size) of bins binNum/Hbox (number/length)
                              Hbox,       &  !length of half the region (including multiple supercells)
                                             !Recall: we can only  compute up to L/2 of the  total region
                              vshell,     &  !Volume of the infinitesimal
                              dens,       &  !Density
                              r, r2          !distance and distance squared 
      !
      !check to make sure the particular case is met
      select case (cal)
         case(1)
            iunit = 160
         case(2) 
            iunit = 161
         case(3)
            iunit = 162
         case default
      end select
      !
      Hbox = box/2.0_DP
      res  = dble(binNum)/Hbox
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
         gofr = ngdr_x(n)/ ( dble(total_count*mic*norm*dens*vshell) )

         !
         write(iunit,*) (r*ao), gofr
      enddo
      !
   End Subroutine final_rdf
   !
END MODULE rdf_module
