 Module rdf_lib

   implicit none

   integer, parameter   :: DP  = SELECTED_REAL_KIND(15,99)

   contains

      subroutine invert(Mi, Mo, det)

         implicit none

         real(DP), intent(in)     :: Mi(3,3)
         real(DP), intent(out)    :: Mo(3,3), det
         real(DP)                 :: tmp(3,3), s

         integer              :: i,j,k,l  !indexes
         integer              :: n,ir     !int counters

         det = 0.0
         s   = 1.0
         i   = 1
         j   = 2
         k   = 3


         do
            do n=1,3,1
               det= det + s*Mi(1,i)*Mi(2,j)*Mi(3,k)
               l = i
               i = j
               j = k
               k = l
            end do

            i = 2
            j = 1
            k = 3
            s = -s
            if (s .GT. 0.0) then
               exit     
            endif
         enddo

         IF(ABS(det) .LT. 1.0e-20) then
            write(*,*) 'Error: Singular Matrix'
            Stop   
         endif

         i = 1
         j = 2
         k = 3

         do ir=1,3
            tmp(ir,1) = (Mi(2,j)*Mi(3,k) - Mi(2,k)*Mi(3,j)) / det
            tmp(ir,2) = (Mi(3,j)*Mi(1,k) - Mi(3,k)*Mi(1,j)) / det
            tmp(ir,3) = (Mi(1,j)*Mi(2,k) - Mi(1,k)*Mi(2,j)) / det
         
            l = i
            i = j
            j = k
            k = l
         enddo
        
         do l=1,3
            do k=1,3
               Mo(k,l) = tmp(k,l)
            enddo
         enddo 

      end subroutine invert
   
      subroutine r_to_s(R,S,M) 
      !===========================================================================
      !Takes a scaled (out) and real (in) vector AND and transform matrix
      !Must be a real matrix
      !
      !===========================================================================
         implicit none

         real(DP), intent(in)     :: R(3), M(3,3)
         real(DP), intent(out)    :: S(3)
      
         integer              :: i, j

         do i=1,3,1
      
            s(i) = 0.0
   
            do j=1,3,1
               S(i) = S(i) + R(j)*M(j,i)
            enddo
         enddo
      
      end subroutine r_to_s

      subroutine s_to_r(S,R,M)
      !===========================================================================
      !Takes a scaled (in) and real (out) vector AND and transform matrix
      !Must be a real matrix
      !
      !===========================================================================
         implicit none   

         real(DP), intent(in)     :: S(3), M(3,3)
         real(DP), intent(out)    :: R(3)

         R(1) = S(1)*M(1,1) + S(2)*M(2,1) + S(3)*M(3,1)
         R(2) = S(1)*M(1,2) + S(2)*M(2,2) + S(3)*M(3,2)
         R(3) = S(1)*M(1,3) + S(2)*M(2,3) + S(3)*M(3,3)


      end subroutine s_to_r

 End Module rdf_lib 
