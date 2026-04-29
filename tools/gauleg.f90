subroutine gauleg(x1,x2,x,w,n)
!
! Computes Gaussian weights and locations for Legendre polynomials
!
 implicit none
 integer, intent(in) :: n
 real*8, intent(in) :: x1, x2
 real*8, dimension(n), intent (out) :: x, w
 real*8, parameter :: eps=3.0D-10
 integer :: i, j, m
 real*8 :: p1, p2, p3, pp, xl, xm, z, z1, pi
!
 pi = -acos(-1.0)
 m = (n+1)/2
 xm = 0.5D0*(x1 + x2)
 xl = 0.5D0*(x2 - x1)
 do i = 1,m
   !z = cos(pi*(i - 0.25D0)/(n + 5.0D0)) ! first guess value (not good) NR
   z = cos(0.5D0*pi*(4.0D0*float(i) -1.0D0)/(2.0D0*n + 1.0D0)) ! first guess value - Durran (p. 205)
   do   
     p1 = 1.0D0
     p2 = 0.0D0
     do j = 1,n
       p3 = p2
       p2 = p1
       p1 = ((2.0D0*j - 1.0D0)*z*p2 - (j - 1.0D0)*p3)/float(j)  ! evaluation Legendre polynomial at z (recurrence)
     enddo
     pp = n*(z*p1 - p2)/(z*z - 1.0D0) ! derivate 
     z1 = z
     z = z1 - p1/pp ! Newton step
     if ((abs(z - z1)) < eps) exit
   enddo 
   x(i) = xm - xl*z
   x(n+1-i) = xm + xl*z
   w(i) = 2.0D0*xl/((1.0D0 - z*z)*pp*pp)
   w(n+1-i) = w(i)
 enddo
 return
end subroutine gauleg  
