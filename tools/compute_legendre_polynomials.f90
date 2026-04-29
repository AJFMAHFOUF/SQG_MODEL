program compute_legendre_polynomials
!
! Requires the following additionals routines :
! - plgndr.f90 : computes a Legendre associated polynomial of degree n and order m for a value x
! - factrl.f90 : computes the factorial of an integer (for normalization factor)
! - gammln.f90 : computes the logarithm of a Gamma function (factorial for large integers)
! - gauleg.f90 : computes the Gaussian latitudes (locations where to evaluate Legendre polynomials) 
!                and the corresponding Gaussian weights
!
! All these routines have been coded/adapted from Numerical Recipes
!
!   Jean-François Mahfouf (05/06/2023)
!
 implicit none
 integer, parameter :: mm = 42 ! maximum wavenumber (truncation) 
 integer, parameter :: nn = (3*mm + 1)/2 + 1 ! number of latitudes to avoid aliasing for quadratic terms
 integer, parameter :: mmax = (mm+1)*(mm+4)/2 ! number of wavenumbers to be stored 
 real*8 :: x1, x2
 real*8 :: pi, z1, z2, z3, factrl
 real*8, dimension(nn) :: x, w
 real*8, dimension(mmax,nn) :: pol_legndr
 real*8 :: plgndr, norm_fact
 character (len=3) :: ctrunc
 character (len=2) :: ctrunc1
 integer :: n, m, k, j
 x1 = -1.0D0
 x2 = 1.0D0
 pi = acos(-1.0)
!
 print *,' For a triangular truncation T',mm,' the number of Gaussian latitudes is ',nn  
 
 if (mm < 100) then
   write(ctrunc1,'(i2)') mm
   ctrunc = '0'//ctrunc1
 else
   write(ctrunc,'(i3)') mm
 endif
!
! Compute Gaussian latitudes and Gaussian weights
! 
 call gauleg(x1,x2,x,w,nn)
!
! Compute normalised Legendre polynomials
! 
 do k = 1,nn ! loop over Gaussian latitudes
   j = 0 
   do m = 0,mm  ! loop over zonal wave number
     do n = abs(m),mm+1 ! loop over total wave number (one wave number added)
       j = j + 1
       z1 = log(2.0*n + 1.0)
       z2 = factrl(n - m)
       z3 = factrl(n + m)
       norm_fact = exp(0.5*(z1 + z2 - z3))
       pol_legndr(j,k) = norm_fact*plgndr(n,m,real(x(k)))
     enddo
   enddo
 enddo   
!
! Store Gaussian latitudes and weights + Legendre polynomials
!
 open (unit=20,file='../data_in/gaulat_legpol_T'//ctrunc//'a.dat',status='unknown')
 write (20,*) (x(j),j=1,nn)
 write (20,*) (w(j),j=1,nn)
 write (20,*) ((pol_legndr(j,k),k=1,nn),j=1,mmax) ! if necessary the storage can be reduced by two : only > 0 or < 0 latitudes
 close (unit=20)
 stop
end program compute_legendre_polynomials
