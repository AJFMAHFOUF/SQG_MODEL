function plgndr(l,m,x)  
!
! Computes the associated Legendre polynomials P_l,m(x)
! Here m and l are integers satisfying 0 <= m <= l while
! x lies in the range [-1, +1]
!
! Taken from Numerical Recipes - pages 247-248
!
! Jean-François Mahfouf (30/05/2023) 
!
 implicit none
 integer, intent(in) :: l,m
 real, intent(in)    :: x
 integer :: i, ll
 real :: plgndr, fact, pll, pmm, pmmp1, somx2
 if (m < 0 .or. m > l .or. abs(x) > 1.0) then
   print *,'bad argument in plgndr',' m=',m,' l=',l, 'x=',x
   stop
 endif
 pmm = 1.0
 if (m > 0) then
   somx2 = sqrt((1.0 - x)*(1.0 + x))
   fact = 1.0
   do i=1,m
     pmm = -pmm*fact*somx2
     fact = fact + 2.0
   enddo
 endif
 if (l == m) then
   plgndr = pmm
 else
   pmmp1 = x*float(2*m + 1)*pmm
   if (l == m+1) then
     plgndr = pmmp1
   else
     do ll=m+2,l
       pll = (x*(2*ll - 1)*pmmp1 - (ll + m -1)*pmm)/float(ll - m)
       pmm = pmmp1
       pmmp1 = pll
     enddo
     plgndr = pll
   endif  
 endif
 return
end function plgndr           
