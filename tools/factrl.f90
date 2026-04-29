function factrl(n)
!
! Computes the log of a factorial 
!
 implicit none
 integer :: n
 real :: factrl
 integer :: j, ntop
 real :: a(33), gammln
 save ntop, a
 data ntop,a(1)/0,0.0/
 if (n < 0) then
   print*, 'negative factorial in factrl'
   return
 elseif (n < ntop) then
   factrl = a(n+1)
 elseif (n < 32) then
   do j=ntop+1,n
     a(j+1) = log(float(j)) + a(j)
   enddo
   ntop = n
   factrl = a(n+1)
 else
   factrl = gammln(n+1.0)
 endif
 return     
end function factrl
