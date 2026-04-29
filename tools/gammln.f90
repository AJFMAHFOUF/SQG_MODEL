function gammln(xx)
 implicit none
 real :: gammln, xx
 integer :: j
 real*8 :: ser, stp, tmp, x, y, cof(6)
 save cof, stp
 data cof,stp/76.18009172947146D0,-86.50532032941677D0, &
 &            24.01409824083091D0,-1.231739572450155D0, &
 &            0.1208650973866179D-2,-0.5395239384953D-5, &
 &            2.5066282746310005D0/
 x = xx
 y = x
 tmp = x + 5.5d0
 tmp = (x + 0.5d0)*log(tmp) - tmp
 ser = 1.0000000001900015D0
 do j = 1,6
   y = y + 1.0D0
   ser = ser + cof(j)/y
 enddo
 gammln = tmp + log(stp*ser/x) 
 return
end function gammln
