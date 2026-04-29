subroutine compute_balance_equation

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 real, dimension (nlon,nlat,nlev)     :: linbal
 complex, dimension(nlat,-mm:mm,nlev) :: linbal_m
 complex, dimension(mmax,nlev)        :: linbal_mn

 integer :: i1,i2,js,jsm,jsp,ms,j1,j_index2,ilev


! Compute u wind component from updated streamfunction (after PV inversion)

 do ilev = 1,nlev
   u_mn(:,ilev) = (0.,0.)
   do i1 = 0,mm
     ms=abs(i1)
     do i2 = ms,mm
       js  = j_index2(mm,ms,i2)
       jsm = j_index2(mm,ms,max0(0,i2-1))
       jsp = j_index2(mm,ms,i2+1)
       u_mn(js,ilev) = (i2-1)*eps(i1,i2)*psi_mn(jsm,ilev) - (i2+2)*eps(i1,i2+1)*psi_mn(jsp,ilev)   
     enddo
     i2 = mm + 1
     js  = j_index2(mm,ms,i2)
     jsm = j_index2(mm,ms,i2-1)
     u_mn(js,ilev) = (i2-1)*eps(i1,i2)*psi_mn(jsm,ilev)    
   enddo  
 enddo

! Back to physical space u wind component

  do ilev = 1,nlev
    call legt_i(u_m(:,:,ilev),u_mn(:,ilev),1)
    call fft_i(u(:,:,ilev),u_m(:,:,ilev)) 
  enddo 
! 
! Physical space - compute div (f x nabla psi) - Eq.(7.5.3) from Daley (1991)!
! 
 do ilev = 1,nlev
   do j1 = 1,nlat
     linbal(:,j1,ilev) = vor(:,j1,ilev)*f(j1)  - 2.0*u(:,j1,ilev)*omega/(a*a)
   enddo  
 enddo
 
! Spectral space

 do ilev = 1,nlev
   call fft_d(linbal(:,:,ilev),linbal_m(:,:,ilev)) 
   call legt_d(linbal_m(:,:,ilev),linbal_mn(:,ilev))
 enddo
 
! Solve linear balance equation to obtain the geopotential (inverse Laplacian)

 do ilev = 1,nlev
   phi_mn(1,ilev) = cmplx(phibar(ilev))  ! mean geopotential
   do i1 = 0,mm
     ms=abs(i1)
     do i2 = ms,mm
       js = j_index2(mm,ms,i2)
       if (i2 > 0) then
         phi_mn(js,ilev) = -a*a/(i2*(i2+1.0))*linbal_mn(js,ilev)     
       endif  
     enddo
   enddo
 enddo   

! Back to physical space to compute temperature from Phi vertical gradient 
 
 do ilev = 1,nlev
   call legt_i(phi_m(:,:,ilev),phi_mn(:,ilev),0)
   call fft_i(phi(:,:,ilev),phi_m(:,:,ilev))
 enddo
 
! Temperature diagnostic - hydrostatic equation - could be done in spectral space (linear transform)
 
 do i1 = 1,nlon
   do j1 = 1,nlat
     t(i1,j1,1) = (phi(i1,j1,1) - phi(i1,j1,2))/(Rv*(log(p2*100.0) - log(p1*100.0)))
     t(i1,j1,2) = (phi(i1,j1,2) - phi(i1,j1,3))/(Rv*(log(p3*100.0) - log(p2*100.0)))  
   enddo
 enddo
 return 

end subroutine compute_balance_equation
