subroutine convert_uv2pvor

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: i1, i2, j1, ms, js, j_index2, ilev, ilat, ilon
 real    :: zvor, ztheta0, zthetaw, zm, za, zlon, zlat
 real    :: d_legpol
 
 real, dimension(nlon,nlat)      :: zf, zf2
 complex, dimension(nlat,-mm:mm) :: zf_m, zf2_m 

! Spectral coefficients for Coriolis parameter + modified scaled orography (1+z/H0)

 do ilat = 1,nlat
   zf(:,ilat) = f(ilat)
 enddo
 call fft_d(zf,zf_m)    
 call legt_d(zf_m,f_mn) 
 
 do ilat = 1,nlat
   do ilon = 1,nlon
     zf2(ilon,ilat) = f(ilat)*(1.0 + alt(ilon,ilat)/H0)
   enddo
 enddo  
 call fft_d(zf2,zf2_m)    
 call legt_d(zf2_m,f2_mn)  

! Spectral coefficients of u and v components (transformed fields)

 do ilev = 1,nlev
   call fft_d(u(:,:,ilev),u_m(:,:,ilev))
   call fft_d(v(:,:,ilev),v_m(:,:,ilev))
 enddo

! Spectral coefficients for vorticity 

 do ilev = 1,nlev
   vor_mn(:,ilev) = (0.,0.)
   do i1 = 0,mm
     ms = abs(i1)
     do i2 = ms,mm
       js = j_index2(mm,ms,i2)
       do j1 = 1,nlat        
         vor_mn(js,ilev) = vor_mn(js,ilev) + w(j1)*(j*i1*v_m(j1,i1,ilev)*pl_legendr(js,j1) + &
       &                   u_m(j1,i1,ilev)*d_legpol(j1,i1,i2))/(a*a*(1.0-x(j1)*x(j1)))
       enddo  
     enddo
   enddo
   vor_mn(:,ilev) = 0.5*vor_mn(:,ilev)
 enddo
 
! Spectral coefficients for streamfunction 

 psi_mn(:,:) = (0.,0.)
 do ilev = 1,nlev
   do i1 = 0,mm
     ms=abs(i1)
     do i2 = ms,mm
       js = j_index2(mm,ms,i2)
       if (i2 > 0) then 
         psi_mn(js,ilev) = -a*a/(i2*(i2+1.0))*vor_mn(js,ilev)
       endif  
     enddo
   enddo  
 enddo  
 
! Spectral coefficients for potential vorticity

 pvor_mn(:,:,:) = (0.,0.)
   
! First level - 200 hPa   

 ilev = 1 
 do i1 = 0,mm
   ms=abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     if (i2 > 0) then 
       pvor_mn(js,ilev,1) = vor_mn(js,ilev) - (psi_mn(js,ilev) - psi_mn(js,ilev+1))*R1 + f_mn(js)
     endif  
   enddo
 enddo  
!
! Second level - 500 hPa
!   
 ilev = 2 
 do i1 = 0,mm
   ms=abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     if (i2 > 0) then 
       pvor_mn(js,ilev,1) = vor_mn(js,ilev) + (psi_mn(js,ilev-1) - psi_mn(js,ilev))*R1 &
      &                   - (psi_mn(js,ilev) - psi_mn(js,ilev+1))*R2 + f_mn(js)
     endif  
   enddo
 enddo  
!   
! Third level - 800 hPa
!   
 ilev = 3 
 do i1 = 0,mm
   ms=abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     if (i2 > 0) then 
       pvor_mn(js,ilev,1) = vor_mn(js,ilev) + (psi_mn(js,ilev-1) - psi_mn(js,ilev))*R2 + f2_mn(js) 
     endif  
   enddo
 enddo       
 
!  For initial step
  
! if (.not.l_real_ic) then 
!   call legt_i(vor_m,vor_mn(:,2),0)
!   call fft_i(vor,vor_m)
 
! Analytical initial conditions proposed by Held and Phillips (1987) JAS
 
!   ztheta0 = 0.25*pi
!   zthetaw = 15.0*pi/180. 
!   zm = 4.0
!   za = 8.0E-5
   
!   do j1=1,nlon
!     do i1=1,nlat
!       zlon = (j1 - 1.0)/float(nlon)*2.0*pi
!       zlat = asin(x(i1))
!       zvor = 0.5*za*sqrt(1.0 - x(i1)*x(i1))*exp(-((zlat -ztheta0)/zthetaw)**2)*cos(zm*zlon)
!       vor(j1,i1) = vor(j1,i1) + zvor
!     enddo
!   enddo  

!  call fft_d(vor,vor_m)    
!  call legt_d(vor_m,vor_mn(:,2))
   
! endif  

 pvor_mn(:,:,2) = pvor_mn(:,:,1)
 pvor_mn(:,:,3) = pvor_mn(:,:,1)
      
 return

end subroutine convert_uv2pvor
