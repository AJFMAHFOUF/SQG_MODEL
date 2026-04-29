subroutine save_output(nstep)

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer, intent(in) :: nstep
 character(len=1)    :: ichst1
 character(len=2)    :: ichst2
 character(len=3)    :: ichst3
 character(len=4)    :: ichst
 character(len=3)    :: tt
 character(len=2)    :: tt1
 character(len=3), dimension(3) :: plev   
 
 character(len=3), dimension(3) :: lev1
 integer :: i1, j1, i2, ms, js, j_index2, ihour, its, ilev
 real :: zlon, zlat, zu, zv, zfac  
  
 ihour = nstep*dt/3600.0
 
! This index change allows the initial fields to be stored after first time step 
 
 if (nstep == 0) then
   its = 1
 else
   its = 2
 endif    
 
 if (ihour < 10) then 
   write(ichst1,'(i1)') ihour
   ichst = '000'//ichst1
 elseif (ihour < 100) then
   write(ichst2,'(i2)') ihour
   ichst = '00'//ichst2
 elseif (ihour < 1000) then
   write(ichst3,'(i3)') ihour
   ichst = '0'//ichst3   
 else
   write(ichst,'(i4)') ihour
 endif      
    
 if (mm < 99) then    
   write(tt1,'(i2)') mm 
   tt = '0'//tt1
 else
   write(tt,'(i3)') mm
 endif    
 
 write (plev(1),'(i3)') p1
 write (plev(2),'(i3)') p2
 write (plev(3),'(i3)') p3 
 
! Back to physical space - vorticity

 do ilev = 1,nlev
   call legt_i(vor_m(:,:,ilev),vor_mn(:,ilev),0)
   call fft_i(vor(:,:,ilev),vor_m(:,:,ilev)) 
 enddo
 
! Back to physical space - geopotential

! do ilev = 1,nlev
!   call legt_i(phi_m(:,:,ilev),phi_mn(:,ilev),0)
!   call fft_i(phi(:,:,ilev),phi_m(:,:,ilev))
! enddo 
  
! Back to physical space - wind components

 do ilev = 1,nlev 
   call legt_i(u_m(:,:,ilev),u_mn(:,ilev),1)
   call fft_i(u(:,:,ilev),u_m(:,:,ilev)) 
  
   call legt_i(v_m(:,:,ilev),v_mn(:,ilev),1)
   call fft_i(v(:,:,ilev),v_m(:,:,ilev)) 
 enddo
 
! Compute streamfunction 

 do ilev = 1,nlev
   psi_mn(:,ilev) = (0.,0.)
   do i1 = 0,mm
     ms=abs(i1)
     do i2 = ms,mm
       js = j_index2(mm,ms,i2)
       if (i2 > 0) then
         psi_mn(js,ilev) = (-a*a/(i2*(i2+1.0)))*vor_mn(js,ilev)
       endif  
     enddo
   enddo
 enddo  
 
! Back to physical space -  psi

 do ilev = 1,nlev
   call legt_i(psi_m,psi_mn,0)
   call fft_i(psi,psi_m)
 enddo 
 
! Write results in ASCII file for plotting purposes
  
 do ilev = 1,nlev 
 open (unit=20+ilev-1,file='../data_out/SQG_T'//tt//'_lev_'//plev(ilev)//'_step_'//ichst//'_expid_'//expid//'.dat',status='unknown')
   do j1=1,nlat 
     zfac = 1.0/(a*sqrt(1.0 - x(j1)*x(j1)))
     do i1=1,nlon
       zlon = float(i1-1)/float(nlon)*360.0
       zlat = asin(x(j1))*180.0/pi
       zu = u(i1,j1,ilev)*zfac   ! real U wind
       zv = v(i1,j1,ilev)*zfac   ! real V wind
       write(20+ilev-1,*) zlon,zlat,vor(i1,j1,ilev),zu,zv,psi(i1,j1,ilev),phi(i1,j1,ilev),t(i1,j1,min(ilev,2))
     enddo
     zu = u(1,j1,ilev)*zfac
     zv = v(1,j1,ilev)*zfac
     write(20+ilev-1,*) 360.0,zlat,vor(1,j1,ilev),zu,zv,psi(1,j1,ilev),phi(1,j1,ilev),t(1,j1,min(ilev,2))
   enddo  
   close (unit=20+ilev-1) 
 enddo  
  
 return  
 
end subroutine save_output
