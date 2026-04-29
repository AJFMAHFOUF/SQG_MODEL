subroutine convert_vor2uv

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: i1, i2, ms, js, jsm, jsp, j_index2, ilev
   
! Physical fields required for non linear terms   
  
 do ilev = 1,nlev
   call legt_i(vor_m(:,:,ilev),vor_mn(:,ilev),0)
   call fft_i(vor(:,:,ilev),vor_m(:,:,ilev))   
 enddo
     
! Stream function (spectral)

 do ilev = 1,nlev
   psi_mn(:,ilev) = (0.,0.)
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
!
!  U and V wind components (spectral)
!
 do ilev = 1,nlev
   u_mn(:,ilev) = (0.,0.)
   v_mn(:,ilev) = (0.,0.)
   do i1 = 0,mm
     ms=abs(i1)
     do i2 = ms,mm
       js  = j_index2(mm,ms,i2)
       jsm = j_index2(mm,ms,max0(0,i2-1))
       jsp = j_index2(mm,ms,i2+1)
       u_mn(js,ilev) = (i2-1)*eps(i1,i2)*psi_mn(jsm,ilev) - (i2+2)*eps(i1,i2+1)*psi_mn(jsp,ilev)   
       v_mn(js,ilev) = j*i1*psi_mn(js,ilev)  
     enddo
     i2 = mm + 1
     js  = j_index2(mm,ms,i2)
     jsm = j_index2(mm,ms,i2-1)
     u_mn(js,ilev) =   (i2-1)*eps(i1,i2)*psi_mn(jsm,ilev)    
     v_mn(js,ilev) =   0.0
   enddo  
 enddo

! Back to physical space (Wind components and potential vorticity)

  do ilev = 1,nlev
    call legt_i(u_m(:,:,ilev),u_mn(:,ilev),1)
    call fft_i(u(:,:,ilev),u_m(:,:,ilev)) 
  
    call legt_i(v_m(:,:,ilev),v_mn(:,ilev),1)
    call fft_i(v(:,:,ilev),v_m(:,:,ilev))
    
    call legt_i(pvor_m(:,:,ilev),pvor_mn(:,ilev,2),1)
    call fft_i(pvor(:,:,ilev),pvor_m(:,:,ilev))   
  enddo 
 
 return

end subroutine convert_vor2uv
