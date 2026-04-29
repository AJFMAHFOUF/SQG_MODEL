subroutine convert_psi2uv

! Compute U and V components @ timestep (1) for dissipative processes
! using streamfunction psi2_mn obtained for PV inversion 

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: i1, i2, ms, js, jsm, jsp, j_index2, ilev 
!
!  U and V wind components (spectral)
!
 do ilev = 1,nlev
   u2_mn(:,ilev) = (0.,0.)
   v2_mn(:,ilev) = (0.,0.)
   do i1 = 0,mm
     ms=abs(i1)
     do i2 = ms,mm
       js  = j_index2(mm,ms,i2)
       jsm = j_index2(mm,ms,max0(0,i2-1))
       jsp = j_index2(mm,ms,i2+1)
       u2_mn(js,ilev) = (i2-1)*eps(i1,i2)*psi2_mn(jsm,ilev) - (i2+2)*eps(i1,i2+1)*psi2_mn(jsp,ilev)   
       v2_mn(js,ilev) = j*i1*psi2_mn(js,ilev)  
     enddo
     i2 = mm + 1
     js  = j_index2(mm,ms,i2)
     jsm = j_index2(mm,ms,i2-1)
     u2_mn(js,ilev) = (i2-1)*eps(i1,i2)*psi2_mn(jsm,ilev)    
     v2_mn(js,ilev) = (0.0,0.0)
   enddo  
 enddo

! Back to physical space wind components 

  do ilev = 1,nlev
    call legt_i(u2_m(:,:,ilev),u2_mn(:,ilev),1)
    call fft_i(u2(:,:,ilev),u2_m(:,:,ilev)) 
  
    call legt_i(v2_m(:,:,ilev),v2_mn(:,ilev),1)
    call fft_i(v2(:,:,ilev),v2_m(:,:,ilev))
  enddo 
 
 return

end subroutine convert_psi2uv
