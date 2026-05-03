subroutine compute_pvorticity_dissipation(tend_pvor_mn)

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 complex, dimension(mmax,nlev), intent(out)  :: tend_pvor_mn
 integer :: i1, i2, j1, ms, js, j_index2
 real    :: d_legpol

! Physical space  
 
 do j1 = 1,nlat
   do i1  = 1,nlon
     uvar(i1,j1,nlev) = u2(i1,j1,nlev)*k_ep(i1,j1)
     vvar(i1,j1,nlev) = v2(i1,j1,nlev)*k_ep(i1,j1)
   enddo
 enddo   
 
! Spectral coefficients for k*U and k*V products

 call fft_d(uvar(:,:,nlev),upvor_m(:,:,nlev))
 call fft_d(vvar(:,:,nlev),vpvor_m(:,:,nlev))  

! Spectral coefficients for Ekman dissipation (level 3 only)

 tend_pvor_mn(:,:) = (0.,0.)
 do i1 = 0,mm
   ms = abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     do j1 = 1,nlat        
       tend_pvor_mn(js,nlev) = tend_pvor_mn(js,nlev) + w(j1)*(j*i1*vpvor_m(j1,i1,nlev)*pl_legendr(js,j1) + &
                               upvor_m(j1,i1,nlev)*d_legpol(j1,i1,i2))/(a*a*(1.0-x(j1)*x(j1)))
     enddo  
   enddo
 enddo
 tend_pvor_mn(:,nlev) = -0.5*tend_pvor_mn(:,nlev)
 
! Include temperature relaxation contribution  
 
 tend_pvor_mn(:,1) = R1*(psi2_mn(:,1) - psi2_mn(:,2))/tau_R
 tend_pvor_mn(:,2) = R2*(psi2_mn(:,2) - psi2_mn(:,3))/tau_R - R1*(psi2_mn(:,1) - psi2_mn(:,2))/tau_R
 tend_pvor_mn(:,3) = tend_pvor_mn(:,3) - R2*(psi2_mn(:,2) - psi2_mn(:,3))/tau_R
 
 tend_pvor_mn(:,:) = (0.0,0.0)
      
 return

end subroutine compute_pvorticity_dissipation
