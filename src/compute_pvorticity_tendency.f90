subroutine compute_pvorticity_tendency(tend_pvor_mn)

 use params
 use model_vars
 use spectral_vars 
 
 implicit none
 
 complex, dimension(mmax,nlev), intent(out)  :: tend_pvor_mn
 real    :: d_legpol
 integer :: i1, i2, j1, j_index2, js, ms, ilev
 
! Physical space  (absolute vorticity)
 
 do ilev = 1,nlev
   do j1 = 1,nlat
     do i1  = 1,nlon
       uvar(i1,j1,ilev) = u(i1,j1,ilev)*pvor(i1,j1,ilev)
       vvar(i1,j1,ilev) = v(i1,j1,ilev)*pvor(i1,j1,ilev)
     enddo
   enddo
 enddo      
 
! Fourier space for non linear terms

 do ilev = 1,nlev
   call fft_d(uvar(:,:,ilev),upvor_m(:,:,ilev)) 
   call fft_d(vvar(:,:,ilev),vpvor_m(:,:,ilev))
 enddo 

! Compute tendencies in spectral space

 do ilev = 1,nlev
   tend_pvor_mn(:,ilev) = (0.,0.)
   do i1 = 0,mm
     ms = abs(i1)
     do i2 = ms,mm 
       js = j_index2(mm,ms,i2)
       do j1 = 1,nlat    
         tend_pvor_mn(js,ilev) = tend_pvor_mn(js,ilev) - (j*i1*upvor_m(j1,i1,ilev)*pl_legendr(js,j1) - &
        &                        vpvor_m(j1,i1,ilev)*d_legpol(j1,i1,i2))*w(j1)/(a*a*(1.0 - x(j1)*x(j1)))                 
       enddo
     enddo
   enddo 
 enddo                          
 tend_pvor_mn(:,:) = 0.5*tend_pvor_mn(:,:) 
 
 return 

end subroutine compute_pvorticity_tendency
