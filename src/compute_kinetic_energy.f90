subroutine compute_kinetic_energy

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: j1, ilev
 
! Physical space - include surface geopotential to kinetic energy (i.e. total energy)
 
 do ilev=1,nlev
   do j1=1,nlat
     ke(:,j1,ilev) = 0.5*(u(:,j1,ilev)*u(:,j1,ilev) + v(:,j1,ilev)*v(:,j1,ilev))/(a*a*(1.0 - x(j1)*x(j1)))
   enddo
 enddo     
 
! Spectral space

 do ilev=1,nlev
   call fft_d(ke(:,:,ilev),ke_m(:,:,ilev)) 
   call legt_d(ke_m(:,:,ilev),ke_mn(:,ilev))
 enddo
 
 return 

end subroutine compute_kinetic_energy
