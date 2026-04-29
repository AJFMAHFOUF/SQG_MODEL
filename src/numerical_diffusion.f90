subroutine numerical_diffusion(sf_mn,dt1,itype)

 use params
 
 implicit none
 
 complex, dimension(mmax,nlev), intent(inout) :: sf_mn
 real, intent(in)                             :: dt1
 integer, intent(in)                          :: itype

! itype = 1 for vorticity and divergence and 0 for geopotential 
 
 integer :: i1, i2, ms, js, j_index2, ilev

! Apply horizontal diffusion to spectral coefficients - Nabla^4 

 do ilev =1,nlev
   do i1 = 0,mm
     ms = abs(i1)
     do i2 = ms,mm
       js = j_index2(mm,ms,i2)        
       sf_mn(js,ilev) = sf_mn(js,ilev)/(1.0 + 2.0*dt1*kdiff*((i2*(i2+1.0))**2 -4.0*itype)/a**4) 
     enddo  
   enddo
 enddo
      
 return

end subroutine numerical_diffusion
