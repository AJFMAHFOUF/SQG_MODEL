subroutine fft_i(b,bm)
 
  use fft99_mod
  use params
  use spectral_vars
  
  implicit none
  
  real, dimension(nlon,nlat), intent(out)      :: b
  complex, dimension(nlat,-mm:mm), intent(in)  :: bm 
  real, dimension(nfft*(nlon+2))               :: zcoef
  integer :: j1, i1
  
  do j1 = 1,nlat
    zcoef(:) = 0.0
    do i1 = 1,2*mm+1,2
       zcoef(i1)   = 0.5*real(bm(j1,(i1-1)/2) + bm(j1,(1-i1)/2))
       zcoef(i1+1) = 0.5*aimag(bm(j1,(i1-1)/2) - bm(j1,(1-i1)/2))    
    enddo
    call fft991(zcoef,work,trigs,ifax,1,0,nlon,nfft,1)
    b(1:nlon,j1) = zcoef(1:nlon)
  enddo
   
  return
 
end subroutine fft_i 
