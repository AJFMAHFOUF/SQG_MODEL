subroutine fft_d(b,bm)
 
  use fft99_mod
  use params
  use spectral_vars
  
  implicit none
  
  real, dimension(nlon,nlat), intent(in)       :: b
  complex, dimension(nlat,-mm:mm), intent(out) :: bm 
  real, dimension(nlon+2)                      :: zcoef
  integer :: j1, i1
  
  zcoef(:) = 0.0
  do j1 = 1,nlat
    zcoef(1:nlon) = b(1:nlon,j1)
    call fft991(zcoef,work,trigs,ifax,1,0,nlon,nfft,-1)
    do i1 = 1,2*mm+1,2
      bm(j1,(i1-1)/2) = zcoef(i1) + j*zcoef(i1+1)
      bm(j1,(1-i1)/2) = zcoef(i1) - j*zcoef(i1+1)  
    enddo
  enddo
   
  return
 
end subroutine fft_d   
        
  
