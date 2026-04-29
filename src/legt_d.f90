subroutine legt_d(am,amn)
 
  use params
  use spectral_vars
  
  implicit none
  
  complex, dimension(nlat,-mm:mm), intent(in) :: am 
  complex, dimension(mmax), intent(out)       :: amn
  integer :: i1, i2, j1, j_index2, js, ms
  
  amn(:) = (0.,0.)
  do i1 = 0,mm
    ms = abs(i1)
    do i2 = ms,mm+1
      js = j_index2(mm,ms,i2)
      do j1 = 1,nlat        
        amn(js) = amn(js) + w(j1)*am(j1,i1)*pl_legendr(js,j1)
      enddo  
    enddo
  enddo
  amn(:) = 0.5*amn(:) 
   
  return
 
end subroutine legt_d   
        
  
