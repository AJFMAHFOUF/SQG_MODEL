subroutine legt_i(am,amn,ival)
 
  use params
  use spectral_vars
  
  implicit none
  
  complex, dimension(nlat,-mm:mm), intent(out)  :: am 
  complex, dimension(mmax), intent(in)          :: amn
  integer, intent(in)                           :: ival
  integer :: i1, i2, j1, j_index2, js, ms, mm_max
  
  mm_max = mm
  if (ival == 1) mm_max = mm + 1
  
  am(:,:) = (0.,0.)
  do j1 = 1,nlat
    do i1 = -mm,mm
      ms = abs(i1)      
      do i2 = ms,mm_max
        js = j_index2(mm,ms,i2)
        if (i1 >= 0) then 
          am(j1,i1) = am(j1,i1) + amn(js)*pl_legendr(js,j1)
        else
          am(j1,i1) = am(j1,i1) + conjg(amn(js))*pl_legendr(js,j1)
        endif  
      enddo  
    enddo
  enddo
   
  return
 
end subroutine legt_i  
        
  
