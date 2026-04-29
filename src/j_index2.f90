function j_index2(mm,ms,ns)
 implicit none
 integer :: mm, ms, ns
 integer :: j_index2, js, i
 js = 0 
 do i=0,abs(ms) - 1
   js = js + (mm + 2 - i)
 enddo  
 j_index2 = js + ns + 1 - abs(ms)  
 return
end function j_index2
