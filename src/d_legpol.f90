function d_legpol(j1,ms,ns)
 
  use params
  use spectral_vars
  
  implicit none
  
  integer, intent(in)  :: j1 ! index of Gaussian latitude
  integer, intent(in)  :: ms ! zonal wave number 
  integer, intent(in)  :: ns ! total wave number 
  integer :: j_index2, js1, js2
  real    :: d_legpol
  
  if (abs(ms) > mm .or. ns > mm) then
    print *,'Cannot compute derivative of Legendre polynomial for',ms,ns
    stop
  endif  
    
  js1 = j_index2(mm,ms,max0(0,ns-1)) 
  js2 = j_index2(mm,ms,ns+1)
  d_legpol = (ns+1)*eps(ms,ns)*pl_legendr(js1,j1) - & 
            float(ns)*eps(ms,ns+1)*pl_legendr(js2,j1) 
   
  return
 
end function d_legpol   
        
  
