subroutine inverse_pvorticity

! Derive streamfunction and vorticity from potential vorticity 

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 real, dimension (nlev,nlev)          :: zb 
 complex, dimension (nlev)            :: zc
 integer :: js,ms,i2,i1,j_index2
 real    :: zzz, delta
 
! Explicit inversion of a (3x3) linear system A*Psi = (Pvor - f) => Psi = B*(Pvor - f)
! for each wave number (m,n) 
 
 psi_mn(:,:)  = (0.0,0.0)
 vor_mn(:,:)  = (0.0,0.0)
 psi2_mn(:,:) = (0.0,0.0)
 
 do i1 = 0,mm
   ms=abs(i1)    
   do i2 = ms,mm
   js = j_index2(mm,ms,i2)
     if (i2 > 0) then
       zzz = -i2*(i2 + 1.0)/a**2
         
       delta = zzz**3 + 3.0*zzz*R2*R1 - 2.0*zzz**2*(R1 + R2)
       zb(1,1) = R2*(R1 - 2.0*zzz) - zzz*(R1 - zzz)
       zb(1,2) = -R1*(zzz - R2)
       zb(1,3) = R1*R2
       zb(2,1) = -R1*(zzz - R2)
       zb(2,2) = (zzz - R1)*(zzz - R2)
       zb(2,3) = -R2*(zzz - R1)
       zb(3,1) = R1*R2
       zb(3,2) = -R2*(zzz - R1)
       zb(3,3) = R1*(R2 - 2.0*zzz) - zzz*(R2 - zzz)
       
       zb(:,:) = zb(:,:)/delta 
 
! Use PV at middle time step for dynamics with leapfrog scheme  
     
       zc(1) = pvor_mn(js,1,2) - f_mn(js)
       zc(2) = pvor_mn(js,2,2) - f_mn(js)
       zc(3) = pvor_mn(js,3,2) - f2_mn(js)
       
       psi_mn(js,1) = zc(1)*zb(1,1) + zc(2)*zb(2,1) + zc(3)*zb(3,1) 
       psi_mn(js,2) = zc(1)*zb(1,2) + zc(2)*zb(2,2) + zc(3)*zb(3,2)
       psi_mn(js,3) = zc(1)*zb(1,3) + zc(2)*zb(2,3) + zc(3)*zb(3,3)
         
       vor_mn(js,:) = psi_mn(js,:)*zzz      

! Stream function at previous time step for dissipative processes      
       
       zc(1) = pvor_mn(js,1,1) - f_mn(js)
       zc(2) = pvor_mn(js,2,1) - f_mn(js)
       zc(3) = pvor_mn(js,3,1) - f2_mn(js)
       
       psi2_mn(js,1) = zc(1)*zb(1,1) + zc(2)*zb(2,1) + zc(3)*zb(3,1) 
       psi2_mn(js,2) = zc(1)*zb(1,2) + zc(2)*zb(2,2) + zc(3)*zb(3,2)
       psi2_mn(js,3) = zc(1)*zb(1,3) + zc(2)*zb(2,3) + zc(3)*zb(3,3)
       
     endif  
   enddo
 enddo  
 
 return 

end subroutine inverse_pvorticity
