subroutine compute_ke_spectrum(nstep)

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer, intent(in)   :: nstep
 integer :: i1, i2, ms, js, j_index2, ihour
 character(len=1)      :: ichst1
 character(len=2)      :: ichst2
 character(len=3)      :: ichst3
 character(len=4)      :: ichst
 character(len=3)      :: tt
 character(len=2)      :: tt1
 real, dimension(0:mm) :: ke_spectrum, zcoef
 
 ihour = nstep*dt/3600.0
 
 if (ihour < 10) then 
   write(ichst1,'(i1)') ihour
   ichst = '000'//ichst1
 elseif (ihour < 100) then
   write(ichst2,'(i2)') ihour
   ichst = '00'//ichst2
 elseif (ihour < 1000) then
   write(ichst3,'(i3)') ihour
   ichst = '0'//ichst3  
 else
   write(ichst,'(i4)') ihour
 endif      
 
 if (mm < 99) then    
   write(tt1,'(i2)') mm 
   tt = '0'//tt1
 else
   write(tt,'(i3)') mm
 endif  
   
! Compute kinetic energy spectrum (formulation from Lambert, 1984)   
 
 zcoef = 0.0
 do i2 = 1,mm
   zcoef(i2) = a*a/float(i2*(i2+1))
 enddo   
 
 ke_spectrum(:) = 0.0
 do i1 = 0,mm
   ms = abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     if (i1 == 0) then
       ke_spectrum(i2) = ke_spectrum(i2) + 0.5*zcoef(i2)*abs(vor_mn(js,2))**2 ! vorticity at level 2 (500 hPa)
     else
       ke_spectrum(i2) = ke_spectrum(i2) + zcoef(i2)*abs(vor_mn(js,2))**2 
     endif 
   enddo   
 enddo
 
! Store results in file 
 
 open (unit=25,file='../data_out/SQG_spectrum_T'//tt//'_lev_500_step_'//ichst//'_expid_'//expid//'.dat',status='unknown')
 
 do i1=1,mm
   write(25,*) i1,ke_spectrum(i1)
 enddo  
 
 close(unit=25)
      
 return

end subroutine compute_ke_spectrum
