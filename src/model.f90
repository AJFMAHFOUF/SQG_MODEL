subroutine model(xin,xout,dt1,npdt_max,loutput)

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 type (prog_var), intent(in)   :: xin
 type (prog_var), intent(out)  :: xout
 real, intent(in)              :: dt1
 integer, intent(in)           :: npdt_max
 logical, intent(in)           :: loutput
 
 complex, dimension(mmax,nlev) :: tend1_pvor_mn, tend2_pvor_mn, filter, pvor2_mn

 integer :: nstep
 
! Put initial conditions in spectral arrays
 
 pvor_mn(:,1,1) = xin%pvormn1
 pvor_mn(:,2,1) = xin%pvormn2
 pvor_mn(:,3,1) = xin%pvormn3

! Start temporal loop

 do nstep = 0,npdt_max-1
       
   call convert_vor2uv
   
!  Tendencies from dynamics (advection of PV by rotational wind)   - use time step 2
 
   call compute_pvorticity_tendency(tend1_pvor_mn)
   
!  Tendencies from dissipative processes - use time step1

   call compute_pvorticity_dissipation(tend2_pvor_mn)   
   
   !tend1_pvor_mn(:,:) = (0.0,0.0)
   !tend2_pvor_mn(:,:) = (0.0,0.0)
   
   if (nstep > 0) then
     
     pvor_mn(:,:,3) = pvor_mn(:,:,1) + 2.0*dt1*(tend1_pvor_mn(:,:) + tend2_pvor_mn(:,:))
     
! Apply horizontal diffusion in spectral space

! 1) Remove stationary components before filtering 

     pvor2_mn(:,1) = pvor_mn(:,1,3) - f_mn(:)
     pvor2_mn(:,2) = pvor_mn(:,2,3) - f_mn(:) 
     pvor2_mn(:,3) = pvor_mn(:,3,3) - f2_mn(:) 
     
     call numerical_diffusion(pvor2_mn(:,:),dt1,0)
     
! 2) Add stationary components after filtering     
     
     pvor_mn(:,1,3) = pvor2_mn(:,1) + f_mn(:)
     pvor_mn(:,2,3) = pvor2_mn(:,2) + f_mn(:) 
     pvor_mn(:,3,3) = pvor2_mn(:,3) + f2_mn(:) 
 
! Apply Robert Asselin Williams filter to remove 2*dt noise
     
     filter(:,:) = pvor_mn(:,:,1) - 2.0*pvor_mn(:,:,2) + pvor_mn(:,:,3)
     pvor_mn(:,:,2) = pvor_mn(:,:,2) + nu*wk*filter(:,:)
     pvor_mn(:,:,3) = pvor_mn(:,:,3) - nu*(1.0-wk)*filter(:,:)
      
! Swap time steps    
   
     pvor_mn(:,:,1) = pvor_mn(:,:,2)
     pvor_mn(:,:,2) = pvor_mn(:,:,3)
        
   else
   
     pvor_mn(:,:,2) = pvor_mn(:,:,1) + dt1*(tend1_pvor_mn(:,:) + tend2_pvor_mn(:,:))
       
   endif     
   
! Derive streamfunction and vorticity from PV inversion

   call inverse_pvorticity   

! Compute U and V wind components at previous time step (for dissipative processes)   
   
   call convert_psi2uv 
  
! Write fields in physical space - spectral transforms in the subroutine
   
   if (mod(nstep,nfreq) == 0 .and. loutput) then
     call compute_balance_equation
     call save_output(nstep)
     call compute_ke_spectrum(nstep)
   endif
   
 enddo  
 
! Put final results in output arrays 

   xout%pvormn1 = pvor_mn(:,1,2)
   xout%pvormn2 = pvor_mn(:,2,2)
   xout%pvormn3 = pvor_mn(:,3,2)
 
 return

end subroutine model
