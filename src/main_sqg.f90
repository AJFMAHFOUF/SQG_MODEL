program main_sqg
!======================================================================================
!
! Global spectral Quasi-Geostrophic model in Fortran 90
! with potential vorticity at 3 levels (200, 500, 800 hPa) as prognostic variables
! as proposed by Marshall and Molteni (1993) J.A.S., vol 50, 1792-1818
!
! Explicit leapfrog scheme for advection - backward scheme for dissipation processes 
! Triangular truncation (defined in module params)
!
! Derived from a spectral non-divergent barotropic model (one level)
!
! Jean-François Mahfouf (16/04/2026)
!
! Use of FFT99 defined by Temperton (ECMWF) 
!
!=======================================================================================
 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: nstep, npdt_max
 real    :: t1, t2, dt1 
 logical :: loutput
 
 type (prog_var) :: xin
 type (prog_var) :: xout

! Required initialisations and read fields in physical space

 call init

 call convert_uv2pvor
 
 print *,'Initial fields for PV stored'
 
 call cpu_time(time=t1)
 
! Set-up initial fields in spectral space

 xin%pvormn1 = pvor_mn(:,1,1)
 xin%pvormn2 = pvor_mn(:,2,1)
 xin%pvormn3 = pvor_mn(:,3,1)     
 
! Model integration with initial fields

 dt1 = dt
 npdt_max = npdt + 1
 loutput = .true. 
    
 call model(xin,xout,dt1,npdt_max,loutput)
 
 call cpu_time(time=t2)
 
 print *,'Model execution CPU time =',t2-t1
 
 stop

end program main_sqg
