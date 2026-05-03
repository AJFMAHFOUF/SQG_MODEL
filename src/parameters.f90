module params

 implicit none

 integer, parameter :: mm = 31                ! maximum wave number
 integer, parameter :: nlat = (3*mm+1)/2 + 1  ! number of latitudes
 integer, parameter :: nlon = 2*nlat          ! number of longitudes
 integer, parameter :: mmax = (mm+1)*(mm+4)/2 ! number of stored wavenumbers
 integer, parameter :: nfft = 1               ! number of FFT to be done
 integer, parameter :: nlev = 3               ! number of vertical levels
 
 complex, parameter :: j = (0,1)              ! square root of -1 
 real, parameter    :: a = 6371.22E3          ! Earth radius
 real, parameter    :: pi = acos(-1.0)        ! Pi constant
 real, parameter    :: Rv = 287.0             ! Perfect gas constant for dry air
 real, parameter    :: g = 9.80616            ! Earth gravitational acceleration
 real, parameter    :: omega = 2.0*pi/86164.1 ! Earth angular speed (stellar day)
 real, parameter    :: nu = 0.02, wk = 0.53   ! tunable parameters for 2*dt filter 
 real, parameter    :: kdiff = 5.00E16        ! Coefficient for horizontal diffusion
 real, parameter    :: dt  = 3600.0           ! model time step
 integer, parameter :: nhtot = 300            ! number of hours of model integration
 integer, parameter :: npdt = nhtot*3600/dt   ! number of model time steps
 integer, parameter :: nfreq = 24*3600/dt     ! hourly output archiving frequency
 character(len=3)   :: expid='105'            ! experiment identifier
 character(len=8)   :: cdate='21121978'       ! DDMMYYY : date of initial conditions  
 character(len=2)   :: chour='00'             ! HH : GMT hour of initial conditions
 logical            :: lreaduv=.true.         ! logical to use u v at initial time
 logical            :: l_real_ic=.false.      ! logical to consider true or analytical initial conditions

! Following parameters taken from Marshall and Molteni (1993)

 integer, parameter    :: p1 = 200, p2 = 500, p3 = 800 ! Vertical grid wil 3 pressure levels in hPa 
 real, parameter       :: R1 = 1.0/(6.0E5)**2          ! Rossby radius of deformation (between layers 1 and 2)
 real, parameter       :: R2 = 1.0/(3.0E5)**2          ! Rossby radius of deformation (between layers 2 and 3)
 real, parameter       :: H0 = 9.0E3                   ! Vertical scale of height (m) 
 real, parameter       :: tau_E = 1.0*86400.0          ! Time scale for Ekman dissipation
 real, parameter       :: tau_R = 25.0*86400.0         ! Time scale for temperature relaxation
 
 
end module params 
 
module model_vars  
    
 use params   
 
 implicit none
 
 real, dimension (nlon,nlat,nlev) :: pvor               ! prognostic variable  in physical space
 real, dimension (nlon,nlat,nlev) :: psi, u, v, vor     ! diagnostic variables in physical space
 real, dimension (nlon,nlat,nlev) :: u2, v2             ! winds at previous time step for dissipative processes
 real, dimension (nlon,nlat,nlev) :: phi                ! geopotential (diagnosed from balance equation)
 real, dimension (nlon,nlat,nlev) :: utr, vtr           ! u and v winds in geographical coordinates
 real, dimension (nlon,nlat,nlev) :: uvar, vvar         ! products for advection in physical space
 real, dimension (nlon,nlat,nlev-1) :: t                ! temperature at mid-levels (diagnostic)
 real, dimension (nlev)           :: phibar             ! mean value of geopotential to solve linear balance equation
 
 real, dimension (nlon,nlat)      :: alt                ! orography in m 
 real, dimension (nlon,nlat)      :: lsm                ! land-sea mask
 real, dimension (nlon,nlat)      :: k_ep               ! surface drag coefficient (Ekman pumping)
 
 complex, dimension(nlat,-mm:mm,nlev) :: pvor_m, vor_m
 complex, dimension(nlat,-mm:mm,nlev) :: upvor_m, vpvor_m
 complex, dimension(nlat,-mm:mm,nlev) :: psi_m, u_m, v_m  
 complex, dimension(nlat,-mm:mm,nlev) :: phi_m
 complex, dimension(nlat,-mm:mm,nlev) :: u2_m, v2_m
 
 complex, dimension(mmax,nlev,3) :: pvor_mn               ! prognostic variable in spectral space (3 time steps)
 complex, dimension(mmax,nlev)   :: psi_mn                ! streamfunction (diagnosed)
 complex, dimension(mmax,nlev)   :: u_mn, v_mn, vor_mn    ! U/V wind components + vorticity (diagnosed)
 complex, dimension(mmax,nlev)   :: phi_mn                ! geopotential (diagnosed)
 complex, dimension(mmax)        :: f_mn, f2_mn           ! Coriolis factor and modified Coriolis factor
 complex, dimension(mmax,nlev)   :: psi2_mn, u2_mn, v2_mn ! additional storage of variables at time (t-dt) 
 
 type prog_var
   complex, dimension(mmax) :: pvormn1
   complex, dimension(mmax) :: pvormn2
   complex, dimension(mmax) :: pvormn3
 end type prog_var
 
end module model_vars

module spectral_vars

 use params
 
 implicit none
  
 real, dimension(mmax,nlat) :: pl_legendr ! Legendre polynomials
 real, dimension(nlat)      :: w, x       ! Gaussian weights and latitudes (sinus value)
 real, dimension(nlat)      :: f          ! Coriolis parameter
 
 real, dimension(-mm:mm,0:mm+1) :: eps    ! array for Legendre polynomials derivatives
 
! arrays for FFT991 

 integer, dimension(13)         :: ifax
 real, dimension(3*nlon/2+1)    :: trigs 
 real, dimension(nfft*(nlon+2)) :: acoef
 real, dimension(nfft*(nlon+1)) :: work
 
end module spectral_vars
