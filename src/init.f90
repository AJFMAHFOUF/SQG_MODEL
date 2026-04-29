subroutine init
 
 use fft99_mod
 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer            :: i1, i2, j1, ilev
 real               :: zlon, zlat, zzz, zweight
 real, dimension(3) :: zfield
 character(len=3)   :: tt
 character(len=2)   :: tt1
 character(len=3), dimension(3) :: plev   
 
 real, dimension(nlat) :: z_umean, z_vmean

 
! Define character for truncation

 if (mm < 99) then    
   write(tt1,'(i2)') mm 
   tt = '0'//tt1
 else
   write(tt,'(i3)') mm
 endif     
 
 write (plev(1),'(i3)') p1
 write (plev(2),'(i3)') p2
 write (plev(3),'(i3)') p3
 
! Initialisation prior FFT

 call set99(trigs,ifax,nlon)
 
! Define array for recurrence formula (Legendre polynomials)

 eps(:,:) = 0.0
 do i1 = -mm,mm
   do i2 = abs(i1),mm+1
     if (i1 /= i2) eps(i1,i2) = sqrt(float(i2*i2 - i1*i1)/float(4*i2*i2 - 1))
   enddo
 enddo 

! Read initial physical fields (geopotential, u and v winds)  + orography
     
 open (unit=11,file='../data_in/U_'//cdate//'_'//chour//'_'//plev(1)//'_T'//tt//'gg.dat',status='old')
 open (unit=12,file='../data_in/V_'//cdate//'_'//chour//'_'//plev(1)//'_T'//tt//'gg.dat',status='old') 
 open (unit=13,file='../data_in/PHI_'//cdate//'_'//chour//'_'//plev(1)//'_T'//tt//'gg.dat',status='old') 
 open (unit=14,file='../data_in/U_'//cdate//'_'//chour//'_'//plev(2)//'_T'//tt//'gg.dat',status='old')
 open (unit=15,file='../data_in/V_'//cdate//'_'//chour//'_'//plev(2)//'_T'//tt//'gg.dat',status='old') 
 open (unit=16,file='../data_in/PHI_'//cdate//'_'//chour//'_'//plev(2)//'_T'//tt//'gg.dat',status='old') 
 open (unit=17,file='../data_in/U_'//cdate//'_'//chour//'_'//plev(3)//'_T'//tt//'gg.dat',status='old')
 open (unit=18,file='../data_in/V_'//cdate//'_'//chour//'_'//plev(3)//'_T'//tt//'gg.dat',status='old') 
 open (unit=19,file='../data_in/PHI_'//cdate//'_'//chour//'_'//plev(3)//'_T'//tt//'gg.dat',status='old')
 open (unit=20,file='../data_in/ALT_'//cdate//'_'//chour//'_T'//tt//'gg.dat',status='old')  
 open (unit=21,file='../data_in/LSM_'//cdate//'_'//chour//'_T'//tt//'gg.dat',status='old')  

 phibar(:) = 0.0
 do ilev = 1,nlev
   zweight = 0.0
   !z_umean(:) = 0.0
   !z_vmean(:) = 0.0
   do j1 = 1,nlat   
     do i1 = 1,nlon+1
       read(11+(ilev-1)*3,*) zlon, zlat, zfield(1)
       read(12+(ilev-1)*3,*) zlon, zlat, zfield(2)
       read(13+(ilev-1)*3,*) zlon, zlat, zfield(3)
       if (i1 /= nlon+1) then
         zzz=cos(pi*zlat/180.)
         !z_umean(j1) = z_umean(j1) + zfield(1)
         !z_vmean(j1) = z_vmean(j1) + zfield(2)
         utr(i1,j1,ilev) = zfield(1)
         vtr(i1,j1,ilev) = zfield(2) 
         phibar(ilev) = phibar(ilev) + zfield(3)*zzz 
         zweight = zweight + zzz 
!         if (.not.l_real_ic) then
!           utr(i1,j1) = 25.0*zzz - 30.0*zzz**3 + 300.0*(1.0-zzz**2)*zzz**6
!           vtr(i1,j1) = 0.0
!         endif
       endif
     enddo 
   enddo  
   !do i1 = 1,nlon
   !  utr(i1,:,ilev) = z_umean(:)/float(nlon)
   !  vtr(i1,:,ilev) = z_vmean(:)/float(nlon)
   !enddo  
 enddo
 phibar(:) = phibar(:)/zweight
 
! Read orography

 do j1 = 1,nlat  
   do i1 = 1,nlon+1
     read(20,*) zlon, zlat, zfield(1)
     if (i1 /= nlon+1) then
       alt(i1,j1) = zfield(1)/g 
     endif
   enddo
 enddo  

! Read land-sea mask 
 
 do j1 = 1,nlat  
   do i1 = 1,nlon+1
     read(21,*) zlon, zlat, zfield(1)
     if (i1 /= nlon+1) then
       lsm(i1,j1) = zfield(1)
     endif
   enddo
 enddo 
 
! Compute surface drag coefficient depending upon land-sea mask and orography 
 
 do j1 = 1,nlat
   do i1 = 1,nlon
     k_ep(i1,j1) = (1.0 + 0.5*lsm(i1,j1) + 0.5*(1.0-exp(-(max(0.0,alt(i1,j1))/1000.0))))/tau_E
   enddo
 enddo   
 
 print *,'Initial fields read from input files - mean Z @ 3 levels (km) ', phibar(:)/g*1E-3
 
 close(unit=11)
 close(unit=12)
 close(unit=13)
 close(unit=14)
 close(unit=15)
 close(unit=16)
 close(unit=17)
 close(unit=18)
 close(unit=19)
 close(unit=20)
 close(unit=21)

! Read Gaussian latitudes (sin) a and Gaussian weights w 
 
 open (30,file='../data_in/gaulat_legpol_T'//tt//'a.dat',status='old')
 
 read(30,*) (x(j1),j1=1,nlat)
 read(30,*) (w(j1),j1=1,nlat)
 
 print *,'Gaussian latitudes and weights have been read'
 
 read(30,*) ((pl_legendr(i1,j1),j1=1,nlat),i1=1,mmax)
 
 print *,'Legendre polynomials have been read'    
       
 close(unit=30)
 
! Compute Coriolis parameter and transformed winds
 
 do j1 = 1,nlat
   f(j1) = 2.*omega*x(j1)
   u(:,j1,:) = utr(:,j1,:)*a*sqrt(1.0 - x(j1)*x(j1))
   v(:,j1,:) = vtr(:,j1,:)*a*sqrt(1.0 - x(j1)*x(j1))
 enddo
 
! Define U/V winds at previous time step

  u2(:,:,:) = u(:,:,:)
  v2(:,:,:) = v(:,:,:) 
 
! Spectral coefficients for vorticity (assumes vorticity field defined)

! call fft_d(vor,vor_m)    
! call legt_d(vor_m,vor_mn(:,2))
                                                                                                                                   
! Fill other arrays
   
! vor_mn(:,1) = vor_mn(:,2)
! vor_mn(:,3) = vor_mn(:,2)   
 
 print *,'Exit from initialisation subroutine'      
      
 return     
end subroutine init
