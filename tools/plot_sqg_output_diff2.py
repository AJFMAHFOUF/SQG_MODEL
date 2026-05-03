import matplotlib.ticker as mticker
import matplotlib.pyplot as plt 
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
nlats=48 ; nlons=97
nstep1='0288'
nstep0='0000'
trunc='031'
expid='105'
var='v'
level1='800'
level2='200'
omega = 2.*np.pi/86400.0
#plt.figure(figsize=(10,10))
plt.figure(figsize=(8,4))
           
lon1,lat1,vor1,u1,v1,psi1,phi1,t1=np.loadtxt('../data_out/SQG_T'+trunc+'_lev_'+level1+'_step_'+nstep1+'_expid_'+expid+'.dat',unpack=True)

lon1,lat1,vor0,u0,v0,psi0,phi0,t0=np.loadtxt('../data_out/SQG_T'+trunc+'_lev_'+level1+'_step_'+nstep0+'_expid_'+expid+'.dat',unpack=True)

lon1,lat1,vor3,u3,v3,psi3,phi3,t3=np.loadtxt('../data_out/SQG_T'+trunc+'_lev_'+level2+'_step_'+nstep1+'_expid_'+expid+'.dat',unpack=True)

lon1,lat1,vor2,u2,v2,psi2,phi2,t2=np.loadtxt('../data_out/SQG_T'+trunc+'_lev_'+level2+'_step_'+nstep0+'_expid_'+expid+'.dat',unpack=True)

lat2 = lat1.reshape((nlats,nlons))
lon2 = lon1.reshape((nlats,nlons))

#proj = ccrs.Orthographic(central_longitude=-180.0,central_latitude=60.0)
proj = ccrs.NorthPolarStereo(central_longitude=20.0)
#proj = ccrs.PlateCarree()
#proj = ccrs.LambertConformal(central_longitude=-85.0,central_latitude=45.0,cutoff=0.0)
#proj = ccrs.LambertConformal(central_longitude=-65.0,central_latitude=45.0,cutoff=0.0)
ax = plt.axes(projection=proj)
ax.coastlines('110m',linewidth=0.5,color='purple')
ax.set_extent([-180, 180, 90, 20], ccrs.PlateCarree()) # North hemisphere
#ax.set_extent([-120, -20, 20, 85]) # North America
gl = ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1,color='black',alpha=0.5)
gl.ylocator = mticker.FixedLocator(np.arange(-90,90,30))
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))

if var == 'psi':
	field0 = psi0.reshape((nlats,nlons))
	field1 = psi1.reshape((nlats,nlons))
	field2 = psi2.reshape((nlats,nlons))
	field3 = psi3.reshape((nlats,nlons))
if var == 't':
	field0 = t0.reshape((nlats,nlons))
	field1 = t1.reshape((nlats,nlons))
	field2 = t2.reshape((nlats,nlons))
	field3 = t3.reshape((nlats,nlons))	
if var == 'phi':
	field0 = phi0.reshape((nlats,nlons))
	field1 = phi1.reshape((nlats,nlons))
	field2 = phi2.reshape((nlats,nlons))
	field3 = phi3.reshape((nlats,nlons))
if var == 'vor':
	field0 = vor0.reshape((nlats,nlons))
	field1 = vor1.reshape((nlats,nlons))
	field2 = vor2.reshape((nlats,nlons))
	field3 = vor3.reshape((nlats,nlons))
	#field2 = lat1.reshape((nlats,nlons))
	#field3 = phi.reshape((nlats,nlons))
	#field = (field + 2.*omega*np.sin(np.pi/180.*field2))/field3
if var == 'u':
	field0 = u0.reshape((nlats,nlons))
	field1 = u1.reshape((nlats,nlons))	
	field2 = u2.reshape((nlats,nlons))
	field3 = u3.reshape((nlats,nlons))	
if var == 'v':
	field0 = v0.reshape((nlats,nlons))
	field1 = v1.reshape((nlats,nlons))
	field2 = v2.reshape((nlats,nlons))
	field3 = v3.reshape((nlats,nlons))
if var == 'uv':
	field1 = u.reshape((nlats,nlons))
	field2 = v.reshape((nlats,nlons))

field_a = field1 - field0 
field_b = field3 - field2 

if var == 'phi':
	#cs = map.contourf(x,y,field/9.81,levels=np.arange(4900,5200,25),cmap='jet')	
	#cs = map.contour(x,y,field/98.1,levels=np.arange(460,590,10),colors='black')
	#cs = ax.contour(lon2,lat2,field/9.81,levels=np.arange(10500,12800,100),linewidths=0.1,transform=ccrs.PlateCarree(),colors='black')
	cs = ax.contourf(lon2,lat2,field_a/9.81,transform=ccrs.PlateCarree(),cmap='jet')
	cs = ax.contour(lon2,lat2,field_a/9.81,linewidths=0.5,transform=ccrs.PlateCarree(),colors='black')
	#cs = ax.contour(lon2,lat2,field_b/9.81,linewidths=1.1,transform=ccrs.PlateCarree(),colors='black')
	#cs = ax.contourf(lon2,lat2,field/9.81,levels=np.arange(4600,5900,100),transform=ccrs.PlateCarree(),cmap='jet')
	#cs = map.contour(x,y,field,10,colors='black')
	#cs = map.contourf(x,y,field,10,cmap='seismic')
else:
	#z = np.sqrt(field1*field1 + field2*field2)
	#cs = ax.contour(lon2,lat2,field,levels=np.arange(-50,60,10),transform=ccrs.PlateCarree(),linewidths=0.1,colors='black')
	#cs = ax.contourf(lon2,lat2,field,levels=np.arange(-50,60,10),transform=ccrs.PlateCarree(),cmap='seismic')
	#cs = ax.contour(lon2,lat2,field,levels=np.arange(-25,60,10),linewidths=0.1,transform=ccrs.PlateCarree(),colors='black')
	#cs = ax.contourf(lon2,lat2,field_b,transform=ccrs.PlateCarree(),cmap='seismic')#'RdPu')
	cs = ax.contour(lon2,lat2,field_a,linewidths=1.1,transform=ccrs.PlateCarree(),colors='blue')
	cs = ax.contour(lon2,lat2,field_b,linewidths=1.1,transform=ccrs.PlateCarree(),colors='black')
	#cs = ax.contourf(lon2,lat2,field,transform=ccrs.PlateCarree(),cmap='seismic')#'RdPu')
	#map.quiver(x[::1,::1],y[::1,::1],field1[::1,::1],field2[::1,::1],z[::1,::1],width=0.002,headwidth=2,scale=0.1,scale_units='xy',cmap='jet')
#plt.colorbar(shrink=0.5)
cbar = plt.colorbar(cs,orientation='vertical',shrink=0.5,pad=0.07)
cbar.ax.tick_params()
#level='650'
plt.title('SQG '+var+' diff2 '+level1+' hpa '+nstep1+'h - T'+trunc+' expid='+expid)
plt.savefig('../plots/SQG_'+var+'_'+level1+'_T'+trunc+'_step_'+nstep1+'_expid_'+expid+'_diff2.pdf')
plt.show()
