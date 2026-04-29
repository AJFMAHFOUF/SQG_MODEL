import matplotlib.ticker as mticker
import matplotlib.pyplot as plt 
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
nlats=48 ; nlons=97
nstep='0024'
trunc='031'
var='PHI'
date='26121978'
hour='00'
level='500'
omega = 2.*np.pi/86400.0
#plt.figure(figsize=(10,10))
plt.figure(figsize=(8,4))

#map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
#            llcrnrlon=0,urcrnrlon=360,resolution='c')
#map = Basemap(projection='ortho',lat_0=40,lon_0=-90,resolution='l')  
#map = Basemap(projection='stere',width=18000000,height=18000000,lat_0=90.0,lon_0=-90.0,lat_ts=60.0,resolution='l')
#map.drawcoastlines(linewidth=0.25)
#map.drawcountries(linewidth=0.25)
#map.drawmeridians(np.arange(0,360,30))
#map.drawparallels(np.arange(-90,90,30)) 
           
lon1,lat1,phi=np.loadtxt('../data_in/'+var+'_'+date+'_'+hour+'_'+level+'_T'+trunc+'gg.dat',unpack=True)
#lon1,lat1,z0=np.loadtxt('../data_in/'+var+'_'+date+'_'+hour+'_T'+trunc+'gg.dat',unpack=True)

lat2 = lat1.reshape((nlats,nlons))
lon2 = lon1.reshape((nlats,nlons))

proj = ccrs.Orthographic(central_longitude=-90.0,central_latitude=40.0)
proj = ccrs.PlateCarree()
proj = ccrs.NorthPolarStereo(central_longitude=-90.0)
#proj = ccrs.LambertConformal(central_longitude=-85.0,central_latitude=45.0,cutoff=0.0)
#proj = ccrs.LambertConformal(central_longitude=-65.0,central_latitude=45.0,cutoff=0.0)
ax = plt.axes(projection=proj)
ax.coastlines('110m',linewidth=1.5,color='purple')
#ax.set_extent([-130, -40, 13, 85]) # North America
#ax.set_extent([-120, -20, 20, 85]) # North America
ax.set_extent([-180, 180, 90, 20], ccrs.PlateCarree()) # North hemisphere
gl = ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1,color='black',alpha=0.5)
gl.ylocator = mticker.FixedLocator(np.arange(-90,90,15))
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,15))


if var == 'T':
	field = t.reshape((nlats,nlons))
if var == 'ALT':
	field = alt.reshape((nlats,nlons))
if var == 'Z0':
	field = z0.reshape((nlats,nlons))
if var == 'LSM':
	field = lsm.reshape((nlats,nlons))
if var == 'DIV':
	field = div.reshape((nlats,nlons))
if var == 'W':
	field = w.reshape((nlats,nlons))
if var == 'PHI':
	field = phi.reshape((nlats,nlons))
if var == 'PV':
	field = pv.reshape((nlats,nlons))
if var == 'VOR':
	field = vor.reshape((nlats,nlons))
	#field2 = lat1.reshape((nlats,nlons))
	#field3 = phi.reshape((nlats,nlons))
	#field = (field + 2.*omega*np.sin(np.pi/180.*field2))/field3
if var == 'U':
	field = u.reshape((nlats,nlons))
if var == 'V':
	field = v.reshape((nlats,nlons))
if var == 'uv':
	field1 = u.reshape((nlats,nlons))
	field2 = v.reshape((nlats,nlons))

if var == 'PHI':
	#cs = map.contourf(x,y,field/9.81,levels=np.arange(4900,5200,25),cmap='jet')	
	#cs = map.contour(x,y,field/98.1,levels=np.arange(460,590,10),colors='black')	
        #cs = ax.contour(lon2,lat2,field/9.81,levels=np.arange(10500,12800,100),linewidths=0.1,transform=ccrs.PlateCarree(),colors='black')
        #cs = ax.contourf(lon2,lat2,field/9.81,levels=np.arange(10500,12800,100),transform=ccrs.PlateCarree(),cmap='jet')
	cs = ax.contour(lon2,lat2,field/9.81,linewidths=0.1,levels=np.arange(4600,5900,100),transform=ccrs.PlateCarree(),colors='black')
	cs = ax.contourf(lon2,lat2,field/9.81,levels=np.arange(4600,5900,100),transform=ccrs.PlateCarree(),cmap='jet')
	#cs = map.contour(x,y,(field)/9.81+4450,levels=np.arange(9150,10300,50),colors='black')
	#cs = map.contourf(x,y,(field)/9.81+4450,levels=np.arange(9150,10300,50),cmap='jet')
	#cs = map.contour(x,y,field,10,colors='black')
	#cs = map.contourf(x,y,field,10,cmap='seismic')
else:
	#z = np.sqrt(field1*field1 + field2*field2)
	#cs = ax.contourf(lon2,lat2,field,10,transform=ccrs.PlateCarree(),cmap='jet') 
	cs = ax.contour(lon2,lat2,field,levels=np.arange(210,265,5),linewidths=0.1,transform=ccrs.PlateCarree(),colors='black')
	cs = ax.contourf(lon2,lat2,field,levels=np.arange(210,265,5),transform=ccrs.PlateCarree(),cmap='jet')#'RdPu')
	#cs = ax.contour(lon2,lat2,field*1E6,levels=np.arange(-200,250,50),linewidths=0.1,transform=ccrs.PlateCarree(),colors='black')
	#cs = ax.contourf(lon2,lat2,field*1E6,levels=np.arange(-200,250,50),transform=ccrs.PlateCarree(),cmap='seismic')#'RdPu')
	#map.quiver(x[::1,::1],y[::1,::1],field1[::1,::1],field2[::1,::1],z[::1,::1],width=0.002,headwidth=2,scale=0.1,scale_units='xy',cmap='jet')
#plt.colorbar(shrink=0.5)
cbar = plt.colorbar(cs,orientation='vertical',shrink=0.5,pad=0.07)
cbar.ax.tick_params()
plt.title('SQG '+var+' @ '+level+' hpa '+hour+' h - T'+trunc+' ERA5')
plt.savefig('../plots/ERA5_'+var+'_'+date+'_'+hour+'_'+level+'_T'+trunc+'_IC.pdf')
#plt.savefig('../plots/ERA5_'+var+'_'+date+'_'+hour+'_T'+trunc+'_IC.pdf')
plt.show()
