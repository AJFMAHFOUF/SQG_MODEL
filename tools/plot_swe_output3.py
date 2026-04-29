from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt 
import numpy as np
nlats=48 ; nlons=97
nstep='0120'
trunc='031'
expid='002'
var='qv'
var='uv'
plt.figure(figsize=(12,6))
map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')
#map = Basemap(projection='ortho',lat_0=40,lon_0=-90,resolution='l')  
#map = Basemap(projection='stere',width=18000000,height=18000000,lat_0=90.0,lon_0=-90.0,lat_ts=60.0,resolution='l')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))            
lon1,lat1,vor,div,phi,u,v,psi,khi,qv=np.loadtxt('../data_out/SSWE_model_T'+trunc+'_step_'+nstep+'_expid_'+expid+'.dat',unpack=True)
if var == 'psi':
	field = psi.reshape((nlats,nlons))
if var == 'phi':
	field = phi.reshape((nlats,nlons))
if var == 'khi':
	field = khi.reshape((nlats,nlons))
if var == 'vor':
	field = vor.reshape((nlats,nlons))
if var == 'div':
	field = div.reshape((nlats,nlons))
if var == 'u':
	field = u.reshape((nlats,nlons))
if var == 'v':
	field = v.reshape((nlats,nlons))
if var == 'uv':
	field1 = u.reshape((nlats,nlons))
	field2 = v.reshape((nlats,nlons))
if var == 'qv':
	field = qv.reshape((nlats,nlons))
lon2 = lon1.reshape((nlats,nlons))
lat2 = lat1.reshape((nlats,nlons))
x,y = map(lon2,lat2)
if var == 'phi':
	#cs = map.contourf(x,y,field/9.81,levels=np.arange(4600,6200,100),cmap='jet')	
	cs = map.contour(x,y,field/98.1,levels=np.arange(460,590,10),colors='black')
	cs = map.contourf(x,y,field/98.1,levels=np.arange(460,590,10),cmap='jet')
	#cs = map.contour(x,y,(field)/9.81+4450,levels=np.arange(9150,10300,50),colors='black')
	#cs = map.contourf(x,y,(field)/9.81+4450,levels=np.arange(9150,10300,50),cmap='jet')
	#cs = map.contourf(x,y,field/9.81,15,cmap='jet')
else:
	z = np.sqrt(field1*field1 + field2*field2)
	#cs = map.contourf(x,y,field,10,cmap='jet')
	map.quiver(x[::1,::1],y[::1,::1],field1[::1,::1],field2[::1,::1],z[::1,::1],width=0.002,headwidth=2,scale=2,scale_units='xy',cmap='jet')
#plt.colorbar(shrink=0.5)
#plt.title(var+' 500 hpa '+nstep+'h - T'+trunc+' expid='+expid)
#plt.savefig('../plots/'+var+'_500_T'+trunc+'_step_'+nstep+'_expid_'+expid+'_NH.pdf')
plt.show()
