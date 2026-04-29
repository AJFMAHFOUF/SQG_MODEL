from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
nlats=96 ; nlons=193
nstep='0000'
trunc='063'
expid1='001'
expid2='002'
var='div'
plt.figure(figsize=(12,6))
map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')
#map = Basemap(projection='ortho',lat_0=40,lon_0=-90,resolution='l')           
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))            
lon1,lat1,vor,div,phi,u,v,psi,khi,qv=np.loadtxt('../data_out/SSWE_model_T'+trunc+'_step_'+nstep+'_expid_'+expid1+'.dat',unpack=True)
if var == 'psi':
	field1 = psi.reshape((nlats,nlons))
if var == 'phi':
	field1 = phi.reshape((nlats,nlons))
if var == 'khi':
	field1 = khi.reshape((nlats,nlons))
if var == 'vor':
	field1 = vor.reshape((nlats,nlons))
if var == 'div':
	field1 = div.reshape((nlats,nlons))
if var == 'u':
	field1 = u.reshape((nlats,nlons))
if var == 'v':
	field1 = v.reshape((nlats,nlons))
if var == 'qv':
	field1 = qv.reshape((nlats,nlons))
lon1,lat1,vor,div,phi,u,v,psi,khi,qv=np.loadtxt('../data_out/SSWE_model_T'+trunc+'_step_'+nstep+'_expid_'+expid2+'.dat',unpack=True)
if var == 'psi':
	field2 = psi.reshape((nlats,nlons))
if var == 'phi':
	field2 = phi.reshape((nlats,nlons))
if var == 'khi':
	field2 = khi.reshape((nlats,nlons))
if var == 'vor':
	field2 = vor.reshape((nlats,nlons))
if var == 'div':
	field2 = div.reshape((nlats,nlons))
if var == 'u':
	field2 = u.reshape((nlats,nlons))
if var == 'v':
	field2 = v.reshape((nlats,nlons))
if var == 'qv':
	field2 = qv.reshape((nlats,nlons))		
lon2 = lon1.reshape((nlats,nlons))
lat2 = lat1.reshape((nlats,nlons))
x,y = map(lon2,lat2)
if var == 'phi':
	cs = map.contourf(x,y,(field2-field1)/9.81,15,cmap='seismic')
	#cs = map.contourf(x,y,field/9.81,15,cmap='jet')
else:
	cs = map.contourf(x,y,field2-field1,15,cmap='seismic')
print ('mean value',np.mean(field2-field1))	
plt.colorbar(shrink=0.5)
plt.title(var+' 500 hpa '+nstep+'h - T'+trunc+' expid='+expid2+'-'+expid1)
plt.savefig('../plots/'+var+'_500_T'+trunc+'_step_'+nstep+'_expid_'+expid2+'-'+expid1+'.png')
plt.show()
