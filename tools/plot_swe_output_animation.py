from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.animation as animation
nlats=48 ; nlons=97
trunc='031'
expid='002'
fig,ax=plt.subplots(figsize=(12,6))
liste = []
i = 0
ech = 121
step = 2
while i < ech:
	if i < 10:
		j='000'+str(i)
	elif i < 100:
		j='00'+str(i)
	else:
		j='0'+str(i)
	liste.append(j)
	i += step
field1 = np.zeros((np.size(liste),nlats,nlons))
field2 = np.zeros((np.size(liste),nlats,nlons))
i=0
# Read output files 
for nstep in liste:
	lon1,lat1,vor,div,phi,u,v,psi,khi,qv=np.loadtxt('../data_out/SSWE_model_T'+trunc+'_step_'+nstep+'_expid_'+expid+'.dat',unpack=True)
	field1[i,:,:] = u.reshape((nlats,nlons))
	field2[i,:,:] = v.reshape((nlats,nlons))
	lon2 = lon1.reshape((nlats,nlons))
	lat2 = lat1.reshape((nlats,nlons))
	i += 1
	
print ('Fields for animations have been read') 

def animate(index):
	ax.clear()
	#map = Basemap(projection='ortho',lat_0=40,lon_0=-90,resolution='l')
	map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution='c')
	#map = Basemap(projection='stere',width=18000000,height=18000000,lat_0=90.0,lon_0=90.0,lat_ts=-60.0,resolution='l')
	map.drawcoastlines(linewidth=0.25)
	map.drawcountries(linewidth=0.25)
	map.drawmeridians(np.arange(0,360,30))
	map.drawparallels(np.arange(-90,90,30))  
	x,y = map(lon2,lat2)
	hour = step*index
	modulus = np.sqrt(field1[index,:,:]**2 + field2[index,:,:]**2)
	ax.set_title('Vectot wind after  '+str(hour)+'  hours')
	#map.contour(x,y,field[index,:,:]/98.1,levels=np.arange(460,590,10),colors='black',linewidths=0.5)
	#map.contourf(x,y,field[index,:,:]/98.1,levels=np.arange(460,590,10),cmap='jet')
	#map.contour(x,y,field[index,:,:],levels=np.arange(-50,50,5),colors='black',linewidths=0.5)
	#map.contourf(x,y,field[index,:,:],levels=np.arange(-50,50,5),cmap='seismic')
	map.quiver(x,y,field1[index,:,:],field2[index,:,:],modulus,width=0.002,headwidth=2,scale=2.5,scale_units='xy',cmap='jet')
	
ani = animation.FuncAnimation(fig=fig,func=animate,frames=range(np.size(liste)))
#plt.show()
ani.save('uv.gif',dpi=300,writer=animation.PillowWriter(fps=3))
