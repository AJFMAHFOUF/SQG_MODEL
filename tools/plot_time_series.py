import matplotlib.pyplot as plt
import numpy as np
x,y=np.loadtxt('../data_out/time_series_dfi.dat',unpack=True)
x,y1=np.loadtxt('../data_out/time_series_nodfi.dat',unpack=True)
plt.figure(figsize=(6,6))
plt.xlim(0.0,120)
plt.xlabel('Hours')
plt.ylabel('h(m)')
plt.title('Geopotential 500 hPa - 21/12/1978')
plt.plot(x/4,y1/9.81+4450)
plt.plot(x/4,y/9.81+4450,linewidth=3,color='red')
plt.savefig('../plots/time_series_phi500_dfi.pdf')
