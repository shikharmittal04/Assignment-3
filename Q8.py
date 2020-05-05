#Question 8: 2D Fourier transform.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax1 = fig.add_subplot(121, projection='3d')

N=64
x_min,x_max=-10,10
y_min,y_max=-10,10
x = np.linspace(x_min, x_max, N)
y = np.linspace(y_min, y_max, N)

Dx=(x_max-x_min)/(N-1)
Dy=(y_max-y_min)/(N-1)
##################################################
x,y = np.meshgrid(x, y)

z=np.exp(-x**2-y**2)        #Given function
ax1.plot_surface(x,y,z,cmap=cm.jet,linewidth=1)
ax1.set_xlabel('$x$',fontsize=14)
ax1.set_ylabel('$y$',fontsize=14)
ax1.set_zlabel(r'$f$',fontsize=14)
##################################################

ax2 = fig.add_subplot(122, projection='3d')

Kx = 2*np.pi*np.fft.fftfreq(N,d=Dx)
Ky = 2*np.pi*np.fft.fftfreq(N,d=Dy)

kx,ky = np.meshgrid(Kx, Ky)
F=Dy*Dx*np.exp(-1j*kx*x_min)*np.exp(-1j*ky*y_min)*(N/2/np.pi)*np.fft.fft2(z,norm="ortho")

odr = Kx.argsort()[::1]
Kx = Kx[odr]
Ky = Ky[odr]
Kx,Ky = np.meshgrid(Kx, Ky)

F=F[odr]
F=F[:,odr]

ax2.plot_surface(Kx,Ky,np.real(F),cmap=cm.hot,linewidth=1)
ax2.set_xlabel('$k_x$',fontsize=14)
ax2.set_ylabel('$k_y$',fontsize=14)
ax2.set_zlabel(r'$\tilde{f}$',fontsize=14)
######################################################
Fk=0.5*np.exp(-(Kx**2+Ky**2)/4) #Analytical version 
ax2.plot_surface(Kx,Ky,Fk,cmap=cm.cool,linewidth=1)

plt.tight_layout()
plt.show()
