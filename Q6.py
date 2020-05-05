#Question 6: Fourier transform of a constant function.
#Analytically it should be a Dirac delta function.
import numpy as np
import matplotlib.pyplot as plt

N=64
x_min=-5
x_max=5
x=np.linspace(x_min,x_max,N)
Dx=x[1]-x[0]
fx=np.ones(N)

FTfx=np.fft.fft(fx,norm='ortho')
k =2*np.pi* np.fft.fftfreq(N, d=Dx)
FTfx=Dx*np.exp(-1j*k*x_min)*FTfx*np.sqrt(N/2/np.pi)

plt.subplot(1,2,1)
plt.plot(x,fx,'b')
plt.xlabel(r'$x$',fontsize=16)
plt.ylabel(r'$f(x)$',fontsize=16)

plt.grid(True)

plt.subplot(1,2,2)
plt.plot(k,FTfx,'b')
plt.xlabel(r'$k$',fontsize=16)
plt.ylabel(r'$\tilde{f}(k)$',fontsize=16)

plt.grid(True)
plt.show()

