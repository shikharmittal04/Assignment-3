#Question 1: Fourier transform of f(x)=sinc(x)
#WARNING: Run this code only after running the Q2 and Q3 C codes!

import numpy as np
import matplotlib.pyplot as plt

def rect(x):
    return np.where(abs(x)<=1,np.sqrt(np.pi/2),0)

N=128
x_min=-10
x_max=10
x=np.linspace(x_min,x_max,N)
Dx=x[1]-x[0]
fx=np.sinc(x/np.pi)

FTfx=np.fft.fft(fx,norm='ortho')
k = 2*np.pi* np.fft.fftfreq(N, d=Dx)
FTfx=Dx*np.exp(-1j*k*x_min)*FTfx*np.sqrt(N/2/np.pi)

plt.subplot(1,2,1)
plt.plot(x,fx,'b')
plt.xlabel(r'$x$',fontsize=16)
plt.ylabel(r'$f(x)$',fontsize=16)

plt.grid(True)


##################################################################

f=open("Q_2.txt","r")           #Read the data obtained from C.
Data=f.readlines()
i=0
K=np.zeros(N) 
fK=np.zeros(N)
for D in Data:
    D1,D2=D.split()
    K[i]=float(D1)
    fK[i]=float(D2)
    i=i+1

plt.subplot(1,2,2)
odr = k.argsort()[::1]
k = k[odr]
FTfx = FTfx[odr]
plt.plot(k,np.real(FTfx),'b')
plt.plot(K,rect(K),'r')
plt.plot(K,fK,'g')
plt.xlabel(r'$k$',fontsize=16)
plt.ylabel(r'$\tilde{f}(k)$',fontsize=16)

plt.grid(True)
plt.legend(['Numpy.fft','Analytical','FFT from C'],fontsize=16)
plt.show()
