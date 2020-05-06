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
k1 = 2*np.pi* np.fft.fftfreq(N, d=Dx)
FTfx=Dx*np.exp(-1j*k1*x_min)*FTfx*np.sqrt(N/2/np.pi)

plt.subplot(1,2,1)
plt.plot(x,fx,'b')
plt.xlabel(r'$x$',fontsize=16)
plt.ylabel(r'$f(x)$',fontsize=16)

plt.grid(True)

##################################################################

f_Q2=open("Q_2.txt","r")           #Read the data obtained from Q2.c
Data=f_Q2.readlines()
i=0
K2=np.zeros(N) 
fK2=np.zeros(N)
for D in Data:
    D1,D2=D.split()
    K2[i]=float(D1)
    fK2[i]=float(D2)
    i=i+1

odr = k1.argsort()[::1]
k1 = k1[odr]
FTfx = FTfx[odr]

plt.subplot(1,2,2)
plt.plot(k1,np.real(FTfx),'b')
plt.plot(K2,rect(K2),'r')
plt.plot(K2,fK2,'g')
plt.xlabel(r'$k$',fontsize=16)
plt.ylabel(r'$\tilde{f}(k)$',fontsize=16)

plt.grid(True)

##################################################################
f_Q3=open("Q_3.txt","r")           #Read the data obtained from Q2.c
data=f_Q3.readlines()
i=0
K3=np.zeros(N) 
fK3=np.zeros(N)
for D in data:
    D1,D2=D.split()
    K3[i]=float(D1)
    fK3[i]=float(D2)
    i=i+1

plt.plot(K3,fK3,'pink')
plt.xlabel(r'$k$',fontsize=16)
plt.ylabel(r'$\tilde{f}(k)$',fontsize=16)

plt.grid(True)
plt.legend(['Numpy.fft','Analytical','FFT from C','GSL'],fontsize=16)
plt.show()
