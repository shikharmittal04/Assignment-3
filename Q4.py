import numpy as np
import matplotlib.pyplot as plt
f=open("Q_4.txt","r")           #Read the data obtained from C.
Data=f.readlines()
i=0
N=64                            #This N should be the same as used in the corresponding C code.
k=np.zeros(N) 
fk=np.zeros(N)
for D in Data:
    D1,D2=D.split()
    k[i]=float(D1)
    fk[i]=float(D2)
    i=i+1

x=np.linspace(-5,5,N)           #The x_min and x_max should be the same as in the corresponding C code.
fx=np.exp(-x**2)                #The function whose Fourier transform is required.
plt.subplot(1,2,1)
plt.plot(x,fx,'b')
plt.xlabel(r'$x$',fontsize=16)
plt.ylabel(r'$f(x)$',fontsize=16)

plt.grid(True)

plt.subplot(1,2,2)
plt.plot(k,fk,'b')                              #From the C code using FFT.
plt.plot(k,np.sqrt(0.5)*np.exp(-k**2/4),'r')    #Analytical result.
plt.xlabel(r'$k$',fontsize=16)
plt.ylabel(r'$\tilde{f}(k)$',fontsize=16)
plt.legend(['From C (using FFT)','Analytical'],fontsize=16)

plt.grid(True)
plt.show()
