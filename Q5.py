#Question 5
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

x_min=-5
x_max=5
i=0
T1=np.zeros(252)
T2=np.zeros(252)

for N in range(4,256):
    
    x=np.linspace(x_min,x_max,N)
    Dx=x[1]-x[0]
    fx=1/(x**2+1)   #This only an example. I have generated the using a Lorentzian.
##################################################################
# Direct computation of Fourier transform
    st1=timer()
    E=np.zeros((N,N))
    for p in range(N):
        for r in range(N):
            q=-N/2+r
            E[r][p]=np.exp(-2*np.pi*1j*p*q/N)

    np.dot(E,fx)       
    en1=timer()
    T1[i]=en1-st1
##################################################################
# Fourier tranform using fft.
    st2=timer()
    np.fft.fft(fx,norm='ortho')
    en2=timer()
    T2[i]=en2-st2
    i=i+1

plt.loglog(range(4,256),T1,'b.')
plt.loglog(range(4,256),T2,'r.')
plt.legend(['Direct computation','Using numpy.fft'])
plt.show()
