#Question 9: Convolution of box function with itself.
import numpy as np
import matplotlib.pyplot as plt

def box_f(x):
    return np.where(abs(x)<=1,1,0)

N=256
x_min=-5
x_max=5
x=np.linspace(x_min,x_max,N)
Dx=x[1]-x[0]

k=2*np.pi* np.fft.fftfreq(N,d=Dx)
odr=k.argsort()[::1]

#fx=np.exp(-x**2)
fx=box_f(x)

fk=np.fft.fft(fx,norm="ortho")
Fk=fk*fk

Fx=Dx*np.sqrt(N)*np.fft.ifft(Fk,norm='ortho')
Fx=Fx[odr]

plt.plot(x,fx,'b--')
plt.plot(x,np.real(Fx),'r')
plt.legend(['$f(x)$',r'$[f\otimes f](x)$'],fontsize=14)
plt.xlabel('$x$')
plt.show()
