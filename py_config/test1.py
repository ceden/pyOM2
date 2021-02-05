import numpy as np
import numpy.fft
import pylab as plt

u_0 = 0.5
z=np.linspace( -1000,0,1000)
dz = z[1]-z[0]
I=complex(0,1)
print dz

u=u_0*z
m=2*np.pi*numpy.fft.fftfreq(u.shape[0],dz)
U=numpy.fft.fft(u)
uz= np.real( numpy.fft.ifft( I*m*U)   )

fig=plt.figure(figsize=(12,8))
fig.clf()
ax=fig.add_subplot(221)
ax.plot(z,uz)
ax=fig.add_subplot(222)
ax.plot(m,np.log10( np.abs(U)**2) )
ax=fig.add_subplot(223)
ax.plot(np.log10( np.abs(U)**2 ))

plt.show()


