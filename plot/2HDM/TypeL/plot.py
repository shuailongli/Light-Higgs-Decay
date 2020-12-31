import matplotlib.pyplot as plt
from pylab import *
import numpy as np
from matplotlib import patches

rc('font',**{'family':'serif','serif':['Helvetica']})
rc('text', usetex=True)

data=np.loadtxt("data.txt")
mA=data[:,0]/1000.
OA1=data[:,1]/1000.
OA2=data[:,2]/1000.
OA3=data[:,3]/1000.
OA4=data[:,4]/1000.
OA5=data[:,5]/1000.

def length2width(x):
    return 1.9732705*1e-16 / x

def width2length(x):
    return 1.9732705*1e-16 / x


plt.plot(mA,length2width(OA5),ls='-',c='c',label=r'$\tan\beta=10^2$')
plt.plot(mA,length2width(OA4),ls='-',c='m',label=r'$\tan\beta=10$')
plt.plot(mA,length2width(OA3),ls='-',c='b',label=r'$\tan\beta=1$')
plt.plot(mA,length2width(OA2),ls='-',c='g',label=r'$\tan\beta=0.1$')
plt.plot(mA,length2width(OA1),ls='-',c='r',label=r'$\tan\beta=10^{-2}$')

plt.legend(loc=1,prop={'size': 14},framealpha=0.0)

plt.xlim(1e-1,10.0)
#plt.ylim(1e-3,1e3)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m_A$ (GeV)',fontsize=16)
plt.ylabel(r'$c\tau$ (m)',fontsize=16)
plt.axes().xaxis.set_ticks_position("both")
#plt.yticks([0.1,0.2,0.5,1,2,5,10,20,50],[0.1,0.2,0.5,1,2,5,10,20,50],fontsize=16)
plt.title(r'Type-L',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

secax = plt.axes().secondary_yaxis('right', functions=(width2length, length2width))
secax.set_ylabel(r'$\Gamma(A)$ (GeV)',fontsize='16')

plt.axes().tick_params(direction='in', which='both',labelsize=16)
secax.tick_params(direction='in', which='both',labelsize=16)

#plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
#plt.axes().yaxis.set_minor_locator(MultipleLocator(1))
plt.tight_layout()
#plt.show()
savefig('DecayWidth.pdf')
