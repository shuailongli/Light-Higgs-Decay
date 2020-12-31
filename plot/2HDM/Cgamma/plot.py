import matplotlib.pyplot as plt
from pylab import *
import numpy as np
from matplotlib import patches

data=np.loadtxt("Cgamma_w_usd.txt")
mA=data[:,0]/1000.
OA1=data[:,2]*(-100000)
OA2=data[:,4]*(-100000)
OA3=data[:,6]*(-100000)

data=np.loadtxt("Cgamma_wo_usd.txt")
OA4=data[:,2]*(-100000)
OA5=data[:,4]*(-100000)
OA6=data[:,6]*(-100000)


plt.plot([],[],ls='-',c='k',label=r'w/ $u,d,s$')
plt.plot([],[],ls='--',c='k',label=r'w/o $u,d,s$')

plt.plot(mA,OA1,ls='-',c='r',label=r'$\tan\beta=0.1$')
plt.plot(mA,OA2,ls='-',c='g',label=r'$\tan\beta=1$')
plt.plot(mA,OA3,ls='-',c='b',label=r'$\tan\beta=10$')
plt.plot(mA,OA4,ls='--',c='r')
plt.plot(mA,OA5,ls='--',c='g')
plt.plot(mA,OA6,ls='--',c='b')

plt.legend(loc=1,prop={'size': 14},framealpha=0.0)

plt.xlim(1e-3,10.0)
plt.ylim(0,10)
plt.xscale('log')
plt.xlabel(r'$m_A$ (GeV)',fontsize=16)
plt.ylabel(r'$-C_\gamma\times 10^5$ (1/MeV)',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.yticks([0.1,0.2,0.5,1,2,5,10,20,50],[0.1,0.2,0.5,1,2,5,10,20,50],fontsize=16)
plt.title(r'Im$(C_\gamma)$ in Type-II',fontsize=16)
plt.axes().yaxis.set_ticks_position("both")
plt.axes().yaxis.set_label_position("left")
#plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
plt.axes().yaxis.set_minor_locator(MultipleLocator(1))
plt.tight_layout()
#plt.show()
savefig('Cgamma_imag.pdf')
