import matplotlib.pyplot as plt
from pylab import *
import numpy as np
from matplotlib import patches

rc('font',**{'family':'serif','serif':['Helvetica']})
rc('text', usetex=True)

data=np.loadtxt("data_width_tanb10.txt")
mA=data[:,0]/1000.
gaga=data[:,1]/1000.
ee=data[:,2]/1000.
mumu=data[:,3]/1000.
tautau=data[:,4]/1000.
hadron=data[:,5]/1000.
quark=data[:,6]/1000.
gluon=data[:,7]/1000.


plt.plot(mA,gaga,ls='-',c='0.5',label=r'$\gamma\gamma$')
plt.plot(mA,ee,ls='-',c='fuchsia',label=r'$ee$')
plt.plot(mA,mumu,ls='-',c='r',label=r'$\mu\mu$')
plt.plot(mA,tautau,ls='-',c='orangered',label=r'$\tau\tau$')
plt.plot(mA,hadron,ls='-',c='g',label=r'hardon')
plt.plot(mA,quark,ls='-',c='brown',label=r'quarks')
plt.plot(mA,gluon,ls='-',c='orange',label=r'$gg$')

plt.legend(loc=2,prop={'size': 14},framealpha=0.0)

plt.xlim(1e-1,10.0)
plt.ylim(1e-18,1e-6)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m_A$ (GeV)',fontsize=16)
plt.ylabel(r'$\Gamma$ (GeV)',fontsize=16)
plt.axes().xaxis.set_ticks_position("both")
plt.axes().yaxis.set_ticks_position("both")
#plt.yticks([0.1,0.2,0.5,1,2,5,10,20,50],[0.1,0.2,0.5,1,2,5,10,20,50],fontsize=16)
plt.title(r'Type-I: $\tan\beta=10$',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.axes().tick_params(direction='in', which='both',labelsize=16)

#plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
#plt.axes().yaxis.set_minor_locator(MultipleLocator(1))
plt.tight_layout()
#plt.show()
savefig('DecayWidth.pdf')
