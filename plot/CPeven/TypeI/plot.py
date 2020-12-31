import matplotlib.pyplot as plt
from pylab import *
import numpy as np

rc('font',**{'family':'serif','serif':['Helvetica']})
rc('text', usetex=True)

data=np.loadtxt("CPeven_TypeI.txt")
mC=data[:,0]
tanb=data[:,1]
length=data[:,2]
gammalength=data[:,3]

mC=np.reshape(mC,(50,50))
tanb=np.reshape(tanb,(50,50))
length=np.reshape(length,(50,50))
gammalength=np.reshape(gammalength,(50,50))

contours1=plt.contour(mC,tanb,length,[0.001,0.01,0.1,10,1000],colors='r')
contours2=plt.contourf(mC,tanb,gammalength,[102,13081.7],colors='0.8')
#plt.clabel(contours1, inline=True, fontsize=16,fmt='%1.3f')

fmt = {}
strs = [r'$10^{-3}$',r'$10^{-2}$',r'0.1',r'$10$',r'$10^3$']
for l, s in zip(contours1.levels, strs):
    fmt[l] = s
manual_locations = [(0.05, 50), (0.2,50), (0.6, 50), (1, 20), (3,20)]
labels1=plt.clabel(contours1, contours1.levels, inline=True,inline_spacing=0, fmt=fmt,manual=manual_locations, fontsize=16)

for l in labels1:
    l.set_rotation(0)

text(1, 1.5, r'Decay Length [m]', color='r', fontsize=16)

plt.xlim(0.01,10)
plt.ylim(1,1000)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$m_{H}$ (GeV)",fontsize=18)
plt.ylabel(r"$\tan\beta$",fontsize=18)
plt.xticks(fontsize=16)
plt.title(r'Type-I: $\cos(\beta-\alpha)=1/\tan\beta$',fontsize=18)
#plt.yticks((1,2,3,4,6,10,20,30,50),(1,2,3,4,6,10,20,30,50),fontsize=16)
plt.yticks(fontsize=16)
#plt.axes().xaxis.set_minor_locator(MultipleLocator(50))

plt.axes().tick_params(direction='in', which='both',labelsize=18)

plt.tight_layout()
#plt.show()
savefig("DecayLength_CPeven.pdf")

