import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve
g = Gaussian1DKernel(stddev=20)

######### Josh Lothringer models
#d= np.genfromtxt("wasp33b_from_Josh/wasp33_self_consistent.csv", delimiter = ',') #Angstroms, log10(erg/s/cm2/cm)
d= np.genfromtxt("wasp33b_from_Josh/WASP33b.h2o=1.0e-04.co=1.0e-03.INV3.7.spec.txt", delimiter = ',') #Angstroms, log10(erg/s/cm2/cm)
#0.0004 at 1 micron
d[:,0] /= 1.e4        #convert wavelength to microns
d[:,1] = 10.**(d[:,1])  #convert to ergs/s/cm2/cm 
#plt.plot(d[:,0], d[:,1])     
wp = d[:,0]
fp = d[:,1]
    

#star: wavelength (microns); flux at 1 AU, in W/m2/micron ; flux at stellar surface in erg/cm2/sec/Hz 
star = np.genfromtxt("W33b_star_NEXTGEN_B.dat")       #W/m2/micron (column 1)
#0.0002 micron at 1 micron

#plt.plot(star[:,0], star[:,1]*2.998e10)
ws = star[:,0]
fs = star[:,1]*2.998e10

fp = np.interp(ws, wp, fp)

fpfs = fp/fs*(0.103)**2*1e6/np.pi/ws

plt.plot(ws, convolve(fpfs, g, boundary = 'extend'))

plt.xlim(0.7, 1.7)
plt.ylim(0, 2000)
#plt.ylim(0,5)
plt.show() 
