import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel
import seaborn as sns
import scipy.stats as st
from astropy.convolution import Gaussian1DKernel, convolve

sns.set_context("talk", font_scale=1.5)
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})

h = 6.626e-34   #J/s
c = 3.0e8       #m/s
kb = 1.4e-23     #J/K


def boxcar_smooth(x, nsmooth):
	n = len(x)
	for i in np.arange(0, n):
		lower = np.max(np.array([0, i - int(nsmooth/2)]))
		upper = np.min(np.array([n-1, i + int(nsmooth/2)]))
		x[i] = 	np.mean(x[lower:upper])
	return x

def bin_at_obs_res(waves, model_waves, model):
    delta_lambda = (waves[1] - waves[0])/2.0
    binned_model = np.zeros_like(waves)
    for i in range(0, len(waves)):
        ind = (model_waves>waves[i]-delta_lambda)&(model_waves<waves[i]+delta_lambda)
        binned_model[i] = np.average(model[ind])
    return binned_model


def get_significance(chisq, dof):
	alpha = (1. - st.chi2.cdf(chisq, dof))/2.
	z = st.norm.ppf(1.-alpha)
	return z


def blackbody(l,T): 
	return 2*h*c**2/(l**5*(np.exp(h*c/(l*kb*T))-1))

def best_fit_bb(w, y, e, rprs):
    Ts = np.linspace(3140, 3170, 50)
    #off1s = np.linspace(-20, 0, 50)/1.e6
    #off2s = np.linspace(130, 150, 50)/1.e6
    off1s = np.linspace(-200, 200, 2)/1.e8
    off2s = np.linspace(-200, 200, 2)/1.e8
    chibest = 10000.
    Tbest = 0.  
    resid = 0.


    #get stellar spectrum
    g = Gaussian1DKernel(stddev=3)
    star = np.genfromtxt("W33b_star_NEXTGEN_B.dat")       #W/m2/micron (column 1)
    #star_bb = np.interp(w, star[:,0], convolve(star[:,1]*22423.,g, boundary = 'extend'))
    star_bb = bin_at_obs_res(w, star[:,0], star[:,1])*7.316156e9
    print "arbitrary rescaling of stellar spectrum!!!!"
    model = blackbody(w*1.0e-6, 7400.)

    chis = []

    for T in Ts:
        for off1 in off1s:
            for off2 in off2s:
                model = blackbody(w*1.0e-6, T)/star_bb*rprs**2
                ytemp = np.copy(y)
                ytemp[0:23] -= off1  
                ytemp[23::] -= off2
                chi2 = np.sum((ytemp - model)**2/e**2)
                chis.append(chi2)
                if chi2 < chibest: 
                    chibest, Tbest, off1best, off2best, resid = chi2, T, off1, off2, (y - model)/e
                    print  chibest, Tbest, off1best, off2best 

    waves_hires = np.linspace(0.7, 5.0, 100)
    star_bb_hires = bin_at_obs_res(waves_hires, star[:,0], star[:,1])*7.316156e9

    #return waves_hires, blackbody(waves_hires*1.0e-6, Tbest)/star_bb_hires*rprs**2, Tbest, off1best, off2best, get_significance(chibest,(len(w)-3))
    return waves_hires, blackbody(waves_hires*1.0e-6, Tbest)/star_bb_hires*rprs**2, Tbest, off1best, off2best, get_significance(chibest,(len(w)-1))


g = Gaussian1DKernel(stddev=1.1)

files  = ["wasp33b_from_Madhu/spec.dat", "wasp33b_from_Madhu/spec_notio.dat"]
labels  = ["best fit model", "best fit with no TiO"]
colors = ['blue', 'red']
linestyles = ['dotted', 'dashed']

plt.figure(figsize = (12,6))

"""for i, f in enumerate(files):
	d = np.genfromtxt(f, skip_header = 2)
	plt.plot(d[:,1]*1.e4, convolve(d[:,4], g, boundary = 'extend'), label = labels[i], color = '0.5', alpha = 0.5, linestyle = linestyles[i])"""

s = np.genfromtxt("w33_data_haynes.txt")
#plt.errorbar(s[:,0], s[:,1]/100., s[:,2]/100., fmt = 'ow', markersize = 4, ecolor = 'k', markeredgecolor = 'k', markeredgewidth = 1., linewidth = 1., linestyle='none', zorder=100, label="G141 data (Haynes)")

s = np.genfromtxt("w33_g141_espec_091018.txt")
d = np.genfromtxt("w33_g102_espec_083118.txt")


x = np.append(d[:,0], s[:,0])
y = np.append(d[:,1], s[:,1])
err = np.append(d[:,2], s[:,2])


rprs = 0.103

xm, ym, Tbest, off1best, off2best, chi2  = best_fit_bb(x, y, err, rprs)
print "best fit T, chi2", Tbest, chi2
#plt.plot(xm, ym, color='0.5',  label = 'blackbody fit', alpha = 0.5, linestyle='dotted', zorder = 0.5)

off1best = 0.
off2best = 0.

#plt.errorbar(d[:,0], d[:,1] - off1best, d[:,2], fmt = '.b', zorder=100, label = "G102 data")
#plt.errorbar(s[:,0], s[:,1] - off2best, s[:,2], marker='.', color='r', linestyle='none', zorder=100, label="G141 data (Kreidberg)")

plt.errorbar(d[:,0], d[:,1] - off1best, d[:,2], fmt = '.k', zorder=100)#, label = "G102 data")
plt.errorbar(s[:,0], s[:,1] - off2best, s[:,2], marker='.', color='k', linestyle='none', zorder=100, label="data")


d = np.genfromtxt("SE-W33b-TiO-NoDrag-AllK-RpoRs-0.1055.dat", skip_header = 1, delimiter = ',')
#plt.plot(d[:,1], d[:,2]*0.95316816783, label = 'GCM', color = '0.1')       #multiplied by correction factor for rp/rs (Vivien assuemd 0.1055)

#plt.ylim(0.0, 1.8e-3)
#plt.xlim(0.75, 1.7)
plt.ylim(0.0, 4.8e-3)
plt.xlim(0.75, 5.7)


#plt.gca().annotate('TiO features', xy=(1.05, 0.0008), xytext=(0.8, 0.0005), arrowprops=dict(facecolor='black', shrink=0.05),)
#plt.gca().annotate('', xy=(0.97, 0.00095), xytext=(0.95, 0.0006), arrowprops=dict(facecolor='black', shrink=0.05),)


g = Gaussian1DKernel(stddev=50)
offset = 0.0000

d = np.genfromtxt("wasp33b_from_Caroline/wasp33b_rfacv1.0_m0.0_co1.0nc.flx")
plt.plot(d[:,0], convolve(d[:,1], g, boundary = 'extend') + offset, label = "forward model, rfacv1.0_m0.0")

d = np.genfromtxt("wasp33b_from_Caroline/wasp33b_rfacv0.85_m0.0_co1.0nc.flx")
plt.plot(d[:,0], convolve(d[:,1], g, boundary = 'extend') + offset, label = "forward model, rfacv0.85_m0.0")

#d = np.genfromtxt("wasp33b_m0.0_co1.0nc.spec")
#plt.plot(d[:,0], convolve(d[:,1], g, boundary = 'extend') + offset, label = "forward model, with TiO")

#d = np.genfromtxt("wasp33b_m0.0_co1.0nc_noTiO.spec")
#plt.plot(d[:,0], convolve(d[:,1], g, boundary = 'extend') + offset,  label = "forward model, no TiO")


plt.errorbar(4.5, 4250e-6, 160e-6, fmt = 'xk')
plt.errorbar(3.6, 3506e-6, 173e-6, fmt = 'xk')

plt.tight_layout()
plt.xlabel("Wavelength (microns)")
plt.ylabel("Planet-to-star flux")
plt.legend(loc = "upper left", frameon=True, fontsize=14)

plt.savefig("w33_models.pdf")
#plt.show() 
