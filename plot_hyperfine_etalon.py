#!/usr/bin/ipython3

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.stats import norm, cauchy
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import scipy.constants as const
from scipy.misc import derivative

args = sys.argv
argnum = 6
etalon_len = 0.2
etalon_freq = 1e-6*const.c/(2*etalon_len)
print(etalon_freq)
usage = '''Usage:
.\{0} probe reference cull_low cull_high etalon
or in ipython console:
%run {0} probe reference cull_low cull_high etalon'''.format(args[0])

def get_chi_squared(func,params,x,y):
    assert(x.shape == y.shape)
    diff = np.abs(x[0]-x)
    diff_nonzero = diff[np.nonzero(diff)]
    error = 0.25*diff_nonzero[diff_nonzero.argmin()]/2
    chi_squared = np.sum((y - func(x,*params))**2)/error**2
    ndof = x.size - params.size
    return chi_squared/ndof

if len(args) < argnum:
    print(usage)
    sys.exit()

etalon = np.loadtxt(args[5],dtype=[('t',np.float64,()),('I',np.float64,())],skiprows=1,delimiter=',')
etalon_mod = etalon
peaks, props = find_peaks(etalon_mod['I'],np.max(etalon_mod['I'])*0.75,distance = 500)
#diff = [etalon_mod['t'][peaks[i]] - etalon_mod['t'][peaks[i-1]] for i in range(1,peaks.size)]
t_loc = etalon_mod['t'][peaks[:]]
nu_loc = [i*etalon_freq for i in range(t_loc.size)]

t_err = 5.5e-6 # Etalon top level half width
map_t_nu = interp1d(t_loc,nu_loc,fill_value='extrapolate')
transform = map_t_nu
inverse = interp1d(nu_loc,t_loc,fill_value='extrapolate')
trans_error = lambda x : derivative(map_t_nu,inverse(x),dx=1e-6)*t_err

fitfunc = lambda x, a, b, A, m, s: a*x + b + A*lorentzian(x,m,s)
params = ['a','b','A','m','s']

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

gaussian = lambda x,m,s : norm.pdf(x,m,s)
lorentzian = lambda x, m, s : cauchy.pdf(x,m,s)

probe = np.loadtxt(args[1],dtype=[('t',np.float64,()),('I',np.float64,())],skiprows=1,delimiter=',')
ref = np.loadtxt(args[2],dtype=[('t',np.float64,()),('I',np.float64,())],skiprows=1,delimiter=',')

#probe['t'] *= scaling
#ref['t'] *= scaling
cull_low = int(args[3])
cull_high = -int(args[4])


plt.plot(t_loc,nu_loc,'x-')
plt.plot(ref['t'][cull_low:cull_high],transform(ref['t'][cull_low:cull_high]))
plt.xlabel("t [s]")
plt.ylabel("$\\nu$ [MHz]")
plt.show()


probe['t'] = transform(probe['t'])

modref = ref['I'][cull_low:cull_high]
modref_t = ref['t'][cull_low:cull_high]
modref_t = transform(modref_t)

reffunc = interp1d(modref_t[4:-5],moving_average(modref,10),kind='linear',fill_value=(modref[0],modref[-1]),bounds_error=False)

#assert(probe['t'].all()==ref['t'].all())
#probe = probe[250:int(len(probe)/3)]
#ref = ref[:int(-250+len(ref)/3)]

#modprobe_t = probe['t'][10:-10]
#modprobe = moving_average(probe['I'],21)
#modref = moving_average(ref['I'],21)

modprobe_t = probe['t'][cull_low:cull_high]
modprobe = probe['I'][cull_low:cull_high]


bg = lambda x, t0, a, b : a*reffunc(x + t0) - b

#def fitfunc(x,B,M1,M2,M3,M4,A1,A2,A3,A4,S,T0):
#    return(B-A1*gaussian(x,T0-M1,S)-A1*gaussian(x,T0-M2,S)-A1*gaussian(x,T0-M3,S)-A1*gaussian(x,T0-M4,S)-A1*gaussian(x,T0+M1,S)-A1*gaussian(x,T0+M2,S)-A1*gaussian(x,T0+M3,S)-A1*gaussian(x,T0+M4,S)).flatten()

#fit, cov = curve_fit(fitfunc,data['t'],data['I'],p0=[48,0.0034,0.0088,0.0217,0.032,1,2,1,0.25,0.001,-0.0034])

bgopt,bgcov = curve_fit(bg,modprobe_t,modprobe,p0=[-70,1,1],maxfev=10000)

etalon_mod['t'] = transform(etalon_mod['t'])

removed_bg = modprobe-bg(modprobe_t,*bgopt)
bgref = bg(modprobe_t,*bgopt)

plt.ion()
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(modprobe_t,modprobe)
ax1.plot(modprobe_t,bg(modprobe_t,*bgopt))
#ax1.plot(etalon_mod['t'],etalon_mod['I'])
#ax1.plot(etalon_mod['t'][peaks],etalon_mod['I'][peaks],'x')
ax1.axvspan(etalon_mod['t'][peaks[0]],etalon_mod['t'][peaks[-1]],alpha=0.25,hatch='/',color='green')
ax2.axvspan(etalon_mod['t'][peaks[0]],etalon_mod['t'][peaks[-1]],alpha=0.25,hatch='/',color='green')
ax2.plot(modprobe_t,removed_bg)
plt.xlabel("Energy [MHz]")
ax1.set_ylabel("Intensity [a.u.]")
ax2.set_ylabel("Intensity [a.u.]")
f.tight_layout()
plt.show()

def fit_line(xmin,xmax):
    x = []
    y = []
    for i in range(modprobe_t.size):
        if(modprobe_t[i] > xmin and modprobe_t[i] < xmax):
            x.append(modprobe_t[i])
            y.append(removed_bg[i])
        elif(modprobe_t[i] >= xmax):
            break
    x = np.array(x)
    y = np.array(y)
    fit, cov = curve_fit(fitfunc,x,y,p0=[0,0,1,(xmin+xmax)/2,10],maxfev = 10000)
    ax2.plot(x,fitfunc(x,*fit))
    for param in zip(params,fit,sp.sqrt(np.diag(cov))):
        print(param[0], ':', param[1], '+-', param[2])
    print("Chi-square/Ndof :", get_chi_squared(fitfunc,fit,x,y))
    return (fit,cov)

def fit_splitting(amin,amax,bmin,bmax):
    print("Params for line A:")
    afit, acov = fit_line(amin, amax)
    print("Params for line B:")
    bfit, bcov = fit_line(bmin, bmax)
    print('')
    print('Splitting:', np.abs(afit[3] - bfit[3]), '+-', np.sqrt(acov[3,3]**2 + bcov[3,3]**2))

def structure(x, bg, a1, a2, a3, a12, a23, a13, m1, m2, m3, s1, s2, s3, s12, s23, s13):
    m12 = (m1+m2)/2
    m23 = (m2+m3)/2
    m13 = (m1+m3)/2
    peak1 = a1*lorentzian(x,m1,s1)
    peak2 = a2*lorentzian(x,m2,s2)
    peak3 = a3*lorentzian(x,m3,s3)
    peak12 = a12*lorentzian(x,m12,s12)
    peak23 = a23*lorentzian(x,m23,s23)
    peak13 = a13*lorentzian(x,m13,s13)
    return bg + peak1 + peak2 + peak3 + peak12 + peak23 + peak13
struct_params = ['bg', 'a1', 'a2', 'a3', 'a12', 'a23', 'a13', 'm1', 'm2', 'm3', 's1', 's2', 's3', 's12', 's23', 's13']

def fit_structure(xmin,xmax,m1=0,m2=0,m3=0,p0=None,plot=True):
    x = []
    y = []
    for i in range(modprobe_t.size):
        if(modprobe_t[i] > xmin and modprobe_t[i] < xmax):
            x.append(modprobe_t[i])
            y.append(removed_bg[i])
        elif(modprobe_t[i] >= xmax):
            break
    x = np.array(x)
    y = np.array(y)
    pos = [0,np.inf]
    scaling = 250000
    HWHM = [0,0.0001*scaling]
    bounds = np.array([[-np.inf,sp.inf],[-np.inf,np.inf],pos,pos,pos,pos,pos,[m1-0.0001*scaling,m1+0.0001*scaling],[m2-0.0001*scaling,m2+0.0001*scaling],[m3-0.0001*scaling,m3+0.0001*scaling],HWHM,HWHM,HWHM,HWHM,HWHM,HWHM]).transpose()
    if(p0 is None):
        fit, cov = curve_fit(structure,x,y,p0=[-1,100,100,100,100,100,100,m1,m2,m3,0.00005*scaling,0.00005*scaling,0.00005*scaling,0.00005*scaling,0.00005*scaling,0.00005*scaling],bounds = bounds,maxfev = 10000)
    else:
        fit, cov = curve_fit(structure,x,y,p0=p0,bounds = bounds,maxfev = 10000)
    ax2.plot(x,structure(x,*fit))
    for i in range(7,10):
        cov[i,i] += trans_error(fit[i])**2
    for param in zip(struct_params,fit,sp.sqrt(np.diag(cov))):
        print(param[0], ':', param[1], '+-', param[2])
    print("Splittings")
    print('1-2 :', np.abs(fit[7] - fit[8]), '+-', np.sqrt(cov[7,7]**2 + cov[8,8]**2))
    print('2-3 :', np.abs(fit[9] - fit[8]), '+-', np.sqrt(cov[9,9]**2 + cov[8,8]**2))
    print('1-3 :', np.abs(fit[7] - fit[9]), '+-', np.sqrt(cov[7,7]**2 + cov[9,9]**2))
    print("Chi-square/Ndof :", get_chi_squared(structure,fit,x,y))
    if(plot):
        f2, (a1,a2,a3) = plt.subplots(3, 1, sharex='col')
        a3.set_xlabel("$\\nu$ [MHz]")
        a1.set_ylabel("Intensity [a.u.]")
        a2.set_ylabel("Intensity [a.u.]")
        a3.set_ylabel("Intensity [a.u.]")
        lim = np.where(np.logical_and(modprobe_t >= xmin, xmax >= modprobe_t))
        a1.plot(modprobe_t[lim],modprobe[lim])
        a1.plot(modprobe_t[lim],bg(modprobe_t,*bgopt)[lim])
        a2.plot(modprobe_t[lim],removed_bg[lim])
        a2.plot(modprobe_t[lim],structure(modprobe_t,*fit)[lim],'r')
        res = removed_bg-structure(modprobe_t,*fit)
        a3.plot(modprobe_t[lim],res[lim])
        f2.tight_layout()
    return (fit, cov)

def fit_total(amin,amax,bmin,bmax,am1,am2,am3,bm1,bm2,bm3):
    assert(am1<am2)
    assert(am2<am3)
    assert(bm1<bm2)
    assert(bm2<bm3)
    assert(am3<bm1)
    print("Fit A parameters:")
    afit, acov = fit_structure(amin,amax,am1,am2,am3,plot=False)
    print("Fit B parameters:")
    bfit, bcov = fit_structure(bmin,bmax,bm1,bm2,bm3,plot=False)
    print("\nS1/2 Splitting:")
    print("a2-b3:", np.abs(afit[8] - bfit[9]), '+-', np.sqrt(acov[8,8]**2 + bcov[9,9]**2))
    print("a1-b2:", np.abs(afit[7] - bfit[8]), '+-', np.sqrt(acov[7,7]**2 + bcov[8,8]**2))
    alim = np.where(np.logical_and(modprobe_t >= amin, amax >= modprobe_t))
    blim = np.where(np.logical_and(modprobe_t >= bmin, bmax >= modprobe_t))
    f2, ((a11, a12), (a21, a22), (a31, a32)) = plt.subplots(3, 2, sharex='col')
    a31.set_xlabel("$\\nu$ [MHz]")
    a32.set_xlabel("$\\nu$ [MHz]")
    a11.set_ylabel("Intensity [a.u.]")
    a12.set_ylabel("Intensity [a.u.]")
    a21.set_ylabel("Intensity [a.u.]")
    a22.set_ylabel("Intensity [a.u.]")
    a31.set_ylabel("Intensity [a.u.]")
    a32.set_ylabel("Intensity [a.u.]")
    a11.set_xlim(amin,amax)
    a12.set_xlim(bmin,bmax)
    f2.show()
    a11.plot(modprobe_t[alim],modprobe[alim])
    a11.plot(modprobe_t[alim],bg(modprobe_t,*bgopt)[alim])
    a12.plot(modprobe_t[blim],modprobe[blim])
    a12.plot(modprobe_t[blim],bg(modprobe_t,*bgopt)[blim])
    a21.plot(modprobe_t[alim],removed_bg[alim])
    a21.plot(modprobe_t[alim],structure(modprobe_t,*afit)[alim],'r')
    a22.plot(modprobe_t[blim],removed_bg[blim])
    a22.plot(modprobe_t[blim],structure(modprobe_t,*bfit)[blim],'r')
    res1 = removed_bg-structure(modprobe_t,*afit)
    a31.plot(modprobe_t[alim],res1[alim])
    res2 = removed_bg-structure(modprobe_t,*bfit)
    a32.plot(modprobe_t[blim],res2[blim])
    f2.tight_layout()
