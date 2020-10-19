#!/usr/bin/ipython3

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import scipy.constants as const

args = sys.argv
argnum = 5
etalon_len = 0.2
etalon_freq = 1e-6*const.c/(2*etalon_len)
print(etalon_freq)
usage = '''Usage:
.\{0} probe reference cull_low cull_high [etalon]
or in ipython console:
%run {0} probe reference cull_low cull_high [etalon]'''.format(args[0])

def get_chi_squared(func,params,x,y):
    assert(x.shape == y.shape)
    error = np.abs(x[0] - x[np.abs(x[0]-x).argmin()])/2
    chi_squared = np.sum((y - func(x,params))**2)/error**2
    ndof = x.size - params.size
    return chi_squared/ndof

if len(args) < argnum:
    print(usage)
    sys.exit()

scaling = 1
if len(args) > argnum:
    etalon = np.loadtxt(args[argnum],dtype=[('t',np.float64,()),('I',np.float64,())],skiprows=1,delimiter=',')
    etalon_mod = etalon[30000:-10000]
    peaks, props = find_peaks(etalon_mod['I'],np.max(etalon_mod['I'])*0.75,distance = 500)
    diff = [etalon_mod['t'][peaks[i]] - etalon_mod['t'][peaks[i-1]] for i in range(1,peaks.size)]
    scaling = etalon_freq/np.abs(np.mean(diff))
    print("Scaling:",scaling,"+-",np.std(diff)*scaling,'MHz/s')


fitfunc = lambda x, a, b, A, m, s: a*x + b + A*gaussian(x,m,s)
params = ['a','b','A','m','s']

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

gaussian = lambda x,m,s : norm.pdf(x,m,s)

probe = np.loadtxt(args[1],dtype=[('t',np.float64,()),('I',np.float64,())],skiprows=1,delimiter=',')
ref = np.loadtxt(args[2],dtype=[('t',np.float64,()),('I',np.float64,())],skiprows=1,delimiter=',')

probe['t'] *= scaling
ref['t'] *= scaling

ref['I'] -= np.mean(ref['I'])

reffunc = interp1d(ref['t'],ref['I'],fill_value="extrapolate")

#assert(probe['t'].all()==ref['t'].all())
#probe = probe[250:int(len(probe)/3)]
#ref = ref[:int(-250+len(ref)/3)]

#modprobe_t = probe['t'][10:-10]
#modprobe = moving_average(probe['I'],21)
#modref = moving_average(ref['I'],21)
cull_low = int(args[3])
cull_high = -int(args[4])

modprobe_t = probe['t'][cull_low:cull_high]
modprobe = probe['I'][cull_low:cull_high] - np.mean(probe['I'][cull_low:cull_high])
modref = ref['I']

bg = lambda x, t0, a, b : (a)*reffunc(x + t0) - b

#def fitfunc(x,B,M1,M2,M3,M4,A1,A2,A3,A4,S,T0):
#    return(B-A1*gaussian(x,T0-M1,S)-A1*gaussian(x,T0-M2,S)-A1*gaussian(x,T0-M3,S)-A1*gaussian(x,T0-M4,S)-A1*gaussian(x,T0+M1,S)-A1*gaussian(x,T0+M2,S)-A1*gaussian(x,T0+M3,S)-A1*gaussian(x,T0+M4,S)).flatten()

#fit, cov = curve_fit(fitfunc,data['t'],data['I'],p0=[48,0.0034,0.0088,0.0217,0.032,1,2,1,0.25,0.001,-0.0034])

bgopt,bgcov = curve_fit(bg,modprobe_t,modprobe,p0=[-0.0001,1,1],maxfev=10000)

if len(args) > argnum:
    etalon_mod['t'] *= scaling

plt.ion()
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(modprobe_t,modprobe)
ax1.plot(modprobe_t,bg(modprobe_t,*bgopt))
if len(args) > argnum:
    ax1.plot(etalon_mod['t'],etalon_mod['I'])
    ax1.plot(etalon_mod['t'][peaks],etalon_mod['I'][peaks],'x')
ax2.plot(modprobe_t,modprobe-bg(modprobe_t,*bgopt))
removed_bg = modprobe-bg(modprobe_t,*bgopt)
if scaling == 1:
    plt.xlabel("Time [s]")
else:
    plt.xlabel("Energy [Hz]")
ax1.set_ylabel("Intensity [a.u.]")
ax2.set_ylabel("Intensity [a.u.]")
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
    fit, cov = curve_fit(fitfunc,x,y,p0=[0,0,1,(xmin+xmax)/2,0.0001*scaling],maxfev = 10000)
    ax2.plot(x,fitfunc(x,*fit))
    for param in zip(params,fit,sp.sqrt(np.diag(cov))):
        print(param[0], ':', param[1], '+-', param[2])
    return (fit,cov)

def fit_splitting(amin,amax,bmin,bmax):
    print("Params for line A:")
    afit, acov = fit_line(amin, amax)
    print("Params for line B:")
    bfit, bcov = fit_line(bmin, bmax)
    print('')
    print('Splitting:', np.abs(afit[3] - bfit[3]), '+-', np.sqrt(acov[3,3]**2 + bcov[3,3]**2))
