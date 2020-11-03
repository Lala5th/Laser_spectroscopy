#!/usr/bin/ipython3

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import scipy.constants as const

args = sys.argv
argnum = 2
etalon_len = 0.2
etalon_freq = 1e-6*const.c/(2*etalon_len)
print(etalon_freq)
usage = '''Usage:
.\{0} dataset [etalon]
or in ipython console:
%run {0} dataset [etalon]'''.format(args[0])

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

scaling = 1
if len(args) > argnum:
    etalon = np.loadtxt(args[argnum],dtype=[('t',np.float64,()),('I',np.float64,())],skiprows=1,delimiter=',')
    etalon_mod = etalon[0:-1]
    peaks, props = find_peaks(etalon_mod['I'],np.max(etalon_mod['I'])*0.75,distance = 500)
    diff = [etalon_mod['t'][peaks[i]] - etalon_mod['t'][peaks[i-1]] for i in range(1,peaks.size)]
    scaling = etalon_freq/np.abs(np.mean(diff))
    print("Scaling:",scaling,"+-",np.std(diff)*scaling,'MHz/s')

fitfunc = lambda x, a, b, A, m, s: a*x + b - A*gaussian(x,m,s)
params = ['a','b','A','m','s']

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def fit_line(xmin,xmax):
    x = []
    y = []
    for i in range(data['t'].size):
        if(data['t'][i] > xmin and data['t'][i] < xmax):
            x.append(data['t'][i])
            y.append(data['I'][i])
        elif(data['t'][i] >= xmax):
            break
    x = np.array(x)
    y = np.array(y)
    fit, cov = curve_fit(fitfunc,x,y,p0=[0,np.max(y),1,(xmax+xmin)/2,0.001*scaling])
    plt.plot(x,fitfunc(x,*fit))
    for param in zip(params,fit,sp.sqrt(np.diag(cov))):
        print(param[0], ':', param[1], '+-', param[2])
    print("Chi-square/Ndof :", get_chi_squared(fitfunc,fit,x,y))
    return (fit,cov)

gaussian = lambda x,m,s : norm.pdf(x,m,s)

data = np.loadtxt(args[1],dtype=[('t',np.float64,()),('I',np.float64,())],skiprows=1,delimiter=',')
#def fitfunc(x,B,M1,M2,M3,M4,A1,A2,A3,A4,S,T0):
#    return(B-A1*gaussian(x,T0-M1,S)-A1*gaussian(x,T0-M2,S)-A1*gaussian(x,T0-M3,S)-A1*gaussian(x,T0-M4,S)-A1*gaussian(x,T0+M1,S)-A1*gaussian(x,T0+M2,S)-A1*gaussian(x,T0+M3,S)-A1*gaussian(x,T0+M4,S)).flatten()

#fit, cov = curve_fit(fitfunc,data['t'],data['I'],p0=[48,0.0034,0.0088,0.0217,0.032,1,2,1,0.25,0.001,-0.0034])
data['t'] *= scaling
if len(args) > argnum:
    etalon_mod['t'] *= scaling
plt.ion()
plt.plot(data['t'][5:-4],moving_average(data['I'],n=10))
if len(args) > argnum:
    plt.plot(etalon_mod['t'],etalon_mod['I'])
    plt.plot(etalon_mod['t'][peaks],etalon_mod['I'][peaks],'x')
if scaling == 1:
    plt.xlabel("Time [s]")
else:
    plt.xlabel("Energy [MHz]")
plt.ylabel("Intensity [a.u.]")
plt.show()

def fit_splitting(amin,amax,bmin,bmax):
    print("Params for line A:")
    afit, acov = fit_line(amin, amax)
    print("Params for line B:")
    bfit, bcov = fit_line(bmin, bmax)
    print('')
    print('Splitting:', np.abs(afit[3] - bfit[3]), '+-', np.sqrt(acov[3,3]**2 + bcov[3,3]**2))
