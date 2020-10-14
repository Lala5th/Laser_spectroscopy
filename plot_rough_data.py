#!/usr/bin/ipython3

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit

args = sys.argv
argnum = 2
usage = '''Usage:
.\{0} dataset [E_scaling]
or in ipython console:
%run {0} dataset [E_scaling]'''.format(args[0])

if len(args) < argnum:
    print(usage)
    sys.exit()

scaling = 1
if len(args) >= 3:
    scaling = float(args[2])

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
    fit, cov = curve_fit(fitfunc,x,y,p0=[50,0,1,(xmin+xmax)/2,0.001])
    plt.plot(x,fitfunc(x,*fit))
    for param in zip(params,fit,sp.sqrt(np.diag(cov))):
        print(param[0], ':', param[1], '+-', param[2])

gaussian = lambda x,m,s : norm.pdf(x,m,s)

data = np.loadtxt(args[1],dtype=[('t',np.float64,()),('I',np.float64,())],skiprows=1,delimiter=',')
#def fitfunc(x,B,M1,M2,M3,M4,A1,A2,A3,A4,S,T0):
#    return(B-A1*gaussian(x,T0-M1,S)-A1*gaussian(x,T0-M2,S)-A1*gaussian(x,T0-M3,S)-A1*gaussian(x,T0-M4,S)-A1*gaussian(x,T0+M1,S)-A1*gaussian(x,T0+M2,S)-A1*gaussian(x,T0+M3,S)-A1*gaussian(x,T0+M4,S)).flatten()

#fit, cov = curve_fit(fitfunc,data['t'],data['I'],p0=[48,0.0034,0.0088,0.0217,0.032,1,2,1,0.25,0.001,-0.0034])

data['t'] *= scaling
plt.ion()
plt.plot(data['t'][5:-4],moving_average(data['I'],n=10))
if scaling == 1:
    plt.xlabel("Time [s]")
else:
    plt.xlabel("Energy [Hz]")
plt.ylabel("Intensity [a.u.]")
plt.show()
