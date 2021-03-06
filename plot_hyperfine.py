#!/usr/bin/ipython3

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.stats import norm,cauchy
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

args = sys.argv
argnum = 5
usage = '''Usage:
.\{0} probe reference cull_low cull_high [E_scaling]
or in ipython console:
%run {0} probe reference cull_low cull_high [E_scaling]'''.format(args[0])

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
if len(args) >= argnum + 1:
    scaling = float(args[5])
lorentzian = lambda x, m, s : cauchy.pdf(x,m,s)

fitfunc = lambda x, a, b, A, m, s: a*x + b + A*lorentzian(x,m,s)
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

plt.ion()
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(modprobe_t,modprobe)
ax1.plot(modprobe_t,bg(modprobe_t,*bgopt))
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

def fit_structure(xmin,xmax,m1=0,m2=0,m3=0,p0=None):
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
    HWHM = [0,0.0001*scaling]
    bounds = np.array([[-np.inf,sp.inf],pos,pos,pos,pos,pos,pos,[m1-0.0001*scaling,m1+0.0001*scaling],[m2-0.0001*scaling,m2+0.0001*scaling],[m3-0.0001*scaling,m3+0.0001*scaling],HWHM,HWHM,HWHM,HWHM,HWHM,HWHM]).transpose()
    if(p0 == None):
        fit, cov = curve_fit(structure,x,y,p0=[-1,100,100,100,100,100,100,m1,m2,m3,0.00005*scaling,0.00005*scaling,0.00005*scaling,0.00005*scaling,0.00005*scaling,0.00005*scaling],bounds = bounds,maxfev = 10000)
    else:
        fit, cov = curve_fit(structure,x,y,p0=p0,bounds = bounds,maxfev = 10000)
    ax2.plot(x,structure(x,*fit))
    for param in zip(struct_params,fit,sp.sqrt(np.diag(cov))):
        print(param[0], ':', param[1], '+-', param[2])
    print("Splittings")
    print('1-2 :', np.abs(fit[7] - fit[8]), '+-', np.sqrt(cov[7,7]**2 + cov[8,8]**2))
    print('2-3 :', np.abs(fit[9] - fit[8]), '+-', np.sqrt(cov[9,9]**2 + cov[8,8]**2))
    print('1-3 :', np.abs(fit[7] - fit[9]), '+-', np.sqrt(cov[7,7]**2 + cov[9,9]**2))
    print("Chi-square/Ndof :", get_chi_squared(structure,fit,x,y))
