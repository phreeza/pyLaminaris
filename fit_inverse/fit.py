import sys
sys.path.append('/home/mccolgan/work/owls/')
import numpy as np
import json
from multielectrode_analysis.lib.clicks import load_clicks
from scipy.optimize import basinhopping,minimize
import scipy
from scipy import signal

params_data = json.load(open(sys.argv[-1]))
filt_freq = 2000.
xlim = (4,8)

z = (np.arange(33)*50.)
stims,t,Fs = load_clicks(params_data)
stims = np.transpose(stims,(1,0,2))
means = stims.mean(axis=1).T

t = t[15200:15800]
means = means[15200:15800]

means = means/np.std(means)

t = np.repeat(t[:,np.newaxis],z.shape[0],1)
z = np.repeat(z[np.newaxis,:],t.shape[0],0)
dists = z[:32,:32]-z[:32,:32].T

#for r in np.arange(10.,200.,10.):
r = 20.
def calc_phi(x):
    r = x[0]
    ivel = x[1]/100.
    vdot = x[2:t.shape[0]+2]
    n = np.abs(x[t.shape[0]+2:t.shape[0]+t.shape[1]+2])
    #n /= n.max()
    #n[n<0] = 0.
    w = 1. / np.sqrt(dists ** 2 + r ** 2)
    w /= w.max()
    vdot_prop = np.zeros_like(t)
    phi = np.zeros_like(z)
    
    #vdot_prop = np.repeat(vdot[:,np.newaxis],t.shape[1],1)
    vdot_prop[:,0] = vdot[:]
    vdot_inter = scipy.interpolate.interp1d(np.arange(600.),vdot,fill_value=0.,bounds_error=False)
    for xn in range(1,phi.shape[1]):
        shift = (z[0,xn]*ivel)
        vdot_prop[:,xn] = vdot_inter(np.arange(600.)+shift)
    phi = np.dot(np.diff(vdot_prop*n,axis=1),w)
    return phi
    #for tn in range(phi.shape[0]):
    #    phi[tn, :] = np.convolve(i[tn,:], w, mode='same')
def func(x):
    phi = calc_phi(x)
    return np.mean((means-phi)**2) 

x0 = np.abs(np.random.randn(t.shape[0]+t.shape[1]+2))
#u,s,v = np.linalg.svd(np.cumsum(means,axis=1))
x0[2:t.shape[0]+2] = -means.mean(axis=1)
x0[t.shape[0]+2:] = [ 0.,   0.35507181,   0.93417839,   1.51885425,
         2.19865669,   3.0813484 ,   4.18111583,   5.36838795,
         6.6368594 ,   7.87631038,   9.09188345,  10.12999729,
        10.95484373,  11.52748123,  11.85101632,  11.82552131,
        11.51326291,  10.83071441,   9.89191209,   8.95150976,
         8.09254976,   7.26675513,   6.51613888,   5.80996425,
         5.13953978,   4.498051  ,   3.84134159,   3.19273031,
         2.53670732,   1.87643104,   1.17558873,   0.44915316,  0.]
x0[1] = -0.2
x0[0] = r

#print np.mean((means)**2)
#print np.mean((means-np.outer(u[:,0],v[0,:])*s[0])**2)
#print func(x0)
#ret = basinhopping(func,x0,disp=True)
print "r",r
ret = minimize(func,x0,options={'disp':True})

phi = calc_phi(ret.x)
print "R^2:",np.corrcoef(phi.flatten(),means.flatten())[0,1]

print "vel:",Fs/(ret.x[1]/100)*1e-6
print "dist:",ret.x[0]

def save_fit(x,suffix,r):
    r = x[0]
    ivel = x[1]/100.
    vdot = x[2:t.shape[0]+2]
    n = np.abs(x[t.shape[0]+2:t.shape[0]+t.shape[1]+2])

    np.savez('fit_'+str(r)+'_'+suffix,r=r,ivel=ivel,vdot=vdot,n=n)

save_fit(ret.x,params_data['prefix'],r)
