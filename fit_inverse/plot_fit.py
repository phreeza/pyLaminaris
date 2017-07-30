import sys
sys.path.append('/home/mccolgan/work/owls/')
import numpy as np
import json
from multielectrode_analysis.lib.clicks import load_clicks
from scipy.optimize import basinhopping,minimize
import scipy
from scipy import signal
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Liberation Sans']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
from matplotlib import pyplot as plt
plt.ion()
#import seaborn as sns
import matplotlib.gridspec as gridspec
import glob

#sns.set_style("whitegrid")

params_data = json.load(open(sys.argv[-1]))

filt_freq_lp = 1000.
filt_freq_hp = 2500.

xlim = (4,7)

z = (np.arange(33)*50.)
stims,t,Fs = load_clicks(params_data)
stims = np.transpose(stims,(1,0,2))
means = stims.mean(axis=1).T

t = t[15200:15800]
means = means[15200:15800]
scale = means.std()
means = means/scale

t = np.repeat(t[:,np.newaxis],z.shape[0],1)
z = np.repeat(z[np.newaxis,:],t.shape[0],0)
dists = z[:32,:32]-z[:32,:32].T

#for r in np.arange(10.,200.,10.):
def calc_phi(r,ivel,vdot,n):
    #n /= n.max()
    #n[n<0] = 0.
    w = 1. / np.sqrt(dists ** 2 + r ** 2)
    w /= w.max()
    vdot_prop = np.zeros_like(t)
    phi = np.zeros_like(z)
    
    #vdot_prop = np.repeat(vdot[:,np.newaxis],t.shape[1],1)
    vdot_prop[:,0] = vdot[:]
    vdot_inter = scipy.interpolate.interp1d(np.arange(600.),vdot,fill_value=(0.,0.),bounds_error=False)
    for xn in range(1,phi.shape[1]):
        shift = (z[0,xn]*ivel)
        vdot_prop[:,xn] = vdot_inter(np.arange(600.)+shift)
    I = np.diff(vdot_prop*n,axis=1)*4*np.pi*0.33/50.*scale

    phi = np.dot(I,w)/(4*np.pi*0.33)*50./scale # mV = (nA/um)/(um*(S/m))*um
    p = (z[0:1,:-1]*(I-I.sum(axis=1,keepdims=True))).sum(axis=1)*50.
    #plt.figure()
    #plt.plot(I.sum(axis=1))
    #plt.show()
    #plt.figure()
    print "Maximal dipole moment:",np.max(np.abs(p)), "nA*um"
    print "                     =",np.max(np.abs(p))*1e-6,"uA*mm"
    return phi
    #for tn in range(phi.shape[0]):
    #    phi[tn, :] = np.convolve(i[tn,:], w, mode='same')


fit_fname = glob.glob('fit_*_'+params_data['prefix']+'*')[0]
fitfile = np.load(fit_fname)
print 'r:',fitfile['r'],z[0,1]
phi = calc_phi(fitfile['r'],fitfile['ivel'],fitfile['vdot'],fitfile['n'])

b,a = signal.butter(5,2*filt_freq_lp/Fs,btype='lowpass')
b_hp,a_hp = signal.butter(5,2*filt_freq_hp/Fs,btype='highpass')


plt.figure()
gs = gridspec.GridSpec(5,2)
ax1 = plt.subplot(gs[1:,0])
plt.xlim(xlim)
gl = plt.plot(t[:,0],-10*signal.lfilter(b,a,means,axis=0) + z[0,:-1], color='orange',label='data', lw=2)
rl = plt.plot(t[:,0],-10*signal.lfilter(b,a,phi,axis=0) + z[0,:-1], color='black',label='model', lw=1)
plt.ylabel('Penetration depth [$\mu m$]')
plt.title('Low frequency (<1kHz)')
plt.xlabel('Time [ms]')
plt.xticks(np.arange(*xlim))
plt.grid(axis='x')
ax1.text(-0.0, 1.05, 'C', transform=ax1.transAxes,
  fontsize=16, fontweight='bold', va='top', ha='right')

ax2 = plt.subplot(gs[1:,1],sharey=ax1)
plt.xlim(xlim)
gl = plt.plot(t[:,0],-signal.lfilter(b_hp,a_hp,means,axis=0)*10 + z[0,:-1], color='orange', lw=2)
rl = plt.plot(t[:,0],-signal.lfilter(b_hp,a_hp,phi,axis=0)*10 + z[0,:-1], color='black', lw=1)
plt.plot([7.2,7.2],np.array([0,1])*10./scale+1400. , color='black', lw=4,solid_capstyle='butt', clip_on=False)
plt.text(7.3,1400+10/scale,'1mV', clip_on=False)
ax2.arrow(3.85, 700, 0., 400, head_width=.1, head_length=40,lw=2, fc='k', ec='k',clip_on=False)
ax2.arrow(3.85, 1000, 0., -400, head_width=.1, head_length=40,lw=2, fc='k', ec='k',clip_on=False)
plt.text(3.75,430,'dorsal', clip_on=False, rotation='vertical')
plt.text(3.75,1200,'ventral', clip_on=False, rotation='vertical')
plt.xlabel('Time [ms]')
plt.title('High frequency (>2.5kHz)')
for label in ax2.get_yticklabels():
    label.set_visible(False)
plt.ylim(-150,1590)
plt.legend((gl[0],rl[0]),('data','model'),loc='upper right')
plt.xticks(np.arange(*xlim))
plt.grid(axis='x')

ax2.text(-0.0, 1.05, 'D', transform=ax2.transAxes,
  fontsize=16, fontweight='bold', va='top', ha='right')
plt.gca().invert_yaxis()

ax3 = plt.subplot(gs[0,0])
plt.plot(t[:,0],-fitfile['vdot'], lw=2)
plt.xlabel('Time [ms]')
plt.ylabel('Membrane\npotential [a.u.]')
plt.xlim(xlim)
plt.xticks(np.arange(*xlim))
plt.yticks(np.arange(3))
plt.grid(axis='x')
ax3.text(-0.0, 1.4, 'A', transform=ax3.transAxes,
  fontsize=16, fontweight='bold', va='top', ha='right')

ax4 = plt.subplot(gs[0,1])
plt.plot(z[0,:], fitfile['n'], lw=2)
plt.xlabel('Penetration depth [$\mu m$]')
plt.ylabel('Fiber number [a.u.]')
plt.xticks(np.arange(3)*500.)
plt.yticks(np.arange(1,5)*2.5)
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax4.text(-0.0, 1.4, 'B', transform=ax4.transAxes,
  fontsize=16, fontweight='bold', va='top', ha='right')

plt.subplots_adjust(wspace=0.1,hspace=1.5, top=0.95)
plt.gcf().set_size_inches((6,9))
for fmt in ['.png', '.pdf']:
    plt.savefig('figs/manuscript_fig4_'+params_data['prefix'] + fmt)
plt.show()
