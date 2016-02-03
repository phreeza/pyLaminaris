import sys
sys.path.append('/home/mccolgan/work/owls/')
import pyXdPhys as xd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from pyPhonic import click
from scipy import signal
import json

params_data = json.load(open(sys.argv[-2]))
fnames = params_data['fnames']
filt_freq = params_data['filt_freq']
name_base = params_data['name_base']

stims = [xd.Stimulation(name_base+fn+'.itd') for fn in fnames]
means = np.array([n.traces.mean(axis=0) for n in stims]).T
depths = np.array([s.params['DV'] for s in stims]) - 10200

b,a = signal.butter(5,2*filt_freq/stims[1].params['daFc'],btype='lowpass')
fted = signal.lfilter(b,a,means,axis=0)
b_hp,a_hp = signal.butter(5,2*filt_freq/stims[1].params['daFc'],btype='highpass')
fted_hp = signal.lfilter(b_hp,a_hp,means,axis=0)
times = stims[0].times - 11.

#ax = plt.subplot(2,3,1)
#plt.plot(times,-means/100 + depths)
#plt.xlim((0,4))
##plt.ylim((0,1600))
#plt.gca().invert_yaxis()
#plt.grid()
#plt.title('30Hz-15kHz')
#plt.ylabel('penetration depth[$\mu m$]')
#plt.xticks(np.arange(1.,4.))
#for label in ax.get_xticklabels():
#    label.set_visible(False)

ax = plt.subplot(2,2,1)
plt.gca().invert_yaxis()
#for label in ax2.get_yticklabels():
#    label.set_visible(False)
plt.plot(times,-fted/100 + depths,color='blue',lw=2)
plt.xlim((0,4))
#plt.ylim((0,1600))
plt.gca().invert_yaxis()
plt.grid()
plt.title('30Hz-2kHz')
plt.ylabel('penetration depth[$\mu m$]')
plt.xticks(np.arange(1.,4.))
for label in ax.get_xticklabels():
    label.set_visible(False)

ax3 = plt.subplot(2,2,2,sharey=ax)
for label in ax3.get_yticklabels():
    label.set_visible(False)
plt.plot(times,-fted_hp/100 + depths,color='green',lw=2)
plt.xlim((0,4))
#plt.ylim((0,1600))
plt.gca().invert_yaxis()
plt.grid()
plt.title('2kHz-15kHz')
plt.xticks(np.arange(1.,4.))
for label in ax3.get_xticklabels():
    label.set_visible(False)

#plt.suptitle(params_data['suptitle'])
#plt.gcf().set_size_inches((11,8))
plt.gcf().set_size_inches((8,6))
plt.tight_layout()
plt.subplots_adjust(top=0.9)

params_sim = json.load(open(sys.argv[-1]))
fname = "data/bundle_pulse_"+params_sim['postfix']+".npz"
filt_freq = params_sim['filt_freq']

if 'raw_scale' in params_sim.keys():
    raw_scale = params_sim['raw_scale']
else:
    raw_scale = 4000.

if 'filt_scale' in params_sim.keys():
    filt_scale = params_sim['filt_scale']
else:
    filt_scale = 2000.
    filt_scale_sim = 1500.

means = np.load(fname)['pot'].T
locs = np.load(fname)['loc'].T


#means = means[500:,3::20]

means = means - means[:3000,:].mean(axis=0)

depths = locs[0,:]# np.arange(4000.,6000.,10.) - 4100.
base_index = ((params_sim['population_params']['root_x_offset']-params_sim['electrode_params']['y_base'])
              /params_sim['electrode_params']['y_d'])
depth_index = ((int(200./params_sim['electrode_params']['y_d'])*np.arange(7)+base_index)*params_sim['electrode_params']['x_N']+4).astype(int)
depths = depths - depths[depth_index[0]] + 200.
print depths[depth_index]
times = np.arange(means.shape[0])*2.5e-6

#b,a = signal.butter(5,filt_freq*times[1],btype='lowpass')

b_30hz,a_30hz = signal.butter(3,2*200.*times[1],btype='highpass')
b,a = signal.butter(3,2*filt_freq*times[1],btype='lowpass')
b_hp,a_hp = signal.butter(3,2*filt_freq*times[1],btype='highpass')

times = times - 0.0102

means = signal.lfilter(b_30hz,a_30hz,means,axis=0)

fted = signal.lfilter(b,a,means,axis=0)
fted_hp = signal.lfilter(b_hp,a,means,axis=0)
#ax4 = plt.subplot(2,3,4,sharex=ax)
#plt.plot(times*1000.,-means[:,depth_index]/raw_scale+ depths[depth_index])
#plt.xlim((0,4))
#plt.gca().invert_yaxis()
#plt.ylim((1600,00))
#plt.grid()
#plt.xlabel('time[ms]')
#plt.ylabel('penetration depth[$\mu m$]')
#plt.xticks(np.arange(1.,4.))

ax5 = plt.subplot(2,2,3,sharex=ax)
plt.gca().invert_yaxis()
#for label in ax5.get_yticklabels():
#    label.set_visible(False)
plt.plot(times*1000.,-fted[:,depth_index]/filt_scale_sim + depths[depth_index],color='blue',lw=2)
plt.xlim((0,4))
plt.gca().invert_yaxis()
plt.ylim((1600,00))
plt.grid()
plt.xlabel('time[ms]')
plt.ylabel('penetration depth[$\mu m$]')
plt.xticks(np.arange(1.,4.))

ax6 = plt.subplot(2,2,4,sharey=ax5,sharex=ax3)
for label in ax6.get_yticklabels():
    label.set_visible(False)

plt.plot(times*1000.,-fted_hp[:,depth_index]/filt_scale_sim + depths[depth_index],color='green',lw=2)
plt.xlim((0,4))
plt.gca().invert_yaxis()
plt.ylim((1600,00))
plt.grid()
plt.xlabel('time[ms]')
plt.xticks(np.arange(1.,4.))
#plt.suptitle(params_sim['suptitle'])
#plt.gcf().set_size_inches((11,8))
plt.gcf().set_size_inches((8,6))
plt.tight_layout()
plt.subplots_adjust(top=0.9,hspace=0.19,wspace=0.05)

plt.figtext(0.08 ,0.92,'A',fontsize=24)
plt.figtext(0.56,0.92,'B',fontsize=24)
#plt.figtext(0.7,0.92,'C',fontsize=24)
plt.figtext(0.08 ,0.48,'C',fontsize=24)
plt.figtext(0.56,0.48,'D',fontsize=24)
#plt.figtext(0.7,0.48,'F',fontsize=24)

for fmt in ['.png','.pdf']:
    plt.savefig('figs/traces_combined_'+params_sim['postfix']+'_'+params_data['postfix']+fmt)
plt.show()
plt.close()
