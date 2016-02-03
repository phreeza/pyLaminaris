import matplotlib as mpl

#mpl.use('Agg')
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

from scipy import signal
import json
import sys


def get_index(row, col, n_rows, n_cols):
    return n_cols * row + col

def run_fig2(pot, n_rows, n_cols, plot_row_lf = 42, plot_row_hf = 59, plot_col = 100, **kwargs):
    dt = 0.0025
    filt_freq = 2000.

    filter_type = ['lowpass', 'highpass']

    scalemax = 1.0

    locs = []
    amps = []
    for freq_n in range(2):
        locs.append([])
        amps.append([])
        b, a = signal.butter(3, 2 * filt_freq * dt / 1000., btype=filter_type[freq_n])
        fted = signal.lfilter(b, a, pot['pot'], axis=1)
        fted = fted-fted.mean(axis=1).reshape((-1, 1))
        ft_abs = np.abs((fted - fted.mean(axis=1).reshape((-1, 1))))

        for m in range(n_cols):
            amps[-1].append([])
            locs[-1].append([])
            for n in range(n_rows):
                index = get_index(n, m, n_rows, n_cols)  # 20 * n + m
                locs[-1][-1].append(pot['loc'][index, :])
                amps[-1][-1].append(fted[index, 4000:6000])
                # amps[-1][-1].append(fted[index,3500:-2499].std())
    #return np.array(amps), np.array(locs)
    amps = np.array(amps)
    amps_psd = np.sum(amps**2,axis=3)
    locs = np.array(locs)

    def norm(a):
        return a/np.max(np.abs(a).squeeze(),axis=1)[:,np.newaxis]
    def norm_first(x):
        return x/x[0]
    
    fig = plt.figure()
    ax1 = plt.subplot(141)
    plot_loc_offset = 0
    plt.plot(locs[0,0,plot_col:,0]-locs[0,0,plot_col-plot_loc_offset,0],norm_first(amps_psd[0,0,plot_col:]),label='low-pass')
    plt.plot(locs[0,0,plot_col:,0]-locs[0,0,plot_col-plot_loc_offset,0],norm_first(amps_psd[1,0,plot_col:]),label='high-pass')

    plt.ylabel('power (normalized)')
    plt.loglog()
    plt.xlabel('axial distance [um]')
    #plt.xlim((locs[0,0,plot_col,0]-locs[0,0,plot_col-plot_loc_offset,0],locs[0,0,-1,0]-locs[0,0,plot_col-20,0]))
    plt.xlim(locs[0,1,0,1],locs[0,-1,0,1])
    plt.ylim(1e-7,1.)
    ax2 = plt.subplot(142,sharey=ax1)
    #ax2 = plt.subplot(222, sharey=ax1)
    plt.plot(locs[0,:,0,1],norm_first(amps_psd[0,:,plot_row_lf]),label='low-pass')
    plt.plot(locs[0,:,0,1],norm_first(amps_psd[1,:,plot_row_hf]),label='high-pass')
    plt.xlim(locs[0,1,0,1],locs[0,-1,0,1])
    plt.ylim(1e-7,1.)
    plt.xlabel('radial distance [um]')
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.legend(loc=3)
    plt.loglog()

    #ax3 = plt.subplot(223, sharex=ax1)
    #plt.plot(locs[0,0,plot_col:,0]-locs[0,0,plot_col,0],norm_first(amps_psd[1,0,plot_col:])/norm_first(amps_psd[0,0,plot_col:]),color='r')
    #plt.xlabel('axial distance [um]')
    #plt.ylabel('power ratio (normalized)')
    ##plt.xlim(0,10000)
    #plt.loglog()
    #ax4 = plt.subplot(224, sharex=ax2)
    #plt.plot(locs[0,:,0,1],norm_first(amps_psd[1,:,plot_row_hf])/norm_first(amps_psd[0,:,plot_row_lf]),color='r')
    #plt.xlabel('radial distance [um]')
    #plt.loglog()
    #plt.ylim(0,1.1)
    #plt.xlim(0,10000)

    #plt.setp(ax1.get_xticklabels(), visible=False)
    #plt.setp(ax2.get_xticklabels(), visible=False)
    #plt.setp(ax2.get_yticklabels(), visible=False)
    #plt.setp(ax4.get_yticklabels(), visible=False)
    #ax3.set_xticks(np.arange(0,15000,5000))
    #ax4.set_xticks(np.arange(0,15000,5000))
    #ax1.set_yticks([0.2, 0.55, 0.76])
    #fig.subplots_adjust(hspace=0.05, wspace = 0.03)
    return fig


def run(params_fname):
    params = json.load(open(params_fname))

    n_rows = params['electrode_params']['y_N']
    n_cols = params['electrode_params']['x_N']

    fname = "data/bundle_pulse_" + params['postfix'] + ".npz"
    potentials = np.load(fname)

    fig = run_fig2(potentials, n_rows, n_cols, 
                      **params)
    for fmt in ['.png', '.pdf']:
        fig.savefig('figs/manuscript_fig2_' + params['postfix'] + fmt)
    plt.close(fig)

if __name__ == '__main__':
    run(sys.argv[-1])
