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

def run_fig3(pot_bif, pot_nobif, n_rows, n_cols, plot_row_lf = 42, plot_row_hf = 59, plot_col = 100, **kwargs):
    dt = 0.0025
    filt_freq = 2000.

    filter_type = ['lowpass', 'highpass']

    scalemax = 1.0

    amps = []
    for pot in [pot_bif,pot_nobif]:
        locs = []
        amps.append([])
        for freq_n in range(2):
            locs.append([])
            amps[-1].append([])
            b, a = signal.butter(3, 2 * filt_freq * dt / 1000., btype=filter_type[freq_n])
            fted = signal.lfilter(b, a, pot['pot'], axis=1)
            fted = fted-fted.mean(axis=1).reshape((-1, 1))
            ft_abs = np.abs((fted - fted.mean(axis=1).reshape((-1, 1))))

            for m in range(n_cols):
                amps[-1][-1].append([])
                locs[-1].append([])
                for n in range(n_rows):
                    index = get_index(n, m, n_rows, n_cols)  # 20 * n + m
                    locs[-1][-1].append(pot['loc'][index, :])
                    amps[-1][-1][-1].append(fted[index, 4000:6000])
                    # amps[-1][-1].append(fted[index,3500:-2499].std())
        #return np.array(amps), np.array(locs)
    amps = np.array(amps)
    amps_psd = np.sum(amps**2,axis=4)
    locs = np.array(locs)

    def norm(a):
        return a/np.max(np.abs(a).squeeze(),axis=1)[:,np.newaxis]
    def norm_first(x):
        return x/x[0]
    
    fig = plt.figure()
    ax1 = plt.subplot(141)
    plot_loc_offset = 0
    plt.plot(locs[0,0,plot_col:,0]-locs[0,0,plot_col-plot_loc_offset,0],norm_first(amps_psd[0,0,0,plot_col:]),label='low-pass',lw=2)
    plt.plot(locs[0,0,plot_col:,0]-locs[0,0,plot_col-plot_loc_offset,0],norm_first(amps_psd[0,1,0,plot_col:]),label='high-pass',lw=2)

    plt.loglog()
    plt.ylabel('power (normalized)')
    plt.xlabel('axial distance [um]')
    
    ax2 = plt.subplot(143,sharex=ax1,sharey=ax1)
    plt.plot(locs[0,:,0,1],norm_first(amps_psd[0,0,:,plot_row_lf]),label='low-pass',lw=2)
    plt.plot(locs[0,:,0,1],norm_first(amps_psd[0,1,:,plot_row_hf]),label='high-pass',lw=2)
    plt.xlabel('radial distance [um]')
    plt.loglog()

    ax3 = plt.subplot(142,sharex=ax1,sharey=ax1)
    plot_loc_offset = 0
    plt.plot(locs[0,0,plot_col:,0]-locs[0,0,plot_col-plot_loc_offset,0],norm_first(amps_psd[1,0,0,plot_col:]),label='low-pass',lw=2)
    plt.plot(locs[0,0,plot_col:,0]-locs[0,0,plot_col-plot_loc_offset,0],norm_first(amps_psd[1,1,0,plot_col:]),label='high-pass',lw=2)
    plt.legend(loc=3,bbox_to_anchor=(-0.95, .02))

    plt.loglog()
    #plt.xlabel('axial distance [um]')
    ax4 = plt.subplot(144,sharex=ax1,sharey=ax1)
    plt.plot(locs[0,:,0,1],norm_first(amps_psd[1,0,:,plot_row_lf]),label='low-pass',lw=2)
    plt.plot(locs[0,:,0,1],norm_first(amps_psd[1,1,:,plot_row_hf]),label='high-pass',lw=2)
    #plt.xlabel('radial distance [um]')
    plt.loglog()


    #plt.xlim((locs[0,0,plot_col,0]-locs[0,0,plot_col-plot_loc_offset,0],locs[0,0,-1,0]-locs[0,0,plot_col-plot_loc_offset,0]))
    plt.xlim(locs[0,1,0,1],locs[0,-1,0,1])
    plt.ylim(1e-4,1.)
    ax1.xaxis.set_label_coords(1.0, -0.05) 
    ax2.xaxis.set_label_coords(1.0, -0.05) 
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    #ax1.set_xticks(np.arange(0,15000,5000))
    #ax2.set_xticks(np.arange(0,15000,5000))
    #ax3.set_xticks(np.arange(0,15000,5000))
    #ax4.set_xticks(np.arange(0,15000,5000))
    #ax1.set_yticks([0.2, 0.55, 0.76])
    fig.subplots_adjust(hspace=0.05, wspace = 0.03)
    plt.figtext(0.125, 0.91, 'A', fontsize=24)
    plt.figtext(0.32, 0.91, 'B', fontsize=24)
    plt.figtext(0.52, 0.91, 'C', fontsize=24)
    plt.figtext(0.71, 0.91, 'D', fontsize=24)
    return fig


def run(params_bif_fname,params_nobif_fname):
    params_bif = json.load(open(params_bif_fname))
    params_nobif = json.load(open(params_nobif_fname))

    n_rows = params_bif['electrode_params']['y_N']
    n_cols = params_bif['electrode_params']['x_N']

    fname = "data/bundle_pulse_" + params_bif['postfix'] + ".npz"
    potentials_bif = np.load(fname)
    fname = "data/bundle_pulse_" + params_nobif['postfix'] + ".npz"
    potentials_nobif = np.load(fname)

    fig = run_fig3(potentials_bif, potentials_nobif, n_rows, n_cols, 
                      **params_bif)
    for fmt in ['.png', '.pdf']:
        fig.savefig('figs/manuscript_fig3_' + params_bif['postfix'] + '_' + params_nobif['postfix'] + fmt)
    plt.close(fig)

if __name__ == '__main__':
    run(sys.argv[-2],sys.argv[-1])
