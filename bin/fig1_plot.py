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

def calc_mua(x,dt=0.0025,f_hp=2000.,f_lp=500.):
    b, a = signal.butter(3, 2 * f_hp * dt / 1000., btype='highpass')
    ret = signal.lfilter(b, a, x, axis=1)
    ret[ret<0.] = 0.
    b, a = signal.butter(3, 2 * f_lp * dt / 1000., btype='lowpass')
    ret = 20*signal.lfilter(b, a, ret, axis=1)
    return ret


def run_fig1(pot, n_rows, n_cols, **params):
    fig = plt.figure()
    gs = [gridspec.GridSpec(10, 14, top=0.9, bottom=0.52),
          gridspec.GridSpec(10, 14, top=0.48, bottom=0.1)]

    dt = 0.0025
    filt_freq = params['filt_freq']

    filter_type = ['lowpass', 'highpass']

    scalemax = 1.0
    norm = mpl.colors.Normalize(vmin=0, vmax=scalemax)

    ax1 = [None, None]
    ax2 = [None, None]
    ax3 = [None, None]
    for freq_n in range(2):
        if freq_n == 0:
            b, a = signal.butter(3, 2 * filt_freq * dt / 1000., btype=filter_type[freq_n])
            fted = signal.lfilter(b, a, pot['pot'], axis=1)
        else:
            fted = calc_mua(pot['pot'])
        ax1[freq_n] = plt.subplot(gs[freq_n][:, :5])

        node_locs = []
        #for ngroup in pot['nodes'][:1]:
        ngroup =  pot['nodes'][0]
        for segment in ngroup[0]:
            segment = np.array(segment)
            node_locs.append(segment)
            ax1[freq_n].plot(segment[:, 1], segment[:, 0], color='k', alpha=1.0)

        # for n in range(70):
        # plt.subplot(gs[n%10,int(n/10)+3])
        # plt.axis('off')
        # m = 900+(n%10)*40+2*int(n/10)+1
        # plt.plot(psds[m,:])
        # #plt.fill([0,100,100,0],[0,0,100000,100000],'k')
        # plt.semilogy()
        #
        #for fmt in ['.png', '.pdf']:
        #    plt.savefig('figs/bundle_pulse_psds' + fmt)

        #plt.close()
        times = np.arange(0., 20., dt)
        locs = []
        for n in range(10):
            for m in range(3):
                ax_t = plt.subplot(gs[freq_n][n, 2 * m + 5:2 * m + 7])
                ax_t.axis('off')
                index = get_index(36 + 3 * n, 2 * m + 4, n_rows, n_cols)  # 760 + 60 * n + 2 * m + 1
                locs.append(pot['loc'][index, :])
                lh = ax_t.plot(times[3500:-2499], fted[index, 3500:-2500] - fted[index, :].mean())[0]
                ymin, ymax = plt.ylim()
                dh = ax1[freq_n].plot(locs[-1][1], locs[-1][0], '.')[0]
                lh.set_color(plt.cm.jet((ymax - ymin) / (scalemax * 1e6)))
                dh.set_color(plt.cm.jet((ymax - ymin) / (scalemax * 1e6)))
                # ax_t.set_ylim(-5e5 / (m + 1), 5e5 / (m + 1))
                #plt.fill([0,100,100,0],[0,0,100000,100000],'k')

        amps = fted[:, 3500:-2500].max(axis=1) - fted[:, 3500:-2500].min(axis=1)

        ft_abs = np.abs((fted - fted.mean(axis=1).reshape((-1, 1))))
        index = [get_index(36 + n, 3, n_rows, n_cols) for n in range(30)]
        if freq_n == 0:
            amp_line = np.array([fted[n, ft_abs[n, :].argmax()] for n in index])
            amp_line = -amp_line + amp_line.mean()
        else:
            amp_line = np.array([np.abs(fted[n, ft_abs[n, :].argmax()]) for n in index])
        amp_line = amp_line / 1e6

        locs_line = pot['loc'][index, 0]

        ax1[freq_n].contour(pot['loc'][:, 1].reshape((n_rows, -1)), pot['loc'][:, 0].reshape((n_rows, -1)),
                            amps.reshape((n_rows, -1)) / (scalemax * 1e6), 0.7 ** np.arange(7), cmap=plt.cm.jet,
                            norm=norm)
        ax1[freq_n].fill([1280, 1280, 1600, 1600], [14000, 15500, 15500, 14000], 'white', edgecolor='white')

        locs = np.array(locs)
        #for segment in pot['nodes'][0][0]:
        #  segment = np.array(segment)
        #  plt.plot(segment[:,1],segment[:,0],color='r',alpha=1.)
        node_locs = np.vstack(node_locs)
        #plt.plot(locs[:,1],locs[:,0],'.')
        #ax1[freq_n].fill([600,600,1100,1100],[16300,16350,16350,16300],'k')
        ax1[freq_n].axis('off')
        #plt.ylim(pot['loc'][-1,0]+100,pot['loc'][480,0]-100)
        ax1[freq_n].set_aspect('equal', 'datalim')
        ax1[freq_n].set_xlim(-100, 1300)
        ax1[freq_n].set_ylim(locs[-1, 0] + 100, locs[0, 0] - 100)

        #ax2[freq_n] = plt.subplot(gs[freq_n][:, 8:11], sharey=ax1[freq_n])
        bins = np.arange(0, 19000, 100.)
        hist, edges = np.histogram(node_locs[:, 0], bins=bins)
        #plt.bar(bottom=bins[:-1], width=hist / 100., height=bins[1], left=0, orientation='horizontal', color='k',
        #        edgecolor='w')
        #ax2[freq_n].tick_params(left="off")
        #plt.xticks(np.arange(0, 18, 4))
        #plt.ylim(locs[-1, 0] + 100, locs[0, 0] - 100)
        #if freq_n == 1:
        #    plt.xlabel(u'nodes\n[1/\u03BCm]')
        #else:
        #    plt.setp(ax2[freq_n].get_xticklabels(), visible=False)

        ax3[freq_n] = plt.subplot(gs[freq_n][:, 11:], sharey=ax1[freq_n])
        plt.subplots_adjust(left=0.05, right=0.95)
        if freq_n == 0:
            plt.bar(bottom=bins[1:-1], width=np.diff(hist) / 100., height=bins[1], left=0, orientation='horizontal',
                    color='k', edgecolor='w')
            plt.xlabel(u'bif.-term.\n[1/\u03BCm]')
            plt.xticks(np.arange(-4, 4.1, 2))
        else:
            plt.bar(bottom=bins[:-1], width=hist / 100., height=bins[1], left=0, orientation='horizontal', color='k',
                    edgecolor='w')
            plt.xlabel(u'nodes\n[1/\u03BCm]')
            ax3[freq_n].tick_params(left="off")
            plt.xticks(np.arange(0, 18, 4))
            plt.ylim(locs[-1, 0] + 100, locs[0, 0] - 100)
        ax3_twin = ax3[freq_n].twiny()
        ax3_twin.plot(amp_line, locs_line, lw=2, color='r')

        plt.setp(ax1[freq_n].get_yticklabels(), visible=False)
        # plt.setp(ax2[freq_n].get_yticklabels(), visible=False)
        plt.setp(ax3[freq_n].get_yticklabels(), visible=False)
        plt.setp(ax1[freq_n].get_xticklabels(), visible=False)
        # plt.setp( ax2[freq_n].get_xticklabels(), visible=False)
        # plt.setp( ax3[freq_n].get_xticklabels(), visible=False)
        plt.ylim(locs[-1, 0] + 100, locs[0, 0] - 100)

    # [left, bottom, width, height]
    cb_ax = fig.add_axes([0.15, .05, .4, .02])

    def form1(x, pos):
        """ This function returns a string with 1 decimal places, given the input x"""
        return '%.1f' % x

    from matplotlib.ticker import FuncFormatter

    cb1 = mpl.colorbar.ColorbarBase(cb_ax, cmap=plt.cm.jet,
                                    norm=norm,
                                    orientation='horizontal', ticks=[0.0, scalemax])
    cb1.set_label('Amplitude [mV]', labelpad=-5)

    formatter = FuncFormatter(form1)
    cb_ax.yaxis.set_major_formatter(FuncFormatter(formatter))

    plt.gcf().set_size_inches((8, 8))

    plt.figtext(0.1, 0.92, 'A', fontsize=24)
    plt.figtext(0.38, 0.92, 'B', fontsize=24)
    # plt.figtext(0.57, 0.92, 'C', fontsize=24)
    plt.figtext(0.77, 0.92, 'C', fontsize=24)

    # sb1_ax = plt.gcf().add_axes([0.05, .0, 0.9 / 14 * 5, .1], sharex=ax1[0])
    # sb1_ax.axis('off')
    # sb1_ax.set_aspect('equal', 'datalim')
    # sb1_ax.set_xlim(-100, 1300)
    # sb1_ax.fill([600, 600, 1100, 1100], [16300, 16350, 16350, 16300], 'k')
    # sb1_ax.text(600, 16370, u'500\u03BCm')
    #
    # sb2_ax = plt.gcf().add_axes([0.5, .0, 0.059, .1], sharex=ax_t, sharey=sb1_ax)
    # sb2_ax.axis('off')
    # t0 = ax_t.get_xlim()[0]
    # sb2_ax.fill([t0, t0, t0 + 5, t0 + 5], [16300, 16350, 16350, 16300], 'k')
    # sb2_ax.text(t0, 16370, '5ms')
    return fig

def run(params_fname):
    params = json.load(open(params_fname))

    n_rows = params['electrode_params']['y_N']
    n_cols = params['electrode_params']['x_N']

    fname = "data/bundle_pulse_" + params['postfix'] + ".npz"
    potentials = np.load(fname)
    fig = run_fig1(potentials, n_rows, n_cols, **params)
    for fmt in ['.png', '.pdf']:
        fig.savefig('figs/manuscript_fig1_' + params['postfix'] + fmt)
    plt.close(fig)

if __name__ == '__main__':
    run(sys.argv[-1])
