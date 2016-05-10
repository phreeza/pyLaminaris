import matplotlib as mpl

# mpl.use('Agg')
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

from scipy import signal
import json
import sys


def get_index(row, col, n_rows, n_cols):
    return n_cols * row + col


def calc_mua(x, dt=0.0025, f_hp=2000., f_lp=500.):
    b, a = signal.butter(3, 2 * f_hp * dt / 1000., btype='highpass')
    ret = signal.lfilter(b, a, x, axis=1)
    ret[ret < 0.] = 0.
    b, a = signal.butter(3, 2 * f_lp * dt / 1000., btype='lowpass')
    ret = signal.lfilter(b, a, ret, axis=1)
    return ret


def run_fig1(pot, n_rows, n_cols, **params):
    fig = plt.figure()
    gs = [gridspec.GridSpec(10, 14, top=0.9, bottom=0.52),
          gridspec.GridSpec(10, 14, top=0.48, bottom=0.1)]
    raw_pot = pot['pot']/1e6
    dt = 0.0025
    filt_freq = params['filt_freq']

    filter_type = ['lowpass', 'highpass']

    scalemax = .075
    norm = mpl.colors.Normalize(vmin=0, vmax=scalemax)

    ax1 = [None, None]
    ax2 = [None, None]
    ax3 = [None, None]

    tree_n = [1, 12]
    for freq_n in range(2):
        if freq_n == 0:
            b, a = signal.butter(3, 2 * filt_freq * dt / 1000., btype=filter_type[freq_n])
            fted = signal.lfilter(b, a, raw_pot, axis=1)
        else:
            fted = calc_mua(raw_pot)
        ax1[freq_n] = plt.subplot(gs[freq_n][:, :5])

        node_locs = []
        # ngroup =  pot['nodes'][0]
        for n, ngroup in enumerate(pot['nodes']):
            for m, segment_group in enumerate(ngroup):
                for segment in segment_group:
                    segment = np.array(segment)
                    node_locs.append(segment)
                    # if n==1 and m==1:
                    # ax1[freq_n].plot(segment[:, 1], segment[:, 0], color='k', alpha=1.0)

        def jittered_line(start, end, sigma=0.5):
            # generate some naturalistic looking jitter on the axon path
            def norm(x):
                return np.sqrt(np.sum(x ** 2))

            line = np.arange(0., norm(end - start), 10.)[np.newaxis, :] * (end - start)[:, np.newaxis] / norm(
                end - start)
            jitter = sigma * np.random.randn(*line.shape)
            jitter[:, 0] = 0.
            jitter = np.cumsum(jitter, axis=1)
            jitter = np.cumsum(jitter, axis=1)
            jitter -= np.linspace(0., 1., jitter.shape[1])[np.newaxis] * jitter[:, -1:]

            return start[:, np.newaxis] + (line + jitter)

        n_neuron = tree_n[freq_n] #np.random.randint(len(pot['segments'][0]))
        print "Neuron number:", n_neuron
        s, e = pot['segments'][0][n_neuron][0]
        l = jittered_line(s, e, 0.01)
        ax1[freq_n].plot(l[1], l[0], color='k', alpha=1.0, lw=0.5)

        for s, e in pot['segments'][0][n_neuron][1:]:
            l = jittered_line(s, e)
            ax1[freq_n].plot(l[1], l[0], color='k', alpha=1.0, lw=0.5)
        times = np.arange(0., 20., dt)
        locs = []
        for n in range(10):
            for m in range(3):
                ax_t = plt.subplot(gs[freq_n][n, 2 * m + 5:2 * m + 7])
                ax_t.axis('off')
                index = get_index(36 + 3 * n, 2 * m + 4, n_rows, n_cols)  # 760 + 60 * n + 2 * m + 1
                locs.append(pot['loc'][index, :])
                lh = ax_t.plot(times[3500:-2499], fted[index, 3500:-2500])[0]
                ymin, ymax = plt.ylim()
                dh = ax1[freq_n].plot(locs[-1][1], locs[-1][0], '.')[0]
                lh.set_color(plt.cm.jet(norm((ymax - ymin))))
                dh.set_color(plt.cm.jet(norm((ymax - ymin))))
                # ax_t.set_ylim(-5e5 / (m + 1), 5e5 / (m + 1))
                #plt.fill([0,100,100,0],[0,0,100000,100000],'k')

        amps = fted[:, 3500:-2500].max(axis=1) - fted[:, 3500:-2500].min(axis=1)

        ft_abs = np.abs((fted - fted.mean(axis=1).reshape((-1, 1))))
        index = [get_index(36 + n, 2, n_rows, n_cols) for n in range(30)]
        if freq_n == 0:
            amp_line = np.array([fted[n, ft_abs[n, 3500:].argmax() + 3500] for n in index])
            amp_line = -amp_line + amp_line.mean()
        else:
            amp_line = np.array([np.abs(fted[n, ft_abs[n, :].argmax()]) for n in index])
        print amp_line

        locs_line = pot['loc'][index, 0]

        ax1[freq_n].contour(pot['loc'][:, 1].reshape((n_rows, -1)), pot['loc'][:, 0].reshape((n_rows, -1)),
                            amps.reshape((n_rows, -1)), scalemax*(0.7 ** np.arange(7)), cmap=plt.cm.jet,
                            norm=norm)
        ax1[freq_n].fill([1280, 1280, 1600, 1600], [14000, 15500, 15500, 14000], 'white', edgecolor='white')

        locs = np.array(locs)
        node_locs = np.vstack(node_locs)
        ax1[freq_n].axis('off')
        ax1[freq_n].set_aspect('equal', 'datalim')
        ax1[freq_n].set_xlim(-100, 1300)
        ax1[freq_n].set_ylim(locs[-1, 0] + 100, locs[0, 0] - 100)

        bins = np.arange(0, 19000, 100.)
        hist, edges = np.histogram(node_locs[:, 0], bins=bins)

        ax3[freq_n] = plt.subplot(gs[freq_n][:, 11:], sharey=ax1[freq_n])
        plt.subplots_adjust(left=0.05, right=0.95)
        if freq_n == 0:
            plt.bar(bottom=bins[1:-1], width=np.diff(hist) / 1., height=bins[1], left=0, orientation='horizontal',
                    color='k', edgecolor='w')
            plt.xlabel(u'bif.-term.\n[1/\u03BCm]')
            plt.xticks(np.arange(-4, 4.1, 2))
        else:
            plt.bar(bottom=bins[:-1], width=np.array(hist) / (float(bins[1])*params['population_size']), height=bins[1], left=0, orientation='horizontal', color='k',
                    edgecolor='w')
            plt.xlabel('nodes [1/$\mu$m]')
            ax3[freq_n].tick_params(left="off")
            plt.xticks(np.arange(0, .31, 0.1))
        ax3_twin = ax3[freq_n].twiny()
        ax3_twin.plot(amp_line, locs_line, lw=2, color='r')

        plt.setp(ax1[freq_n].get_yticklabels(), visible=False)
        plt.setp(ax3[freq_n].get_yticklabels(), visible=False)
        plt.setp(ax1[freq_n].get_xticklabels(), visible=False)
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
    plt.figtext(0.77, 0.92, 'C', fontsize=24)

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
