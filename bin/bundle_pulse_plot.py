import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

import matplotlib as mpl
from scipy import signal


def run_fig1(pot):
    gs = [gridspec.GridSpec(10, 14, top=0.9, bottom=0.52),
          gridspec.GridSpec(10, 14, top=0.48, bottom=0.1)]

    dt = 0.0025
    filt_freq = 2000.

    filter_type = ['lowpass', 'highpass']

    scalemax = 1.0
    norm = mpl.colors.Normalize(vmin=0, vmax=scalemax)

    ax1 = [None, None]
    ax2 = [None, None]
    ax3 = [None, None]
    for freq_n in range(2):
        b, a = signal.butter(3, 2 * filt_freq * dt / 1000., btype=filter_type[freq_n])
        fted = signal.lfilter(b, a, pot['pot'], axis=1)
        ax1[freq_n] = plt.subplot(gs[freq_n][:, :5])

        node_locs = []
        for ngroup in pot['nodes']:
            for segment in ngroup[0]:
                segment = np.array(segment)
                node_locs.append(segment)
                ax1[freq_n].plot(segment[:, 1], segment[:, 0], color='k', alpha=0.1)

        # for n in range(70):
        # plt.subplot(gs[n%10,int(n/10)+3])
        # plt.axis('off')
        # m = 900+(n%10)*40+2*int(n/10)+1
        # plt.plot(psds[m,:])
        # #plt.fill([0,100,100,0],[0,0,100000,100000],'k')
        #    plt.semilogy()
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
                index = 760 + 60 * n + 2 * m + 1
                locs.append(pot['loc'][index, :])

                lh = ax_t.plot(times[3500:-2499], fted[index, 3500:-2500] - fted[index, :].mean())[0]
                ymin, ymax = plt.ylim()
                dh = ax1[freq_n].plot(locs[-1][1], locs[-1][0], '.')[0]
                lh.set_color(plt.cm.jet((ymax - ymin) / (scalemax * 1000000.)))
                dh.set_color(plt.cm.jet((ymax - ymin) / (scalemax * 1000000.)))
                #plt.fill([0,100,100,0],[0,0,100000,100000],'k')

        amps = fted[:, 3500:-2500].max(axis=1) - fted[:, 3500:-2500].min(axis=1)

        ft_abs = np.abs((fted - fted.mean(axis=1).reshape((-1,1))))
        index = [760 + 20 * n + 1 + 1 for n in range(30)]
        if freq_n == 0:
            amp_line = np.array([fted[n, ft_abs[n, :].argmax()] for n in index])
            amp_line = -amp_line + amp_line.mean()
            amp_line = 3*amp_line/amp_line.max()
        else:
            amp_line = np.array([np.abs(fted[n, ft_abs[n, :].argmax()]) for n in index])
            amp_line = 15*amp_line/amp_line.max()

        locs_line = pot['loc'][index, 0]

        ax1[freq_n].contour(pot['loc'][:, 1].reshape((70, -1)), pot['loc'][:, 0].reshape((70, -1)),
                            amps.reshape((70, -1)) / (scalemax * 1000000.), 0.7 ** np.arange(7), cmap=plt.cm.jet,
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

        plt.plot(amp_line, locs_line,lw=2,color='r')

        plt.setp(ax1[freq_n].get_yticklabels(), visible=False)
        #plt.setp(ax2[freq_n].get_yticklabels(), visible=False)
        plt.setp(ax3[freq_n].get_yticklabels(), visible=False)
        plt.setp(ax1[freq_n].get_xticklabels(), visible=False)
        #plt.setp( ax2[freq_n].get_xticklabels(), visible=False)
        #plt.setp( ax3[freq_n].get_xticklabels(), visible=False)
        plt.ylim(locs[-1, 0] + 100, locs[0, 0] - 100)

    # [left, bottom, width, height]
    cb_ax = plt.gcf().add_axes([0.15, .05, .4, .02])

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

    for fmt in ['.png', '.pdf']:
        plt.savefig('figs/bundle_pulse_potentials' + fmt)
    plt.close()


def run_fig2(pot):
    def psd(x):
        return np.abs(np.fft.rfft(x, axis=1)) ** 2

    ind = [0, 2, 20, 60, 100]
    colors_ind = ['r', 'b', 'g', 'purple', 'y']

    psds = psd(pot['pot'])
    psd_freqs = np.fft.fftfreq(pot['pot'].shape[1], 0.0025)[:psds.shape[1]]

    # f,axarr = (10,5)
    axarr = gridspec.GridSpec(1, 30, bottom=0.15)
    # cb_ax= plt.gcf().add_axes([0.44,.94,.46,.03])
    # norm = mpl.colors.LogNorm(vmin=psd_freqs[1], vmax=psd_freqs[100])
    # cb1 = mpl.colorbar.ColorbarBase(cb_ax, cmap=plt.cm.autumn,
    # norm=norm,
    # orientation='horizontal')
    #cb_ax.xlim(psd_freqs[1],psd_freqs[100])
    #cb_ax.semilogx()
    #cb1.set_label('frequency [kHz]')
    #plt.setp( cb_ax.get_xticklabels(), visible=False)
    plt.subplot(axarr[0, :9])
    plt.plot([1e1, 1e5], 3 * np.array([5e0, 5e-4]), color='k', linestyle='dashed')
    plt.plot([1e1, 1e5], 11 * np.array([5e0, 5e-8]), color='k', linestyle='dotted')

    #for n in np.arange(-20,20,0.5):
    #  plt.plot([1e0,1e5],10**n*np.array([1e0,1e-5]),color='k',alpha=0.1)
    #  plt.plot([1e0,1e5],10**n*np.array([1e0,1e-10]),color='k',alpha=0.1)

    n = 5

    k = 480 + (n % 10) * 100

    def normalize(x):
        return x / x[0]

    for m in range(len(ind) - 1):
        plt.plot(pot['loc'][k:k + 20, 1], normalize(psds[k:k + 20, ind[m]:ind[m + 1]].mean(axis=1)), '-',
                 color=colors_ind[m], label=str(round(psd_freqs[10 * m], 1)) + 'kHz')

    plt.loglog()
    #if n == 0:
    #  plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #     ncol=2, mode="expand", borderaxespad=0.)
    if n == 5:
        plt.ylabel('Power [mV^2]')

    plt.xlim((50, 2e4))
    plt.ylim((1e-6, 1e0))
    plt.xlabel(u'distance from trunk [\u03BCm]')
    #for fmt in ['.png', '.pdf']:
    #    plt.savefig('figs/bundle_pulse_scalings' + fmt)

    #plt.close()

    slopes = np.array([[np.polyfit(np.log(pot['loc'][n * 20 + 1:n * 20 + 20, 1]),
                                   np.log(psds[n * 20 + 1:n * 20 + 20, m:m + 20].mean(axis=1)), 1)[0] for m in
                        range(100)] for n in range(30, 60)])
    plt.subplot(axarr[:, 12:])
    plt.plot(psd_freqs[1:101], slopes.mean(axis=0))
    plt.fill_between(psd_freqs[1:101], slopes.mean(axis=0) + slopes.std(axis=0),
                     slopes.mean(axis=0) - slopes.std(axis=0), alpha=0.2)
    plt.semilogx()
    plt.xlim(psd_freqs[1], psd_freqs[100])

    ind.append(-1)
    for m in range(len(ind) - 1):
        plt.fill([max(0.001, psd_freqs[ind[m]]), psd_freqs[ind[m + 1]], psd_freqs[ind[m + 1]],
                  max(0.001, psd_freqs[ind[m]])], [-2.5, -2.5, -2.6, -2.6], colors_ind[m], edgecolor=None)
        if ind[m] > 0:
            plt.axvline(psd_freqs[ind[m]], color='k')

    plt.xlabel('f[kHz]')
    plt.ylabel('Distance scaling coefficient', labelpad=0)

    plt.gcf().set_size_inches((8, 4))

    plt.figtext(0.08, 0.92, 'A', fontsize=24)
    plt.figtext(0.4, 0.92, 'B', fontsize=24)

    for fmt in ['.png', '.pdf']:
        plt.savefig('figs/bundle_pulse_freq_slope' + fmt)

    plt.close()


if __name__ == '__main__':
    potentials = np.load('data/parallel_bundle_potential.npz')
    run_fig1(potentials)
    run_fig2(potentials)
