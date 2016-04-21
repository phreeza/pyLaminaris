import matplotlib
#matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
from scalebars import add_scalebar

lw=1.5
subfig_fontsize = 12
snip_spacing = 40

def draw_snip(x0,x1,y,slant):
    plt.plot([x0,x1],[y+snip_spacing,y+snip_spacing+slant],color='white',lw=lw)
    plt.plot([x0,x1],[y,y+slant],color='black',lw=1)
    plt.plot([x0,x1],[y+2*snip_spacing,y+2*snip_spacing+slant],color='black',lw=1)

pots =  np.load('pots.npz')['pots']
pots2 = np.load('pots_2.npz')['pots']
#pots3 = np.load('pots_3.npz')['pots'].sum(axis=1)
pots3 = np.load('pots_bifterm.npz')['pots']
pots4 = np.load('pots_3_sum.npz')['pots']
scale = 50.
scale2 = scale*8
scale3 = scale*300
scale4 = scale3*10
window = np.exp(-(np.arange(1000)-500.)**2/2.**2)
window /= window.sum()
window = window*5

window2 = np.exp(-(np.arange(1000)-500.)**2/2.**2)
window2 /= window2.sum()
window2 = window2*5

axes_w = 0.1
axes_h = 0.52
axes_ws = 0.08
axes_hs = 0.46

ax1 = plt.axes([2*axes_w,axes_h,axes_ws,axes_hs]) #plt.subplot(2,8,2*order[n_p]+1,sharey=ax1)
for n in range(80,120,4):
    plt.plot(np.arange(pots.shape[1])*0.0025,pots[n,:].T/scale-(n-100)*100.,color='black')
plt.plot(np.arange(pots.shape[1])*0.0025,pots[100,:].T/scale,color='black',lw=lw)
plt.xlim(1.4,1.9)
plt.axis('off')

add_scalebar(ax1,sizex=.2,sizey=0.02e6/scale,matchx=False,matchy=False,labelx='0.2 ms',
        labely='20 $\mu$V',loc=8,bbox_to_anchor=(120.5,170.5))

plt.axes([4*axes_w,axes_h,axes_ws,axes_hs],sharey=ax1) #plt.subplot(2,8,2*order[n_p]+1,sharey=ax1)
for n in range(180,220,4):
    plt.plot(np.arange(pots.shape[1])*0.0025,pots[n,:].T/scale-(n-200)*100.,color='black')
plt.plot(np.arange(pots.shape[1])*0.0025,pots[200,:].T/scale,color='black',lw=lw)
plt.xlim(2.0,2.5)
plt.axis('off')

ax2 = plt.axes([6*axes_w,axes_h,axes_ws,axes_hs],sharey=ax1) #plt.subplot(2,8,2*order[n_p]+1,sharey=ax1)
for n in range(80,120,4):
    plt.plot(np.arange(pots.shape[1])*0.0025,pots2[n,:].T/scale2-(n-100)*100.,color='black')
plt.plot(np.arange(pots.shape[1])*0.0025,pots2[100,:].T/scale2,color='black',lw=lw)
plt.xlim(1.4,1.9)
plt.axis('off')

add_scalebar(ax2,sizex=.2,sizey=0.2e6/scale2,matchx=False,matchy=False,labelx='0.2 ms',
        labely='0.2 mV',loc=8,bbox_to_anchor=(250.5,170.5))

plt.axes([8*axes_w,axes_h,axes_ws,axes_hs],sharey=ax1) #plt.subplot(2,8,2*order[n_p]+1,sharey=ax1)
for n in range(80,120,4):
    plt.plot(np.arange(pots3.shape[1])*0.0025,pots3[n,:].T/scale2-(n-100)*100.,color='black')
plt.plot(np.arange(pots3.shape[1])*0.0025,pots3[96,:].T/scale2+400.,color='black',lw=lw)
plt.plot(np.arange(pots3.shape[1])*0.0025,pots3[104,:].T/scale2-400.,color='black',lw=lw)
plt.xlim(1.4,1.9)
plt.axis('off')

order = [1,0,2,3]
c = ['r','g','b']
d_sep = 20
for n_p in range(4):
    ax = plt.axes([axes_w*(2*order[n_p]+1)+0.017,axes_h,axes_ws,axes_hs],sharey=ax1) #plt.subplot(2,8,2*order[n_p]+1,sharey=ax1)
    plt.arrow(-15,1800,0,-800,width=1,head_width=8,facecolor='black',head_length=100)
    if n_p == 0:
        plt.plot([0,0],[2000,0],color='black',lw=lw)
        plt.plot([0],[0],color='black',lw=lw,marker='o',mec='none',ms=3.5)
        draw_snip(-5,5,1820,10)
    elif n_p == 1:
        plt.plot([0,0],[-2000,2000],color='black',lw=lw)

        draw_snip(-5,5,1820,10)
        draw_snip(-5,5,-1860,10)
    elif n_p == 2:
        y_sep = [0.,100.,100.,1900.]
        x_sep = [20,10,5]
        p = [[[0.,0.] for m in range(2**n)] for n in range(4)]
        p[0][0][1] = 100.
        plt.plot([0,0],[2000,100],color='black',lw=lw)
        for n in range(1,4):
            for m in range(2**n):
                p[n][m][1] = p[n-1][0][1] - y_sep[n-1]
                if m%2 == 0:
                    p[n][m][0] = p[n-1][m/2][0] + x_sep[n-1]
                else:
                    p[n][m][0] = p[n-1][(m-1)/2][0] - x_sep[n-1]
                    plt.plot([p[n][m][0],p[n][m-1][0]],[p[n][m][1],p[n][m-1][1]],color='black',lw=lw)
                plt.plot([p[n][m][0],p[n][m][0]],[p[n][m][1],p[n][m][1]-y_sep[n]],color='black',lw=lw)

        draw_snip(-5,5,1820,10)
        draw_snip(p[-1][-1][0]-5,p[-1][0][0]+5,-1880,30)

        #plt.plot([5,15],[-1835,-1845],color='white',lw=lw)
        #plt.plot([5,15],[-1860,-1870],color='black',lw=lw)
        #plt.plot([5,15],[-1810,-1820],color='black',lw=lw)
    elif n_p == 3:
        y_sep = [0.,100.,100.,700.]
        x_sep = [20,10,5]
        p = [[[0.,0.] for m in range(2**n)] for n in range(4)]
        p[0][0][1] = 500.
        plt.plot([0,0],[2000,p[0][0][1]],color='black',lw=lw)
        for n in range(1,4):
            for m in range(2**n):
                p[n][m][1] = p[n-1][0][1] - y_sep[n-1]
                if m%2 == 0:
                    p[n][m][0] = p[n-1][m/2][0] + x_sep[n-1]
                else:
                    p[n][m][0] = p[n-1][(m-1)/2][0] - x_sep[n-1]
                    plt.plot([p[n][m][0],p[n][m-1][0]],[p[n][m][1],p[n][m-1][1]],color='black',lw=lw)
                plt.plot([p[n][m][0],p[n][m][0]],[p[n][m][1],p[n][m][1]-y_sep[n]],color='black',lw=lw)
                if n == 3:
                    plt.plot([p[n][m][0]],[p[n][m][1]-y_sep[n]],color='black',lw=lw,marker='o',mec='none',ms=3.5)

        draw_snip(-5,5,1820,10)

        #plt.plot([p[-1][0][0]+5,p[-1][-1][0]-5],[-1855,-1835],color='white',lw=lw)
        #plt.plot([p[-1][0][0]+5,p[-1][-1][0]-5],[-1880,-1860],color='black',lw=lw)
        #plt.plot([p[-1][0][0]+5,p[-1][-1][0]-5],[-1830,-1810],color='black',lw=lw)

        #for k in range(-1,2):
        #    p1 = 1000*np.random.random()-500.
        #    p2 = 1000*np.random.random()
        #    p3 = 1000*np.random.random()
        #    plt.plot([5*k,5*k],[2000,p1],color=c[k+1],lw=lw)
        #    plt.plot([5*k+10,5*k+10],[p1-p2,p1],color=c[k+1],lw=lw)
        #    plt.plot([5*k-10,5*k-10],[p1-p3,p1],color=c[k+1],lw=lw)
        #    plt.plot([5*k-10,5*k+10],[p1,p1],color=c[k+1],lw=lw)
        #    plt.plot([5*k+10],[p1-p2],marker='o',color=c[k+1],lw=lw,mec='none')
        #    plt.plot([5*k-10],[p1-p3],marker='o',color=c[k+1],lw=lw,mec='none')
        #plt.plot([-10,10],[1855,1825],color='white',lw=lw)
        #plt.plot([-10,10],[1880,1850],color='black',lw=lw)
        #plt.plot([-10,10],[1830,1800],color='black',lw=lw)
    plt.xlim(-40,40)
    plt.axis('off')

ax = plt.subplot(256,sharey=ax1)

y_sep = [0.,500.,500.,1000.]
x_sep = [20,10,5]
for k in range(3):
    p = [[[5*(k-1),0.] for m in range(2**n)] for n in range(4)]
    p[0][0][1] = 1000.
    plt.plot([5*(k-1),5*(k-1)],[2000,p[0][0][1]],color=c[k],lw=lw)
    for n in range(1,4):
        for m in range(2**n):
            if m%2 == 0:
                p[n][m+1][1] = p[n][m][1] = p[n-1][0][1] - y_sep[n-1] + 200*np.random.randn()
                p[n][m][0] = p[n-1][m/2][0] + x_sep[n-1]
                plt.plot([p[n-1][m/2][0],p[n-1][m/2][0]],[p[n][m][1],p[n-1][m/2][1]],color=c[k],lw=lw)
            else:
                p[n][m][0] = p[n-1][(m-1)/2][0] - x_sep[n-1]
                plt.plot([p[n][m][0],p[n][m-1][0]],[p[n][m][1],p[n][m-1][1]],color=c[k],lw=lw)
            if n == 3:
                lm = y_sep[n] + 200*np.random.randn()
                plt.plot([p[n][m][0],p[n][m][0]],[p[n][m][1],p[n][m][1]-lm],color=c[k],lw=lw)
                plt.plot([p[n][m][0]],[p[n][m][1]-lm],color=c[k],lw=lw,marker='o',mec='none',ms=3.5)

draw_snip(-10,10,1820,15)
plt.xlim(-45,45)
plt.axis('off')

plt.subplot(258,sharey=ax1)
csd = -np.diff(np.diff(pots4,axis=0),axis=0)
m = np.max(np.abs(csd))
plt.imshow(csd,vmin=-m,vmax=m,interpolation='nearest',origin='upper',aspect='auto',extent=(0,5,-12000,10000))
plt.xlim(1.4,1.9)
plt.axis('off')

ax3 = plt.subplot(257,sharey=ax1)
ax3.eventplot(np.array([[1.5,1.5,1.5]]).T, colors=c, lineoffsets=(2000.-160.*np.arange(1,4)),
              linelengths=[150,150,150])
for n in range(88,120,4):
    plt.plot(np.arange(pots4.shape[1])*0.0025,pots4[n,:].T/(scale3)-(n-100)*100.,color='black')
plt.xlim(1.4,1.9)
plt.ylim(-2000,2100)
plt.axis('off')

add_scalebar(ax3,sizex=.2,sizey=2e6/(scale3),matchx=False,matchy=False,labelx='0.2 ms',
        labely='2 mV',loc=8,bbox_to_anchor=(120.5,20.0))

from pyLaminaris import helper
import glob
files = glob.glob('pots_3_[1-9]*')
pots = np.sum([np.load(f)['pots'] for f in files],axis=0)
csd_i = np.sum([np.load(f)['csd'] for f in files],axis=0)
#csd_i = np.load('pots_3.npz')['csd']
#pots = pots-pots[:,:,100:101]
pots -= pots[:,:,1000:].mean(axis=2,keepdims=True)
def pulse(t):
    return 2*(0.1
            + 0.9 * np.exp(-((t - 10.)/1.) ** 2))
ret = np.zeros((pots.shape[0],int(20./0.0025)))
ret_csd = np.zeros((pots.shape[0],int(20./0.0025)))
for n in range(pots.shape[0]):
    ret[n,:] = np.convolve(pots[n,:,500:].sum(axis=0),pulse(np.arange(ret.shape[1]+pots.shape[2]-500-1)*0.0025)*0.0025*1000,mode='valid')
    ret_csd[n,:] = np.convolve(csd_i[n,:,500:].sum(axis=0),pulse(np.arange(ret.shape[1]+pots.shape[2]-500-1)*0.0025)*0.0025*1000,mode='valid')
#n_reps = 1000
#for rep in range(n_reps):
#    if rep%100 == 0:
#        print rep
#    for n_neuron in range(pots.shape[1]):
#        for t in helper.inhom_poisson(pulse,20.,0.,6.):
#            t_ind = int(t/0.0025)
#            l = min(600,ret.shape[1]-t_ind)
#            if l>0: ret[:,t_ind:t_ind+l] += pots[:,n_neuron,500:500+l]/n_reps
print pots.shape

ax4 = plt.subplot(2,5,9,sharey=ax1)
for n in range(88,120,4):
    plt.plot(np.arange(ret.shape[1])*0.0025-10.,ret[n,:]/(scale4)-(n-100)*100.,color='black')
plt.xlim(-10,10)
plt.axis('off')

times = [helper.inhom_poisson(pulse,20.,0.,6.)-10. for n in range(3)]
ax4.eventplot(times, colors=c, lineoffsets=(2000.-160.*np.arange(1,4)),
              linelengths=[150,150,150])

ttimes = np.arange(ret.shape[1])*0.0025-10.
plt.plot(ttimes,160.*3./2.*pulse(ttimes+10)+2000.-4*160,color='gray')

add_scalebar(ax4,sizex=10,sizey=2e6/(scale4),matchx=False,matchy=False,labelx='10 ms',
        labely='2 mV',loc=8,bbox_to_anchor=(220.5,20.0))

plt.subplot(2,5,10,sharey=ax1)
csd = -np.diff(np.diff(ret,axis=0),axis=0)
m = np.max(np.abs(csd[80:,3000:]))
plt.imshow(csd,vmin=-m,vmax=m,interpolation='nearest',origin='upper',aspect='auto',extent=(-10,10,-12000,10000))
#plt.xlim(-2.4,3.6)
plt.xlim(-10,10)
plt.axis('off')

plt.ylim(-2000,2100)

plt.figtext(0.10, 0.95, 'A', fontsize=subfig_fontsize)
plt.figtext(0.30, 0.95, 'B', fontsize=subfig_fontsize)
plt.figtext(0.50, 0.95, 'C', fontsize=subfig_fontsize)
plt.figtext(0.70, 0.95, 'D', fontsize=subfig_fontsize)

plt.figtext(0.14, 0.48, 'E', fontsize=subfig_fontsize)
plt.figtext(0.30, 0.48, 'F', fontsize=subfig_fontsize)
plt.figtext(0.45, 0.48, 'G', fontsize=subfig_fontsize)
plt.figtext(0.61, 0.48, 'H', fontsize=subfig_fontsize)
plt.figtext(0.76, 0.48, 'I', fontsize=subfig_fontsize)
plt.gcf().set_size_inches(4.48,4.45)
plt.savefig('figs/simple_axon_comparison_unfiltered.pdf')
plt.savefig('figs/simple_axon_comparison_unfiltered.png')
plt.show()
