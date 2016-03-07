import numpy as np
from matplotlib import pyplot as plt
from scalebars import add_scalebar

pots = np.load('/Users/tom/Downloads/pots.npz')['pots']
pots2 = np.load('/Users/tom/Downloads/pots_2.npz')['pots']
pots3 = np.load('/Users/tom/Downloads/pots_3.npz')['pots'].sum(axis=1)
scale = 300.
scale2 = scale*100
window = np.exp(-(np.arange(1000)-500.)**2/2.**2)
window /= window.sum()
window = window*5

window2 = np.exp(-(np.arange(1000)-500.)**2/2.**2)
window2 /= window2.sum()
window2 = window2*5

ax1 = plt.subplot(182)
for n in range(90,110,2):
    plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots[n,:],mode='same').T/scale-(n-90)*200.+2000,color='black')
plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots[100,:],mode='same').T/scale,color='black',lw=2)
plt.xlim(0.6,1.2)
plt.axis('off')
add_scalebar(ax1,sizex=.2,sizey=scale,matchx=False,matchy=False,labelx='0.2 ms',
        labely='1 $\mu$V',loc=1,borderpad=0.0,bbox_to_anchor=(100.,100.))

plt.subplot(186,sharey=ax1)
for n in range(190,210,2):
    plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots[n,:],mode='same').T/scale-(n-190)*200.+2000,color='black')
plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots[200,:],mode='same').T/scale,color='black',lw=2)
plt.xlim(1.6,2.2)
plt.axis('off')

plt.subplot(184,sharey=ax1)
for n in range(90,110,2):
    plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots2[n,:],mode='same').T/scale-(n-90)*200.+2000,color='black')
plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots2[100,:],mode='same').T/scale,color='black',lw=2)
plt.xlim(0.6,1.2)
plt.axis('off')

plt.subplot(188,sharey=ax1)
for n in range(90,110,2):
    plt.plot(np.arange(pots3.shape[1])*0.0025,np.convolve(window2,pots3[n,:],mode='same').T/scale2-(n-90)*200.+2000,color='black')
plt.plot(np.arange(pots3.shape[1])*0.0025,np.convolve(window2,pots3[100,:],mode='same').T/scale2,color='black',lw=2)
plt.xlim(1.6,2.2)
plt.axis('off')

order = [2,0,1,3]
c = ['r','g','b']
d_sep = 20
for n_p in range(4):
    ax = plt.subplot(1,8,2*order[n_p]+1,sharey=ax1)
    plt.arrow(-15,1800,0,-800,width=1,head_width=8,facecolor='black',head_length=100)
    if n_p == 0:
        plt.plot([0,0],[2000,0],color='black',lw=2)
        plt.plot([0],[0],color='black',lw=2,marker='o',mec='none')
        plt.plot([-5,5],[1845,1835],color='white',lw=2)
        plt.plot([-5,5],[1870,1860],color='black',lw=2)
        plt.plot([-5,5],[1820,1810],color='black',lw=2)
    elif n_p == 1:
        plt.plot([0,0],[-2000,2000],color='black',lw=2)
        plt.plot([-5,5],[1845,1835],color='white',lw=2)
        plt.plot([-5,5],[1870,1860],color='black',lw=2)
        plt.plot([-5,5],[1820,1810],color='black',lw=2)
        plt.plot([-5,5],[-1835,-1845],color='white',lw=2)
        plt.plot([-5,5],[-1860,-1870],color='black',lw=2)
        plt.plot([-5,5],[-1810,-1820],color='black',lw=2)
    elif n_p == 2:
        plt.plot([0,0],[2000,0],color='black',lw=2)
        plt.plot([10,10],[-2000,0],color='black',lw=2)
        plt.plot([-10,-10],[-2000,0],color='black',lw=2)
        plt.plot([-10,10],[0,0],color='black',lw=2)

        plt.plot([-5,5],[1845,1835],color='white',lw=2)
        plt.plot([-5,5],[1870,1860],color='black',lw=2)
        plt.plot([-5,5],[1820,1810],color='black',lw=2)

        plt.plot([-15,-5],[-1835,-1845],color='white',lw=2)
        plt.plot([-15,-5],[-1860,-1870],color='black',lw=2)
        plt.plot([-15,-5],[-1810,-1820],color='black',lw=2)
        plt.plot([5,15],[-1835,-1845],color='white',lw=2)
        plt.plot([5,15],[-1860,-1870],color='black',lw=2)
        plt.plot([5,15],[-1810,-1820],color='black',lw=2)
    elif n_p == 3:
        for k in range(-1,2):
            p1 = 1000*np.random.random()-500.
            p2 = 1000*np.random.random()
            p3 = 1000*np.random.random()
            plt.plot([5*k,5*k],[2000,p1],color=c[k+1],lw=2)
            plt.plot([5*k+10,5*k+10],[p1-p2,p1],color=c[k+1],lw=2)
            plt.plot([5*k-10,5*k-10],[p1-p3,p1],color=c[k+1],lw=2)
            plt.plot([5*k-10,5*k+10],[p1,p1],color=c[k+1],lw=2)
            plt.plot([5*k+10],[p1-p2],marker='o',color=c[k+1],lw=2,mec='none')
            plt.plot([5*k-10],[p1-p3],marker='o',color=c[k+1],lw=2,mec='none')
        plt.plot([-10,10],[1855,1825],color='white',lw=2)
        plt.plot([-10,10],[1880,1850],color='black',lw=2)
        plt.plot([-10,10],[1830,1800],color='black',lw=2)
    plt.xlim(-40,40)
    plt.axis('off')

plt.figtext(0.14, 0.84, 'A', fontsize=24)
plt.figtext(0.34, 0.84, 'B', fontsize=24)
plt.figtext(0.54, 0.84, 'C', fontsize=24)
plt.figtext(0.74, 0.84, 'D', fontsize=24)
plt.savefig('figs/simple_axon_comparison_unfiltered.pdf')
plt.savefig('figs/simple_axon_comparison_unfiltered.png')
