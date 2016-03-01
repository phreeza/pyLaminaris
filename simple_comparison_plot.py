import numpy as np
from matplotlib import pyplot as plt

pots = np.load('/Users/tom/Downloads/pots.npz')['pots']
pots2 = np.load('/Users/tom/Downloads/pots_2.npz')['pots']
pots3 = np.load('/Users/tom/Downloads/pots_3.npz')['pots']
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
    plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots[n,:],mode='same').T/scale+(n-90)*200.-2000,color='black')
plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots[100,:],mode='same').T/scale,color='black',lw=2)
plt.xlim(0.6,1.2)
plt.axis('off')

plt.subplot(186,sharey=ax1)
for n in range(190,210,2):
    plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots[n,:],mode='same').T/scale+(n-190)*200.-2000,color='black')
plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots[200,:],mode='same').T/scale,color='black',lw=2)
plt.xlim(1.6,2.2)
plt.axis('off')

plt.subplot(184,sharey=ax1)
for n in range(90,110,2):
    plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots2[n,:],mode='same').T/scale+(n-90)*200.-2000,color='black')
plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window,pots2[100,:],mode='same').T/scale,color='black',lw=2)
plt.xlim(0.6,1.2)
plt.axis('off')

plt.subplot(188,sharey=ax1)
for n in range(90,110,2):
    plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window2,pots3[n,:],mode='same').T/scale2+(n-90)*200.-2000,color='black')
plt.plot(np.arange(pots.shape[1])*0.0025-10.,np.convolve(window2,pots3[100,:],mode='same').T/scale2,color='black',lw=2)
plt.xlim(0.6,1.2)
plt.axis('off')

order = [2,0,1,3]
c = ['r','g','b']
for n_p in range(4):
    ax = plt.subplot(1,8,2*order[n_p]+1,sharey=ax1)
    if n_p == 0:
        plt.plot([0,0],[-2000,0],color='black',lw=2)
    elif n_p == 1:
        plt.plot([0,0],[-2000,2000],color='black',lw=2)
    elif n_p == 2:
        plt.plot([0,0],[-2000,0],color='black',lw=2)
        plt.plot([10,10],[2000,0],color='black',lw=2)
        plt.plot([-10,-10],[2000,0],color='black',lw=2)
        plt.plot([-10,10],[0,0],color='black',lw=2)
        plt.xlim(-40,40)
    elif n_p == 3:
        
        for k in range(-1,2):
            p1 = -500.+1000*np.random.random()
            p2 = 1000*np.random.random()
            p3 = 1000*np.random.random()
            plt.plot([5*k,5*k],[-2000,p1],color=c[k+1],lw=2)
            plt.plot([5*k+10,5*k+10],[p1+p2,p1],color=c[k+1],lw=2)
            plt.plot([5*k-10,5*k-10],[p1+p3,p1],color=c[k+1],lw=2)
            plt.plot([5*k-10,5*k+10],[p1,p1],color=c[k+1],lw=2)
        plt.xlim(-40,40)
    plt.axis('off')

plt.figtext(0.14, 0.92, 'A', fontsize=24)
plt.figtext(0.34, 0.92, 'B', fontsize=24)
plt.figtext(0.54, 0.92, 'C', fontsize=24)
plt.figtext(0.74, 0.92, 'D', fontsize=24)
plt.savefig('figs/simple_axon_comparison_unfiltered.pdf')
plt.savefig('figs/simple_axon_comparison_unfiltered.png')
