# IPython log file
import numpy as np
from pyLaminaris import helper
pots = np.load('/Users/tom/Downloads/pots_3.npz')['pots']
pots = pots-pots[:,:,100:101]
def pulse(t):
    return (0.1
            + 0.5 * np.exp(-1.* (t - 10.) ** 2))
ret = np.zeros((pots.shape[0],int(20./0.0025)))
for rep in range(3000):
    print rep
    for n_neuron in range(100):
        for t in helper.inhom_poisson(pulse,20.,0.,6.):
            t_ind = int(t/0.0025)
            l = min(600,ret.shape[1]-t_ind)
            if l>0: ret[:,t_ind:t_ind+l] += pots[:,n_neuron,500:500+l]
