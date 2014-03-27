__author__ = 'mccolgan'
import neurons
import numpy as np


class NMNeuronPopulation:
    def __init__(self, size=100, side='ipsi'):
        self.side = side
        self.size = size
        x_offset = 4500.
        x_spread = 200.
        self.neurons = [neurons.NMNeuron(root_point=x_offset + x_spread * np.random.rand()) for n in range(size)]

    def set_stimulation(self, stimtype='mod_click', freq=4000.):
        #TODO: make this take more than a single stimtype
        import helper

        def mod_click(t):
            return (0.1
                    + 2.9 * np.exp(np.sin(4 * 2 * np.pi * t) - 1. - 16 * (t - 10) ** 2)
                    + 1.5 * np.exp(np.sin(4 * 2 * np.pi * t) - 1. - 16 * (t - 11) ** 2)
                    - 0.1 * np.exp(-16 * (t - 10.5) ** 2))

        for n in self.neurons:
            n.set_spiketimes(helper.inhom_poisson(mod_click, 20., 0., 6.))
        pass