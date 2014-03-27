__author__ = 'mccolgan'

import axons
import numpy as np
from neuron import h


class NMNeuron:
    #TODO: fully parametrize root point and growth direction
    def __init__(self, root_point):
        self.axon = axons.ProbTree(root_point=root_point,
                                   mode='mod_click', structure='logistic')
        self.spiketimes = []
        self.stims = []

    def set_spiketimes(self, spiketimes):
        self.spiketimes = spiketimes
        for d in self.spiketimes:
            stim = h.IClamp(self.axon.root.nodes[0](0.5))
            stim.delay = d
            stim.amp = 1.
            stim.dur = 0.1
            self.stims.append(stim)