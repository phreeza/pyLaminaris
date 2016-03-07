__author__ = 'mccolgan'

import axons
import numpy as np
from neuron import h


class NMNeuron:
    #TODO: fully parametrize root point and growth direction
    def __init__(self, root_point, structure='logistic', record=True, **params):
        self.record = record
        self.axon = axons.ProbTree(root_point=root_point,
                                   mode='mod_click', structure=structure, record=self.record, **params)
        self.spiketimes = []
        self.stims = []

    def set_spiketimes(self, spiketimes):
        self.spiketimes = spiketimes
        for d in self.spiketimes:
            stim = h.IClamp(self.axon.root.nodes[0](0.5))
            stim.delay = d
            stim.amp = 0.5
            stim.dur = 0.1
            self.stims.append(stim)

    def nodes_imem_loc(self):
        return self.axon.nodes_imem_loc()

class SimpleNeuron(NMNeuron):
    def __init__(self, mtype, record=True, **params):
        self.record = record
        self.mtype = mtype
        if mtype == 'bif':
            self.axon = axons.Tree(depth=1,root_point=np.array([10000.,0.,0.]), record=self.record, **params)
        elif mtype == 'long':
            self.axon = axons.Tree(depth=0,root_point=np.array([20000.,0.,0.]), record=self.record, **params)
        elif mtype == 'pop_sync':
            self.axon = axons.Tree(depth=1,root_point=np.array([9500.+1000.*np.random.random(),0.,0.]), record=self.record,mtype=mtype, **params)
        self.spiketimes = []
        self.stims = []
