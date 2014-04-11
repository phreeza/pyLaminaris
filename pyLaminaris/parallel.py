__author__ = 'mccolgan'

import numpy as np
import populations as pops
import experiment
import recording
from mpi4py import MPI
from neuron import h


class ParallelExperiment:
    def __init__(self, n=2):
        self.pc = h.ParallelContext()
        self.nhost = int(self.pc.nhost())
        self.rank = int(self.pc.id())

        s = "mpi4py thinks I am %d of %d, NEURON thinks I am %d of %d\n"
        self.cw = MPI.COMM_WORLD
        print s % (self.cw.rank, self.cw.size,
                   self.pc.id(), self.pc.nhost())

        print "Minion %i of %i reporting for duty" % (self.rank, self.nhost)

        self.sizes = np.array([int(n / self.nhost) for i in range(n)])
        self.sizes[:(n % self.nhost)] += 1

        self.mysize = self.sizes[self.rank]

    def setup(self):
        self.exp = experiment.Experiment()
        self.ipsi_pop = pops.NMNeuronPopulation(side='ipsi', size=self.mysize, record=False)
        self.contra_pop = pops.NMNeuronPopulation(side='contra', size=self.mysize, record=False)
        self.electrode = recording.Electrode(location=np.array([0., 0., 0.]))
        self.ipsi_pop.set_stimulation(stimtype='click', freq=4000.)
        self.contra_pop.set_stimulation(stimtype='click', freq=4000.)
        self.exp.add_population(self.ipsi_pop)
        self.exp.add_population(self.contra_pop)
        self.exp.add_electrode(self.electrode)

    def run(self, t=20., mode='step'):
        self.exp.run(t=t, mode=mode)
        pot = self.cw.gather(self.electrode.recorded_potential, root=0)
        if self.rank == 0:
            return pot
