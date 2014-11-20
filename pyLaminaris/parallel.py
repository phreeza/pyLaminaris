__author__ = 'mccolgan'

import numpy as np
import populations as pops
import experiment
import recording
from mpi4py import MPI
from neuron import h


class ParallelExperiment:
    """Parallel Experiments can be run with mpirun -np 4 python ..."""

    def __init__(self, n=2):
        self.cw = MPI.COMM_WORLD
        self.nhost = int(self.cw.size)
        self.rank = int(self.cw.rank)

        if n is None:
            n = self.nhost

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
        MPI.Finalize()
