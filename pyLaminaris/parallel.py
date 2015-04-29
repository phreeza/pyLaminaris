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

class ParallelBundleExperiment:
    """Parallel Experiments can be run with mpirun -np 4 python ..."""

    def __init__(self, n=200):
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
        self.pop = pops.NMNeuronPopulation(size=self.mysize, record=True)
        self.electrodes = []

        for n in range(30):
            for m in range(10):
                self.electrodes.append(recording.Electrode(location=np.array([3000.+100.*n, 200.*m, 0.])))
        self.pop.set_stimulation(stimtype='pulse')
        self.exp.add_population(self.pop)
        for e in self.electrodes:
            self.exp.add_electrode(e)

    def run(self, t=20., mode='batch'):
        self.exp.run(t=t, mode=mode)
        pot = self.cw.gather([e.recorded_potential for e in self.electrodes], root=0)
        if self.rank == 0:
            loc = [e.location for e in self.electrodes]
            pot = np.sum(pot,axis=0)[:,:,0]
            np.savez('potential.npz', pot=pot,loc=loc)
        MPI.Finalize()
