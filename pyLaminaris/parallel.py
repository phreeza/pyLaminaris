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

    def __init__(self, n=120, stimtype='pulse'):
        self.cw = MPI.COMM_WORLD
        self.nhost = int(self.cw.size)
        self.rank = int(self.cw.rank)
        self.stimtype = stimtype

        if n is None:
            n = self.nhost

        self.sizes = np.array([int(n / self.nhost) for i in range(n)])
        self.sizes[:(n % self.nhost)] += 1

        self.mysize = self.sizes[self.rank]
        print self.mysize

    def setup(self):
        self.exp = experiment.Experiment()
        self.pop = pops.NMNeuronPopulation(size=self.mysize, record=True)
        self.electrodes = []

        for n in range(200):
            for m in range(30):
                self.electrodes.append(recording.Electrode(location=np.array([10000.+100.*n,50.*m*m, 0.])))
        self.pop.set_stimulation(stimtype=self.stimtype)
        self.exp.add_population(self.pop)
        for e in self.electrodes:
            self.exp.add_electrode(e)

    def run(self, t=20., mode='batch', fname='data/parallel_bundle_potential.npz'):
        self.exp.run(t=t, mode=mode)
        print "gathering",self.rank
        nodes = self.cw.gather(self.pop.nodes_imem_loc()[1],root=0)
        pot = [self.cw.gather(e.recorded_potential,root=0) for e in self.electrodes]
        print "gathering complete",self.rank
        if self.rank == 0:
            loc = [e.location for e in self.electrodes]
            pot = np.sum(pot,axis=1)
            np.savez(fname, pot=pot,loc=loc,nodes=nodes)
        MPI.Finalize()
