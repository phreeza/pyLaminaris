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

    def __init__(self, population_size=120, stimtype='pulse', simulation_time=20., **params):
        self.cw = MPI.COMM_WORLD
        self.nhost = int(self.cw.size)
        self.rank = int(self.cw.rank)
        self.stimtype = stimtype
        self.params = params
        self.simulation_time = simulation_time
        if "postfix" in params.keys():
            self.fname = "data/bundle_pulse_"+params['postfix']+".npz"
        else:
            self.fname = "data/bundle_pulse.npz"
        self.electrode_params = dict(x_N=30, x_d=20, x_base=0, x_quadscale=True, y_N=300, y_d=20, y_base=10000,
                                     y_quadscale=False)
        if "electrode_params" in params.keys():
            self.electrode_params.update(params["electrode_params"])

        if population_size is None:
            population_size = self.nhost

        self.sizes = np.array([int(population_size / self.nhost) for i in range(population_size)])
        self.sizes[:(population_size % self.nhost)] += 1

        self.mysize = self.sizes[self.rank]

    def setup(self):
        self.exp = experiment.Experiment()
        if 'population_params' in self.params.keys():
            pop_params = self.params['population_params']
        else:
            pop_params = {}
        self.pop = pops.NMNeuronPopulation(size=self.mysize, record=True, **pop_params)
        self.electrodes = []

        for n in range(self.electrode_params['y_N']):
            if self.electrode_params['y_quadscale']:
                nn = n * n
            else:
                nn = n
            for m in range(self.electrode_params['x_N']):
                if self.electrode_params['x_quadscale']:
                    mm = m * m
                else:
                    mm = m
                self.electrodes.append(recording.Electrode(location=np.array([
                    self.electrode_params['y_base'] + self.electrode_params['y_d'] * nn,
                    self.electrode_params['x_base'] + self.electrode_params['x_d'] * mm,
                    0.])))
        self.pop.set_stimulation(stimtype=self.stimtype)
        self.exp.add_population(self.pop)
        for e in self.electrodes:
            self.exp.add_electrode(e)

    def run(self, mode='batch'):
        self.exp.run(t=self.simulation_time, mode=mode)
        print "gathering", self.rank
        nodes = self.cw.gather(self.pop.get_node_locations(), root=0)
        segments = self.cw.gather(self.pop.get_segment_locations(), root=0)
        pot = [self.cw.gather(e.recorded_potential, root=0) for e in self.electrodes]
        print "gathering complete", self.rank
        if self.rank == 0:
            loc = [e.location for e in self.electrodes]
            pot = np.sum(pot, axis=1)
            np.savez(self.fname, pot=pot, loc=loc, nodes=nodes, segments=segments)
        MPI.Finalize()

