__author__ = 'mccolgan'


class Experiment:
    def __init__(self):
        self.populations = []
        self.electrodes = []

    def add_electrode(self, electrode):
        self.electrodes.append(electrode)

    def add_electrodes(self, electrodes):
        for e in electrodes.get_electrode_list():
            self.add_electrode(e)

    def add_population(self, population):
        self.populations.append(population)

    def run(self):
        from neuron import h
        import neuron
        import numpy as np

        h.celsius = 40.
        h.dt = 0.0025
        h.finitialize(-75)
        neuron.init()
        neuron.run(5)
        imem = []
        iloc = []
        for p in self.populations:
            imem_n, iloc_n = p.nodes_imem_loc()
            imem.append(imem_n)
            iloc.append(iloc_n)

        imem, iloc = np.hstack(imem), np.hstack(iloc)
        for e in self.electrodes:
            e.calc_fields(imem, iloc)

