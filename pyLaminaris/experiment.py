__author__ = 'mccolgan'


class Experiment:
    def __init__(self):
        self.populations = []
        self.neurons = []
        self.electrodes = []

    def add_electrode(self, electrode):
        self.electrodes.append(electrode)

    def add_electrodes(self, electrodes):
        for e in electrodes.get_electrode_list():
            self.add_electrode(e)

    def add_population(self, population):
        self.populations.append(population)

    def add_neuron(self, neuron):
        self.neurons.append(neuron)

    def run(self, t=20, mode='batch'):
        from neuron import h
        import neuron
        import numpy as np

        h.celsius = 40.
        h.dt = 0.0025
        h.finitialize(-75)
        neuron.init()

        assert mode in {'batch', 'step'}
        if mode == 'batch':
            for p in self.populations:
                assert p.record
            neuron.run(t)
            imem = []
            iloc = []
            for p in self.populations:
                imem_n, iloc_n = p.nodes_imem_loc()
                imem.append(imem_n)
                iloc.append(iloc_n)

            for n in self.neurons:
                imem_n, iloc_n = n.nodes_imem_loc()
                imem.append(imem_n)
                iloc.append(iloc_n)

            imem, iloc = np.vstack(imem), np.vstack(iloc)
            for e in self.electrodes:
                e.calc_fields(iloc, imem)

        if mode == 'step':
            from time import time

            counter = 0
            interval = 10

            t0 = time()
            ti = h.t
            print "starting stepwise run"
            for p in self.populations:
                assert not p.record
            #TODO: stepwise initialisation
            # do not record currents over time
            imem = []
            iloc = []

            while h.t < t:
                h.fadvance()
                counter += 1.
                if np.mod(counter, interval) == 0:
                    iloc = []
                    imem.append(np.array([[]]))
                    rtfactor = (h.t - ti) * 1E-3 / (time() - t0)
                    print 't = %.0f, realtime factor: %.3f' % (h.t, rtfactor)
                    t0 = time()
                    ti = h.t

                    for p in self.populations:
                        imem_n, iloc_n = p.nodes_imem_loc()
                        imem[-1] = np.append(imem[-1], imem_n)
                        iloc.append(iloc_n)

                    for n in self.neurons:
                        imem_n, iloc_n = n.nodes_imem_loc()
                        imem[-1] = np.append(imem[-1], imem_n)
                        iloc.append(iloc_n)

            imem, iloc = np.vstack(imem), np.vstack(iloc)
            for e in self.electrodes:
                e.calc_fields(iloc, imem.T)
