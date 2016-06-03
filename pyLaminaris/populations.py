__author__ = 'mccolgan'
import neurons as _neurons
import numpy as np

class NMNeuronPopulation:
    def __init__(self, size=100, side='ipsi', record=True, root_x_offset=14500., root_x_spread=200., root_y_offset=-10., y_spread=20.,
                 structure='logistic', **params):
        self.side = side
        self.size = size

        self.x_offset = root_x_offset
        self.x_spread = root_x_spread
        self.y_offset = root_y_offset
        self.y_spread = y_spread
        self.record = record
        self.neurons = [
            _neurons.NMNeuron(
                root_point=[
                    self.x_offset + self.x_spread * np.random.rand(),
                    self.y_offset + self.y_spread * np.random.rand(),
                    0.0
                ], record=self.record, structure=structure, **params) for n in
            range(size)]

    def set_stimulation(self, stimtype='mod_click', freq=4000.):
        import helper

        def mod_click(t):
            return (0.1
                    + 2.9 * np.exp(np.sin(4 * 2 * np.pi * t) - 1. - 16 * (t - 10) ** 2)
                    + 1.5 * np.exp(np.sin(4 * 2 * np.pi * t) - 1. - 16 * (t - 11) ** 2)
                    - 0.1 * np.exp(-16 * (t - 10.5) ** 2))

        def pulse(t):
            return (0.1
                    + 2.9 * np.exp(-16 * (t - 10.) ** 2))

        if stimtype == 'mod_click':
            fun = mod_click
        elif stimtype == 'pulse':
            fun = pulse
        elif stimtype == 'single':
            for n in self.neurons:
                n.set_spiketimes([1.])
            return
        else:
            raise NameError, "Unknown stimulation type"

        for n in self.neurons:
            n.set_spiketimes(helper.inhom_poisson(fun, 20., 0., 6.))

    def nodes_imem_loc(self):
        imem = []
        iloc = []
        for n in self.neurons:
            imem_n, iloc_n = n.nodes_imem_loc()
            imem.append(imem_n)
            iloc.append(iloc_n)
        imem = np.vstack(imem)
        iloc = np.vstack(iloc)
        if hasattr(self, 'side') and self.side == 'contra':
            iloc[:, 0] = 2 * (self.x_offset + self.x_spread) - iloc[:, 0]
        return imem, iloc

    def get_node_locations(self):
        ret = []
        for ax in self.neurons:
            ret.append([])
            for seg in ax.axon.segments():
                ret[-1].append(seg.node_locations)

        return ret

    def get_segment_locations(self):
        ret = []
        for ax in self.neurons:
            ret.append([])
            for seg in ax.axon.segments():
                l = np.vstack((seg.start, seg.end))
                ret[-1].append(l)
        return ret

    def draw_2d(self, thin=10):
        for n in self.neurons[::thin]:
            n.axon.draw_2d(alpha=0.1)

class SimpleNeuronPopulation(NMNeuronPopulation):
    def __init__(self, size=1, record=True, mtype='long', **params):
        self.size = size
        assert mtype in ['long','bif','bif_term','pop_sync']
        self.mtype = mtype
        self.record   = record
        self.neurons  = [
            _neurons.SimpleNeuron(record=self.record, mtype=self.mtype, **params) for n in
            range(size)]

    def set_stimulation(self):
        import helper

        for n in self.neurons:
            #n.set_spiketimes([10.])
            n.set_spiketimes([1.])
