__author__ = 'mccolgan'
import numpy as np

from brian2.units import siemens,meter,um,mamp,mvolt
from brian2.units import dot,get_dimensions


class Electrode:
    def __init__(self, location, node_locs=None):
        self.location = location*um
        self.recorded_potential = np.array([])
        if node_locs is not None:
            self.build_dist_coeffs(node_locs)
        else:
            self.dist_coeffs = None

    def calc_fields(self,
                    node_locs, node_imem,
                    conductivity=0.33*siemens/meter):
        node_imem, node_locs = np.vstack(node_imem/mamp)*mamp, np.vstack(node_locs/um)*um
        if self.dist_coeffs is None:
            self.build_dist_coeffs(node_locs)

        ret = dot(node_imem.T, self.dist_coeffs) / (conductivity)
        self.recorded_potential = ret

    def calc_csd(self,
                 node_locs, node_imem,
                 csd_range=75.*um):
        node_imem, node_locs = np.vstack(node_imem/mamp)*mamp, np.vstack(node_locs/um)*um
        self.csd = np.sum(node_imem[np.abs(self.location - node_locs)[:, 0] < csd_range, :], axis=0) / (csd_range * 2.)


    def build_dist_coeffs(self, node_locs):
        self.dist_coeffs = (
            1. / (
                4. * np.pi * np.sqrt(np.sum((node_locs - self.location) ** 2, axis=1))
            )
        )

        # self.dist_coeffs = np.zeros(node_locs.shape)
        # for j, xl in enumerate(node_locs):
        # self.dist_coeffs[j] = (
        # 1. / (
        #             4. * np.pi * np.sqrt(np.sum((xl - self.location) * (xl - self.location)))
        #         )
        #     )


class SingleUnitElectrode:
    """ Magical Electrode that records potentials from every unit separately. Useful for
    sped up simulations in which only a single spike per fiber is simulated."""

    def __init__(self, location, node_locs=None):
        self.location = location*um
        self.recorded_potential = np.array([])
        if node_locs is not None:
            self.build_dist_coeffs(node_locs)
        else:
            self.dist_coeffs = None

    def calc_fields(self,
                    node_locs, node_imem,
                    conductivity=0.33*siemens/meter):
        self.recorded_potential = []
        if self.dist_coeffs is None:
            self.build_dist_coeffs(node_locs)
        for imem, dcoef in zip(node_imem, self.dist_coeffs):
            ret = dot(imem.T, dcoef) / conductivity
            self.recorded_potential.append(ret[:, 0])

    def calc_csd(self,
                 node_locs, node_imem,
                 csd_range=50.*um):
        self.csd = []
        self.n_nodes = []
        for imem, iloc in zip(node_imem, node_locs):
            self.csd.append(np.sum(imem[np.abs(self.location - iloc)[:, 0] < csd_range, :], axis=0) / (csd_range * 2.))
            self.n_nodes.append(np.sum(np.abs(self.location - iloc)[:, 0] < csd_range, axis=0))

    def build_dist_coeffs(self, node_locs):
        self.dist_coeffs = []
        for n in node_locs:
            self.dist_coeffs.append(np.zeros(n.shape))
            for j, xl in enumerate(n):
                self.dist_coeffs[-1][j] = (
                    1. / (
                        4. * np.pi * np.sqrt(np.sum((xl - self.location) * (xl - self.location)))
                    )
                )

class IntracellularElectrode:

    def __init__(self, compartment):
        import neuron
        neuron.load_mechanisms('neuron')
        from neuron import h
        self.recorded_potential = h.Vector()
        self.recorded_potential.record(compartment(0.5)._ref_v)

