__author__ = 'mccolgan'
import numpy as np


class Electrode:
    def __init__(self, location, node_locs=None):
        self.location = location
        self.recorded_potential = np.array([])
        if node_locs is not None:
            self.build_dist_coeffs(node_locs)
        else:
            self.dist_coeffs = None

    def calc_fields(self,
                    node_locs, node_imem,
                    conductivity=1. / (330. * 1e4)):  #conductivity in ohm*um
        imem, iloc = np.vstack(imem), np.vstack(iloc)
        if self.dist_coeffs is None:
            self.build_dist_coeffs(node_locs)

        ret = np.dot(node_imem.T, self.dist_coeffs) / conductivity
        self.recorded_potential = ret[:,0]

    def build_dist_coeffs(self, node_locs):
        self.dist_coeffs = np.zeros(node_locs.shape)
        for j, xl in enumerate(node_locs):
            self.dist_coeffs[j] = (
                1. / (
                    4. * np.pi * np.sqrt(np.sum((xl - self.location) * (xl - self.location)))
                )
            )


class SingleUnitElectrode:
    """ Magical Electrode that records potentials from every unit separately. Useful for
    sped up simulations in which only a single spike is simulated."""
    def __init__(self, location, node_locs=None):
        self.location = location
        self.recorded_potential = np.array([])
        if node_locs is not None:
            self.build_dist_coeffs(node_locs)
        else:
            self.dist_coeffs = None

    def calc_fields(self,
                    node_locs, node_imem,
                    conductivity=1. / (330. * 1e4)):  #conductivity in ohm*um
        self.recorded_potential = []
        if self.dist_coeffs is None:
            self.build_dist_coeffs(node_locs)
        for imem,dcoef in zip(node_imem,self.dist_coeffs):
            ret = np.dot(imem.T, dcoef) / conductivity
            self.recorded_potential.append(ret[:,0])

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
