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
        #TODO make node_locs argument optional
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


class ElectrodeArray:
    def __init__(self):
        pass

    def get_electrode_list(self):
        return []
