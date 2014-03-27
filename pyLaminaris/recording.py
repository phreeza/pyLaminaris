__author__ = 'mccolgan'
import numpy as np


class Electrode:
    def __init__(self, location):
        self.location = location
        self.recorded_potential = np.array()


    def calc_kernels_3d(self,
                        node_locs, node_imem,
                        conductivity=1. / (330. * 1e4)):  #conductivity in ohm*um
        ret = np.zeros(node_imem.shape[1])
        for j, xl in enumerate(node_locs):
            ret[:] += (
                node_imem[j, :] / (
                    4. * np.pi * conductivity * np.sqrt(np.sum((xl - self.location) * (xl - self.location)))
                )
            )
        self.recorded_potential = ret


class ElectrodeArray:
    def __init__(self):
        pass

    def get_electrode_list(self):
        return []
