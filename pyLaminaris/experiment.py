__author__ = 'mccolgan'


class Experiment:
    def __init__(self):
        pass

    def add_electrode(self, electrode):
        pass

    def add_electrodes(self, electrodes):
        for e in electrodes.get_electrode_list():
            self.add_electrode(e)

    def add_population(self, population):
        pass

    def run(self):
        pass
