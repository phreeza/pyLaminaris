__author__ = 'mccolgan'


def run():
    import pyLaminaris.populations as pops
    import pyLaminaris.experiment
    import pyLaminaris.recording
    import numpy as np

    exp = pyLaminaris.experiment.Experiment()
    ipsi_pop = pops.NMNeuronPopulation(side='ipsi', size=200)
    contra_pop = pops.NMNeuronPopulation(side='contra', size=200)
    electrode = pyLaminaris.recording.Electrode(location=np.array([4500., 0., 0.]))
    ipsi_pop.set_stimulation(stimtype='click', freq=4000.)
    contra_pop.set_stimulation(stimtype='click', freq=4000.)

    exp.add_population(ipsi_pop)
    exp.add_population(contra_pop)
    exp.add_electrode(electrode)

    exp.run(t=20)

    np.savez('data/simple_binaural_click_fields', recorded_potentials=electrode.recorded_potential)