def test_axons_instances():
    import pyLaminaris.axons

    s = pyLaminaris.axons.Segment()
    t = pyLaminaris.axons.Tree()
    p = pyLaminaris.axons.ProbTree()


def test_populations_instances():
    import pyLaminaris.populations as pops

    mypop = pops.NMNeuronPopulation(size=2)


def test_experiment_setup():
    import pyLaminaris.populations as pops
    import pyLaminaris.experiment
    import pyLaminaris.recording

    exp = pyLaminaris.experiment.Experiment()
    ipsi_pop = pops.NMNeuronPopulation(side='ipsi', size=5)
    contra_pop = pops.NMNeuronPopulation(side='contra', size=5)
    electrodes = pyLaminaris.recording.ElectrodeArray()
    ipsi_pop.set_stimulation(stimtype='click', freq=4000.)
    contra_pop.set_stimulation(stimtype='click', freq=4000.)
    exp.add_population(ipsi_pop)
    exp.add_population(contra_pop)
    exp.add_electrodes(electrodes)
    exp.run()
