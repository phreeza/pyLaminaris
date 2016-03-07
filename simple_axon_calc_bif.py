import numpy as np
import pyLaminaris
pop = pyLaminaris.populations.SimpleNeuronPopulation(mtype='bif')
pop.set_stimulation()
expt = pyLaminaris.experiment.Experiment()
expt.add_population(pop)
for n in range(220):
    expt.add_electrode(pyLaminaris.recording.Electrode(location=np.array([n*100,150.,0])))
expt.run()
pots = np.array([e.recorded_potential for e in expt.electrodes])
np.savez('pots_2',pots=pots)
