import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import pyLaminaris
pop = pyLaminaris.populations.SimpleNeuronPopulation(mtype='pop_sync',size=100)
pop.set_stimulation()
expt = pyLaminaris.experiment.Experiment()
expt.add_population(pop)
for n in range(220):
    expt.add_electrode(pyLaminaris.recording.SingleUnitElectrode(location=np.array([n*100,150.,0])))
expt.run(t=5)
pots = np.array([e.recorded_potential for e in expt.electrodes])
np.savez('pots_3',pots=pots)
