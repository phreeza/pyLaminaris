import numpy as np
import pyLaminaris
from brian2.units import uvolt, get_dimensions
pop = pyLaminaris.populations.SimpleNeuronPopulation(mtype='long')
pop.set_stimulation()
expt = pyLaminaris.experiment.Experiment()
expt.add_population(pop)
for n in range(220): expt.add_electrode(pyLaminaris.recording.Electrode(location=np.array([n*100,150.,0])))
expt.run()
pots = np.array([e.recorded_potential/uvolt for e in expt.electrodes])
csd = np.array([e.csd for e in expt.electrodes])
assert get_dimensions(pots).is_dimensionless
np.savez('pots',pots=pots,csd=csd)
