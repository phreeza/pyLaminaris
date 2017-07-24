import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import pyLaminaris
from brian2.units import uvolt, get_dimensions
pop = pyLaminaris.populations.SimpleNeuronPopulation(mtype='pop_sync',size=100)
pop.set_stimulation()
expt = pyLaminaris.experiment.Experiment()
expt.add_population(pop)
for n in range(220):
    expt.add_electrode(pyLaminaris.recording.SingleUnitElectrode(location=np.array([n*100,150.,0])))
expt.run(t=5)
pots = np.array([e.recorded_potential/uvolt for e in expt.electrodes])
csd = np.array([e.csd for e in expt.electrodes])
n_nodes = np.array([e.n_nodes for e in expt.electrodes])
uid = str(int(np.random.random()*1000))
assert get_dimensions(pots).is_dimensionless
np.savez('pots_pop_'+uid,pots=pots,csd=csd,n_nodes=n_nodes)
np.savez('pots_pop_sum_'+uid,pots=pots.sum(axis=1),csd=csd.sum(axis=1))
