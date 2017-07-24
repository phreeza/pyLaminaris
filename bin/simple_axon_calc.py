import numpy as np
import pyLaminaris
from brian2.units import uvolt, get_dimensions
pop = pyLaminaris.populations.SimpleNeuronPopulation(mtype='long')
pop.set_stimulation()
expt = pyLaminaris.experiment.Experiment()
expt.add_population(pop)
for n in range(220): expt.add_electrode(pyLaminaris.recording.Electrode(location=np.array([n*100,150.,0])))

root_nodes = pop.neurons[0].axon.root.nodes
intracellular_target_comp = root_nodes[len(root_nodes)/2]
intracellular_electrode = pyLaminaris.recording.IntracellularElectrode(intracellular_target_comp)
expt.add_intracellular_electrode(intracellular_electrode)

intracellular_target_comp = root_nodes[-1]
intracellular_electrode = pyLaminaris.recording.IntracellularElectrode(intracellular_target_comp)
expt.add_intracellular_electrode(intracellular_electrode)
expt.run()
pots = np.array([e.recorded_potential/uvolt for e in expt.electrodes])
csd = np.array([e.csd for e in expt.electrodes])
intracellular_pot = np.array([e.recorded_potential for e in expt.intracellular_electrodes])

assert get_dimensions(pots).is_dimensionless
np.savez('pots',pots=pots,csd=csd,intracellular_pot=intracellular_pot)
