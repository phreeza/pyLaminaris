from pyLaminaris import parallel
import numpy as np
from neuron import h

e = parallel.ParallelExperiment(n=20)

e.setup()
ret = e.run(t=2.)
if ret is not None:
    print np.sum(ret, axis=0)
    print np.sum(ret, axis=0).shape
e.pc.barrier()
e.pc.done()
h.quit()
