from pyLaminaris import parallel
import numpy as np

e = parallel.ParallelExperiment(n=None)
e.setup()

ret = e.run(t=2.)
if ret is not None:
    print np.sum(ret, axis=0)
    print np.sum(ret, axis=0).shape
