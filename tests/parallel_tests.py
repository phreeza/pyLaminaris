__author__ = 'mccolgan'


def test_stepwise_run():
    import pyLaminaris.parallel
    import numpy as np

    parallel_expt = pyLaminaris.parallel.ParallelExperiment(n=1)
    parallel_expt.setup()
    data = np.sum(parallel_expt.run(t=1., mode='step'), axis=0)
    print data
    print data.shape
