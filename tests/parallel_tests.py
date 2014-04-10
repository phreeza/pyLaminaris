__author__ = 'mccolgan'


def test_stepwise_run():
    import pyLaminaris.parallel

    parallel_expt = pyLaminaris.parallel.ParallelExperiment(n=200)
    data = parallel_expt.run(t=1., mode='step')
    print data.shape
