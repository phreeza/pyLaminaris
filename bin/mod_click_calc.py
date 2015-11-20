import pyLaminaris.parallel as p

expt = p.ParallelBundleExperiment(n=40*100, stimtype='mod_click')
expt.setup()
expt.run(fname='data/mod_click_bundle.npz')
