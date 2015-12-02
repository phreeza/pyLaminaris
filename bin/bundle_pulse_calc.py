import pyLaminaris.parallel as p
import sys

params = json.load(open(sys.argv[-1]))

exp = p.ParallelBundleExperiment(**params)

exp.setup()
exp.run()
