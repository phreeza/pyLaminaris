
def run(fname):
    import pyLaminaris.parallel as p
    import json
    params = json.load(open(fname))

    exp = p.ParallelBundleExperiment(**params)

    exp.setup()
    exp.run()

if __name__ == '__main__':
    import sys
    run(sys.argv[-1])
