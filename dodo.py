__author__ = 'mccolgan'
import matplotlib as mpl
mpl.use('Agg')

def task_simple_binaural_click():
    '''Timecourses of features over long time courses'''

    from bin import simple_binaural_click_calc as c
    from bin import simple_binaural_click_plot as p

    yield {
        'name': 'calc',
        'actions': [(c.run,)],
        'targets': ['data/simple_binaural_click_fields.npz'],
        'file_dep': ['bin/simple_binaural_click_calc.py']
    }

    yield {
        'name': 'plot',
        'actions': [(p.run,)],
        'targets': ['figs/simple_binaural_click.png', 'figs/simple_binaural_click.pdf'],
        'file_dep': ['bin/simple_binaural_click_plot.py', 'data/simple_binaural_click_fields.npz']
    }

def task_parallel_bundle():
    '''Generic bundle of axons for figures 1&2'''
    from bin import bundle_pulse_calc as c
    from bin import bundle_pulse_plot as p
    import glob
    import json
    for fname in glob.glob('params/*.params'):
        params = json.load(open(fname))
        postfix = params['postfix'] 
        if "n_jobs" in params.keys():
            n_jobs = params['n_jobs']
        else:
            n_jobs = 10
        yield {
            'name': 'calc_'+postfix,
            'actions': ["mpirun -n "+str(n_jobs)+' python bin/bundle_pulse_calc.py '+fname],
            'targets': ['data/bundle_pulse_'+postfix+'.npz'],
            'file_dep': ['bin/bundle_pulse_calc.py']
        }

        yield {
            'name': 'plot_'+postfix,
            'actions': [(p.run,[fname])],
            'targets': ['figs/bundle_pulse_potentials_'+postfix+'.png', 
                        'figs/bundle_pulse_potentials_'+postfix+'.pdf'],
            'file_dep': ['bin/bundle_pulse_plot.py', 'data/bundle_pulse_'+postfix+'.npz']
        }
