__author__ = 'mccolgan'


def task_simple_binaural_click():
    """Timecourses of features over long time courses"""

    from bin import simple_binaural_click_calc as c
    from bin import simple_binaural_click_plot as p
    import glob

    yield {
        'name': 'calc',
        'actions': [(c.run,)],
        'targets': ['data/simple_binaural_click_fields.npz'],
        'file_dep': ["bin/simple_binaural_click_calc.py"]
    }

    yield {
        'name': 'plot',
        'actions': [(p.run,)],
        'targets': ["figs/simple_binaural_click.png", "figs/simple_binaural_click.pdf"],
        'file_dep': ["bin/simple_binaural_click_plot.py", "data/simple_binaural_click_fields.npz"]
    }
