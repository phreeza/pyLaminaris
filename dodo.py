__author__ = 'mccolgan'


def task_simple_binaural_click():
    """Timecourses of features over long time courses"""

    from bin import simple_binaural_click_calc as c
    from bin import simple_binaural_click_plot as p
    import glob

    yield {
        'name': 'calc',
        'actions': [(c.run,)],
        'targets': [""],
        'file_dep': ["script/psth_isih_from_txt.py"]
    }

    yield {
        'name': 'plot',
        'actions': [(p.run,)],
        'targets': [""],
        'file_dep': ["script/psth_isih_from_txt.py"]
    }
