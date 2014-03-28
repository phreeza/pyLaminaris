__author__ = 'mccolgan'

import numpy as np
from pyplot import matplotlib as plt


def run():
    saved_data = np.load('data/simple_binaural_click_fields')
    potentials = saved_data['recorded_potentials']
    plt.plot(potentials)
    for fmt in ['.png', '.pdf']:
        plt.savefig('figs/simple_binaural_click' + fmt)
