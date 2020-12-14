from classes import Model
import numpy as np

from bead_params import bead_params_init, bead_params_fair

# define models
modelN = 2
models = [Model() for i in range(modelN)]

models[0].name = 'inducer_SMAD'
models[0].label = 'Without\nnbhd'
models[0].index_string = 'A'
models[0].plot_color = 'C8'
models[0].bead_params = bead_params_fair
models[0].threshold = 0.298905377484657
models[0].nbhd_size = None
models[0].inducer_scaling = 10.0 ** (1.07852964487069)
models[0].inhibitor_scaling = 10.0 ** (0.82009652169455)

models[1].name = 'inducer_SMAD_nbhd'
models[1].label = 'With\nnbhd'
models[1].index_string = 'B'
models[1].plot_color = 'C1'
models[1].bead_params = bead_params_fair
models[1].threshold = 0.489438862449442
models[1].nbhd_size = 2*np.floor(52.7675295356746) - 1
models[1].inducer_scaling = 10.0 ** (-2.47545875046506)
models[1].inhibitor_scaling = 10.0 ** 1.80154516772199


