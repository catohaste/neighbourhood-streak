from classes import Model
import numpy as np

from bead_params import bead_params_init, bead_params_A, bead_params_B, bead_params_fair

# define models
modelN = 2
models = [Model() for i in range(modelN)]

models[0].name = 'inducer_SMAD_nbhd'
models[0].label = 'With\nnbhd'
models[0].index_string = 'A'
models[0].plot_color = 'C1'
models[0].bead_params = bead_params_fair
models[0].threshold = 0.498959407670072
models[0].nbhd_size = 2*np.floor(70) - 1
models[0].inducer_scaling = 10.0 ** (-0.369861982125379)
models[0].inhibitor_scaling = 10.0 ** 2.03049314837895


models[1].name = 'inducer_SMAD'
models[1].label = 'Without\nnbhd'
models[1].index_string = 'B'
models[1].plot_color = 'C8'
models[1].bead_params = bead_params_fair
models[1].threshold = 0.283033459349818
models[1].nbhd_size = None
models[1].inducer_scaling = 10.0 ** (0.522844014349102)
models[1].inhibitor_scaling = 10.0 ** (0.83050960833543)

# models[2].name = 'inducer_SMAD_nbhd'
# models[2].index_string = 'C'
# models[2].plot_color = 'C8'
# models[2].bead_params = bead_params_A
# models[2].threshold = 0.1