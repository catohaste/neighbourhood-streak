from classes import Model

from bead_params import bead_params_init, bead_params_A, bead_params_B, bead_params_fair

# define models
modelN = 2
models = [Model() for i in range(modelN)]

models[0].name = 'inducer_SMAD_nbhd'
models[0].label = 'With\nnbhd'
models[0].index_string = 'A'
models[0].plot_color = 'C1'
models[0].bead_params = bead_params_fair
models[0].threshold = 0.497
models[0].nbhd_size = 2*70 - 1
models[0].inducer_scaling = 10.0 ** (-2.63)
models[0].inhibitor_scaling = 10.0 ** 2.60


models[1].name = 'inducer_SMAD'
models[1].label = 'Without\nnbhd'
models[1].index_string = 'B'
models[1].plot_color = 'C8'
models[1].bead_params = bead_params_fair
models[1].threshold = 0.234
models[1].nbhd_size = None
models[1].inducer_scaling = 10.0 ** (-0.110)
models[1].inhibitor_scaling = 10.0 ** (-0.435)

# models[2].name = 'inducer_SMAD_nbhd'
# models[2].index_string = 'C'
# models[2].plot_color = 'C8'
# models[2].bead_params = bead_params_A
# models[2].threshold = 0.1