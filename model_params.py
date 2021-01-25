from classes import Model
import numpy as np

from bead_params import bead_params_init, bead_params_fair, bead_params_play, bead_params_all_exps, bead_params_cell_pellet, bead_params_activin_ant, bead_params_bmp4_ant, bead_params_threshold

# define models
modelN = 2
models = [Model() for i in range(modelN)]

##### before ##############################################################################

# models[0].name = 'inducer_SMAD'
# models[0].label = 'Without\nnbhd'
# models[0].index_string = 'A'
# models[0].plot_color = 'C8'
# models[0].bead_params = bead_params_fair
# models[0].threshold = 0.298905377484657
# models[0].nbhd_size = None
# models[0].inducer_scaling = 10.0 ** (1.07852964487069)
# models[0].inhibitor_scaling = 10.0 ** (0.82009652169455)
#
# models[1].name = 'inducer_SMAD_nbhd'
# models[1].label = 'With\nnbhd'
# models[1].index_string = 'B'
# models[1].plot_color = 'C1'
# models[1].bead_params = bead_params_fair
# models[1].threshold = 0.489438862449442
# models[1].nbhd_size = 2*np.floor(52.7675295356746) - 1
# models[1].inducer_scaling = 10.0 ** (-2.47545875046506)
# models[1].inhibitor_scaling = 10.0 ** 1.80154516772199

############################################################################################
##### testing ##############################################################################

# models[0].name = 'inducer_SMAD'
# models[0].label = 'Without\nnbhd'
# models[0].index_string = 'A'
# models[0].plot_color = 'C8'
# models[0].bead_params = bead_params_play
# models[0].threshold = 0.473131763864675
# models[0].nbhd_size = None
# models[0].inducer_scaling = 10.0 ** (0.313337594353507)
# models[0].inhibitor_scaling = 10.0 ** (1.66978719065307)
#
# models[1].name = 'inducer_SMAD_nbhd'
# models[1].label = 'With\nnbhd'
# models[1].index_string = 'B'
# models[1].plot_color = 'C1'
# models[1].bead_params = bead_params_play
# models[1].threshold = 0.473131763864675
# models[1].nbhd_size = 2*np.floor(55.8743571703637) - 1
# models[1].inducer_scaling = 10.0 ** (0.313337594353507)
# models[1].inhibitor_scaling = 10.0 ** 1.66978719065307

#############################################################################################
##### all exps ##############################################################################

models[0].name = 'inducer_SMAD'
models[0].label = 'Without\nnbhd'
models[0].index_string = 'A'
models[0].plot_color = 'C8'
models[0].bead_params = bead_params_all_exps
models[0].threshold = 0.232585001821855
models[0].nbhd_size = None
models[0].inducer_scaling = 10.0 ** (2.22201734127452)
models[0].inhibitor_scaling = 10.0 ** (2.83382643185952)

models[1].name = 'inducer_SMAD_nbhd'
models[1].label = 'With\nnbhd'
models[1].index_string = 'B'
models[1].plot_color = 'C1'
models[1].bead_params = bead_params_all_exps
models[1].threshold = 0.44950416545415
models[1].nbhd_size = 2*np.floor(59.433636026949) - 1
models[1].inducer_scaling = 10.0 ** (-1.58937622237702)
models[1].inhibitor_scaling = 10.0 ** 2.58202368084432

################################################################################################
##### cell pellet ##############################################################################

# models[0].name = 'inducer_SMAD'
# models[0].label = 'Without\nnbhd'
# models[0].index_string = 'A'
# models[0].plot_color = 'C8'
# models[0].bead_params = bead_params_fair
# models[0].threshold = 0.298905377484657
# models[0].nbhd_size = None
# models[0].inducer_scaling = 10.0 ** (1.07852964487069)
# models[0].inhibitor_scaling = 10.0 ** (0.82009652169455)
#
# models[1].name = 'inducer_SMAD_nbhd'
# models[1].label = 'With\nnbhd'
# models[1].index_string = 'B'
# models[1].plot_color = 'C1'
# models[1].bead_params = bead_params_fair
# models[1].threshold = 0.489438862449442
# models[1].nbhd_size = 2*np.floor(52.7675295356746) - 1
# models[1].inducer_scaling = 10.0 ** (-2.47545875046506)
# models[1].inhibitor_scaling = 10.0 ** 1.80154516772199

################################################################################################
##### activin ant ##############################################################################

# models[0].name = 'inducer_SMAD'
# models[0].label = 'Without\nnbhd'
# models[0].index_string = 'A'
# models[0].plot_color = 'C8'
# models[0].bead_params = bead_params_fair
# models[0].threshold = 0.298905377484657
# models[0].nbhd_size = None
# models[0].inducer_scaling = 10.0 ** (1.07852964487069)
# models[0].inhibitor_scaling = 10.0 ** (0.82009652169455)
#
# models[1].name = 'inducer_SMAD_nbhd'
# models[1].label = 'With\nnbhd'
# models[1].index_string = 'B'
# models[1].plot_color = 'C1'
# models[1].bead_params = bead_params_fair
# models[1].threshold = 0.489438862449442
# models[1].nbhd_size = 2*np.floor(52.7675295356746) - 1
# models[1].inducer_scaling = 10.0 ** (-2.47545875046506)
# models[1].inhibitor_scaling = 10.0 ** 1.80154516772199

################################################################################################
##### bmp4 ant #################################################################################

# models[0].name = 'inducer_SMAD'
# models[0].label = 'Without\nnbhd'
# models[0].index_string = 'A'
# models[0].plot_color = 'C8'
# models[0].bead_params = bead_params_fair
# models[0].threshold = 0.298905377484657
# models[0].nbhd_size = None
# models[0].inducer_scaling = 10.0 ** (1.07852964487069)
# models[0].inhibitor_scaling = 10.0 ** (0.82009652169455)
#
# models[1].name = 'inducer_SMAD_nbhd'
# models[1].label = 'With\nnbhd'
# models[1].index_string = 'B'
# models[1].plot_color = 'C1'
# models[1].bead_params = bead_params_fair
# models[1].threshold = 0.489438862449442
# models[1].nbhd_size = 2*np.floor(52.7675295356746) - 1
# models[1].inducer_scaling = 10.0 ** (-2.47545875046506)
# models[1].inhibitor_scaling = 10.0 ** 1.80154516772199

################################################################################################
##### threshold ################################################################################

# models[0].name = 'inducer_SMAD'
# models[0].label = 'Without\nnbhd'
# models[0].index_string = 'A'
# models[0].plot_color = 'C8'
# models[0].bead_params = bead_params_fair
# models[0].threshold = 0.298905377484657
# models[0].nbhd_size = None
# models[0].inducer_scaling = 10.0 ** (1.07852964487069)
# models[0].inhibitor_scaling = 10.0 ** (0.82009652169455)
#
# models[1].name = 'inducer_SMAD_nbhd'
# models[1].label = 'With\nnbhd'
# models[1].index_string = 'B'
# models[1].plot_color = 'C1'
# models[1].bead_params = bead_params_fair
# models[1].threshold = 0.489438862449442
# models[1].nbhd_size = 2*np.floor(52.7675295356746) - 1
# models[1].inducer_scaling = 10.0 ** (-2.47545875046506)
# models[1].inhibitor_scaling = 10.0 ** 1.80154516772199

################################################################################################