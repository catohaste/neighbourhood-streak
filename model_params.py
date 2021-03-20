from classes import Model
import numpy as np

from bead_params import bead_params_play, bead_params_all_exps, bead_params_cell_pellet, bead_params_activin_ant, bead_params_bmp4_ant, bead_params_threshold


def load_models(select_exp):
    """ load models with correct parameters based upon experiment """
    
    # define models
    modelN = 2
    models = [Model() for i in range(modelN)]

    # parameters which stay the same regardless of experiment
    models[0].name = 'inducer_SMAD'
    models[0].label = 'Without\nnbhd'
    models[0].index_string = 'A'
    models[0].plot_color = 'C8'

    models[1].name = 'inducer_SMAD_nbhd'
    models[1].label = 'With\nnbhd'
    models[1].index_string = 'B'
    models[1].plot_color = 'C1'
            
    if select_exp is 'testing':
        
        models[0].bead_params = bead_params_play
        models[0].threshold = 0.497687138515397
        models[0].nbhd_size = None
        models[0].inhibitor_scaling = 10.0 ** (1.45030857599371)
        models[0].inducer_scaling = 10.0 ** (1.29838186230979)

        models[1].bead_params = bead_params_play
        models[1].threshold = 0.421573076609568
        models[1].nbhd_size = 2*np.floor(56.5325206131255) + 1
        models[1].inhibitor_scaling = 10.0 ** -0.329085126426299
        models[1].inducer_scaling = 10.0 ** (-2.61174561686766)
        
    elif select_exp is 'all_exps':
        
        models[0].bead_params = bead_params_all_exps
        models[0].threshold = 0.427922233370511
        models[0].nbhd_size = None
        models[0].inhibitor_scaling = 10.0 ** (2.28635614897028)
        models[0].inducer_scaling = 10.0 ** (2.40585923229449)

        models[1].bead_params = bead_params_all_exps
        models[1].threshold = 0.454828043213113
        models[1].nbhd_size = 2*np.floor(62.8098261951711) + 1
        models[1].inhibitor_scaling = 10.0 ** 1.57491873987476
        models[1].inducer_scaling = 10.0 ** (-1.79339961627464)
        
    elif select_exp is 'cell_pellet':
        
        models[0].bead_params = bead_params_cell_pellet
        models[0].threshold = 0.497687138515397
        models[0].nbhd_size = None
        models[0].inhibitor_scaling = 10.0 ** (1.45030857599371)
        models[0].inducer_scaling = 10.0 ** (1.29838186230979)

        models[1].bead_params = bead_params_cell_pellet
        models[1].threshold = 0.421573076609568
        models[1].nbhd_size = 2*np.floor(56.5325206131255) + 1
        models[1].inhibitor_scaling = 10.0 ** -0.329085126426299
        models[1].inducer_scaling = 10.0 ** (-2.61174561686766)
        
    elif select_exp is 'activin_ant':
        
        models[0].bead_params = bead_params_activin_ant
        models[0].threshold = 0.496121824802409
        models[0].nbhd_size = None
        models[0].inhibitor_scaling = 10.0 ** (0.826017091134685)
        models[0].inducer_scaling = 10.0 ** (1.57477776175289)

        models[1].bead_params = bead_params_activin_ant
        models[1].threshold = 0.431647005869817
        models[1].nbhd_size = 2*np.floor(130.477438526608) + 1
        models[1].inhibitor_scaling = 10.0 ** 2.22004529783025
        models[1].inducer_scaling = 10.0 ** (-1.62424826505605)
        
    elif select_exp is 'bmp4_ant':
        
        models[0].bead_params = bead_params_bmp4_ant
        models[0].threshold = 0.485598910878967
        models[0].nbhd_size = None
        models[0].inhibitor_scaling = 10.0 ** (2.22676213806672)
        models[0].inducer_scaling = 10.0 ** (2.67031533897687)

        models[1].bead_params = bead_params_bmp4_ant
        models[1].threshold = 0.245747297043042
        models[1].nbhd_size = 2*np.floor(33.4112433863738) + 1
        models[1].inhibitor_scaling = 10.0 ** 2.39830585243343
        models[1].inducer_scaling = 10.0 ** (-2.70035400851161)
        
    elif select_exp is 'threshold':
        
        models[0].bead_params = bead_params_threshold
        models[0].threshold = 0.456084933585625
        models[0].nbhd_size = None
        models[0].inhibitor_scaling = 10.0 ** (-2.87675720588008)
        models[0].inducer_scaling = 10.0 ** (0.912230470757642)

        models[1].bead_params = bead_params_threshold
        models[1].threshold = 0.450426293687721
        models[1].nbhd_size = 2*np.floor(62.0953642599517) + 1
        models[1].inhibitor_scaling = 10.0 ** 1.43501518320614
        models[1].inducer_scaling = 10.0 ** (1.30491072217825)
        
    else:
        print("Unexpected input for 'select_exp'.\nPlease choose valid choice from 'experiment_options'.\nLoaded 'testing' parameter values.")
        
        models[0].bead_params = bead_params_play
        models[0].threshold = 0.497687138515397
        models[0].nbhd_size = None
        models[0].inhibitor_scaling = 10.0 ** (1.45030857599371)
        models[0].inducer_scaling = 10.0 ** (1.29838186230979)

        models[1].bead_params = bead_params_play
        models[1].threshold = 0.421573076609568
        models[1].nbhd_size = 2*np.floor(56.5325206131255) + 1
        models[1].inhibitor_scaling = 10.0 ** -0.329085126426299
        models[1].inducer_scaling = 10.0 ** (-2.61174561686766)
    
    return models
