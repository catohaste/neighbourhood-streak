from shutil import copy2
import os
import glob
from copy import deepcopy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from functions_pyDREAM import run_pyDREAM, check_success_rate
from plot_pyDREAM import save_pyDREAM_out_dataframe, create_pyDREAM_figs
from bayes_factor import set_bayes_factor_threshold, calculate_model_score

from classes import Embryo
from functions import define_initial_protein_concentrations, setup_embryos, run_model, check_embryos_success, define_experiment_groups, set_params_from_df

from dicts import param_priors_dict, param_lims_dict, axes_labels_dict

from initial_params import initial_params
from model_params import load_models

anterior = range(150,451)
whole_embryo = range(600)

# select_embryos = [1,2,3,4] # testing
# sub_directory = 'results/dream/testing_bead_lim_005_cell_pellet_500_stage_XII_disp_50/'
# likelihood_region = [anterior for i in select_embryos]
# param_names_with_nbhd = ['threshold', '$b_B$', '$b_V$', 'n', 'vg1_cell_conc', 'bmp4_cell_conc', 'cell_pellet_spread']
# param_names_no_nbhd = ['threshold', '$b_B$', '$b_V$', 'vg1_cell_conc', 'bmp4_cell_conc', 'cell_pellet_spread']

# select_embryos = [1,2,3,4] # cell pellet
# sub_directory = 'results/dream/bead_lim_005_cell_pellet_5000_stage_XII_disp_50/'
# likelihood_region = [anterior for i in select_embryos]
# param_names_with_nbhd = ['threshold', '$b_B$', '$b_V$', 'n', 'vg1_cell_conc', 'bmp4_cell_conc', 'cell_pellet_spread']
# param_names_no_nbhd = ['threshold', '$b_B$', '$b_V$', 'vg1_cell_conc', 'bmp4_cell_conc', 'cell_pellet_spread']

select_embryos = [5,6,7,8,9] # microbead repeat of cell pellet
sub_directory = 'results/dream/bead_lim_005_activin_ant_5000_stage_XII_disp_50/'
likelihood_region = [anterior for i in select_embryos]
param_names_with_nbhd = ['threshold', '$b_B$', '$b_V$', 'n', 'activin_10_conc', 'activin_10_spread', 'bmp4_12_conc', 'bmp4_25_conc', 'bmp4_12_spread', 'bmp4_25_spread' ]
param_names_no_nbhd = ['threshold', '$b_B$', '$b_V$', 'activin_10_conc', 'activin_10_spread', 'bmp4_12_conc', 'bmp4_25_conc', 'bmp4_12_spread', 'bmp4_25_spread']

# select_embryos = [10,11] # dm / bmp4 anterior
# sub_directory = 'results/dream/bead_lim_005_bmp4_ant_5000_stage_XII_disp_50/'
# likelihood_region = [anterior for i in select_embryos]
# param_names_with_nbhd = ['threshold', '$b_B$', '$b_V$', 'n', 'bmp4_12_conc', 'bmp4_12_spread', 'DM_conc', 'AG1X2_spread']
# param_names_no_nbhd = ['threshold', '$b_B$', '$b_V$', 'bmp4_12_conc', 'bmp4_12_spread', 'DM_conc', 'AG1X2_spread']

# select_embryos = [12,13,14] # threshold
# sub_directory = 'results/dream/bead_lim_005_threshold_reduced_5000_stage_XII_disp_50/'
# likelihood_region = [anterior for i in select_embryos]
# param_names_with_nbhd = ['threshold', '$b_B$', '$b_V$', 'n', 'activin_2_conc', 'activin_2_spread',  'bmp4_6_conc', 'bmp4_6_spread']
# param_names_no_nbhd = ['threshold', '$b_B$', '$b_V$', 'activin_2_conc', 'activin_2_spread',  'bmp4_6_conc', 'bmp4_6_spread']


# select_embryos = list(range(1,15)) # all exps
# sub_directory = 'results/dream/bead_lim_005_cell_pellet_activin_ant_bmp4_ant_threshold_reduced_10000_stage_XII_disp_50/'
# likelihood_region = [anterior for i in select_embryos]
# param_names_with_nbhd = ['threshold', '$b_B$', '$b_V$', 'n', 'activin_2_conc', 'activin_10_conc', 'activin_2_spread', 'activin_10_spread','bmp4_6_conc', 'bmp4_12_conc', 'bmp4_25_conc', 'bmp4_6_spread', 'bmp4_12_spread', 'bmp4_25_spread', 'DM_conc', 'AG1X2_spread', 'vg1_cell_conc', 'bmp4_cell_conc', 'cell_pellet_spread']
# param_names_no_nbhd = ['threshold', '$b_B$', '$b_V$', 'activin_2_conc', 'activin_10_conc', 'activin_2_spread', 'activin_10_spread','bmp4_6_conc', 'bmp4_12_conc', 'bmp4_25_conc', 'bmp4_6_spread', 'bmp4_12_spread', 'bmp4_25_spread', 'DM_conc', 'AG1X2_spread', 'vg1_cell_conc', 'bmp4_cell_conc', 'cell_pellet_spread']


if not os.path.isdir(sub_directory):
    os.mkdir(sub_directory)

dream_params = {
    'nchains' : 5,
    'niterations' : 5000,
    'GRlim' : 1.2 # GR Convergence limit
}

models = load_models('all_exps') # it doesn't matter which experiment is loaded as parameter values vary

#### with_nbhd ###########################################################################################
save_directory = sub_directory + 'with_nbhd/'
if not os.path.isdir(save_directory):
    os.mkdir(save_directory)


param_N = len(param_names_with_nbhd)
axes_labels = [axes_labels_dict[key] for key in param_names_with_nbhd]
param_lims = [param_lims_dict[key] for key in param_names_with_nbhd]
parameters_to_sample = [param_priors_dict[key] for key in param_names_with_nbhd]

like_model_with_nbhd = deepcopy(models[1])
def likelihood_with_nbhd(param_vector):

    Delta = 0.05

    embryoN = 15
    embryos = [Embryo('title', initial_params['number_of_cells']) for i in range(embryoN)]
    initial_concentrations = define_initial_protein_concentrations(initial_params)

    param_dict = dict(zip(param_names_with_nbhd, param_vector))
    param_df = pd.DataFrame([param_dict])

    set_params_from_df(param_df, like_model_with_nbhd)

    total_logp = np.ndarray((len(select_embryos)), dtype=float)
    embryos = setup_embryos(embryos, like_model_with_nbhd, initial_concentrations)
    for idx, emb_idx in enumerate(select_embryos):
        embryo = embryos[emb_idx]
        run_model(embryo, like_model_with_nbhd)

        likelihood = [0.5 * (1 + np.tanh( embryo.target[cell_idx] * (embryo.model_value[cell_idx] - like_model_with_nbhd.threshold) * (1.0 / Delta) ) ) for cell_idx in likelihood_region[idx]]
        total_logp[idx] = np.sum(np.log10(likelihood))

    return np.sum(total_logp)


dream_out_suffix = 'dream_out/'
dream_out_directory = save_directory + dream_out_suffix
if not os.path.isdir(dream_out_directory):
    os.mkdir(dream_out_directory)
else:
    files = glob.glob(dream_out_directory + '*')
    for f in files:
        os.remove(f)
run_pyDREAM(parameters_to_sample, likelihood_with_nbhd, dream_params, dream_out_directory)

save_pyDREAM_out_dataframe(param_names_with_nbhd, dream_params, save_directory, dream_out_suffix)

verify_model_with_nbhd = deepcopy(models[1])
check_success_rate(select_embryos, verify_model_with_nbhd, save_directory)

create_pyDREAM_figs(dream_params, param_names_with_nbhd, param_lims, axes_labels, verify_model_with_nbhd.plot_color, save_directory)


#### no_nbhd ###########################################################################################

save_directory = sub_directory + 'no_nbhd/'
if not os.path.isdir(save_directory):
    os.mkdir(save_directory)

param_N = len(param_names_no_nbhd)
axes_labels = [axes_labels_dict[key] for key in param_names_no_nbhd]
param_lims = [param_lims_dict[key] for key in param_names_no_nbhd]
parameters_to_sample = [param_priors_dict[key] for key in param_names_no_nbhd]

like_model_no_nbhd = deepcopy(models[0])
def likelihood_no_nbhd(param_vector):

    Delta = 0.05

    embryoN = 15
    embryos = [Embryo('title', initial_params['number_of_cells']) for i in range(embryoN)]
    initial_concentrations = define_initial_protein_concentrations(initial_params)

    param_dict = dict(zip(param_names_no_nbhd, param_vector))
    param_df = pd.DataFrame([param_dict])
    
    set_params_from_df(param_df, like_model_no_nbhd)

    total_logp = np.ndarray((len(select_embryos)), dtype=float)
    embryos = setup_embryos(embryos, like_model_no_nbhd, initial_concentrations)
    for idx, emb_idx in enumerate(select_embryos):
        embryo = embryos[emb_idx]
        run_model(embryo, like_model_no_nbhd)

        likelihood = [0.5 * (1 + np.tanh( embryo.target[cell_idx] * (embryo.model_value[cell_idx] - like_model_no_nbhd.threshold) * (1.0 / Delta) ) ) for cell_idx in likelihood_region[idx]]
        total_logp[idx] = np.sum(np.log10(likelihood))

    return np.sum(total_logp)


dream_out_suffix = 'dream_out/'
dream_out_directory = save_directory + dream_out_suffix
if not os.path.isdir(dream_out_directory):
    os.mkdir(dream_out_directory)
else:
    files = glob.glob(dream_out_directory + '*')
    for f in files:
        os.remove(f)
run_pyDREAM(parameters_to_sample, likelihood_no_nbhd, dream_params, dream_out_directory)

save_pyDREAM_out_dataframe(param_names_no_nbhd, dream_params, save_directory, dream_out_suffix)

verify_model_no_nbhd = deepcopy(models[0])
check_success_rate(select_embryos, verify_model_no_nbhd, save_directory)

create_pyDREAM_figs(dream_params, param_names_no_nbhd, param_lims, axes_labels, verify_model_no_nbhd.plot_color, save_directory)

##########################################################################################################

code_directory = sub_directory + 'code/'
if not os.path.isdir(code_directory):
    os.mkdir(code_directory)

filenames = ['main_pyDREAM.py', 'classes.py', 'functions.py', 'dicts.py', 'plot_functions.py', 'model_params.py', 'bead_params.py', 'initial_params.py', 'functions_pyDREAM.py', 'plot_pyDREAM.py', 'bayes_factor.py']
for filename in filenames:
    copy2(filename, code_directory + filename)
