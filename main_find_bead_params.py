from shutil import copy2
import os
import glob
from copy import deepcopy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from functions_pyDREAM import run_pyDREAM
from functions_find_bead_params import set_params_from_df_2models, check_success_rate_2models
from plot_pyDREAM import save_pyDREAM_out_dataframe, create_pyDREAM_figs_2models

from classes import Embryo
from functions import define_initial_protein_concentrations, setup_embryos, run_model, check_embryos_success, define_experiment_groups

from dicts import param_priors_dict, param_lims_dict, axes_labels_dict

from initial_params import initial_params
from model_params import load_models

anterior = range(150,451)
whole_embryo = range(600)

# select_embryos = [1,2,3,4] # cell pellet
# sub_directory = 'results/dream/bead_lim_entire_cell_pellet_500_stage_XII_disp_50_find_bead/'
# likelihood_region = [anterior for i in select_embryos]
# param_names = ['threshold$^A$', '$b_B^A$', '$b_V^A$',  'n', 'threshold$^B$', '$b_B^B$', '$b_V^B$', 'vg1_cell_conc', 'bmp4_cell_conc', 'cell_pellet_spread']

select_embryos = [5,6,7,8,9] # microbead repeat of cell pellet
sub_directory = 'results/dream/bead_lim_entire_activin_ant_5000_stage_XII_disp_50_find_bead/'
likelihood_region = [anterior for i in select_embryos]
param_names = ['threshold$^A$', '$b_B^A$', '$b_V^A$',  'n', 'threshold$^B$', '$b_B^B$', '$b_V^B$', 'activin_10_conc', 'activin_10_spread', 'bmp4_12_conc', 'bmp4_25_conc', 'bmp4_12_spread', 'bmp4_25_spread']

# select_embryos = [10,11] # dm / bmp4 anterior
# sub_directory = 'results/dream/bead_lim_005_bmp4_ant_5000_stage_XII_disp_50_find_bead/'
# likelihood_region = [anterior for i in select_embryos]
# param_names = ['threshold$^A$', '$b_B^A$', '$b_V^A$',  'n', 'threshold$^B$', '$b_B^B$', '$b_V^B$', 'bmp4_12_conc', 'bmp4_12_spread', 'DM_conc', 'AG1X2_spread']

# select_embryos = [12,13,14] # threshold
# sub_directory = 'results/dream/bead_lim_005_threshold_reduced_5000_stage_XII_disp_50_find_bead/'
# likelihood_region = [anterior for i in select_embryos]
# param_names = ['threshold$^A$', '$b_B^A$', '$b_V^A$',  'n', 'threshold$^B$', '$b_B^B$', '$b_V^B$', 'activin_2_conc', 'activin_2_spread', 'bmp4_6_conc', 'bmp4_6_spread']

# select_embryos = list(range(1,15)) # all exps
# sub_directory = 'results/dream/bead_lim_005_cell_pellet_activin_ant_bmp4_ant_threshold_reduced_10000_stage_XII_disp_50_find_bead/'
# likelihood_region = [anterior for i in select_embryos]
# param_names = ['threshold$^A$', '$b_B^A$', '$b_V^A$',  'n', 'threshold$^B$', '$b_B^B$', '$b_V^B$', 'activin_2_conc', 'activin_2_spread', 'activin_10_conc', 'activin_10_spread', 'bmp4_6_conc','bmp4_12_conc', 'bmp4_25_conc', 'bmp4_6_spread', 'bmp4_12_spread', 'bmp4_25_spread', 'DM_conc', 'AG1X2_spread', 'vg1_cell_conc', 'bmp4_cell_conc', 'cell_pellet_spread']

##########################################################################################################
# save code
if not os.path.isdir(sub_directory):
    os.mkdir(sub_directory)

code_directory = sub_directory + 'code/'
if not os.path.isdir(code_directory):
    os.mkdir(code_directory)

filenames = ['main_find_bead_params.py', 'classes.py', 'functions.py', 'dicts.py', 'plot_functions.py', 'model_params.py', 'bead_params.py', 'initial_params.py', 'functions_pyDREAM.py', 'plot_pyDREAM.py', 'functions_find_bead_params.py']
for filename in filenames:
    copy2(filename, code_directory + filename)
    
##########################################################################################################

dream_params = {
    'nchains' : 5,
    'niterations' : 5000,
    'GRlim' : 1.2 # GR Convergence limit
}

save_directory = sub_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory)

param_N = len(param_names)
axes_labels = [axes_labels_dict[key] for key in param_names]
param_lims = [param_lims_dict[key] for key in param_names]
parameters_to_sample = [param_priors_dict[key] for key in param_names]

models = load_models('all_exps') # it doesn't matter which experiment is loaded as parameter values vary
like_model_with_nbhd = deepcopy(models[1])
like_model_no_nbhd = deepcopy(models[0])
def likelihood(param_vector):

    Delta = 0.05
    
    param_dict = dict(zip(param_names, param_vector))
    param_df = pd.DataFrame([param_dict])
    set_params_from_df_2models(param_df, like_model_with_nbhd, like_model_no_nbhd)
    
    total_logp = np.ndarray((len(select_embryos)*2), dtype=float)
    total_logp_with_nbhd = np.ndarray((len(select_embryos)), dtype=float)
    total_logp_no_nbhd = np.ndarray((len(select_embryos)), dtype=float)
    
    # with nbhd
    embryoN = 15
    embryos = [Embryo('title', initial_params['number_of_cells']) for i in range(embryoN)]
    initial_concentrations = define_initial_protein_concentrations(initial_params)
    
    embryos = setup_embryos(embryos, like_model_with_nbhd, initial_concentrations)
    for idx, emb_idx in enumerate(select_embryos):
        embryo = embryos[emb_idx]
        run_model(embryo, like_model_with_nbhd)

        likelihood_with_nbhd = [0.5 * (1 + np.tanh( embryo.target[cell_idx] * (embryo.model_value[cell_idx] - like_model_with_nbhd.threshold) * (1.0 / Delta) ) ) for cell_idx in likelihood_region[idx]]
        total_logp_with_nbhd[idx] = np.sum(np.log10(likelihood_with_nbhd))
        
    # no nbhd
    embryoN = 15
    embryos = [Embryo('title', initial_params['number_of_cells']) for i in range(embryoN)]
    initial_concentrations = define_initial_protein_concentrations(initial_params)
    
    embryos = setup_embryos(embryos, like_model_no_nbhd, initial_concentrations)
    for idx, emb_idx in enumerate(select_embryos):
        embryo = embryos[emb_idx]
        run_model(embryo, like_model_no_nbhd)

        likelihood_no_nbhd = [0.5 * (1 + np.tanh( embryo.target[cell_idx] * (embryo.model_value[cell_idx] - like_model_no_nbhd.threshold) * (1.0 / Delta) ) ) for cell_idx in likelihood_region[idx]]
        total_logp_no_nbhd[idx] = np.sum(np.log10(likelihood_no_nbhd))
        
    total_logp = np.sum(total_logp_with_nbhd) + np.sum(total_logp_no_nbhd)

    return total_logp

dream_out_suffix = 'dream_out/'
dream_out_directory = save_directory + dream_out_suffix
if not os.path.isdir(dream_out_directory):
    os.mkdir(dream_out_directory)
else:
    files = glob.glob(dream_out_directory + '*')
    for f in files:
        os.remove(f)
run_pyDREAM(parameters_to_sample, likelihood, dream_params, dream_out_directory)

save_pyDREAM_out_dataframe(param_names, dream_params, save_directory, dream_out_suffix)

verify_model_with_nbhd = deepcopy(models[1])
verify_model_no_nbhd = deepcopy(models[0])
dream_success_df = check_success_rate_2models(select_embryos, verify_model_with_nbhd, verify_model_no_nbhd, save_directory)

create_pyDREAM_figs_2models(dream_params, param_names, param_lims, axes_labels, verify_model_no_nbhd.plot_color, verify_model_with_nbhd.plot_color, save_directory)



