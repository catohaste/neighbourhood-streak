import copy
import os.path
from shutil import copy2

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from initial_params import initial_params
from model_params import load_models

from classes import Embryo
from functions import define_initial_protein_concentrations, setup_embryos, run_model, check_embryos_success, define_experiment_groups, set_params_from_df
from plot_functions import create_presentation_fig_arrays, save_presentation_figs, save_method_figs, save_method_figs_poster, save_results_figs , save_results_figs_poster, set_up_protein_fig, set_up_fig_trio

# initialize embryos
embryoN = 15
embryos = [Embryo('title', initial_params['number_of_cells']) for i in range(embryoN)]

# select experiment
experiment_options = ['testing', 'all_exps', 'cell_pellet', 'activin_ant', 'bmp4_ant', 'threshold', 'run2_A_B']
select_exp = 'activin_ant'
models = load_models(select_exp)

# set up save directory
results_directory = 'output/'
if not os.path.isdir(results_directory):
    os.mkdir(results_directory)
save_directory = 'output/' + select_exp + '/'
if not os.path.isdir(save_directory):
    os.mkdir(save_directory)
    
# save code, to ensure reproducibility
code_directory = save_directory + 'code/'
# poster_code_directory = poster_directory + 'code/'
# paper_code_directory = paper_directory + 'code/'
# report_code_directory = report_directory + 'code/'
if not os.path.isdir(code_directory):
    os.mkdir(code_directory)
# if not os.path.isdir(poster_code_directory):
#     os.mkdir(poster_code_directory)
# if not os.path.isdir(paper_code_directory):
#     os.mkdir(paper_code_directory)
# if not os.path.isdir(report_code_directory):
#     os.mkdir(report_code_directory)

# copy2(param_filename, code_directory + 'best_params.csv')
filenames = ['main_signal_slope.py', 'classes.py', 'functions.py', 'plot_functions.py', 'model_params.py', 'bead_params.py', 'initial_params.py']
for filename in filenames:
    copy2(filename, code_directory + filename)
    # copy2(filename, poster_code_directory + filename)
    # copy2(filename, paper_code_directory + filename)
    # copy2(filename, report_code_directory + filename)

# initialize arrays for plots
modelN = len(models)
model_values = np.ndarray((modelN, embryoN, initial_params['number_of_cells']), dtype=float)
model_ylim = np.ndarray((modelN, embryoN, 2), dtype=float)

for model_idx, model in enumerate(models):
    
    initial_concentrations = define_initial_protein_concentrations(initial_params)
    embryos = setup_embryos(embryos, model, initial_concentrations)
    
    for embryo in embryos:
        run_model(embryo, model)
        embryo.find_streaks()
        
    successN, failureN = check_embryos_success(embryos)
    experiments = define_experiment_groups(embryos)
    for exp in experiments:
        exp.find_plot_model_ylim()
        
    model_values[model_idx,:,:], model_ylim[model_idx,:,:] = create_presentation_fig_arrays(embryos)
        
# save_presentation_figs(models, embryos, model_values, model_ylim, 'results/presentation_figs/')

save_method_figs( models, embryos, model_values, model_ylim, 'Arial', save_directory + 'method/' )
save_results_figs( models, embryos, model_values, model_ylim, 'Arial', save_directory + 'results/' )

# poster_directory = 'results/poster_figures/' + select_exp + '/'
# if not os.path.isdir(poster_directory):
#     os.mkdir(poster_directory)
# save_method_figs_poster( models, embryos, model_values, model_ylim, 'Arial', poster_directory + 'method/' )
# save_results_figs_poster( models, embryos, model_values, model_ylim, 'Arial', poster_directory + 'results/' )

# paper_directory = 'results/paper_figures/' + select_exp + '/'
# if not os.path.isdir(paper_directory):
#     os.mkdir(paper_directory)
# save_method_figs( models, embryos, model_values, model_ylim, 'Arial', paper_directory + 'method/' )
# save_results_figs( models, embryos, model_values, model_ylim, 'Arial', paper_directory + 'results/' )

# report_directory = 'results/report/' + select_exp + '/'
# if not os.path.isdir(report_directory):
#     os.mkdir(report_directory)
# save_method_figs( models, embryos, model_values, model_ylim, 'Clear Sans', report_directory + 'method/' )
# save_results_figs( models, embryos, model_values, model_ylim, 'Clear Sans', report_directory + 'results/' )

