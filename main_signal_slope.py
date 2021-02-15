import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import copy
import datetime
from pprint import pprint
from shutil import copy2
import os.path

from initial_params import initial_params
# from bead_params import bead_params_init, bead_params_A, bead_params_B
from model_params import modelN, models

from classes import Embryo
from functions import define_initial_protein_concentrations, setup_embryos, run_model, check_embryos_success, define_experiment_groups, set_params_from_df
from plot_functions import create_presentation_fig_arrays, save_presentation_figs, save_method_figs, save_results_figs

# set up save directory
# x = datetime.datetime.now()
# sub_directory = x.strftime('%Y') + '_' + x.strftime('%m') + '_' + x.strftime('%d') + '_' + x.strftime('%H') + x.strftime('%M') + x.strftime('%S') + '/'
sub_directory = 'testing/'

if sub_directory[-1] != '/':
    sub_directory = sub_directory + '/'

save_directory = 'results/' + sub_directory

# initialize embryos
embryoN = 30
embryos = [Embryo('title', initial_params['number_of_cells']) for i in range(embryoN)]

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
    temp_model_values, temp_model_ylim = create_presentation_fig_arrays(embryos)
        
# save_presentation_figs(models, embryos, model_values, model_ylim, 'results/presentation_figs/')

save_method_figs( models, embryos, model_values, model_ylim, 'Arial', save_directory + 'method/' )
save_results_figs( models, embryos, model_values, model_ylim, 'Arial', save_directory + 'results/' )

exp_name = 'all_exps'

paper_directory = 'results/paper_figures/' + exp_name + '/'
if not os.path.isdir(paper_directory):
    os.mkdir(paper_directory)
save_method_figs( models, embryos, model_values, model_ylim, 'Arial', paper_directory + 'method/' )
save_results_figs( models, embryos, model_values, model_ylim, 'Arial', paper_directory + 'results/' )

report_directory = 'results/report/' + exp_name + '/'
if not os.path.isdir(report_directory):
    os.mkdir(report_directory)
save_method_figs( models, embryos, model_values, model_ylim, 'Clear Sans', report_directory + 'method/' )
save_results_figs( models, embryos, model_values, model_ylim, 'Clear Sans', report_directory + 'results/' )

code_directory = save_directory + 'code/'
paper_code_directory = paper_directory + 'code/'
report_code_directory = report_directory + 'code/'
if not os.path.isdir(code_directory):
    os.mkdir(code_directory)
if not os.path.isdir(paper_code_directory):
    os.mkdir(paper_code_directory)
if not os.path.isdir(report_code_directory):
 os.mkdir(report_code_directory)

# copy2(param_filename, code_directory + 'best_params.csv')
filenames = ['main_signal_slope.py', 'classes.py', 'functions.py', 'plot_functions.py', 'model_params.py', 'bead_params.py', 'initial_params.py']
for filename in filenames:
    copy2(filename, code_directory + filename)
    copy2(filename, paper_code_directory + filename)
    copy2(filename, report_code_directory + filename)

