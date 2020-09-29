import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import copy
import datetime
from pprint import pprint
from shutil import copy2
import os.path
import pandas as pd

from initial_params import initial_params
from bead_params import bead_params_init, bead_params_A, bead_params_B
from model_params import models

from classes import Embryo
from functions import define_initial_protein_concentrations, setup_embryos, run_model, check_embryos_success, define_experiment_groups
from plot_functions import save_standard_figs, save_model_figs

# set up save directory
sub_directory = 'report/'
x = datetime.datetime.now()
sub_directory = x.strftime('%Y') + '_' + x.strftime('%m') + '_' + x.strftime('%d') + '_' + x.strftime('%H') + x.strftime('%M') + x.strftime('%S') + '/'
sub_directory = 'testing/'
sub_directory = 'vary_params/inducer_inhibitor/'

if sub_directory[-1] != '/':
    sub_directory = sub_directory + '/'

save_directory = 'results/' + sub_directory

# initialize embryos
embryoN = 23
embryos = [Embryo('title', initial_params['number_of_cells']) for i in range(embryoN)]
select_indices = [0,1,2,3]

# # vary certain params
# mins = [0.1, 0.2, 0.3, 0.4, 0.5]
# maxs = [0.4, 0.5, 0.6, 0.7, 0.8]
# mins = [0.1]
# maxs = [0.4, 0.7]
# lims = []
# for min_val in mins:
#     for max_val in maxs:
#         if min_val < max_val:
#             lims.append((min_val, max_val))

b_V_list = np.arange(-3,3,0.5)
b_B_list = np.arange(-3,3,0.5)
            
results = pd.DataFrame(columns=['model', '$b_V$', '$b_B$', 'success', 'streakN', 'positions'])
counter = 0
for model in models:
    for b_V in b_V_list:
        for b_B in b_B_list:
            model.inducer_scaling = 10.0 ** b_V
            model.inhibitor_scaling = 10.0 ** b_B
            initial_concentrations = define_initial_protein_concentrations(initial_params)
            embryos = setup_embryos(embryos, model, initial_concentrations)
            for embryo_idx in select_indices:
                run_model(embryos[embryo_idx], model)
                embryos[embryo_idx].find_streaks()
                save_standard_figs(embryos[embryo_idx], model, save_directory)
            successN, failureN = check_embryos_success(embryos)
            # print(successN, failureN)
            experiments = define_experiment_groups(embryos)
            for exp in [experiments[7]]:
                exp.find_plot_model_ylim()
            for embryo_idx in select_indices:
                results.at[counter,'$b_V$'] = b_V
                results.at[counter,'$b_B$'] = b_B
                results.at[counter,'embryo_idx'] = embryo_idx
                results.at[counter, 'success'] = embryos[embryo_idx].success
                results.at[counter, 'streakN'] = embryos[embryo_idx].streakN
                results.at[counter, 'positions'] = embryos[embryo_idx].streak_positions
                results.at[counter,'model'] = model.name
                save_model_figs(embryos[embryo_idx], model, save_directory, '_' + str(b_B) + '_' + str(b_V))
                counter += 1
                
print('success percentage', (len(results[results.success])*100) / len(results) )

print('success percentage', (len(results[(results.success) & (results.model == 'inducer_SMAD_nbhd')])*100) / len(results[(results.model == 'inducer_SMAD_nbhd')]) )

code_directory = save_directory + 'code/'
if not os.path.isdir(code_directory):
    os.mkdir(code_directory)

results.to_csv(save_directory + 'results.tsv', sep='\t')

filenames = ['main_vary_params.py', 'classes.py', 'functions.py', 'plot_functions.py', 'model_params.py', 'bead_params.py', 'initial_params.py']
for filename in filenames:
    copy2(filename, code_directory + filename)

