import os.path
import numpy as np
import pandas as pd
from copy import deepcopy
from pydream.core import run_dream
from pydream.convergence import Gelman_Rubin

from classes import Embryo
from functions import define_initial_protein_concentrations, setup_embryos, run_model, check_embryos_success, define_experiment_groups, set_params_from_df

from initial_params import initial_params
    

def run_pyDREAM(parameters_to_sample, likelihood_function, dream_params, save_directory):
    
    nchains = dream_params['nchains']
    niterations = dream_params['niterations']
    GRlim = dream_params['GRlim']
    
    converged = False
    total_iterations = niterations

	#Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    sampled_params, log_ps = run_dream(parameters_to_sample, likelihood_function, niterations=niterations, nchains=nchains, multitry=False, parallel=False, verbose=False, model_name=save_directory)
    
	#Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save(save_directory + 'sampled_params_chain' + str(chain) + '_' + str(total_iterations), sampled_params[chain])
        np.save(save_directory + 'logps_chain' + str(chain) + '_' + str(total_iterations), log_ps[chain])

	# Check convergence and continue sampling if not converged
    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ' + str(total_iterations) + ' GR = ', GR)
    np.savetxt(save_directory + 'GelmanRubin_iteration_' + str(total_iterations) + '.out', GR)

    old_samples = sampled_params
    if np.any(GR > GRlim):
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
        while not converged:
            total_iterations += niterations

            sampled_params, log_ps = run_dream(parameters_to_sample, likelihood_function, niterations=niterations, nchains=nchains, multitry=False, parallel = False, verbose=False, model_name=save_directory)
            
            for chain in range(len(sampled_params)):
                np.save(save_directory + 'sampled_params_chain' + str(chain) + '_' + str(total_iterations), sampled_params[chain])
                np.save(save_directory + 'logps_chain' + str(chain) + '_' + str(total_iterations), log_ps[chain])
                
            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt(save_directory + 'GelmanRubin.txt', GR)
            
            if np.all(GR < 1.2):
                converged = True
                
    return
    
    
def check_success_rate(select_embryos, current_model, save_directory):
        
    df = pd.read_csv(save_directory + 'dream_out.tsv', sep='\t', index_col=0, header=0)
    df = df.drop(columns=['chainID'])
    df = df.drop_duplicates()
    df = df.sort_values(by='logp', ascending=False)
    
    top_params_N = int(np.ceil(len(df) * 1))
    top_params = df.iloc[:top_params_N,:]

    for idx, emb_idx in enumerate(select_embryos):
        colN = len(top_params.columns)
        top_params.insert(colN, str(emb_idx) , False)
    colN = len(top_params.columns)
    top_params.insert(colN, 'success_proportion' , np.nan)
    
    for index, row in top_params.iterrows():
    
        embryoN = 15
        embryos = [Embryo('title', initial_params['number_of_cells']) for i in range(embryoN)]
        initial_concentrations = define_initial_protein_concentrations(initial_params)

        row_df = row.to_frame().T
        current_model = set_params_from_df(row_df, current_model)
    
        embryos = setup_embryos(embryos, current_model, initial_concentrations)
        for idx, emb_idx in enumerate(select_embryos):
            embryo = embryos[emb_idx]
            run_model(embryo, current_model)
            embryo.find_streaks()
    
        successN, failureN = check_embryos_success(embryos)
        for idx, emb_idx in enumerate(select_embryos):
            top_params.at[index, str(emb_idx)] = embryos[emb_idx].success
        
        top_params.at[index,'success_proportion'] = successN / len(select_embryos)
    
    top_params.to_csv(save_directory + 'top_params.tsv', sep='\t')
    
    return

