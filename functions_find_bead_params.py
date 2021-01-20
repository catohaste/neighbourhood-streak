import os.path
import numpy as np
import pandas as pd
from copy import deepcopy
from pydream.core import run_dream
from pydream.convergence import Gelman_Rubin

from classes import Embryo
from functions import define_initial_protein_concentrations, setup_embryos, run_model, check_embryos_success, define_experiment_groups, set_params_from_df

from initial_params import initial_params
from model_params import models

def set_params_from_df_2models(df, model_with_nbhd, model_no_nbhd):
    
    if '$b_B^A$' in df.columns:
        model_with_nbhd.inhibitor_scaling = 10.0 ** df.iloc[0]['$b_B^A$']
    if '$b_B^B$' in df.columns:
        model_no_nbhd.inhibitor_scaling = 10.0 ** df.iloc[0]['$b_B^B$']
    if '$b_V^A$' in df.columns:
        model_with_nbhd.inducer_scaling = 10.0 ** df.iloc[0]['$b_V^A$']
    if '$b_V^B$' in df.columns:
        model_no_nbhd.inducer_scaling = 10.0 ** df.iloc[0]['$b_V^B$']
    if 'threshold$^A$' in df.columns:
        model_with_nbhd.threshold = df.iloc[0]['threshold$^A$']
    if 'threshold$^B$' in df.columns:
        model_no_nbhd.threshold = df.iloc[0]['threshold$^B$']
    if 'n' in df.columns:
        model_with_nbhd.nbhd_size = 2*np.floor(df.iloc[0]['n']) - 1
    if 'activin_conc' in df.columns:
        model_with_nbhd.bead_params['activin_2_conc'] = df.iloc[0]['activin_conc']
        model_no_nbhd.bead_params['activin_2_conc'] = df.iloc[0]['activin_conc']
        model_with_nbhd.bead_params['activin_10_conc'] = df.iloc[0]['activin_conc']
        model_no_nbhd.bead_params['activin_10_conc'] = df.iloc[0]['activin_conc']
    if 'activin_2_conc' in df.columns:
        model_with_nbhd.bead_params['activin_2_conc'] = df.iloc[0]['activin_2_conc']
        model_no_nbhd.bead_params['activin_2_conc'] = df.iloc[0]['activin_2_conc']
    if 'activin_10_conc' in df.columns:
        model_with_nbhd.bead_params['activin_10_conc'] = df.iloc[0]['activin_10_conc']
        model_no_nbhd.bead_params['activin_10_conc'] = df.iloc[0]['activin_10_conc']
    if 'activin_spread' in df.columns:
        model_with_nbhd.bead_params['heparin_2_spread'] = df.iloc[0]['activin_spread']
        model_with_nbhd.bead_params['heparin_10_spread'] = df.iloc[0]['activin_spread']
        model_no_nbhd.bead_params['heparin_2_spread'] = df.iloc[0]['activin_spread']
        model_no_nbhd.bead_params['heparin_10_spread'] = df.iloc[0]['activin_spread']
    if 'activin_2_spread' in df.columns:
        model_with_nbhd.bead_params['heparin_2_spread'] = df.iloc[0]['activin_2_spread']
        model_no_nbhd.bead_params['heparin_2_spread'] = df.iloc[0]['activin_2_spread']
    if 'activin_10_spread' in df.columns:
        model_with_nbhd.bead_params['heparin_10_spread'] = df.iloc[0]['activin_10_spread']
        model_no_nbhd.bead_params['heparin_10_spread'] = df.iloc[0]['activin_10_spread']
    if 'bmp4_conc' in df.columns:
        model_with_nbhd.bead_params['bmp4_6_conc'] = df.iloc[0]['bmp4_conc']
        model_no_nbhd.bead_params['bmp4_6_conc'] = df.iloc[0]['bmp4_conc']
        model_with_nbhd.bead_params['bmp4_12_conc'] = df.iloc[0]['bmp4_conc']
        model_no_nbhd.bead_params['bmp4_12_conc'] = df.iloc[0]['bmp4_conc']
        model_with_nbhd.bead_params['bmp4_25_conc'] = df.iloc[0]['bmp4_conc']
        model_no_nbhd.bead_params['bmp4_25_conc'] = df.iloc[0]['bmp4_conc']
        model_with_nbhd.bead_params['bmp4_50_conc'] = df.iloc[0]['bmp4_conc']
        model_no_nbhd.bead_params['bmp4_50_conc'] = df.iloc[0]['bmp4_conc']
    if 'bmp4_50_conc' in df.columns:
        model_with_nbhd.bead_params['bmp4_50_conc'] = df.iloc[0]['bmp4_50_conc']
        model_no_nbhd.bead_params['bmp4_50_conc'] = df.iloc[0]['bmp4_50_conc']
    if 'bmp4_25_conc' in df.columns:
        model_with_nbhd.bead_params['bmp4_25_conc'] = df.iloc[0]['bmp4_25_conc']
        model_no_nbhd.bead_params['bmp4_25_conc'] = df.iloc[0]['bmp4_25_conc']
    if 'bmp4_12_conc' in df.columns:
        model_with_nbhd.bead_params['bmp4_12_conc'] = df.iloc[0]['bmp4_12_conc']
        model_no_nbhd.bead_params['bmp4_12_conc'] = df.iloc[0]['bmp4_12_conc']
    if 'bmp4_6_conc' in df.columns:
        model_with_nbhd.bead_params['bmp4_6_conc'] = df.iloc[0]['bmp4_6_conc']
        model_no_nbhd.bead_params['bmp4_6_conc'] = df.iloc[0]['bmp4_6_conc']
    if 'bmp4_spread' in df.columns:
        model_with_nbhd.bead_params['afigel_50_spread'] = df.iloc[0]['bmp4_spread']
        model_with_nbhd.bead_params['afigel_25_spread'] = df.iloc[0]['bmp4_spread']
        model_with_nbhd.bead_params['afigel_12_spread'] = df.iloc[0]['bmp4_spread']
        model_with_nbhd.bead_params['afigel_6_spread'] = df.iloc[0]['bmp4_spread']
        model_no_nbhd.bead_params['afigel_50_spread'] = df.iloc[0]['bmp4_spread']
        model_no_nbhd.bead_params['afigel_25_spread'] = df.iloc[0]['bmp4_spread']
        model_no_nbhd.bead_params['afigel_12_spread'] = df.iloc[0]['bmp4_spread']
        model_no_nbhd.bead_params['afigel_6_spread'] = df.iloc[0]['bmp4_spread']
    if 'bmp4_50_spread' in df.columns:
        model_with_nbhd.bead_params['afigel_50_spread'] = df.iloc[0]['bmp4_50_spread']
        model_no_nbhd.bead_params['afigel_50_spread'] = df.iloc[0]['bmp4_50_spread']
    if 'bmp4_25_spread' in df.columns:
        model_with_nbhd.bead_params['afigel_25_spread'] = df.iloc[0]['bmp4_25_spread']
        model_no_nbhd.bead_params['afigel_25_spread'] = df.iloc[0]['bmp4_25_spread']
    if 'bmp4_12_spread' in df.columns:
        model_with_nbhd.bead_params['afigel_12_spread'] = df.iloc[0]['bmp4_12_spread']
        model_no_nbhd.bead_params['afigel_12_spread'] = df.iloc[0]['bmp4_12_spread']
    if 'bmp4_6_spread' in df.columns:
        model_with_nbhd.bead_params['afigel_6_spread'] = df.iloc[0]['bmp4_6_spread']
        model_no_nbhd.bead_params['afigel_6_spread'] = df.iloc[0]['bmp4_6_spread']
    if 'DM_conc' in df.columns:
        model_with_nbhd.bead_params['DM_conc'] = df.iloc[0]['DM_conc']
        model_no_nbhd.bead_params['DM_conc'] = df.iloc[0]['DM_conc']
    if 'AG1X2_spread' in df.columns:
        model_with_nbhd.bead_params['AG1X2_spread'] = df.iloc[0]['AG1X2_spread']
        model_no_nbhd.bead_params['AG1X2_spread'] = df.iloc[0]['AG1X2_spread']
        
    return model_with_nbhd, model_no_nbhd

    
def check_success_rate_2models(select_embryos, model_with_nbhd, model_no_nbhd, save_directory):
    
    df = pd.read_csv(save_directory + 'dream_out.tsv', sep='\t', index_col=0, header=0)
    df = df.drop(columns=['chainID'])
    df = df.drop_duplicates()
    df = df.sort_values(by='logp', ascending=False)
    
    top_params_N = int(np.ceil(len(df) * 1))
    top_params = df.iloc[:top_params_N,:]
    
    colN = len(top_params.columns)
    for idx, emb_idx in enumerate(select_embryos):
        colN = len(top_params.columns)
        top_params.insert(colN, 'A_' + str(emb_idx) , False)
    for idx, emb_idx in enumerate(select_embryos):
        colN = len(top_params.columns)
        top_params.insert(colN + 1, 'B_' + str(emb_idx) , False)
    colN = len(top_params.columns)
    top_params.insert(colN, 'A_success_proportion' , np.nan)
    top_params.insert(colN + 1, 'B_success_proportion' , np.nan)
    
    for index, row in top_params.iterrows():
    
        embryoN = 30
        row_df = row.to_frame().T
        model_with_nbhd, model_no_nbhd = set_params_from_df_2models(row_df, model_with_nbhd, model_no_nbhd)
        
        # with nbhd
        embryos = [Embryo('title', initial_params['number_of_cells']) for i in range(embryoN)]
        initial_concentrations = define_initial_protein_concentrations(initial_params)
    
        embryos = setup_embryos(embryos, model_with_nbhd, initial_concentrations)
        for idx, emb_idx in enumerate(select_embryos):
            embryo = embryos[emb_idx]
            run_model(embryo, model_with_nbhd)
            embryo.find_streaks()
    
        successN, failureN = check_embryos_success(embryos)
        for idx, emb_idx in enumerate(select_embryos):
            top_params.at[index, 'A_' + str(emb_idx)] = embryos[emb_idx].success
            
        top_params.at[index,'A_success_proportion'] = successN / len(select_embryos)
            
        # no nbhd
        embryos = [Embryo('title', initial_params['number_of_cells']) for i in range(embryoN)]
        initial_concentrations = define_initial_protein_concentrations(initial_params)
    
        embryos = setup_embryos(embryos, model_no_nbhd, initial_concentrations)
        for idx, emb_idx in enumerate(select_embryos):
            embryo = embryos[emb_idx]
            run_model(embryo, model_no_nbhd)
            embryo.find_streaks()
    
        successN, failureN = check_embryos_success(embryos)
        for idx, emb_idx in enumerate(select_embryos):
            top_params.at[index, 'B_' + str(emb_idx)] = embryos[emb_idx].success
        
        top_params.at[index,'B_success_proportion'] = successN / len(select_embryos)
    
    top_params.to_csv(save_directory + 'top_params.tsv', sep='\t')
    
    return top_params
