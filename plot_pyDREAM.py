import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import pandas as pd
import seaborn as sns
import glob
import os.path

from itertools import combinations

from plot_functions import set_figs_font_settings

def save_pyDREAM_out_dataframe(param_names, dream_params, save_directory, dream_out_suffix):
    
    dream_out_directory = save_directory + dream_out_suffix    
    param_N = len(param_names)

    col_names = np.concatenate((['chainID'],param_names,['logp']))
    
    nchains = dream_params['nchains']
    niterations = dream_params['niterations']

    chainID = np.ndarray((0,1))
    logp = np.ndarray((0,1))
    params = np.ndarray((0,param_N))

    for chain in range(nchains):
        counter = 0
        for file in glob.glob(dream_out_directory + 'logps_chain' + str(chain) + '_*.npy'):

            temp_logp = np.load(file)
            logp = np.append(logp, temp_logp)
        
            iteration = file.split(sep='.',maxsplit=1)[0].split(sep='_')[-1]
            params_file = dream_out_directory + 'sampled_params_chain' + str(chain) + '_' + iteration + '.npy'
            temp_params = np.load(params_file)
            params = np.append(params, temp_params, axis=0)
        
            temp_chainID = [chain for i in range(niterations)]
            chainID = np.append(chainID, temp_chainID)

    chainID = np.array([[val] for val in chainID])
    logp = np.array([[val] for val in logp])

    data = np.concatenate((chainID,params,logp),axis=1)

    df = pd.DataFrame(columns=col_names, data=data)
    
    df.to_csv(save_directory + 'dream_out.tsv', sep='\t')
    
    pyDREAM_out_dataframe = df
    
    return pyDREAM_out_dataframe
    

    
def plot_logp_over_time(dream_out_df, dream_params, figs_directory):
    
    font_size = set_figs_font_settings()
    
    df = dream_out_df    
    nchains = dream_params['nchains']
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for chain in range(nchains):
        temp_df = df[df['chainID']==chain]
        temp_df = temp_df.reset_index()
        temp_df.plot(ax=ax, y='logp', kind='line', label='chain ' + str(chain))
    ax.set_xlabel('Iteration')
    ax.set_ylabel('log$_{10}$(Likelihood)')
    plt.tight_layout()
    plt.savefig(figs_directory + 'logp_over_time.png')
    plt.close()
    
    return
    
    
def plot_param_pair_grid_logp(dream_out_df, param_names, param_lims, axes_labels, figs_directory):
    
    font_size = set_figs_font_settings()
    
    plt.rc('xtick', labelsize=font_size['SMALLEST_SIZE'])    # fontsize of the tick labels
    plt.rc('ytick', labelsize=font_size['SMALLEST_SIZE'])    # fontsize of the tick labels
    plt.rc('axes', labelsize=font_size['SMALL_SIZE'])    # fontsize of the x and y labels
    
    param_N = len(param_names)
    
    df = dream_out_df
    df = df.drop_duplicates()
    df = df.sort_values(by='logp', ascending=True)
   
    """ find settings for the colorbar """
    logp_min = df['logp'].min()
    logp_max = df['logp'].max()

    cmin = np.floor(logp_min / 10) * 10
    cmax = np.ceil(logp_max / 10) * 10
    
    df = df.reset_index()
    
    color_cutoff = df.at[int(np.ceil(0.1 * len(df))),'logp']
    
    fig, axes = plt.subplots(nrows=param_N, ncols=param_N, figsize=(14,9))
    
    param_pair_indices = combinations(range(param_N), 2)
    for param_pair_index in param_pair_indices:

        par1_idx = param_pair_index[0]
        par2_idx = param_pair_index[1]

        par1 = param_names[param_pair_index[0]]
        par2 = param_names[param_pair_index[1]]

        lims = (param_lims[par1_idx], param_lims[par2_idx])
        xmin = lims[1][0]
        xmax = lims[1][1]
        ymin = lims[0][0]
        ymax = lims[0][1]

        xheight = xmax - xmin
        xmin = xmin - 0.05*xheight
        xmax = xmax + 0.05*xheight
        yheight = ymax - ymin
        ymin = ymin - 0.05*yheight
        ymax = ymax + 0.05*yheight

        scatter_cmap = sns.cubehelix_palette(8, start=.5, rot=-.75, reverse=True, as_cmap=True)
        axes[par1_idx, par2_idx].scatter(x=df[par2],y=df[par1], c=df['logp'], cmap=scatter_cmap, vmin=color_cutoff, vmax=cmax, s=2)
        axes[par1_idx, par2_idx].set_xlim([xmin, xmax])
        axes[par1_idx, par2_idx].set_ylim([ymin, ymax])

        axes[par1_idx, par2_idx].set_xticklabels([])
        axes[par1_idx, par2_idx].set_yticklabels([])

        kde_cmap =  sns.cubehelix_palette(8, start=.5, rot=-.75, reverse=True, as_cmap=True)
        sns.kdeplot(ax=axes[par2_idx, par1_idx], data=df[par1], data2=df[par2], shade=False, cmap=kde_cmap)
        axes[par2_idx, par1_idx].set_xlim([ymin, ymax])
        axes[par2_idx, par1_idx].set_ylim([xmin, xmax])
        
        axes[par2_idx, par1_idx].set_xlabel('')
        axes[par2_idx, par1_idx].set_ylabel('')
        
        if par2_idx is not (param_N - 1):
            axes[par2_idx, par1_idx].set_xticklabels([])
        if par1_idx is not 0:
            axes[par2_idx, par1_idx].set_yticklabels([])

    distplot_palette =  sns.cubehelix_palette(8, start=.5, rot=-.75, light=0.5, reverse=True, as_cmap=False)
    # distplot_palette = sns.cubehelix_palette(param_N, light=0.5, start=0.5, reverse=True, as_cmap=False)
    for par_idx, par in enumerate(param_names):

        sns.distplot( df[par], ax=axes[par_idx, par_idx], color=distplot_palette[5])
        xmin, xmax = param_lims[par_idx]
        xheight = xmax - xmin
        xmin = xmin - 0.05 * xheight
        xmax = xmax + 0.05 * xheight
        axes[par_idx, par_idx].set_xlim([xmin,xmax])

        axes[par_idx, par_idx].set_xlabel('')
        axes[par_idx, par_idx].set_ylabel('')
        
        if par_idx is 0:
            old_ymin, old_ymax = axes[par_idx, par_idx].get_ylim()
            old_yheight = old_ymax - old_ymin
            axes[par_idx, par_idx].set_yticks([old_ymin + 0.05*old_yheight, old_ymax - 0.05*old_yheight ])
            axes[par_idx, par_idx].set_yticklabels([str(param_lims[par_idx][0]), str(param_lims[par_idx][1])])
        
        if par_idx is not (param_N - 1):
            axes[par_idx, par_idx].set_xticklabels([])
        if par_idx is not 0: 
            axes[par_idx, par_idx].set_yticklabels([])
        
    for par_idx, axis_label in enumerate(axes_labels):

        axes[-1, par_idx].set_xlabel(axis_label)
        axes[par_idx, 0].set_ylabel(axis_label)
        
    cbar_cmap = sns.cubehelix_palette(8, start=.5, rot=-.75, reverse=True, as_cmap=True)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.6])
    norm = colors.Normalize(vmin=color_cutoff, vmax=cmax, clip=False)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cbar_cmap), cax=cbar_ax, extend='min')
    cbar_ax.get_yaxis().labelpad = -100
    cbar_ax.set_ylabel('log$_{10}$(Likelihood)', rotation=90)
    # cbar.cmap.set_under('k')
        
    # plt.tight_layout()
    plt.savefig(figs_directory + 'param_pair_grid_logp.png')
    plt.close()
    
    return
    

def plot_param_pair_grid_success(dream_success_df, param_names, param_lims, axes_labels, figs_directory):
    
    font_size = set_figs_font_settings()
    
    plt.rc('xtick', labelsize=font_size['SMALLEST_SIZE'])    # fontsize of the tick labels
    plt.rc('ytick', labelsize=font_size['SMALLEST_SIZE'])    # fontsize of the tick labels
    plt.rc('axes', labelsize=font_size['SMALL_SIZE'])    # fontsize of the x and y labels
    
    param_N = len(param_names)
    
    df = dream_success_df
   
    """ find settings for the colorbar """
    success_min = df['success_proportion'].min()
    success_max = df['success_proportion'].max()
    
    cmin = 0
    cmax = 1
    
    
    fig, axes = plt.subplots(nrows=param_N, ncols=param_N, figsize=(14,9))
    
    param_pair_indices = combinations(range(param_N), 2)
    for param_pair_index in param_pair_indices:

        par1_idx = param_pair_index[0]
        par2_idx = param_pair_index[1]

        par1 = param_names[param_pair_index[0]]
        par2 = param_names[param_pair_index[1]]

        lims = (param_lims[par1_idx], param_lims[par2_idx])
        xmin = lims[1][0]
        xmax = lims[1][1]
        ymin = lims[0][0]
        ymax = lims[0][1]

        xheight = xmax - xmin
        xmin = xmin - 0.05*xheight
        xmax = xmax + 0.05*xheight
        yheight = ymax - ymin
        ymin = ymin - 0.05*yheight
        ymax = ymax + 0.05*yheight

        cmap = sns.cubehelix_palette(8, start=2.8, rot=-.1, dark=0.8, light=0.2, reverse=False, as_cmap=True)
        kde_cmap = sns.cubehelix_palette(8, start=2.8, rot=-.1, reverse=True, as_cmap=True)
        distplot_palette = sns.cubehelix_palette(8, start=2.8, rot=-.1, light=0.5, reverse=True, as_cmap=False)
        axes[par1_idx, par2_idx].scatter(x=df[par2],y=df[par1], c=df['success_proportion'], cmap=cmap, vmin=cmin, vmax=cmax, s=2)
        axes[par1_idx, par2_idx].set_xlim([xmin, xmax])
        axes[par1_idx, par2_idx].set_ylim([ymin, ymax])

        axes[par1_idx, par2_idx].set_xticklabels([])
        axes[par1_idx, par2_idx].set_yticklabels([])

        sns.kdeplot(ax=axes[par2_idx, par1_idx], data=df[par1], data2=df[par2], shade=False, cmap=kde_cmap)
        axes[par2_idx, par1_idx].set_xlim([ymin, ymax])
        axes[par2_idx, par1_idx].set_ylim([xmin, xmax])
        
        axes[par2_idx, par1_idx].set_xlabel('')
        axes[par2_idx, par1_idx].set_ylabel('')
        
        if par2_idx is not (param_N - 1):
            axes[par2_idx, par1_idx].set_xticklabels([])
        if par1_idx is not 0:
            axes[par2_idx, par1_idx].set_yticklabels([])

    for par_idx, par in enumerate(param_names):

        sns.distplot( df[par], ax=axes[par_idx, par_idx], color=distplot_palette[5])
        xmin, xmax = param_lims[par_idx]
        xheight = xmax - xmin
        xmin = xmin - 0.05 * xheight
        xmax = xmax + 0.05 * xheight
        axes[par_idx, par_idx].set_xlim([xmin,xmax])

        axes[par_idx, par_idx].set_xlabel('')
        axes[par_idx, par_idx].set_ylabel('')
        
        if par_idx is 0:
            old_ymin, old_ymax = axes[par_idx, par_idx].get_ylim()
            old_yheight = old_ymax - old_ymin
            axes[par_idx, par_idx].set_yticks([old_ymin + 0.05*old_yheight, old_ymax - 0.05*old_yheight ])
            axes[par_idx, par_idx].set_yticklabels([str(param_lims[par_idx][0]), str(param_lims[par_idx][1])])
        
        if par_idx is not (param_N - 1):
            axes[par_idx, par_idx].set_xticklabels([])
        if par_idx is not 0: 
            axes[par_idx, par_idx].set_yticklabels([])
        
    for par_idx, axis_label in enumerate(axes_labels):

        axes[-1, par_idx].set_xlabel(axis_label)
        axes[par_idx, 0].set_ylabel(axis_label)
        
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.6])
    norm = colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax)
    cbar_ax.get_yaxis().labelpad = -80
    cbar_ax.set_ylabel('success rate', rotation=90)
        
    # plt.tight_layout()
    plt.savefig(figs_directory + 'param_pair_grid_success.png')
    plt.close()
    
    return
    
def logp_success_correlation(dream_success_df, param_names, model_color, figs_directory):
    
    font_size = set_figs_font_settings()
    
    df = dream_success_df
    
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logp"])
    df.sort_values(by='logp', ascending=False)
    
    logp_min = df.logp.min()
    
    fig = plt.figure(figsize=(8,5))
    ax1 = fig.add_subplot(211)
    ax1.scatter(x=df['logp'], y=df['success_proportion'], c=model_color, s=25, marker='|')
    ax1.set_xlim([logp_min * 1.05 , 0])
    ax1.set_ylim([-0.05,1.05])
    ax1.set_ylabel('success rate')
    ax1.set_xlabel('log$_{10}$(Likelihood)')
    
    df = df.reset_index()
    logp_cutoff = df.at[int(np.ceil(0.7 * len(df))),'logp']
    
    df_filt = df[df['logp'] > logp_cutoff]
    ax2 = fig.add_subplot(212)
    ax2.scatter(x=df_filt['logp'], y=df_filt['success_proportion'], c=model_color, s=25,  marker='|')
    ax2.set_xlim([logp_cutoff * 1.05 , 0])
    ax2.set_ylim([-0.05,1.05])
    ax2.set_ylabel('success rate')
    ax2.set_xlabel('log$_{10}$(Likelihood)')
    
    plt.tight_layout()
    plt.savefig(figs_directory + 'logp_vs_success.png')
    plt.close()
    
    
def plot_dist_from_all_chains(dream_out_df, param_names, param_lims, axes_labels, figs_directory):
    
    font_size = set_figs_font_settings()
    
    plt.rc('xtick', labelsize=font_size['SMALLEST_SIZE'])    # fontsize of the tick labels
    plt.rc('ytick', labelsize=font_size['SMALLEST_SIZE'])    # fontsize of the tick labels
    plt.rc('axes', labelsize=font_size['SMALL_SIZE'])    # fontsize of the x and y labels
    
    param_N = len(param_names)
    
    df = dream_out_df
    chain_N = int(df.chainID.max()) + 1
   
    """ find settings for the colorbar """
    logp_min = df['logp'].min()
    logp_max = df['logp'].max()

    cmin = np.floor(logp_min / 10) * 10
    cmax = np.ceil(logp_max / 10) * 10
    
    df = df.reset_index()
    
    color_cutoff = df.at[int(np.ceil(0.1 * len(df))),'logp']
    
    fig, axes = plt.subplots(nrows=chain_N, ncols=param_N, figsize=(14,9))

    distplot_palette =  sns.cubehelix_palette(8, start=.5, rot=-.75, light=0.5, reverse=True, as_cmap=False)
    # distplot_palette = sns.cubehelix_palette(param_N, light=0.5, start=0.5, reverse=True, as_cmap=False)
    for chain_idx in range(chain_N):
        df_temp = df[df['chainID'] == chain_idx]
        titleN = str(len(df_temp))
        for par_idx, par in enumerate(param_names):
            sns.distplot( df_temp[par], ax=axes[chain_idx, par_idx], color=distplot_palette[5])
            
            xmin, xmax = param_lims[par_idx]
            xheight = xmax - xmin
            xmin = xmin - 0.05 * xheight
            xmax = xmax + 0.05 * xheight
            axes[chain_idx, par_idx].set_xlim([xmin,xmax])

            axes[chain_idx, par_idx].set_xlabel('')
            axes[chain_idx, par_idx].set_ylabel('')
            axes[chain_idx, par_idx].set_yticklabels([])
        
            if chain_idx is not (chain_N - 1):
                axes[chain_idx, par_idx].set_xticklabels([])
            
            if par_idx is 0:
                axes[chain_idx, par_idx].set_ylabel('chain ' + str(chain_idx))
            if chain_idx is (chain_N - 1):
                axes[chain_idx, par_idx].set_xlabel(axes_labels[par_idx])
        
    
    fig.suptitle('Iterations per chain = ' + titleN)
    # plt.tight_layout()
    plt.savefig(figs_directory + 'compare_dist_all_chains.png')
    plt.close()
    
    return
    


def create_pyDREAM_figs(dream_out_df, dream_success_df, dream_params, param_names, param_lims, axes_labels, model_color, save_directory):
    
    figs_directory = save_directory + 'figs/'
    if not os.path.isdir(figs_directory):
        os.mkdir(figs_directory)
        
    # desired output
    plot_param_pair_grid_logp(dream_out_df, param_names, param_lims, axes_labels, figs_directory)
    plot_param_pair_grid_success(dream_success_df, param_names, param_lims, axes_labels, figs_directory)
    
    # checks
    plot_logp_over_time(dream_out_df, dream_params, figs_directory)
    plot_dist_from_all_chains(dream_out_df, param_names, param_lims, axes_labels, figs_directory)
    logp_success_correlation(dream_success_df, param_names, model_color, figs_directory)
    
def create_pyDREAM_figs_2models(dream_out_df, dream_success_df, dream_params, param_names, param_lims, axes_labels, model_A_color, model_B_color, save_directory):
    
    figs_directory = save_directory + 'figs/'
    if not os.path.isdir(figs_directory):
        os.mkdir(figs_directory)
        
    # desired output
    plot_param_pair_grid_logp(dream_out_df, param_names, param_lims, axes_labels, figs_directory)
    plot_param_pair_grid_success(dream_success_df, param_names, param_lims, axes_labels, figs_directory)
    
    # checks
    plot_logp_over_time(dream_out_df, dream_params, figs_directory)
    plot_dist_from_all_chains(dream_out_df, param_names, param_lims, axes_labels, figs_directory)
    # logp_success_correlation(dream_success_df, param_names, model_A_color, figs_directory)
    # logp_success_correlation(dream_success_df, param_names, model_B_color, figs_directory)
    
    