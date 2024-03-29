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
    
    return
    

    
def plot_logp_over_time(dream_out_df, dream_params, figs_directory):
    
    font_size = set_figs_font_settings()
    
    df = dream_out_df    
    nchains = dream_params['nchains']
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(211)
    for chain in range(nchains):
        temp_df = df[df['chainID']==chain]
        temp_df = temp_df.reset_index()
        temp_df.plot(ax=ax, y='logp', kind='line', label='chain ' + str(chain))
    ax.set_xlabel('Iteration')
    ax.set_ylabel('log$_{10}$(Likelihood)')
    
    df = df.reset_index()
    min_logp = df['logp'].min()
    max_logp = df['logp'].max()
    
    if not np.isinf(min_logp):
        logp_cutoff = ( 0.9 * (max_logp - min_logp) ) + min_logp
    else:
        logp_cutoff = 0.9 * max_logp 
    # logp_cutoff = df.at[int(np.ceil(0.7 * len(df))),'logp']
    
    ax_cutoff = fig.add_subplot(212)
    for chain in range(nchains):
        temp_df = df[df['chainID']==chain]
        temp_df = temp_df.reset_index()
        df_filt = temp_df[temp_df['logp'] > logp_cutoff]
        df_filt.plot(ax=ax_cutoff, y='logp', kind='line', label='chain ' + str(chain))
    ax_cutoff.set_xlabel('Iteration')
    ax_cutoff.set_ylabel('log$_{10}$(Likelihood)')
    ax_cutoff.get_legend().remove()
    
    plt.tight_layout()
    plt.savefig(figs_directory + 'logp_over_time.png', dpi=600)
    plt.close()
    
    return
    
    
def plot_param_pair_grid_logp(dream_out_df, param_names, param_lims, axes_labels, model_color, figs_directory):
    
    font_size = set_figs_font_settings()
    
    plt.rc('xtick', labelsize=font_size['SMALLEST_SIZE'])    # fontsize of the tick labels
    plt.rc('ytick', labelsize=font_size['SMALLEST_SIZE'])    # fontsize of the tick labels
    plt.rc('axes', labelsize=font_size['SMALL_SIZE'])    # fontsize of the x and y labels
    
    param_N = len(param_names)

    df = dream_out_df
    # df = df.drop_duplicates()
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
        
        cmap_length = 25
        
        kde_cmap = sns.cubehelix_palette(cmap_length, start=.5, rot=-.75, reverse=True, as_cmap=False)

        scatter_cmap = sns.cubehelix_palette(cmap_length, start=.5, rot=-.75, reverse=True, as_cmap=True)
        axes[par1_idx, par2_idx].scatter(x=df[par2],y=df[par1], c=df['logp'], cmap=scatter_cmap, vmin=color_cutoff, vmax=cmax, s=2)
        
        axes[par1_idx, par2_idx].set_facecolor(kde_cmap[0])
        
        axes[par1_idx, par2_idx].set_xlim([xmin, xmax])
        axes[par1_idx, par2_idx].set_ylim([ymin, ymax])

        axes[par1_idx, par2_idx].set_xticklabels([])
        axes[par1_idx, par2_idx].set_yticklabels([])

        # kde_cmap =  sns.cubehelix_palette(8, start=.5, rot=-.75, reverse=True, as_cmap=True)
        # kde_cmap = sns.cubehelix_palette(cmap_length, start=.5, rot=-.75, reverse=True, as_cmap=False)
        
        kde_idx = -3
        hist_idx = -3
        hist_idx = 5
        hist_line_idx = 5
        
        sns.kdeplot(ax=axes[par2_idx, par1_idx], x=df[par1], y=df[par2], shade=True, color=kde_cmap[hist_idx])
        axes[par2_idx, par1_idx].set_xlim([ymin, ymax])
        axes[par2_idx, par1_idx].set_ylim([xmin, xmax])
        
        axes[par2_idx, par1_idx].set_xlabel('')
        axes[par2_idx, par1_idx].set_ylabel('')
        
        if par2_idx is not (param_N - 1):
            axes[par2_idx, par1_idx].set_xticklabels([])
        if par1_idx is not 0:
            axes[par2_idx, par1_idx].set_yticklabels([])

    # distplot_palette =  sns.cubehelix_palette(8, start=.5, rot=-.75, light=0.5, reverse=True, as_cmap=False)
    # distplot_palette = sns.cubehelix_palette(param_N, light=0.5, start=0.5, reverse=True, as_cmap=False)
    for par_idx, par in enumerate(param_names):

        sns.histplot( data=df[par], ax=axes[par_idx, par_idx], kde=True, color=model_color, element='step')
        # sns.kdeplot( data=df[par], ax=axes[par_idx, par_idx], color=kde_cmap[hist_line_idx])
        
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
    plt.savefig(figs_directory + 'param_pair_grid_logp.png', dpi=600)
    plt.close()
    
    return
    

def plot_param_pair_grid_success(dream_success_df, param_names, param_lims, axes_labels, figs_directory):
    
    font_size = set_figs_font_settings()
    
    plt.rc('xtick', labelsize=font_size['SMALLEST_SIZE'])    # fontsize of the tick labels
    plt.rc('ytick', labelsize=font_size['SMALLEST_SIZE'])    # fontsize of the tick labels
    plt.rc('axes', labelsize=font_size['SMALL_SIZE'])    # fontsize of the x and y labels
    
    param_N = len(param_names)
    
    df = dream_success_df
    df = df.dropna(subset=['success_proportion'])
    max_success = np.max(df['success_proportion'])
    df_max_success = df[df['success_proportion']==max_success]
    
   
    """ find settings for the colorbar """
    success_min = df['success_proportion'].min()
    success_max = df['success_proportion'].max()
    
    cmin = 0
    cmax = 1
    
    """ find settings for the logp colorbar """
    logp_min = df_max_success['logp'].min()
    logp_max = df_max_success['logp'].max()

    cmin_logp = np.floor(logp_min / 10) * 10
    cmax_logp = np.ceil(logp_max / 10) * 10
    
    cbar_cmap_logp = sns.cubehelix_palette(8, start=.5, rot=-.75, reverse=True, as_cmap=True)
    cmap = sns.cubehelix_palette(8, start=2.8, rot=-.1, dark=0.8, light=0.2, reverse=False, as_cmap=True)
    kde_cmap = sns.cubehelix_palette(8, start=2.8, rot=-.1, reverse=True, as_cmap=True)
    distplot_palette = sns.cubehelix_palette(8, start=2.8, rot=-.1, light=0.5, reverse=True, as_cmap=False)
    
    fig, axes = plt.subplots(nrows=param_N, ncols=param_N, figsize=(14,9))
    fig_max_success, axes_max_success = plt.subplots(nrows=param_N, ncols=param_N, figsize=(14,9))
    
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

        
        axes[par1_idx, par2_idx].scatter(x=df[par2],y=df[par1], c=df['success_proportion'], cmap=cmap, vmin=cmin, vmax=cmax, s=2)
        axes_max_success[par1_idx, par2_idx].scatter(x=df_max_success[par2],y=df_max_success[par1], c=df_max_success['logp'], cmap=cbar_cmap_logp, vmin=cmin_logp, vmax=cmax_logp, s=2)
        
        axes[par1_idx, par2_idx].set_xlim([xmin, xmax])
        axes[par1_idx, par2_idx].set_ylim([ymin, ymax])
        axes_max_success[par1_idx, par2_idx].set_xlim([xmin, xmax])
        axes_max_success[par1_idx, par2_idx].set_ylim([ymin, ymax])

        axes[par1_idx, par2_idx].set_xticklabels([])
        axes[par1_idx, par2_idx].set_yticklabels([])
        axes_max_success[par1_idx, par2_idx].set_xticklabels([])
        axes_max_success[par1_idx, par2_idx].set_yticklabels([])

        sns.kdeplot(ax=axes[par2_idx, par1_idx], data=df[par1], data2=df[par2], shade=False, cmap=kde_cmap)
        sns.kdeplot(ax=axes_max_success[par2_idx, par1_idx], data=df_max_success[par1], data2=df_max_success[par2], shade=False)
        axes[par2_idx, par1_idx].set_xlim([ymin, ymax])
        axes[par2_idx, par1_idx].set_ylim([xmin, xmax])
        axes_max_success[par2_idx, par1_idx].set_xlim([ymin, ymax])
        axes_max_success[par2_idx, par1_idx].set_ylim([xmin, xmax])
        
        axes[par2_idx, par1_idx].set_xlabel('')
        axes[par2_idx, par1_idx].set_ylabel('')
        axes_max_success[par2_idx, par1_idx].set_xlabel('')
        axes_max_success[par2_idx, par1_idx].set_ylabel('')
        
        if par2_idx is not (param_N - 1):
            axes[par2_idx, par1_idx].set_xticklabels([])
            axes_max_success[par2_idx, par1_idx].set_xticklabels([])
        if par1_idx is not 0:
            axes[par2_idx, par1_idx].set_yticklabels([])
            axes_max_success[par2_idx, par1_idx].set_yticklabels([])

    for par_idx, par in enumerate(param_names):

        sns.distplot( df[par], ax=axes[par_idx, par_idx], color=distplot_palette[5])
        sns.distplot( df_max_success[par], ax=axes_max_success[par_idx, par_idx])
        xmin, xmax = param_lims[par_idx]
        xheight = xmax - xmin
        xmin = xmin - 0.05 * xheight
        xmax = xmax + 0.05 * xheight
        axes[par_idx, par_idx].set_xlim([xmin,xmax])
        axes_max_success[par_idx, par_idx].set_xlim([xmin,xmax])

        axes[par_idx, par_idx].set_xlabel('')
        axes[par_idx, par_idx].set_ylabel('')
        axes_max_success[par_idx, par_idx].set_xlabel('')
        axes_max_success[par_idx, par_idx].set_ylabel('')
        
        if par_idx is 0:
            old_ymin, old_ymax = axes[par_idx, par_idx].get_ylim()
            old_yheight = old_ymax - old_ymin
            axes[par_idx, par_idx].set_yticks([old_ymin + 0.05*old_yheight, old_ymax - 0.05*old_yheight ])
            axes[par_idx, par_idx].set_yticklabels([str(param_lims[par_idx][0]), str(param_lims[par_idx][1])])
            axes_max_success[par_idx, par_idx].set_yticks([old_ymin + 0.05*old_yheight, old_ymax - 0.05*old_yheight ])
            axes_max_success[par_idx, par_idx].set_yticklabels([str(param_lims[par_idx][0]), str(param_lims[par_idx][1])])
        
        if par_idx is not (param_N - 1):
            axes[par_idx, par_idx].set_xticklabels([])
            axes_max_success[par_idx, par_idx].set_xticklabels([])
        if par_idx is not 0: 
            axes[par_idx, par_idx].set_yticklabels([])
            axes_max_success[par_idx, par_idx].set_yticklabels([])
        
    for par_idx, axis_label in enumerate(axes_labels):

        axes[-1, par_idx].set_xlabel(axis_label)
        axes[par_idx, 0].set_ylabel(axis_label)
        axes_max_success[-1, par_idx].set_xlabel(axis_label)
        axes_max_success[par_idx, 0].set_ylabel(axis_label)
        
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.6])
    norm = colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax)
    cbar_ax.get_yaxis().labelpad = -80
    cbar_ax.set_ylabel('success rate', rotation=90)
    
    
    fig_max_success.subplots_adjust(right=0.8)
    cbar_ax_max_success = fig_max_success.add_axes([0.85, 0.25, 0.03, 0.6])
    norm = colors.Normalize(vmin=cmin_logp, vmax=cmax_logp, clip=False)
    cbar_max_success = fig_max_success.colorbar(cm.ScalarMappable(norm=norm, cmap=cbar_cmap_logp), cax=cbar_ax_max_success, extend='min')
    cbar_ax_max_success.get_yaxis().labelpad = -100
    cbar_ax_max_success.set_ylabel('log$_{10}$(Likelihood)', rotation=90)
        
    # plt.tight_layout()
    fig.savefig(figs_directory + 'param_pair_grid_success.png', dpi=600)
    fig_max_success.savefig(figs_directory + 'param_pair_grid_max_success.png', dpi=600)
    
    plt.close(fig)
    plt.close(fig_max_success)
    
    return
    
def logp_success_correlation(dream_success_df, success_prop_col_name_prefix, model_color, figs_directory):
    
    print(model_color)
    if model_color == 'C8':
        model_letter = 'A'
    elif model_color == 'C1':
        model_letter = 'B'
    else:
        model_letter = ''
        print('Model color value unknown.')
    
    font_size = set_figs_font_settings()
    
    df_raw = dream_success_df
    
    len_raw = len(df_raw)
    
    df = df_raw.replace([np.inf, -np.inf], np.nan).dropna(subset=["logp"])
    df.sort_values(by='logp', ascending=False)
    
    len_no_inf = len(df)
    
    # print(len_raw, len_no_inf)
    
    col_name = success_prop_col_name_prefix + 'success_proportion'
    
    logp_min = df.logp.min()
    
    print(df[col_name].value_counts())
    
    fig = plt.figure(figsize=(8,6))
    
    logp_shade_min = -60
    logp_shade_max = -35
    
    ax1 = fig.add_subplot(211)
    ax1.scatter(x=df['logp'], y=df[col_name], c=model_color, s=25, marker='|')
    ax1.set_xlim([logp_min * 1.05 , 0])
    ax1.set_ylim([-0.05,1.05])
    x = np.linspace(logp_shade_min , logp_shade_max, 101)
    ax1.fill_between(x, -0.05, y2=1.05, color='lightgray', alpha=0.4)
    ax1.set_ylabel('success rate')
    ax1.set_xlabel('log$_{10}$(Likelihood)')
    ax1.set_title('Model ' + model_letter, fontsize=24, fontweight='bold')
    # ax1.set_title('Correlation between likelihood and success rate\nof parameter values found\nwith Bayesian parameter inference\n')
    
    df = df.reset_index()
    # logp_cutoff = df.at[int(np.ceil(0.7 * len(df))),'logp']
    logp_cutoff = logp_shade_min
    
    df_filt = df[df['logp'] > logp_cutoff]
    ax2 = fig.add_subplot(212)
    ax2.scatter(x=df_filt['logp'], y=df_filt[col_name], c=model_color, s=25,  marker='|', label='one set of parameters')
    ax2.set_xlim([logp_cutoff - 1 , logp_shade_max + 1])
    ax2.fill_between(x, -0.05, y2=1.05, color='lightgray', alpha=0.4)
    ax2.set_ylim([-0.05,1.05])
    ax2.set_ylabel('success rate')
    ax2.set_xlabel('log$_{10}$(Likelihood)')
    ax2.legend(loc='lower right')
    
    plt.tight_layout()
    plt.savefig(figs_directory + success_prop_col_name_prefix + 'logp_vs_success.png', dpi=600)
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
    
    fig, axes = plt.subplots(nrows=chain_N + 1, ncols=param_N, figsize=(15,7))

    plot_palette =  sns.cubehelix_palette(8, start=.5, rot=-.75, light=0.5, reverse=True, as_cmap=False)
    # plot_palette = sns.cubehelix_palette(param_N, light=0.5, start=0.5, reverse=True, as_cmap=False)
    ymin = 10
    ymax = 0
    
    for chain_idx in range(chain_N):
        df_temp = df[df['chainID'] == chain_idx]
        titleN = str(len(df_temp))
        for par_idx, par in enumerate(param_names):
            # sns.histplot( df_temp[par], ax=axes[chain_idx, par_idx], color=plot_palette[5], kde=True, fill=True)
            sns.histplot( df_temp[par], ax=axes[chain_idx, par_idx], kde=True)
            # sns.kdeplot( df_temp[par], ax=axes[chain_idx, par_idx], color=plot_palette[5])
            
            temp_ymin, temp_ymax = axes[chain_idx, par_idx].get_ylim()
            if temp_ymin < ymin:
                ymin = temp_ymin
            if temp_ymax > ymax:
                ymax = temp_ymax
                
    for chain_idx in range(chain_N):
        df_temp = df[df['chainID'] == chain_idx]
        titleN = str(len(df_temp))
        for par_idx, par in enumerate(param_names):
            # sns.histplot( df_temp[par], ax=axes[chain_idx, par_idx], color=plot_palette[5], kde=True)
            # sns.kdeplot( df_temp[par], ax=axes[chain_idx, par_idx], color=plot_palette[5])
            sns.histplot( df_temp[par], ax=axes[chain_idx, par_idx], kde=True)
            xmin, xmax = param_lims[par_idx]
            xheight = xmax - xmin
            xmin = xmin - 0.05 * xheight
            xmax = xmax + 0.05 * xheight
            axes[chain_idx, par_idx].set_xlim([xmin,xmax])

            axes[chain_idx, par_idx].set_xlabel('')
            axes[chain_idx, par_idx].set_ylabel('')
            # axes[chain_idx, par_idx].set_ylim([ymin,ymax])
            # axes[chain_idx, par_idx].set_yticklabels([])
            
            
        
            # if chain_idx is not (chain_N-1):
            axes[chain_idx, par_idx].set_xticklabels([])
            
            if par_idx is 0:
                axes[chain_idx, par_idx].set_ylabel('chain ' + str(chain_idx))
            # if chain_idx is not (chain_N - 1):
                # axes[chain_idx, par_idx].set_xlabel(axes_labels[par_idx])
            
            axes[chain_idx, par_idx].tick_params(axis='y', direction='in', labelsize=6)
                
    for par_idx, par in enumerate(param_names):
        # sns.histplot( df_temp[par], ax=axes[chain_idx, par_idx], color=plot_palette[5], kde=True)
        # sns.kdeplot( df_temp[par], ax=axes[chain_idx, par_idx], color=plot_palette[5])
        sns.histplot( df[par], ax=axes[chain_N, par_idx], kde=True)
        xmin, xmax = param_lims[par_idx]
        xheight = xmax - xmin
        xmin = xmin - 0.05 * xheight
        xmax = xmax + 0.05 * xheight
        axes[chain_N, par_idx].set_xlim([xmin,xmax])

        axes[chain_N, par_idx].set_xlabel('')
        axes[chain_N, par_idx].set_ylabel('')
        # axes[chain_N, par_idx].set_ylim([ymin,ymax])
        # axes[chain_N, par_idx].set_yticklabels([])
        axes[chain_N, 0].set_ylabel('All chains')
        axes[chain_N, par_idx].set_xlabel(axes_labels[par_idx])
        axes[chain_N, par_idx].tick_params(axis='y', direction='in', labelsize=6)
        
    
    fig.suptitle('Iterations per chain = ' + titleN)
    # plt.tight_layout()
    plt.savefig(figs_directory + 'compare_dist_all_chains.png', dpi=600)
    plt.close()
    
    return
    
    
def switch_models_A_and_B_find_bead(param_names, save_directory):

    dream_out_df = pd.read_csv(save_directory + 'dream_out.tsv', sep='\t', header=0, index_col=0)
    dream_success_df = pd.read_csv(save_directory + 'top_params.tsv', sep='\t', header=0, index_col=0)
        
    df_out_switch_AB_and_lower = {
        "threshold$^A$": "threshold$^b$",
        "$b_B^A$": "$b_B^b$",
        "$b_V^A$": "$b_V^b$",
        "threshold$^B$": "threshold$^a$",
        "$b_B^B$": "$b_B^a$",
        "$b_V^B$": "$b_V^a$"
    }
    dream_out_df = dream_out_df.rename(columns=df_out_switch_AB_and_lower)
    
    df_out_switch_upper = {
        "threshold$^b$": "threshold$^B$",
        "$b_B^b$": "$b_B^B$",
        "$b_V^b$": "$b_V^B$",
        "threshold$^a$": "threshold$^A$",
        "$b_B^a$": "$b_B^A$",
        "$b_V^a$": "$b_V^A$"
    }
    dream_out_df = dream_out_df.rename(columns=df_out_switch_upper)
    
    col_N = len(dream_success_df.columns)
    param_N = len(param_names)
    
    embryo_N = 15
    
    df_success_switch_AB_and_lower = {
        "threshold$^A$": "threshold$^b$",
        "$b_B^A$": "$b_B^b$",
        "$b_V^A$": "$b_V^b$",
        "threshold$^B$": "threshold$^a$",
        "$b_B^B$": "$b_B^a$",
        "$b_V^B$": "$b_V^a$",
        "A_success_proportion": "b_success_proportion",
        "B_success_proportion": "a_success_proportion"
    }
    for idx in range(embryo_N):
        if 'A_' + str(idx + 1) in dream_success_df.columns:
            df_success_switch_AB_and_lower['A_' + str(idx + 1)] = 'b_' + str(idx + 1)
            df_success_switch_AB_and_lower['B_' + str(idx + 1)] = 'a_' + str(idx + 1)
    dream_success_df = dream_success_df.rename(columns=df_success_switch_AB_and_lower)
    
    df_success_switch_upper = {
        "threshold$^b$": "threshold$^B$",
        "$b_B^b$": "$b_B^B$",
        "$b_V^b$": "$b_V^B$",
        "threshold$^a$": "threshold$^A$",
        "$b_B^a$": "$b_B^A$",
        "$b_V^a$": "$b_V^A$",
        "b_success_proportion": "B_success_proportion",
        "a_success_proportion": "A_success_proportion"
    }
    for idx in range(embryo_N):
        if 'b_' + str(idx + 1) in dream_success_df.columns:
            df_success_switch_upper['b_' + str(idx + 1)] = 'B_' + str(idx + 1)
            df_success_switch_upper['a_' + str(idx + 1)] = 'A_' + str(idx + 1)
    dream_success_df = dream_success_df.rename(columns=df_success_switch_upper)
    
    dream_out_df.to_csv(save_directory + 'dream_out_switched_AB.tsv', sep='\t')
    dream_success_df.to_csv(save_directory + 'top_params_switched_AB.tsv', sep='\t')
    
    # param_names = ['threshold$^A$', '$b_B^A$', '$b_V^A$', 'threshold$^B$', '$b_B^B$', '$b_V^B$', 'n'] + param_names[7:]
    
    return param_names
    

def create_pyDREAM_figs(dream_params, param_names, param_lims, axes_labels, model_color, save_directory):
    
    dream_out_df = pd.read_csv(save_directory + 'dream_out.tsv', sep='\t', index_col=0, header=0)
    dream_success_df = pd.read_csv(save_directory + 'top_params.tsv', sep='\t', index_col=0, header=0)
    
    dream_success_df = dream_success_df.sort_values(by='logp', ascending=False)
    
    figs_directory = save_directory + 'figs/'
    if not os.path.isdir(figs_directory):
        os.mkdir(figs_directory)
        
    # desired output
    # plot_param_pair_grid_logp(dream_success_df, param_names, param_lims, axes_labels, model_color, figs_directory)
    # plot_param_pair_grid_success(dream_success_df, param_names, param_lims, axes_labels, figs_directory)
    
    # checks
    # plot_logp_over_time(dream_out_df, dream_params, figs_directory)
    # plot_dist_from_all_chains(dream_out_df, param_names, param_lims, axes_labels, figs_directory)
    logp_success_correlation(dream_success_df, '', model_color, figs_directory)
    
def create_pyDREAM_figs_2models(dream_params, param_names, param_lims, axes_labels, without_nbhd_color, with_nbhd_color, save_directory):
    
    param_names = switch_models_A_and_B_find_bead(param_names, save_directory)
    
    dream_out_df = pd.read_csv(save_directory + 'dream_out_switched_AB.tsv', sep='\t')
    dream_success_df = pd.read_csv(save_directory + 'top_params_switched_AB.tsv', sep='\t')
    
    figs_directory = save_directory + 'figs/'
    if not os.path.isdir(figs_directory):
        os.mkdir(figs_directory)
        
    param_names_with_nbhd = param_names[:5]
    param_names_without_nbhd = param_names[5:8]
    param_names_bead = param_names[8:]
    
    param_names_with_nbhd_bead = param_names_with_nbhd + param_names_bead
    param_names_without_nbhd_bead = param_names_without_nbhd + param_names_bead
        
    # desired output
    # plot_param_pair_grid_logp(dream_out_df, param_names, param_lims, axes_labels, figs_directory)
    # plot_param_pair_grid_success(dream_success_df, param_names, param_lims, axes_labels, figs_directory)
    
    # posterior distributions
    
    
    # checks
    plot_logp_over_time(dream_out_df, dream_params, figs_directory)
    plot_dist_from_all_chains(dream_out_df, param_names, param_lims, axes_labels, figs_directory)
    logp_success_correlation(dream_success_df, 'A_', without_nbhd_color, figs_directory)
    logp_success_correlation(dream_success_df, 'B_', with_nbhd_color, figs_directory)
    
    
    