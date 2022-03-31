import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.font_manager as font_manager
from matplotlib import rcParams
from scipy import stats
import copy
import os.path


def set_figs_font_settings():

    font_sizes = {
        'SMALLEST_SIZE' : 12,
        "SMALLER_SIZE": 14,
        'SMALL_SIZE' : 16,
        'MEDIUM_SIZE' : 20,
        'BIGGER_SIZE' : 24
    }

    # figure font parameters
    # rcParams['text.usetex'] = True
    # rcParams["text.latex.preamble"] = r"\usepackage{helvet}"
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    # rcParams['font.sans-serif'] = ['Clear Sans']
    plt.rc('font', size=font_sizes['SMALL_SIZE'])          # controls default text sizes
    plt.rc('axes', titlesize=font_sizes['MEDIUM_SIZE'])     # fontsize of the axes title
    plt.rc('axes', labelsize=font_sizes['MEDIUM_SIZE'])    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=font_sizes['SMALL_SIZE'])    # fontsize of the tick labels
    plt.rc('ytick', labelsize=font_sizes['SMALL_SIZE'])    # fontsize of the tick labels
    plt.rc('legend', fontsize=font_sizes['SMALL_SIZE'])    # legend fontsize

    return font_sizes
    
def set_figs_font_settings_poster():

    font_sizes = {
        'SMALLEST_SIZE' : 8,
        "SMALLER_SIZE": 10,
        'SMALL_SIZE' : 12,
        'MEDIUM_SIZE' : 14,
        'BIGGER_SIZE' : 16
    }

    # figure font parameters
    # rcParams['text.usetex'] = True
    # rcParams["text.latex.preamble"] = r"\usepackage{helvet}"
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    # rcParams['font.sans-serif'] = ['Clear Sans']
    plt.rc('font', size=font_sizes['SMALL_SIZE'])          # controls default text sizes
    plt.rc('axes', titlesize=font_sizes['MEDIUM_SIZE'])     # fontsize of the axes title
    plt.rc('axes', labelsize=font_sizes['MEDIUM_SIZE'])    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=font_sizes['SMALL_SIZE'])    # fontsize of the tick labels
    plt.rc('ytick', labelsize=font_sizes['SMALL_SIZE'])    # fontsize of the tick labels
    plt.rc('legend', fontsize=font_sizes['SMALL_SIZE'])    # legend fontsize

    return font_sizes
    
            
def create_presentation_fig_arrays(list_of_embryos, **kwargs):
    
    embryoN = len(list_of_embryos)
    noc = list_of_embryos[0].number_of_cells
    
    model_values = np.ndarray((len(list_of_embryos), noc), dtype=float)
    model_ylim = np.ndarray((len(list_of_embryos), 2), dtype=float)
    
    for emb_idx in range(embryoN):
        model_values[emb_idx,:] = list_of_embryos[emb_idx].model_value
        model_ylim[emb_idx,0] = list_of_embryos[emb_idx].plot_model_ylim[0]
        model_ylim[emb_idx,1] = list_of_embryos[emb_idx].plot_model_ylim[1]
        
    return model_values, model_ylim
    
    

def plot_proteins(ax_pro, embryo, protein_settings):
    
    noc = embryo.number_of_cells
    
    font_sizes = protein_settings['font_sizes']
    rcParams['font.sans-serif'] = [protein_settings['font']]
    if protein_settings['view'] is 'anterior':
        xticklabel = ['post.','ant.','post.']
        roll_idx = 0
    elif protein_settings['view'] is 'posterior':
        xticklabel = ['ant.','post.','ant.']
        roll_idx = int(noc/2)
    else:
        print('protein_settings["view"] must be "anterior" or "posterior".' )
    
    ax_pro.set_xlim([-1, noc])
    ax_pro.set_xticks(np.linspace(0, noc, num=3))
    
    if protein_settings["x_label"]:
        ax_pro.set_xticklabels(xticklabel)
        ax_pro.set_xlabel('Cell location')
    else:
        ax_pro.set_xticklabels([])

    ''' proteins '''
    plot_inducer = np.roll(embryo.inducer.conc, roll_idx)
    plot_inhibitor = np.roll(embryo.inhibitor.conc,roll_idx)
    plot_target = np.roll(embryo.target,roll_idx)
    
    ax_pro.set_title(embryo.name, fontsize=font_sizes['BIGGER_SIZE'])
    ax_pro.plot(range(0,noc), plot_inducer, linewidth=2, marker=None, color='C0', markersize = 1, label='Inducer')
    ax_pro.plot(range(0,noc), plot_inhibitor, linewidth=2, marker=None, color='C3', markersize = 1, label="Inhibitor")
    ax_pro.set_ylim([0,1.3])
    ax_pro.set_ylabel('Protein conc.')
    
    if protein_settings["legend"]:
        ax_pro.legend(loc='upper right')
        
    if protein_settings["include_target"]:
    
        # # With bar and text below
        # [protein_ymin, protein_ymax] = axs[0].get_ylim()
        # fig_height =  protein_ymax - protein_ymin
        # fig_height_large = fig_height*1.3
        # new_protein_ymin = protein_ymax - fig_height_large
        # axs[0].set_ylim([new_protein_ymin, protein_ymax])
        # bar_height = fig_height * 0.1
        # bar_bottom = new_protein_ymin + 0.18 * fig_height
        # bar_top = bar_bottom + bar_height
        # text_yloc = new_protein_ymin + 0.05*fig_height_large
    
        target_color = 'C4'  # purple
        gray_where = np.array([i == -1 for i in plot_target])
        streak_where = np.array([i == 1 for i in plot_target])
        
        # With bar and text below
        [protein_ymin, protein_ymax] = ax_pro.get_ylim()
        fig_height =  protein_ymax - protein_ymin
        fig_height_large = fig_height*1.4
        new_protein_ymin = protein_ymax - fig_height_large
        ax_pro.set_ylim([new_protein_ymin, protein_ymax])
        bar_height = fig_height_large * 0.11
        bar_bottom = new_protein_ymin + 0.15*fig_height_large
        bar_top = bar_bottom + bar_height
        text_yloc = new_protein_ymin + 0.045*fig_height_large
        ax_pro.fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
        ax_pro.fill_between(np.arange(noc), bar_bottom, bar_top, where=streak_where, facecolor=target_color, edgecolor=target_color, alpha=1, step='mid', linewidth=1)
        ax_pro.text(400, text_yloc, 'Target streak', backgroundcolor='lightgray', color=target_color, fontsize=font_sizes['SMALL_SIZE'], fontweight='bold')
        
        # get rid of the frame
        ax_pro.spines['top'].set_visible(False)
        ax_pro.spines['right'].set_visible(False)
        ax_pro.spines['bottom'].set_visible(False)
        ax_pro.spines['left'].set_linewidth(1.5)
        ax_pro.spines['left'].set_bounds([0,protein_ymax])
        ax_pro.axhline(linewidth=1.5, color='k')
        
        return ax_pro

def plot_model(ax_model, embryo, model, model_values, model_ylim, model_settings):
    
    noc = embryo.number_of_cells
    
    font_sizes = model_settings['font_sizes']
    rcParams['font.sans-serif'] = [model_settings['font']]
    if model_settings['view'] is 'anterior':
        xticklabel = ['post.','ant.','post.']
        roll_idx = 0
    elif model_settings['view'] is 'posterior':
        xticklabel = ['ant.','post.','ant.']
        roll_idx = int(noc/2)
    else:
        print('model_settings["view"] must be "anterior" or "posterior".' )
    
    ax_model.set_xlim([-1, noc])
    ax_model.set_xticks(np.linspace(0, noc, num=3))
    
    if model_settings["x_label"]:
        ax_model.set_xticklabels(xticklabel)
        ax_model.set_xlabel('Cell location')
    else:
        ax_model.set_xticklabels([])

    plot_model = np.roll(model_values[:], roll_idx)
    ax_model.set_ylabel(model_settings['y_label_str'])
            
    ax_model.plot(range(0,noc), plot_model, linewidth=2, marker=None, color=model.plot_color, markersize = 1, label=model.label)
    ax_model.plot(range(0,noc), [model.threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)

    # get rid of the frame
    ax_model.spines['top'].set_visible(False)
    ax_model.spines['right'].set_visible(False)
    ax_model.spines['bottom'].set_visible(False)
    ax_model.spines['left'].set_linewidth(1.5)

    # model bar params
    fig_height_increase_multiplier = 1.4
    bar_height_proportion = 0.16
    bar_bottom_proportion = 0.19

    # With bar and text below
    [model_ymin, model_ymax] = [model_ylim[0], model_ylim[1]]
    fig_height =  model_ymax - model_ymin
    fig_height_large = fig_height * fig_height_increase_multiplier
    new_model_ymin = model_ymax - fig_height_large
    ax_model.set_ylim([new_model_ymin, model_ymax])
    bar_height = fig_height * bar_height_proportion
    bar_bottom = new_model_ymin + bar_bottom_proportion * fig_height
    bar_top = bar_bottom + bar_height
    text_yloc = new_model_ymin + 0.04 * fig_height_large

    brachyury_color = model.plot_color  # 'C4'  # purple
    ax_model.fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < model.threshold , facecolor='lightgray', alpha=1, step='mid')
    ax_model.fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= model.threshold ,facecolor=brachyury_color, edgecolor=brachyury_color, alpha=1, step='mid', linewidth=1)
    font_props = font_manager.FontProperties(size=font_sizes['SMALL_SIZE'], weight='semibold')
    predicted_A = ax_model.text(345, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontproperties=font_props)
            
    # if emb_idx in inset_axes:
    #
    #     # inset axes....
    #     axins = axs[1].inset_axes([0.02, 0.6, 0.35, 0.35])
    #     axins.plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[0].plot_color, markersize = 1)
    #     axins.plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
    #     # axins.imshow(Z2, extent=extent, interpolation="nearest", origin="lower")
    #     # sub region of the original image
    #     x1, x2, y1, y2 = 230, 370, 0.2, 0.35
    #     axins.set_xlim(x1, x2)
    #     axins.set_ylim(y1, y2)
    #     axins.set_xticklabels('')
    #     axins.set_yticklabels('')
    #
    #     axs[1].indicate_inset_zoom(axins)
    
    if model_settings["legend"]:
        ax_model.legend(loc='upper right')
        
    return ax_model

def set_up_protein_fig(embryo):
    
    font_sizes = set_figs_font_settings() 
    
    fig_pro = plt.figure()
    ax_pro = fig_pro.add_subplot(111)
    
    protein_settings = {
        "font" : 'Arial',
        "font_sizes" : font_sizes,
        "view" : 'anterior',
        "legend" : True,
        "include_target" : True,
        "x_label" : True
    }
    
    ax_pro = plot_proteins(ax_pro, embryo, protein_settings)
    
    fig_pro.tight_layout()
    
    return fig_pro
    
def set_up_fig_trio(embryo, models, model_values, model_ylim):
    
    font_sizes = set_figs_font_settings() 
    
    fig = plt.figure(figsize=(6,10))
    ax_pro = fig.add_subplot(311)
    ax_model_A = fig.add_subplot(312)
    ax_model_B = fig.add_subplot(313)
    
    global_settings = {
        "font" : 'Arial',
        "font_sizes" : font_sizes,
        "view" : 'anterior'
    }
    
    protein_settings = {
        "font" : 'Arial',
        "font_sizes" : font_sizes,
        "view" : 'anterior',
        "legend" : True,
        "include_target" : True,
        "x_label" : False
    }
    
    model_A_settings = {
        "font" : 'Arial',
        "font_sizes" : font_sizes,
        "view" : 'anterior',
        "legend" : False,
        "x_label" : False,
        "y_label_str" : 'Model A value\n(absolute)'
    }
    
    model_B_settings = {
        "font" : 'Arial',
        "font_sizes" : font_sizes,
        "view" : 'anterior',
        "legend" : False,
        "x_label" : True,
        "y_label_str" : 'Model B value\n(relative)'
    }
    
    ax_pro = plot_proteins(ax_pro, embryo, protein_settings)
    ax_model_A = plot_model(ax_model_A, embryo, models[0], model_values[0], model_ylim[0], model_A_settings)
    ax_model_B = plot_model(ax_model_B, embryo, models[1], model_values[1], model_ylim[1], model_B_settings)
    
    fig.tight_layout()

    return fig   

def save_presentation_figs( models, list_of_embryos, model_values, model_ylim, directory_name ):

    font_sizes = set_figs_font_settings()
    rcParams['font.sans-serif'] = ['Arial']

    noc = list_of_embryos[0].number_of_cells
    half_noc = int(np.floor(noc / 2))

    # create directory
    if directory_name[-1] != '/':
        directory_name = directory_name + '/'
    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)        # I have no idea why this is necessary

    for pos_idx, position in enumerate(['anterior_view', 'posterior_view']):

        full_directory_name = directory_name + position + '/'
        if not os.path.isdir(full_directory_name):
            os.mkdir(full_directory_name)

        # position options
        roll_idx = [0,int(noc/2)]
        xticklabel_options = [['post.','ant.','post.'], ['ant.','post.','ant.']]

        for emb_idx in range(len(list_of_embryos)):

            embryo = list_of_embryos[emb_idx]

            fig_pro = plt.figure()
            ax_pro = fig_pro.add_subplot(111)

            fig_pro_des = plt.figure()
            ax_pro_des = fig_pro_des.add_subplot(111)

            fig_full = plt.figure(figsize=(5,9))
            axs={}
            axs[0] = fig_full.add_subplot(311)
            axs[1] = fig_full.add_subplot(312)
            axs[2] = fig_full.add_subplot(313)

            ax_pro.set_xlim([-1, noc])
            ax_pro.set_xticks(np.linspace(0, noc, num=3))
            ax_pro.set_xticklabels(xticklabel_options[pos_idx])
            ax_pro.set_xlabel('Cell location')

            ax_pro_des.set_xlim([-1, noc])
            ax_pro_des.set_xticks(np.linspace(0, noc, num=3))
            ax_pro_des.set_xticklabels(xticklabel_options[pos_idx])
            ax_pro_des.set_xlabel('Cell location')

            for ax_idx in range(3):
                axs[ax_idx].set_xlim([-1, noc])
                axs[ax_idx].set_xticks(np.linspace(0, noc, num=3))
                axs[ax_idx].set_xticklabels(xticklabel_options[pos_idx])
            axs[2].set_xlabel('Cell location')

            ''' proteins '''
            plot_inducer = np.roll(embryo.inducer.conc, roll_idx[pos_idx])
            plot_inhibitor = np.roll(embryo.inhibitor.conc,roll_idx[pos_idx])
            plot_target = np.roll(embryo.target,roll_idx[pos_idx])
            axs[0].set_title(embryo.name, fontsize=font_sizes['BIGGER_SIZE'])
            axs[0].plot(range(0,noc), plot_inducer, linewidth=1, marker=None, color='C0', markersize = 1, label='Inducer')
            axs[0].plot(range(0,noc), plot_inhibitor, linewidth=1, marker=None, color='C3', markersize = 1, label="Inhibitor")
            axs[0].set_ylim([0,1.3])
            axs[0].set_ylabel('Protein conc.')

            ax_pro.set_title(list_of_embryos[emb_idx].name, fontsize=font_sizes['BIGGER_SIZE'])
            ax_pro.plot(range(0,noc), plot_inducer, linewidth=1, marker=None, color='C0', markersize = 1, label='Inducer')
            ax_pro.plot(range(0,noc), plot_inhibitor, linewidth=1, marker=None, color='C3', markersize = 1, label="Inhibitor")
            ax_pro.set_ylim([0,1.3])
            ax_pro.set_ylabel('Protein conc.')
            ax_pro.legend(loc='upper right')

            ax_pro_des.set_title(list_of_embryos[emb_idx].name, fontsize=font_sizes['BIGGER_SIZE'])
            ax_pro_des.plot(range(0,noc), plot_inducer, linewidth=1, marker=None, color='C0', markersize = 1, label='Inducer')
            ax_pro_des.plot(range(0,noc), plot_inhibitor, linewidth=1, marker=None, color='C3', markersize = 1, label="Inhibitor")
            ax_pro_des.set_ylim([0,1.3])
            ax_pro_des.set_ylabel('Protein conc.')

            # [ymin, ymax] = axs[0].get_ylim()
            # fig_height =  ymax - ymin
            # fig_height_large = fig_height*1.2
            # new_ymin = ymax - fig_height_large
            # axs[0].set_ylim([new_ymin, ymax])
            # bar_height = fig_height_large * 0.1
            # bar_bottom = new_ymin + 0.05*fig_height_large
            # bar_top = bar_bottom + bar_height

            # With bar and text below
            [protein_ymin, protein_ymax] = axs[0].get_ylim()
            fig_height =  protein_ymax - protein_ymin
            fig_height_large = fig_height*1.3
            new_protein_ymin = protein_ymax - fig_height_large
            axs[0].set_ylim([new_protein_ymin, protein_ymax])
            bar_height = fig_height * 0.1
            bar_bottom = new_protein_ymin + 0.18 * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc = new_protein_ymin + 0.05*fig_height_large

            target_color = 'C4'  # purple
            gray_where = np.array([i == -1 for i in plot_target])
            certain_where = np.array([i == 1 for i in plot_target])
            # uncertain_where = np.array([i == 0 for i in embryo.target])
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=target_color, edgecolor=target_color, alpha=1, step='mid', linewidth=1)
            # axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=uncertain_where, facecolor=target_color, edgecolor=target_color, alpha=0.4, step='mid', linewidth=1)
            axs[0].text(393, text_yloc, 'Target streak', backgroundcolor='lightgray', color=target_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')
            axs[0].legend(loc='upper right', fontsize=font_sizes['SMALLEST_SIZE'])

            # [ymin, ymax] = ax_pro_des.get_ylim()
            # fig_height =  ymax - ymin
            # fig_height_large = fig_height*1.2
            # new_ymin = ymax - fig_height_large
            # ax_pro_des.set_ylim([new_ymin, ymax])
            # bar_height = fig_height_large * 0.1
            # bar_bottom = new_ymin + 0.05*fig_height_large
            # bar_top = bar_bottom + bar_height

            # With bar and text below
            [protein_ymin, protein_ymax] = ax_pro_des.get_ylim()
            fig_height =  protein_ymax - protein_ymin
            fig_height_large = fig_height*1.3
            new_protein_ymin = protein_ymax - fig_height_large
            ax_pro_des.set_ylim([new_protein_ymin, protein_ymax])
            bar_height = fig_height_large * 0.08
            bar_bottom = new_protein_ymin + 0.12*fig_height_large
            bar_top = bar_bottom + bar_height
            text_yloc = new_protein_ymin + 0.045*fig_height_large
            ax_pro_des.fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
            ax_pro_des.fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=target_color, edgecolor=target_color, alpha=1, step='mid', linewidth=1)
            ax_pro_des.text(397, text_yloc, 'Target streak', backgroundcolor='lightgray', color=target_color, fontsize=font_sizes['SMALL_SIZE'], fontweight='bold')
            ax_pro_des.legend(loc='upper right')

            ''' model A '''
            plot_model = np.roll(model_values[0, emb_idx, :], roll_idx[pos_idx])
            axs[1].set_ylabel('Model A value')
            axs[1].plot(range(0,noc), plot_model, linewidth=1, marker=None, color=models[0].plot_color, markersize = 1, label=models[0].label)
            axs[1].plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            axs[1].legend(loc='upper right', fontsize=font_sizes['SMALLEST_SIZE'])

            # [ymin, ymax] = [model_ylim[0, emb_idx, 0], model_ylim[0, emb_idx, 1]]
            # fig_height =  ymax - ymin
            # fig_height_large = fig_height*1.2
            # new_ymin = ymax - fig_height_large
            # axs[1].set_ylim([new_ymin, ymax])
            # bar_height = fig_height_large * 0.1
            # bar_bottom = new_ymin + 0.05*fig_height_large
            # bar_top = bar_bottom + bar_height

            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[0, emb_idx, 0], model_ylim[0, emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * 1.3
            new_model_ymin = model_ymax - fig_height_large
            axs[1].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * 0.1
            bar_bottom = new_model_ymin + 0.18 * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc = new_model_ymin + 0.05 * fig_height_large

            brachyury_color = models[0].plot_color  # 'C4'  # purple
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < models[0].threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= models[0].threshold ,facecolor=brachyury_color, edgecolor=brachyury_color, alpha=1, step='mid', linewidth=1)
            axs[1].text(368, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')


            ''' model B '''
            plot_model = np.roll(model_values[1, emb_idx, :], roll_idx[pos_idx])
            axs[2].set_ylabel('Model B value')
            axs[2].plot(range(0,noc), plot_model, linewidth=1, marker=None, color=models[1].plot_color, markersize = 1, label=models[1].label)
            axs[2].plot(range(0,noc), [models[1].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            axs[2].legend(loc='upper right', fontsize=font_sizes['SMALLEST_SIZE'])

            # [ymin, ymax] = [model_ylim[1, emb_idx, 0], model_ylim[1, emb_idx, 1]]
            # fig_height =  ymax - ymin
            # fig_height_large = fig_height*1.2
            # new_ymin = ymax - fig_height_large
            # axs[2].set_ylim([new_ymin, ymax])
            # bar_height = fig_height_large * 0.1
            # bar_bottom = new_ymin + 0.05*fig_height_large
            # bar_top = bar_bottom + bar_height

            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[1, emb_idx, 0], model_ylim[1, emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * 1.3
            new_model_ymin = model_ymax - fig_height_large
            axs[2].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * 0.1
            bar_bottom = new_model_ymin + 0.18 * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc = new_model_ymin + 0.05 * fig_height_large

            brachyury_color = models[1].plot_color # 'C4'  # purple
            axs[2].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < models[1].threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[2].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= models[1].threshold ,facecolor=brachyury_color, edgecolor=brachyury_color, alpha=1, step='mid', linewidth=1)
            axs[2].text(368, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')


            fig_full.tight_layout()
            fig_pro.tight_layout()
            fig_pro_des.tight_layout()

            filename = full_directory_name + embryo.index + embryo.name + '.png'
            filename_pro = full_directory_name + embryo.index + embryo.name + '_pro.png'
            filename_pro_des = full_directory_name + embryo.index + embryo.name + '_pro_des.png'
            fig_full.savefig(filename, format='png')
            fig_pro.savefig(filename_pro, format='png')
            fig_pro_des.savefig(filename_pro_des, format='png')

            plt.close(fig_full)
            plt.close(fig_pro)
            plt.close(fig_pro_des)

def save_method_figs( models, list_of_embryos, model_values, model_ylim, font_string, directory_name ):

    font_sizes = set_figs_font_settings()
    rcParams['font.sans-serif'] = [font_string]

    noc = list_of_embryos[0].number_of_cells
    half_noc = int(np.floor(noc / 2))

    # create directory
    if directory_name[-1] != '/':
        directory_name = directory_name + '/'
    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)        # I have no idea why this is necessary

    for pos_idx, position in enumerate(['anterior_view', 'posterior_view']):

        full_directory_name = directory_name + position + '/'
        if not os.path.isdir(full_directory_name):
            os.mkdir(full_directory_name)

        # position options
        roll_idx = [0,int(noc/2)]
        xticklabel_options = [['post.','ant.','post.'], ['ant.','post.','ant.']]

        for emb_idx in range(len(list_of_embryos)):

            embryo = list_of_embryos[emb_idx]

            fig_pro_des = plt.figure(figsize=(5,4))
            ax_pro_des = fig_pro_des.add_subplot(111)

            fig_full = plt.figure(figsize=(5,5))
            axs={}
            axs[0] = fig_full.add_subplot(211)
            axs[1] = fig_full.add_subplot(212)

            ax_pro_des.set_xlim([-1, noc])
            ax_pro_des.set_xticks(np.linspace(0, noc, num=3))
            ax_pro_des.set_xticklabels(xticklabel_options[pos_idx])
            ax_pro_des.tick_params(axis='x',width=0)
            ax_pro_des.set_xlabel('Cell location')

            for ax_idx in range(2):
                axs[ax_idx].set_xlim([-1, noc])
                axs[ax_idx].set_xticks(np.linspace(0, noc, num=3))
                axs[ax_idx].set_xticklabels(xticklabel_options[pos_idx])
                axs[ax_idx].tick_params(axis='x',width=0)
            axs[0].set_xticklabels([])
            axs[1].set_xlabel('Cell location')

            ''' proteins '''
            plot_inducer = np.roll(embryo.inducer.conc, roll_idx[pos_idx])
            plot_inhibitor = np.roll(embryo.inhibitor.conc,roll_idx[pos_idx])
            plot_target = np.roll(embryo.target,roll_idx[pos_idx])
            target_color = 'C4'  # purple
            gray_where = np.array([i == -1 for i in plot_target])
            certain_where = np.array([i == 1 for i in plot_target])

            protein_ymax = 1.4
            ax_pro_des.plot(range(0,noc), plot_inducer, linewidth=2, linestyle='solid', marker=None, color='C0', markersize = 1, label='Inducer, V')
            ax_pro_des.plot(range(0,noc), plot_inhibitor, linewidth=2, linestyle='dashdot', marker=None, color='C3', markersize = 1, label="Inhibitor, B")
            ax_pro_des.set_ylim([0,protein_ymax])
            ax_pro_des.set_ylabel('Protein conc.')

            # get rid of the frame
            ax_pro_des.spines['top'].set_visible(False)
            ax_pro_des.spines['right'].set_visible(False)
            ax_pro_des.spines['bottom'].set_visible(False)
            ax_pro_des.spines['left'].set_linewidth(1.5)
            ax_pro_des.spines['left'].set_bounds([0,protein_ymax])
            ax_pro_des.axhline(linewidth=1.5, color='k')

            # With bar and text below
            [protein_ymin, protein_ymax] = ax_pro_des.get_ylim()
            fig_height =  protein_ymax - protein_ymin
            fig_height_large = fig_height*1.35
            new_protein_ymin = protein_ymax - fig_height_large
            ax_pro_des.set_ylim([new_protein_ymin, protein_ymax])
            bar_height = fig_height_large * 0.11
            bar_bottom = new_protein_ymin + 0.133*fig_height_large
            bar_top = bar_bottom + bar_height
            text_yloc = new_protein_ymin + 0.045*fig_height_large
            ax_pro_des.fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
            ax_pro_des.fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=target_color, edgecolor=target_color, alpha=1, step='mid', linewidth=1)
            ax_pro_des.text(355, text_yloc, 'Target streak', backgroundcolor='lightgray', color=target_color, fontsize=font_sizes['SMALL_SIZE'], fontweight='bold')
            # ax_pro_des.text(418, text_yloc, 'Target streak', backgroundcolor='lightgray', color=target_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')

            ax_pro_des.legend(loc='upper center', fontsize=font_sizes["SMALLER_SIZE"], ncol=2)

            # model legend params
            legend_height = 0.57

            # model bar params
            fig_height_increase_multiplier = 1.5
            bar_height_proportion = 0.24
            bar_bottom_proportion = 0.27

            predicted_s_loc = 255
            ''' model A '''
            plot_model = np.roll(model_values[0, emb_idx, :], roll_idx[pos_idx])
            axs[0].set_ylabel('Model A value\n(absolute)')
            axs[0].plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[0].plot_color, markersize = 1, label=models[0].label)
            axs[0].plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            # axs[0].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])

            # get rid of the frame
            axs[0].spines['top'].set_visible(False)
            axs[0].spines['right'].set_visible(False)
            axs[0].spines['bottom'].set_visible(False)
            axs[0].spines['left'].set_linewidth(1.5)

            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[0, emb_idx, 0], model_ylim[0, emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * fig_height_increase_multiplier
            new_model_ymin = model_ymax - fig_height_large
            axs[0].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * bar_height_proportion
            bar_bottom = new_model_ymin + bar_bottom_proportion * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc = new_model_ymin + 0.04 * fig_height_large

            brachyury_color = models[0].plot_color  # 'C4'  # purple
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < models[0].threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= models[0].threshold ,facecolor=brachyury_color, edgecolor=brachyury_color, alpha=1, step='mid', linewidth=1)
            axs[0].text(predicted_s_loc, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontsize=font_sizes['SMALL_SIZE'], fontweight='bold')


            ''' model B '''
            plot_model = np.roll(model_values[1, emb_idx, :], roll_idx[pos_idx])
            axs[1].set_ylabel('Model B value\n(relative)')
            axs[1].plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[1].plot_color, markersize = 1, label=models[1].label)
            axs[1].plot(range(0,noc), [models[1].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            # axs[1].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])

            # get rid of the frame
            axs[1].spines['top'].set_visible(False)
            axs[1].spines['right'].set_visible(False)
            axs[1].spines['bottom'].set_visible(False)
            axs[1].spines['left'].set_linewidth(1.5)

            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[1, emb_idx, 0], model_ylim[1, emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * fig_height_increase_multiplier
            new_model_ymin = model_ymax - fig_height_large
            axs[1].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * bar_height_proportion
            bar_bottom = new_model_ymin + bar_bottom_proportion * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc = new_model_ymin + 0.04 * fig_height_large

            brachyury_color = models[1].plot_color # 'C4'  # purple
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < models[1].threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= models[1].threshold ,facecolor=brachyury_color, edgecolor=brachyury_color, alpha=1, step='mid', linewidth=1)
            axs[1].text(predicted_s_loc, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontsize=font_sizes['SMALL_SIZE'], fontweight='bold')

            fig_full.tight_layout()
            fig_pro_des.tight_layout()

            filename = full_directory_name + embryo.index + embryo.name + '.jpg'
            filename_pro_des = full_directory_name + embryo.index + embryo.name + '_pro_des.jpg'
            fig_full.savefig(filename, format='jpg', dpi=600)
            fig_pro_des.savefig(filename_pro_des, format='jpg', dpi=600)

            plt.close(fig_full)
            plt.close(fig_pro_des)
            
def save_method_figs_poster( models, list_of_embryos, model_values, model_ylim, font_string, directory_name ):

    font_sizes = set_figs_font_settings()
    rcParams['font.sans-serif'] = [font_string]

    noc = list_of_embryos[0].number_of_cells
    half_noc = int(np.floor(noc / 2))

    # create directory
    if directory_name[-1] != '/':
        directory_name = directory_name + '/'
    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)        # I have no idea why this is necessary

    for pos_idx, position in enumerate(['anterior_view', 'posterior_view']):

        full_directory_name = directory_name + position + '/'
        if not os.path.isdir(full_directory_name):
            os.mkdir(full_directory_name)

        # position options
        roll_idx = [0,int(noc/2)]
        xticklabel_options = [['post.','ant.','post.'], ['ant.','post.','ant.']]

        for emb_idx in range(len(list_of_embryos)):

            embryo = list_of_embryos[emb_idx]

            fig_pro_des = plt.figure(figsize=(5,4))
            ax_pro_des = fig_pro_des.add_subplot(111)

            fig_full = plt.figure(figsize=(5,5))
            axs={}
            axs[0] = fig_full.add_subplot(211)
            axs[1] = fig_full.add_subplot(212)

            ax_pro_des.set_xlim([-1, noc])
            ax_pro_des.set_xticks(np.linspace(0, noc, num=3))
            ax_pro_des.set_xticklabels(xticklabel_options[pos_idx])
            ax_pro_des.tick_params(axis='x',width=0)
            ax_pro_des.set_xlabel('Cell location')

            for ax_idx in range(2):
                axs[ax_idx].set_xlim([-1, noc])
                axs[ax_idx].set_xticks(np.linspace(0, noc, num=3))
                axs[ax_idx].set_xticklabels(xticklabel_options[pos_idx])
                axs[ax_idx].tick_params(axis='x',width=0)
            axs[0].set_xticklabels([])
            axs[1].set_xlabel('Cell location')

            ''' proteins '''
            plot_inducer = np.roll(embryo.inducer.conc, roll_idx[pos_idx])
            plot_inhibitor = np.roll(embryo.inhibitor.conc,roll_idx[pos_idx])
            plot_target = np.roll(embryo.target,roll_idx[pos_idx])
            target_color = 'C4'  # purple
            gray_where = np.array([i == -1 for i in plot_target])
            certain_where = np.array([i == 1 for i in plot_target])

            protein_ymax = 1.4
            ax_pro_des.plot(range(0,noc), plot_inducer, linewidth=2, linestyle='solid', marker=None, color='C0', markersize = 1, label='Inducer')
            ax_pro_des.plot(range(0,noc), plot_inhibitor, linewidth=2, linestyle='dashdot', marker=None, color='C3', markersize = 1, label="Inhibitor")
            ax_pro_des.set_ylim([0,protein_ymax])
            ax_pro_des.set_ylabel('Protein conc.')

            # get rid of the frame
            ax_pro_des.spines['top'].set_visible(False)
            ax_pro_des.spines['right'].set_visible(False)
            ax_pro_des.spines['bottom'].set_visible(False)
            ax_pro_des.spines['left'].set_linewidth(1.5)
            ax_pro_des.spines['left'].set_bounds([0,protein_ymax])
            ax_pro_des.axhline(linewidth=1.5, color='k')

            # With bar and text below
            [protein_ymin, protein_ymax] = ax_pro_des.get_ylim()
            fig_height =  protein_ymax - protein_ymin
            fig_height_large = fig_height*1.4
            new_protein_ymin = protein_ymax - fig_height_large
            ax_pro_des.set_ylim([new_protein_ymin, protein_ymax])
            bar_height = fig_height_large * 0.12
            bar_bottom = new_protein_ymin + 0.15*fig_height_large
            bar_top = bar_bottom + bar_height
            text_yloc = new_protein_ymin + 0.045*fig_height_large
            ax_pro_des.fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
            ax_pro_des.fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=target_color, edgecolor=target_color, alpha=1, step='mid', linewidth=1)
            ax_pro_des.text(378, text_yloc, 'Target streak', backgroundcolor='lightgray', color=target_color, fontsize=font_sizes['SMALLER_SIZE'], fontweight='bold')
            # ax_pro_des.text(418, text_yloc, 'Target streak', backgroundcolor='lightgray', color=target_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')

            ax_pro_des.legend(loc='upper center', fontsize=font_sizes['SMALLER_SIZE'], ncol=2)
            
            ax_pro_des.set_title(embryo.fig_title, fontweight='bold')

            # model legend params
            legend_height = 0.57

            # model bar params
            fig_height_increase_multiplier = 1.5
            bar_height_proportion = 0.18
            bar_bottom_proportion = 0.24

            predicted_fontsize = font_sizes["SMALLEST_SIZE"]
            predicted_s_loc = 341
            ''' model A '''
            axs[0].set_title("Model predictions", fontweight='bold')
            
            plot_model = np.roll(model_values[0, emb_idx, :], roll_idx[pos_idx])
            axs[0].set_ylabel('Model A value\n(absolute)')
            axs[0].plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[0].plot_color, markersize = 1, label=models[0].label)
            axs[0].plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            # axs[0].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])

            # get rid of the frame
            axs[0].spines['top'].set_visible(False)
            axs[0].spines['right'].set_visible(False)
            axs[0].spines['bottom'].set_visible(False)
            axs[0].spines['left'].set_linewidth(1.5)

            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[0, emb_idx, 0], model_ylim[0, emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * fig_height_increase_multiplier
            new_model_ymin = model_ymax - fig_height_large
            axs[0].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * bar_height_proportion
            bar_bottom = new_model_ymin + bar_bottom_proportion * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc = new_model_ymin + 0.04 * fig_height_large

            brachyury_color = models[0].plot_color  # 'C4'  # purple
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < models[0].threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= models[0].threshold ,facecolor=brachyury_color, edgecolor=brachyury_color, alpha=1, step='mid', linewidth=1)
            axs[0].text(predicted_s_loc, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontsize=predicted_fontsize, fontweight='bold')

            
            ''' model B '''
            plot_model = np.roll(model_values[1, emb_idx, :], roll_idx[pos_idx])
            axs[1].set_ylabel('Model B value\n(relative)')
            axs[1].plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[1].plot_color, markersize = 1, label=models[1].label)
            axs[1].plot(range(0,noc), [models[1].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            # axs[1].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])

            # get rid of the frame
            axs[1].spines['top'].set_visible(False)
            axs[1].spines['right'].set_visible(False)
            axs[1].spines['bottom'].set_visible(False)
            axs[1].spines['left'].set_linewidth(1.5)

            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[1, emb_idx, 0], model_ylim[1, emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * fig_height_increase_multiplier
            new_model_ymin = model_ymax - fig_height_large
            axs[1].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * bar_height_proportion
            bar_bottom = new_model_ymin + bar_bottom_proportion * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc = new_model_ymin + 0.04 * fig_height_large

            brachyury_color = models[1].plot_color # 'C4'  # purple
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < models[1].threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= models[1].threshold ,facecolor=brachyury_color, edgecolor=brachyury_color, alpha=1, step='mid', linewidth=1)
            axs[1].text(predicted_s_loc, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontsize=predicted_fontsize, fontweight='bold')

            fig_full.tight_layout()
            fig_pro_des.tight_layout()

            filename = full_directory_name + embryo.index + embryo.name + '.jpg'
            filename_pro_des = full_directory_name + embryo.index + embryo.name + '_pro_des.jpg'
            fig_full.savefig(filename, format='jpg', dpi=300)
            fig_pro_des.savefig(filename_pro_des, format='jpg', dpi=300)

            plt.close(fig_full)
            plt.close(fig_pro_des)

def save_results_figs( models, list_of_embryos, model_values, model_ylim, font_string, directory_name ):

    font_sizes = set_figs_font_settings()
    rcParams['font.sans-serif'] = [font_string]
    
    font_sizes['SMALLER_SIZE'] = 16
    font_sizes['SMALL_SIZE'] = 20
    font_sizes['MEDIUM_SIZE'] = 24
    
    plt.rc('font', size=font_sizes['SMALL_SIZE'])          # controls default text sizes
    plt.rc('axes', titlesize=font_sizes['MEDIUM_SIZE'])     # fontsize of the axes title
    plt.rc('axes', labelsize=font_sizes['MEDIUM_SIZE'])    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=font_sizes['SMALL_SIZE'])    # fontsize of the tick labels
    plt.rc('ytick', labelsize=font_sizes['SMALL_SIZE'])    # fontsize of the tick labels
    plt.rc('legend', fontsize=font_sizes['SMALL_SIZE'])    # legend fontsize

    noc = list_of_embryos[0].number_of_cells
    half_noc = int(np.floor(noc / 2))
    
    plot_linewidth = 3
    inset_linewidth = 2
    thresh_linewidth = 1.3
    inset_thresh_linewidth = 1.3

    # create directory
    if directory_name[-1] != '/':
        directory_name = directory_name + '/'
    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)        # I have no idea why this is necessary

    inset_axes_A = [5,6,7,8]
    inset_axes_B = [5,6,7,8,9]

    for pos_idx, position in enumerate(['anterior_view', 'posterior_view']):

        full_directory_name = directory_name + position + '/'
        if not os.path.isdir(full_directory_name):
            os.mkdir(full_directory_name)

        # position options
        roll_idx = [0,int(noc/2)]
        xticklabel_options = [['post.','ant.','post.'], ['ant.','post.','ant.']]


        for emb_idx in range(len(list_of_embryos)):

            embryo = list_of_embryos[emb_idx]

            fig_full = plt.figure(figsize=(5,9))
            axs={}
            axs[0] = fig_full.add_subplot(311)
            axs[1] = fig_full.add_subplot(312)
            axs[2] = fig_full.add_subplot(313)

            for ax_idx in range(3):
                axs[ax_idx].set_xlim([-1, noc])
                axs[ax_idx].set_xticks(np.linspace(0, noc, num=3))
                axs[ax_idx].set_xticklabels([])
                axs[ax_idx].tick_params(axis='x',width=0)
            axs[2].set_xlabel('Cell location')
            axs[2].set_xticklabels(xticklabel_options[pos_idx])

            ''' proteins '''
            plot_inducer = np.roll(embryo.inducer.conc, roll_idx[pos_idx])
            plot_inhibitor = np.roll(embryo.inhibitor.conc,roll_idx[pos_idx])
            plot_target = np.roll(embryo.target,roll_idx[pos_idx])
            target_color = 'C4'  # purple
            gray_where = np.array([i == -1 for i in plot_target])
            certain_where = np.array([i == 1 for i in plot_target])

            protein_ymax = 1.55
            if emb_idx in [11]:
                protein_ymax = 3
            axs[0].plot(range(0,noc), plot_inducer, linewidth=plot_linewidth, marker=None, color='C0', markersize = 1, label='Inducer')
            axs[0].plot(range(0,noc), plot_inhibitor, linewidth=plot_linewidth, linestyle='dashdot', marker=None, color='C3', markersize = 1, label="Inhibitor")
            axs[0].set_ylim([0,protein_ymax])
            axs[0].set_ylabel('Protein conc.\n')

            # get rid of the frame
            axs[0].spines['top'].set_visible(False)
            axs[0].spines['right'].set_visible(False)
            axs[0].spines['bottom'].set_visible(False)
            axs[0].spines['left'].set_linewidth(1.5)
            axs[0].spines['left'].set_bounds([0,protein_ymax])
            axs[0].axhline(linewidth=1.5, color='k')

            # With bar and text below
            [protein_ymin, protein_ymax] = axs[0].get_ylim()
            fig_height =  protein_ymax - protein_ymin
            fig_height_large = fig_height*1.5
            new_protein_ymin = protein_ymax - fig_height_large
            axs[0].set_ylim([new_protein_ymin, protein_ymax])
            bar_height = fig_height_large * 0.14
            bar_bottom = new_protein_ymin + 0.17*fig_height_large
            bar_top = bar_bottom + bar_height
            text_yloc_pro = new_protein_ymin + 0.04*fig_height_large
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=target_color, edgecolor=target_color, alpha=1, step='mid', linewidth=1)

            # axs[0].legend(loc='upper center', fontsize=font_sizes['SMALLER_SIZE'], ncol=2)
            axs[0].legend(loc='upper right', fontsize=font_sizes['SMALLER_SIZE'], ncol=1)

            # model legend params
            legend_height = 0.46

            # model bar params
            fig_height_increase_multiplier = 1.5
            bar_height_proportion = 0.20
            bar_bottom_proportion = 0.25

            axs[0].set_title(embryo.fig_title, fontweight='bold')
            # fig_full.suptitle('C-A-C')

            ''' model A '''
            plot_model = np.roll(model_values[0, emb_idx, :], roll_idx[pos_idx])

            # plt.rc('text', usetex=True)
            # axs[1].set_ylabel(r'{\fontsize{20pt}{3em}\selectfont{Arial}Model A value\n}{\fontsize{16pt}{3em}\selectfont{Arial}(absolute)}')
            # plt.rc('text', usetex=False)
            axs[1].set_ylabel('Model A value\n(absolute)')

            axs[1].plot(range(0,noc), plot_model, linewidth=plot_linewidth, marker=None, color=models[0].plot_color, markersize = 1, label=models[0].label)
            axs[1].plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=thresh_linewidth, marker=None, color='black', markersize = 1)
            # axs[1].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])

            # get rid of the frame
            axs[1].spines['top'].set_visible(False)
            axs[1].spines['right'].set_visible(False)
            axs[1].spines['bottom'].set_visible(False)
            axs[1].spines['left'].set_linewidth(1.5)

            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[0, emb_idx, 0], model_ylim[0, emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * fig_height_increase_multiplier
            new_model_ymin = model_ymax - fig_height_large
            axs[1].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * bar_height_proportion
            bar_bottom = new_model_ymin + bar_bottom_proportion * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc_A = new_model_ymin + 0.04 * fig_height_large

            brachyury_color_A = models[0].plot_color  # 'C4'  # purple
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < models[0].threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= models[0].threshold ,facecolor=brachyury_color_A, edgecolor=brachyury_color_A, alpha=1, step='mid', linewidth=1)

            
            inset_connector_alpha=0.6
            if emb_idx in inset_axes_A:

                # inset axes....
                
                axinsA = axs[1].inset_axes([0.02, 0.68, 0.35, 0.3])
                axinsA.plot(range(0,noc), plot_model, linewidth=inset_linewidth, marker=None, color=models[0].plot_color, markersize = 1)
                axinsA.plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=inset_thresh_linewidth, marker=None, color='black', markersize = 1)
                # axinsA.imshow(Z2, extent=extent, interpolation="nearest", origin="lower")
                # sub region of the original image
                x1, x2, y1, y2 = 230, 370, 0.45, 0.55
                axinsA.set_xlim(x1, x2)
                axinsA.set_ylim(y1, y2)
                axinsA.set_xticklabels('')
                axinsA.set_yticklabels('')
                for pos in ['top','bottom','left','right']:
                    axinsA.spines[pos].set_linewidth(inset_thresh_linewidth)

                axs[1].indicate_inset_zoom(axinsA, edgecolor='black', alpha=inset_connector_alpha)


            ''' model B '''
            plot_model = np.roll(model_values[1, emb_idx, :], roll_idx[pos_idx])
            axs[2].set_ylabel('Model B value\n(relative)')
            axs[2].plot(range(0,noc), plot_model, linewidth=plot_linewidth, marker=None, color=models[1].plot_color, markersize = 1, label=models[1].label)
            axs[2].plot(range(0,noc), [models[1].threshold for i in range(noc)], '--', linewidth=thresh_linewidth, marker=None, color='black', markersize = 1)
            # axs[2].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])

            # get rid of the frame
            axs[2].spines['top'].set_visible(False)
            axs[2].spines['right'].set_visible(False)
            axs[2].spines['bottom'].set_visible(False)
            axs[2].spines['left'].set_linewidth(1.5)

            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[1, emb_idx, 0], model_ylim[1, emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * fig_height_increase_multiplier
            new_model_ymin = model_ymax - fig_height_large
            axs[2].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * bar_height_proportion
            bar_bottom = new_model_ymin + bar_bottom_proportion * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc_B = new_model_ymin + 0.04 * fig_height_large

            brachyury_color_B = models[1].plot_color # 'C4'  # purple
            axs[2].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < models[1].threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[2].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= models[1].threshold ,facecolor=brachyury_color_B, edgecolor=brachyury_color_B, alpha=1, step='mid', linewidth=1)

            if emb_idx in inset_axes_B:

                # inset axes....
                
                axinsB = axs[2].inset_axes([0.02, 0.33, 0.35, 0.35])
                axinsB.plot(range(0,noc), plot_model, linewidth=inset_linewidth, marker=None, color=models[1].plot_color, markersize = 1)
                axinsB.plot(range(0,noc), [models[1].threshold for i in range(noc)], '--', linewidth=inset_thresh_linewidth, marker=None, color='black', markersize = 1)
                # axinsB.imshow(Z2, extent=extent, interpolation="nearest", origin="lower")
                # sub region of the original image
                x1, x2, y1, y2 = 230, 370, 0.25, 0.60
                axinsB.set_xlim(x1, x2)
                axinsB.set_ylim(y1, y2)
                axinsB.set_xticklabels('')
                axinsB.set_yticklabels('')
                for pos in ['top','bottom','left','right']:
                    axinsB.spines[pos].set_linewidth(inset_thresh_linewidth)

                axs[2].indicate_inset_zoom(axinsB, edgecolor='black', alpha=inset_connector_alpha)


            # clear yaxis for certain embryos
            yaxis_clear_embryos = [2,3,4,6,7,8,9,13,14]
            if emb_idx in yaxis_clear_embryos:
                for ax_idx in range(3):
                    axs[ax_idx].yaxis.label.set_visible(False)
                    axs[ax_idx].set_yticklabels([])

            # clear legend for certain embryos
            # legend_clear_embryos = [0,8,9,10,11,13,14,15,17,19,20,21]
            legend_clear_embryos = [1,2,3,5,6,7,8,12,13]
            if emb_idx in legend_clear_embryos:
                for ax_idx in range(1):
                    axs[ax_idx].get_legend().remove()
                # for ax_idx in range(3):
                #     axs[ax_idx].get_legend().remove()

            # only add text for certain embryos

            add_text_embryos = [item for item in range(len(list_of_embryos)) if item not in legend_clear_embryos]
            if emb_idx in add_text_embryos:
                if emb_idx in [10,11]:
                    text_xloc_pro = 230             # for bmp ant plots
                    text_xloc_model = 138
                else:
                    text_xloc_pro = 317           # for most plots
                    text_xloc_model = 250
                font_props = font_manager.FontProperties(size=font_sizes['SMALL_SIZE'], weight='semibold')
                target = axs[0].text(text_xloc_pro, text_yloc_pro, 'Target streak', backgroundcolor='lightgray', color=target_color, fontproperties=font_props)
                predicted_A = axs[1].text(text_xloc_model, text_yloc_A, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color_A, fontproperties=font_props)
                predicted_B = axs[2].text(text_xloc_model, text_yloc_B, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color_B, fontproperties=font_props)

                outline_width = 0.1
                target.set_path_effects([path_effects.Stroke(linewidth=outline_width, foreground='black'),
                                       path_effects.Normal()])
                predicted_A.set_path_effects([path_effects.Stroke(linewidth=outline_width, foreground='black'),
                                       path_effects.Normal()])
                predicted_B.set_path_effects([path_effects.Stroke(linewidth=outline_width, foreground='black'),
                                           path_effects.Normal()])


            fig_full.tight_layout()

            filename = full_directory_name + embryo.index + embryo.name + '.jpg'
            fig_full.savefig(filename, format='jpg', dpi=600)

            plt.close(fig_full)
                    
def save_results_figs_poster( models, list_of_embryos, model_values, model_ylim, font_string, directory_name ):

    font_sizes = set_figs_font_settings_poster()
    rcParams['font.sans-serif'] = [font_string]

    noc = list_of_embryos[0].number_of_cells
    half_noc = int(np.floor(noc / 2))

    # create directory
    if directory_name[-1] != '/':
        directory_name = directory_name + '/'
    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)        # I have no idea why this is necessary

    inset_axes_A = [5,6,7,8]
    inset_axes_B = [5,6,7,8,9]

    for pos_idx, position in enumerate(['anterior_view', 'posterior_view']):

        full_directory_name = directory_name + position + '/'
        if not os.path.isdir(full_directory_name):
            os.mkdir(full_directory_name)

        # position options
        roll_idx = [0,int(noc/2)]
        xticklabel_options = [['post.','ant.','post.'], ['ant.','post.','ant.']]


        for emb_idx in range(len(list_of_embryos)):

            embryo = list_of_embryos[emb_idx]

            fig_full = plt.figure(figsize=(3.94,3.23))
            axs={}
            axs[0] = fig_full.add_subplot(211)
            axs[1] = fig_full.add_subplot(413)
            axs[2] = fig_full.add_subplot(414)

            for ax_idx in range(3):
                axs[ax_idx].set_xlim([-1, noc])
                axs[ax_idx].set_xticks(np.linspace(0, noc, num=3))
                axs[ax_idx].set_xticklabels([])
                axs[ax_idx].tick_params(axis='x',width=0)
            axs[2].set_xlabel('Cell location')
            axs[2].set_xticklabels(xticklabel_options[pos_idx])

            ''' proteins '''
            plot_inducer = np.roll(embryo.inducer.conc, roll_idx[pos_idx])
            plot_inhibitor = np.roll(embryo.inhibitor.conc,roll_idx[pos_idx])
            plot_target = np.roll(embryo.target,roll_idx[pos_idx])
            target_color = 'C4'  # purple
            gray_where = np.array([i == -1 for i in plot_target])
            certain_where = np.array([i == 1 for i in plot_target])

            protein_ymax = 1.55
            axs[0].plot(range(0,noc), plot_inducer, linewidth=2, marker=None, color='C0', markersize = 1, label='Inducer')
            axs[0].plot(range(0,noc), plot_inhibitor, linewidth=2, linestyle='dashdot', marker=None, color='C3', markersize = 1, label="Inhibitor")
            axs[0].set_ylim([0,protein_ymax])
            axs[0].set_ylabel('Protein\nconc.\n')

            # get rid of the frame
            axs[0].spines['top'].set_visible(False)
            axs[0].spines['right'].set_visible(False)
            axs[0].spines['bottom'].set_visible(False)
            axs[0].spines['left'].set_linewidth(1.5)
            axs[0].spines['left'].set_bounds([0,protein_ymax])
            axs[0].axhline(linewidth=1.5, color='k')

            # With bar and text below
            [protein_ymin, protein_ymax] = axs[0].get_ylim()
            fig_height =  protein_ymax - protein_ymin
            fig_height_large = fig_height*1.32
            new_protein_ymin = protein_ymax - fig_height_large
            axs[0].set_ylim([new_protein_ymin, protein_ymax])
            bar_height = fig_height_large * 0.1
            bar_bottom = new_protein_ymin + 0.12*fig_height_large
            bar_top = bar_bottom + bar_height
            text_yloc_pro = new_protein_ymin + 0.04*fig_height_large
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=target_color, edgecolor=target_color, alpha=1, step='mid', linewidth=1)

            axs[0].legend(fontsize=font_sizes['SMALLEST_SIZE'], bbox_to_anchor=(1.1, 0.9))

            # model legend params
            legend_height = 0.46

            # model bar params
            fig_height_increase_multiplier = 1.3
            bar_height_proportion = 0.12
            bar_bottom_proportion = 0.155

            axs[0].set_title(embryo.fig_title, fontweight='bold')
            # fig_full.suptitle('C-A-C')

            ''' model A '''
            plot_model = np.roll(model_values[0, emb_idx, :], roll_idx[pos_idx])

            # plt.rc('text', usetex=True)
            # axs[1].set_ylabel(r'{\fontsize{20pt}{3em}\selectfont{Arial}Model A value\n}{\fontsize{16pt}{3em}\selectfont{Arial}(absolute)}')
            # plt.rc('text', usetex=False)
            axs[1].set_ylabel('Model\nA')

            axs[1].plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[0].plot_color, markersize = 1, label=models[0].label)
            axs[1].plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            # axs[1].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])

            # get rid of the frame
            axs[1].spines['top'].set_visible(False)
            axs[1].spines['right'].set_visible(False)
            axs[1].spines['bottom'].set_visible(False)
            axs[1].spines['left'].set_linewidth(1.5)

            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[0, emb_idx, 0], model_ylim[0, emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * fig_height_increase_multiplier
            new_model_ymin = model_ymax - fig_height_large
            axs[1].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * bar_height_proportion
            bar_bottom = new_model_ymin + bar_bottom_proportion * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc_A = new_model_ymin + 0.04 * fig_height_large

            brachyury_color_A = models[0].plot_color  # 'C4'  # purple
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < models[0].threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= models[0].threshold ,facecolor=brachyury_color_A, edgecolor=brachyury_color_A, alpha=1, step='mid', linewidth=1)

            if emb_idx in inset_axes_A:

                # inset axes....
                
                axinsA = axs[1].inset_axes([0.02, 0.65, 0.35, 0.3])
                axinsA.plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[0].plot_color, markersize = 1)
                axinsA.plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
                # axinsA.imshow(Z2, extent=extent, interpolation="nearest", origin="lower")
                # sub region of the original image
                x1, x2, y1, y2 = 230, 370, 0.45, 0.55
                axinsA.set_xlim(x1, x2)
                axinsA.set_ylim(y1, y2)
                axinsA.set_xticklabels('')
                axinsA.set_yticklabels('')

                axs[1].indicate_inset_zoom(axinsA)


            ''' model B '''
            plot_model = np.roll(model_values[1, emb_idx, :], roll_idx[pos_idx])
            axs[2].set_ylabel('Model\nB')
            axs[2].plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[1].plot_color, markersize = 1, label=models[1].label)
            axs[2].plot(range(0,noc), [models[1].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            # axs[2].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])

            # get rid of the frame
            axs[2].spines['top'].set_visible(False)
            axs[2].spines['right'].set_visible(False)
            axs[2].spines['bottom'].set_visible(False)
            axs[2].spines['left'].set_linewidth(1.5)

            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[1, emb_idx, 0], model_ylim[1, emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * fig_height_increase_multiplier
            new_model_ymin = model_ymax - fig_height_large
            axs[2].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * bar_height_proportion
            bar_bottom = new_model_ymin + bar_bottom_proportion * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc_B = new_model_ymin + 0.04 * fig_height_large

            brachyury_color_B = models[1].plot_color # 'C4'  # purple
            axs[2].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < models[1].threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[2].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= models[1].threshold ,facecolor=brachyury_color_B, edgecolor=brachyury_color_B, alpha=1, step='mid', linewidth=1)

            if emb_idx in inset_axes_B:

                # inset axes....
                
                axinsB = axs[2].inset_axes([0.02, 0.25, 0.35, 0.35])
                axinsB.plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[1].plot_color, markersize = 1)
                axinsB.plot(range(0,noc), [models[1].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
                # axinsB.imshow(Z2, extent=extent, interpolation="nearest", origin="lower")
                # sub region of the original image
                x1, x2, y1, y2 = 230, 370, 0.2, 0.65
                axinsB.set_xlim(x1, x2)
                axinsB.set_ylim(y1, y2)
                axinsB.set_xticklabels('')
                axinsB.set_yticklabels('')

                axs[2].indicate_inset_zoom(axinsB)


            # clear yaxis for certain embryos
            yaxis_clear_embryos = [2,3,4,6,7,8,9,13,14]
            if emb_idx in yaxis_clear_embryos:
                for ax_idx in range(3):
                    axs[ax_idx].yaxis.label.set_visible(False)
                    axs[ax_idx].set_yticklabels([])

            # clear legend for certain embryos
            # legend_clear_embryos = [0,8,9,10,11,13,14,15,17,19,20,21]
            legend_clear_embryos = [1,2,3,5,6,7,8,12,13]
            if emb_idx in legend_clear_embryos:
                for ax_idx in range(1):
                    axs[ax_idx].get_legend().remove()
                # for ax_idx in range(3):
                     # axs[ax_idx].get_legend().remove()

            # only add text for certain embryos

            add_text_embryos = [item for item in range(len(list_of_embryos)) if item not in legend_clear_embryos]
            if emb_idx in add_text_embryos:
                if emb_idx in [10,11]:
                    text_xloc_pro = 600             # for bmp ant plots
                    text_xloc_model = 550
                else:
                    text_xloc_pro = 600           # for most plots
                    text_xloc_model = 550
                font_props = font_manager.FontProperties(size=font_sizes['SMALLEST_SIZE'], weight='semibold')
                target = axs[0].text(text_xloc_pro, text_yloc_pro, 'Target streak', backgroundcolor='lightgray', color=target_color, fontproperties=font_props)
                predicted_A = axs[1].text(text_xloc_model, text_yloc_A, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color_A, fontproperties=font_props)
                predicted_B = axs[2].text(text_xloc_model, text_yloc_B, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color_B, fontproperties=font_props)

                outline_width = 0.1
                target.set_path_effects([path_effects.Stroke(linewidth=outline_width, foreground='black'),
                                       path_effects.Normal()])
                predicted_A.set_path_effects([path_effects.Stroke(linewidth=outline_width, foreground='black'),
                                       path_effects.Normal()])
                predicted_B.set_path_effects([path_effects.Stroke(linewidth=outline_width, foreground='black'),
                                           path_effects.Normal()])


            fig_full.tight_layout()

            filename = full_directory_name + embryo.index + embryo.name + '.jpg'
            fig_full.savefig(filename, format='jpg', dpi=300)

            plt.close(fig_full)