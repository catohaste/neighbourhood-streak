import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats
import copy
import os.path


def set_figs_font_settings():

    font_sizes = {
        'SMALLEST_SIZE' : 12,
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

def save_standard_figs(Embryo, Model, directory_name):
    
    font_sizes = set_figs_font_settings()
    
    noc = Embryo.number_of_cells
    half_noc = int(np.floor(Embryo.number_of_cells / 2))
    
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

        plot_strings = ['proteins', 'intermediate']
        figs={}
        axs={}
        for idx, plot_string in enumerate(plot_strings):
        
            figs[idx] = plt.figure()
            axs[idx] = figs[idx].add_subplot(111)

            # axs[idx].set_title(Embryo.name, fontsize=font_sizes['BIGGER_SIZE'])
            axs[idx].set_xlim([-1, noc])
            axs[idx].set_xlabel('Cell location')
            axs[idx].set_xticks(np.linspace(0, noc, num=3))
            axs[idx].set_xticklabels(xticklabel_options[pos_idx])
            axs[idx].tick_params(axis='x',width=0)
            # axs[idx].set_xlabel('Cell index')
            # axs[idx].set_xticks(np.linspace(0, noc, num=7))
            # axs[idx].set_xticklabels(np.linspace(-half_noc, half_noc, num=7, dtype=int))
        
        ''' proteins '''
        # print(Embryo.inducer.conc)
        # print(roll_idx)
        # print(pos_idx)
        # print(roll_idx[pos_idx])
        # print(np.roll(Embryo.inducer.conc,roll_idx[pos_idx]))
        plot_inducer = np.roll(Embryo.inducer.conc,roll_idx[pos_idx])
        plot_inhibitor = np.roll(Embryo.inhibitor.conc,roll_idx[pos_idx])
        axs[0].plot(range(0,noc), plot_inducer, linewidth=1, marker=None, color='C0', markersize = 1, label='Inducer')
        axs[0].plot(range(0,noc), plot_inhibitor, linewidth=1, marker=None, color='C3', markersize = 1, label="Inhibitor")
        
        # get rid of the frame
        axs[0].spines['top'].set_visible(False)
        axs[0].spines['right'].set_visible(False)
        axs[0].spines['bottom'].set_visible(False)
        axs[0].spines['left'].set_linewidth(1.5)
        
        axs[0].spines['left'].set_bounds([0,1.3])
        
        axs[0].axhline(linewidth=1.5, color='k')
          
        axs[0].set_ylim([0,1.3])
        axs[0].set_yticks([0.0,0.5,1.0])
        axs[0].set_ylabel('Protein conc.')
        
        # With bar and text below
        [protein_ymin, protein_ymax] = axs[0].get_ylim()
        fig_height =  protein_ymax - protein_ymin
        fig_height_large = fig_height*1.3
        new_protein_ymin = protein_ymax - fig_height_large
        axs[0].set_ylim([new_protein_ymin, protein_ymax])
        bar_height = fig_height_large * 0.08
        bar_bottom = new_protein_ymin + 0.12*fig_height_large
        bar_top = bar_bottom + bar_height
        text_yloc = new_protein_ymin + 0.05*fig_height_large
        
        desired_color = 'C4'  # purple
        gray_where = np.array([i == -1 for i in Embryo.desired])
        certain_where = np.array([i == 1 for i in Embryo.desired])
        uncertain_where = np.array([i == 0 for i in Embryo.desired])
        gray_where = np.roll(gray_where,roll_idx[pos_idx])
        certain_where = np.roll(certain_where,roll_idx[pos_idx])
        uncertain_where = np.roll(uncertain_where,roll_idx[pos_idx])
        
        axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
        axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where,facecolor=desired_color, edgecolor=desired_color, alpha=1, step='mid', linewidth=1)
        axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=uncertain_where, facecolor=desired_color, edgecolor=desired_color, alpha=0.4, step='mid', linewidth=1)
        
        axs[0].text(401, text_yloc, 'Desired streak', backgroundcolor='lightgray', color=desired_color, fontsize=font_sizes['SMALL_SIZE'], fontweight='bold')
        
        axs[0].legend(loc='upper right')
    
        ''' intermediate '''
        plot_inducer_fraction = np.roll(Embryo.inducer_fraction,roll_idx[pos_idx])
        plot_inhibitor_fraction = np.roll(Embryo.inhibitor_fraction,roll_idx[pos_idx])
        
        nbhd_center = 100
        axs[1].plot(range(0,noc), plot_inducer_fraction, linewidth=1, marker=None, color='C9', markersize = 1,label='Inducer-bound fraction')
        axs[1].plot(range(0,noc), plot_inhibitor_fraction, linewidth=1, marker=None, color='C6', markersize = 1, label='Inhibitor-bound fraction')
        ymin, ymax = axs[1].get_ylim()
        if Model.name == 'inducer_SMAD_nbhd':
            half_nbhd = int(np.floor(Model.nbhd_size / 2))
            axs[1].vlines(nbhd_center, ymin, ymax, colors='k', label='nbhd centre')
            axs[1].vlines([nbhd_center - half_nbhd, nbhd_center + half_nbhd], ymin, ymax, colors='k', linestyles='dashed', label='nbhd bounds')
        axs[1].set_ylabel('Intermediate model value')
         # Shrink current axis's height by 10% on the bottom
        l, b, w, h = axs[1].get_position().bounds
        new_box = [l + w*0.15, b + h*0.25, w*0.75, h*0.75]
        axs[1].set_position(new_box, which='both')
        axs[1].legend(loc='upper center', bbox_to_anchor=(0.45, -0.2), fancybox=True, ncol=2)
        axs[1].set_title('Intermediate values', fontsize=font_sizes['BIGGER_SIZE'])
        
        for idx in [0]:
            figs[idx].tight_layout()
            
        for idx, plot_string in enumerate(plot_strings):
            filename = full_directory_name + Embryo.index + Embryo.name + '_' + plot_string + '.png'
            figs[idx].savefig(filename, format='png')

        for idx, plot_string in enumerate(plot_strings):
            plt.close(figs[idx])

def save_model_figs(Embryo, Model, directory_name, suffix):
    
    font_sizes = set_figs_font_settings()
    
    noc = Embryo.number_of_cells
    half_noc = int(np.floor(Embryo.number_of_cells / 2))
    
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

        plot_strings = ['model_' +  Model.index_string]
        figs={}
        axs={}
        for idx, plot_string in enumerate(plot_strings):
        
            figs[idx] = plt.figure()
            axs[idx] = figs[idx].add_subplot(111)

            axs[idx].set_xlim([-1, noc])
            axs[idx].set_xlabel('Cell location')
            axs[idx].set_xticks(np.linspace(0, noc, num=3))
            axs[idx].set_xticklabels(xticklabel_options[pos_idx])
            axs[idx].tick_params(axis='x',width=0)
            # axs[idx].set_xlabel('Cell index')
            # axs[idx].set_xticks(np.linspace(0, noc, num=7))
            # axs[idx].set_xticklabels(np.linspace(-half_noc, half_noc, num=7, dtype=int))
            
        ''' model '''
        plot_model = np.roll(Embryo.model_value,roll_idx[pos_idx])
        axs[0].set_ylabel('Model ' + Model.index_string + ' value')
        axs[0].plot(range(0,noc), plot_model, linewidth=1, marker=None, color=Model.plot_color, markersize = 1, label=Model.label)
        axs[0].plot(range(0,noc), [Model.threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
        
        # Change color of title based on success
        # if Embryo.success:
        #     title_color = 'C2' # green
        # else:
        #     title_color = 'C3' # red
        # axs[0].set_title(Model.index_string + '_' + Model.name + suffix, fontsize=font_sizes['BIGGER_SIZE'], color=title_color)
        
        # get rid of the frame
        axs[0].spines['top'].set_visible(False)
        axs[0].spines['right'].set_visible(False)
        axs[0].spines['bottom'].set_visible(False)
        axs[0].spines['left'].set_linewidth(1.5)
        
        # axs[0].set_ylim(self.model_ylim)
        # fig_height =  self.model_ylim[1] - self.model_ylim[0]
        # bar_height = fig_height * 0.1
        # bar_bottom = self.model_ylim[0] + 0.05*fig_height
        # bar_top = bar_bottom + bar_height
        
        # With bar and text below
        [model_ymin, model_ymax] = Embryo.plot_model_ylim
        fig_height =  model_ymax - model_ymin
        fig_height_large = fig_height*1.3
        new_model_ymin = model_ymax - fig_height_large
        axs[0].set_ylim([new_model_ymin, model_ymax])
        bar_height = fig_height_large * 0.08
        bar_bottom = new_model_ymin + 0.12*fig_height_large
        bar_top = bar_bottom + bar_height
        text_yloc = new_model_ymin + 0.045*fig_height_large
        
        brachyury_color =  Model.plot_color #  'C4'  # purple
        axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < Model.threshold , facecolor='lightgray', alpha=1, step='mid')
        axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= Model.threshold ,facecolor=brachyury_color, edgecolor=brachyury_color, alpha=1, step='mid', linewidth=1)
        
        axs[0].text(375, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontsize=font_sizes['SMALL_SIZE'], fontweight='bold')
        
        axs[0].legend(loc='lower right')
        
        for idx in [0]:
            figs[idx].tight_layout()
            
        for idx, plot_string in enumerate(plot_strings):
            filename = full_directory_name + Embryo.index + Embryo.name + '_' + plot_string + suffix + '.png'
            figs[idx].savefig(filename, format='png')

        for idx, plot_string in enumerate(plot_strings):
            plt.close(figs[idx])
            
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
            plot_desired = np.roll(embryo.desired,roll_idx[pos_idx])
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
            
            desired_color = 'C4'  # purple
            gray_where = np.array([i == -1 for i in plot_desired])
            certain_where = np.array([i == 1 for i in plot_desired])
            # uncertain_where = np.array([i == 0 for i in embryo.desired])
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=desired_color, edgecolor=desired_color, alpha=1, step='mid', linewidth=1)
            # axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=uncertain_where, facecolor=desired_color, edgecolor=desired_color, alpha=0.4, step='mid', linewidth=1)
            axs[0].text(393, text_yloc, 'Desired streak', backgroundcolor='lightgray', color=desired_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')
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
            ax_pro_des.fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=desired_color, edgecolor=desired_color, alpha=1, step='mid', linewidth=1)
            ax_pro_des.text(397, text_yloc, 'Desired streak', backgroundcolor='lightgray', color=desired_color, fontsize=font_sizes['SMALL_SIZE'], fontweight='bold')
            ax_pro_des.legend(loc='upper right')
            
            ''' model A '''
            plot_model = np.roll(model_values[0, emb_idx, :], roll_idx[pos_idx])
            axs[1].set_ylabel('Model A value')
            axs[1].plot(range(0,noc), plot_model, linewidth=1, marker=None, color=models[0].plot_color, markersize = 1, label='With\nnbhd')
            axs[1].plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            # axs[1].set_title('Model ' + models[0].index_string + ', with nbhd', fontsize=font_sizes['BIGGER_SIZE'])
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
            axs[2].plot(range(0,noc), plot_model, linewidth=1, marker=None, color=models[1].plot_color, markersize = 1, label='Without\nnbhd')
            axs[2].plot(range(0,noc), [models[1].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            # axs[2].set_title('Model ' + models[1].index_string + ', without nbhd', fontsize=font_sizes['BIGGER_SIZE'])
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



def save_presentation_figs_duo( model, list_of_embryos, model_values, model_ylim, directory_name ):
    
    font_sizes = set_figs_font_settings()
    rcParams['font.sans-serif'] = ['Arial']
    
    noc = list_of_embryos[0].number_of_cells
    half_noc = int(np.floor(noc / 2))
    
    # create directory
    if directory_name[-1] != '/':
        directory_name = directory_name + '/'
    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)        # I have no idea why this is necessary
        
    model_directory_name = directory_name + 'model_' + model.index_string + '/'
    if not os.path.isdir(model_directory_name):
        os.mkdir(directory_name + 'model_' + model.index_string + '/')

    for pos_idx, position in enumerate(['anterior_view', 'posterior_view']):
    
        full_directory_name = model_directory_name + position + '/'            
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

            fig_full = plt.figure(figsize=(5,6))
            axs={}
            axs[0] = fig_full.add_subplot(211)
            axs[1] = fig_full.add_subplot(212)
        
            ax_pro.set_xlim([-1, noc])
            ax_pro.set_xticks(np.linspace(0, noc, num=3))
            ax_pro.set_xticklabels(xticklabel_options[pos_idx])
            ax_pro.set_xlabel('Cell location')
        
            ax_pro_des.set_xlim([-1, noc])
            ax_pro_des.set_xticks(np.linspace(0, noc, num=3))
            ax_pro_des.set_xticklabels(xticklabel_options[pos_idx])
            ax_pro_des.set_xlabel('Cell location')

            for ax_idx in range(2):
                axs[ax_idx].set_xlim([-1, noc])
                axs[ax_idx].set_xticks(np.linspace(0, noc, num=3))
                axs[ax_idx].set_xticklabels(xticklabel_options[pos_idx])
            axs[1].set_xlabel('Cell location')
    
            ''' proteins '''
            plot_inducer = np.roll(embryo.inducer.conc, roll_idx[pos_idx])
            plot_inhibitor = np.roll(embryo.inhibitor.conc,roll_idx[pos_idx])
            plot_desired = np.roll(embryo.desired,roll_idx[pos_idx])
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
        
            desired_color = 'C4'  # purple
            gray_where = np.array([i == -1 for i in plot_desired])
            certain_where = np.array([i == 1 for i in plot_desired])
            # uncertain_where = np.array([i == 0 for i in embryo.desired])
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=desired_color, edgecolor=desired_color, alpha=1, step='mid', linewidth=1)
            # axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=uncertain_where, facecolor=desired_color, edgecolor=desired_color, alpha=0.4, step='mid', linewidth=1)
            axs[0].text(393, text_yloc, 'Desired streak', backgroundcolor='lightgray', color=desired_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')
            axs[0].legend(loc='upper right', fontsize=font_sizes['SMALLEST_SIZE'])
        
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
            ax_pro_des.fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=desired_color, edgecolor=desired_color, alpha=1, step='mid', linewidth=1)
            ax_pro_des.text(397, text_yloc, 'Desired streak', backgroundcolor='lightgray', color=desired_color, fontsize=font_sizes['SMALL_SIZE'], fontweight='bold')
            ax_pro_des.legend(loc='upper right')
                
            ''' model '''
            plot_model = np.roll(model_values[emb_idx, :], roll_idx[pos_idx])
            axs[1].set_ylabel('Model value')
            axs[1].plot(range(0,noc), plot_model, linewidth=1, marker=None, color=model.plot_color, markersize = 1, label=model.label)
            axs[1].plot(range(0,noc), [model.threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            # axs[1].set_title('Model ' + models[0].index_string + ', with nbhd', fontsize=font_sizes['BIGGER_SIZE'])
            axs[1].legend(loc='upper right', fontsize=font_sizes['SMALLEST_SIZE'])
        
            # With bar and text below
            [model_ymin, model_ymax] = [model_ylim[ emb_idx, 0], model_ylim[ emb_idx, 1]]
            fig_height =  model_ymax - model_ymin
            fig_height_large = fig_height * 1.3
            new_model_ymin = model_ymax - fig_height_large
            axs[1].set_ylim([new_model_ymin, model_ymax])
            bar_height = fig_height * 0.1
            bar_bottom = new_model_ymin + 0.18 * fig_height
            bar_top = bar_bottom + bar_height
            text_yloc = new_model_ymin + 0.05 * fig_height_large
    
            brachyury_color = model.plot_color  # 'C4'  # purple
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model < model.threshold , facecolor='lightgray', alpha=1, step='mid')
            axs[1].fill_between(np.arange(noc), bar_bottom, bar_top, where=plot_model >= model.threshold ,facecolor=brachyury_color, edgecolor=brachyury_color, alpha=1, step='mid', linewidth=1)
            axs[1].text(368, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')        

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
            plot_desired = np.roll(embryo.desired,roll_idx[pos_idx])
            desired_color = 'C4'  # purple
            gray_where = np.array([i == -1 for i in plot_desired])
            certain_where = np.array([i == 1 for i in plot_desired])
            
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
            fig_height_large = fig_height*1.28
            new_protein_ymin = protein_ymax - fig_height_large
            ax_pro_des.set_ylim([new_protein_ymin, protein_ymax])
            bar_height = fig_height_large * 0.075
            bar_bottom = new_protein_ymin + 0.12*fig_height_large
            bar_top = bar_bottom + bar_height
            text_yloc = new_protein_ymin + 0.045*fig_height_large
            ax_pro_des.fill_between(np.arange(noc), bar_bottom, bar_top, where=gray_where , facecolor='lightgray', alpha=1, step='mid')
            ax_pro_des.fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=desired_color, edgecolor=desired_color, alpha=1, step='mid', linewidth=1)
            ax_pro_des.text(388, text_yloc, 'Desired streak', backgroundcolor='lightgray', color=desired_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')
            
            ax_pro_des.legend(loc='upper right', fontsize=font_sizes['SMALLEST_SIZE'])
            
            # model legend params
            legend_height = 0.57
            
            # model bar params
            fig_height_increase_multiplier = 1.4
            bar_height_proportion = 0.15
            bar_bottom_proportion = 0.2
            
            ''' model A '''
            plot_model = np.roll(model_values[0, emb_idx, :], roll_idx[pos_idx])
            axs[0].set_ylabel('Model A value')
            axs[0].plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[0].plot_color, markersize = 1, label='With\nnbhd')
            axs[0].plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            axs[0].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])
            
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
            axs[0].text(368, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')

        
            ''' model B '''
            plot_model = np.roll(model_values[1, emb_idx, :], roll_idx[pos_idx])
            axs[1].set_ylabel('Model B value')
            axs[1].plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[1].plot_color, markersize = 1, label='Without\nnbhd')
            axs[1].plot(range(0,noc), [models[1].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            # axs[2].set_title('Model ' + models[1].index_string + ', without nbhd', fontsize=font_sizes['BIGGER_SIZE'])
            axs[1].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])
            
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
            axs[1].text(368, text_yloc, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')
            
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
    
    noc = list_of_embryos[0].number_of_cells
    half_noc = int(np.floor(noc / 2))
    
    # create directory
    if directory_name[-1] != '/':
        directory_name = directory_name + '/'
    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)        # I have no idea why this is necessary
    
    inset_axes = [10, 11]
        
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
            plot_desired = np.roll(embryo.desired,roll_idx[pos_idx])
            desired_color = 'C4'  # purple
            gray_where = np.array([i == -1 for i in plot_desired])
            certain_where = np.array([i == 1 for i in plot_desired])
            
            protein_ymax = 1.4
            axs[0].plot(range(0,noc), plot_inducer, linewidth=2, marker=None, color='C0', markersize = 1, label='Inducer')
            axs[0].plot(range(0,noc), plot_inhibitor, linewidth=2, linestyle='dashdot', marker=None, color='C3', markersize = 1, label="Inhibitor")
            axs[0].set_ylim([0,protein_ymax])
            axs[0].set_ylabel('Protein conc.')
            
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
            axs[0].fill_between(np.arange(noc), bar_bottom, bar_top, where=certain_where, facecolor=desired_color, edgecolor=desired_color, alpha=1, step='mid', linewidth=1)
            
            axs[0].legend(loc='upper right', fontsize=font_sizes['SMALLEST_SIZE'])
            
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
            axs[1].set_ylabel('Model A value')
            axs[1].plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[0].plot_color, markersize = 1, label='With\nnbhd')
            axs[1].plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            axs[1].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])
            
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
            
            if emb_idx in inset_axes:
                
                # inset axes....
                axins = axs[1].inset_axes([0.02, 0.25, 0.35, 0.3])
                axins.plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[0].plot_color, markersize = 1)
                axins.plot(range(0,noc), [models[0].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
                # axins.imshow(Z2, extent=extent, interpolation="nearest", origin="lower")
                # sub region of the original image
                x1, x2, y1, y2 = 230, 370, 0, model_ymax
                axins.set_xlim(x1, x2)
                axins.set_ylim(y1, y2)
                axins.set_xticklabels('')
                axins.set_yticklabels('')

                axs[1].indicate_inset_zoom(axins)

        
            ''' model B '''
            plot_model = np.roll(model_values[1, emb_idx, :], roll_idx[pos_idx])
            axs[2].set_ylabel('Model B value')
            axs[2].plot(range(0,noc), plot_model, linewidth=2, marker=None, color=models[1].plot_color, markersize = 1, label='Without\nnbhd')
            axs[2].plot(range(0,noc), [models[1].threshold for i in range(noc)], '--', linewidth=0.8, marker=None, color='black', markersize = 1)
            axs[2].legend(loc='upper right', bbox_to_anchor=(1,legend_height), fontsize=font_sizes['SMALLEST_SIZE'])
            
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
            
            
            # clear yaxis for certain embryos
            yaxis_clear_embryos = [9,10,11,12,14,15,16,18,19,20,21,22]
            if emb_idx in yaxis_clear_embryos:
                for ax_idx in range(3):
                    axs[ax_idx].yaxis.label.set_visible(False)
                    axs[ax_idx].set_yticklabels([])
                    
            # clear legend for certain embryos
            legend_clear_embryos = [8,9,10,11,13,14,15,17,0,19,20,21]
            if emb_idx in legend_clear_embryos:
                for ax_idx in range(3):
                    axs[ax_idx].get_legend().remove()
                    
            # only add text for certain embryos
            
            add_text_embryos = [item for item in range(len(list_of_embryos)) if item not in legend_clear_embryos]
            if emb_idx in add_text_embryos:
                axs[0].text(408, text_yloc_pro, 'Desired streak', backgroundcolor='lightgray', color=desired_color, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')
                axs[1].text(388, text_yloc_A, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color_A, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')
                axs[2].text(388, text_yloc_B, 'Predicted streak', backgroundcolor='lightgray', color=brachyury_color_B, fontsize=font_sizes['SMALLEST_SIZE'], fontweight='bold')
            
            fig_full.tight_layout()
            
            filename = full_directory_name + embryo.index + embryo.name + '.jpg'
            fig_full.savefig(filename, format='jpg', dpi=300)

            plt.close(fig_full)