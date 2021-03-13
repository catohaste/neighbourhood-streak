import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats
import copy
import os.path

from classes import Experiment

def parabola( protein_conc, stationary_loc, min_conc, max_conc , positive ):

    if not positive:
        temp_min = copy.deepcopy(min_conc)
        min_conc = copy.deepcopy(max_conc)
        max_conc = temp_min

    anterior_cell = int(len(protein_conc)/2)
    b = stationary_loc
    c = min_conc
    a = (max_conc - min_conc) / (anterior_cell**2)

    cells = range(0 - anterior_cell + stationary_loc, len(protein_conc) - anterior_cell + stationary_loc)
    for cell_id in cells:
        new_cell_id = cell_id % len(protein_conc)
        protein_conc[new_cell_id] = a * ((cell_id - b)**2) + c

    return protein_conc


def define_initial_protein_concentrations(initial_params):
    
    IP = initial_params
    
    ''' define RVs '''
    inducer_stageXII_rv = stats.norm(IP['anterior_cell'], IP['inducer_stageXII_std'])
    inducer_stage2_right_rv = stats.expon(IP['posterior_cell'], IP['inducer_stage2_std'])


    '''initialize'''
    inhibitor_postFocus_stageX = np.ndarray((IP['number_of_cells']), dtype=float)
    inhibitor_postFocus_stageXII = np.ndarray((IP['number_of_cells']), dtype=float)
    inhibitor_postFocus_stage2 = np.ndarray((IP['number_of_cells']), dtype=float)
    inhibitor_antFocus_stageX = np.ndarray((IP['number_of_cells']), dtype=float)
    inhibitor_antFocus_stageXII = np.ndarray((IP['number_of_cells']), dtype=float)
    inhibitor_antFocus_stage2 = np.ndarray((IP['number_of_cells']), dtype=float)

    inducer_postFocus_stageX = np.ndarray((IP['number_of_cells']), dtype=float)
    inducer_postFocus_stageXII = np.ndarray((IP['number_of_cells']), dtype=float)
    inducer_postFocus_stage2 = np.ndarray((IP['number_of_cells']), dtype=float)
    inducer_antFocus_stageX = np.ndarray((IP['number_of_cells']), dtype=float)

    ''' use RV to define inducer stage XII '''
    for i in range(IP['number_of_cells']):
        inducer_postFocus_stageXII[i] = IP['inducer_stageXII_baseline'] + inducer_stageXII_rv.pdf(i) * IP['inducer_stageXII_conc_coeff']
    inducer_postFocus_stageXII = np.roll(inducer_postFocus_stageXII, IP['anterior_cell'])

    ''' use RV to define inducer stage 2 '''
    for i in range(1, IP['anterior_cell']):
        pos_idx = i
        neg_idx = IP['number_of_cells'] - i

        inducer_postFocus_stage2[pos_idx] = IP['inducer_stage2_baseline'] + inducer_stage2_right_rv.pdf(i) * IP['inducer_stage2_conc_coeff']
        inducer_postFocus_stage2[neg_idx] = IP['inducer_stage2_baseline'] + inducer_stage2_right_rv.pdf(i) * IP['inducer_stage2_conc_coeff']

    inducer_postFocus_stage2[IP['posterior_cell']] = IP['inducer_stage2_baseline'] + inducer_stage2_right_rv.pdf(IP['posterior_cell']) * IP['inducer_stage2_conc_coeff']
    inducer_postFocus_stage2[IP['anterior_cell']] = IP['inducer_stage2_baseline'] + inducer_stage2_right_rv.pdf(IP['anterior_cell']) * IP['inducer_stage2_conc_coeff']

    
    ''' use parabola to define inducer stage X and inhibitor '''
    inducer_postFocus_stageX = parabola( inducer_postFocus_stageX, IP['posterior_cell'], IP['inducer_stageX_min'], IP['inducer_stageX_max'], False)
    inducer_antFocus_stageX = parabola( inducer_antFocus_stageX, IP['anterior_cell'], IP['inducer_stageX_min'], IP['inducer_stageX_max'], True)

    inhibitor_antFocus_stageX = parabola( inhibitor_antFocus_stageX, IP['anterior_cell'], IP['inhibitor_stageX_min'] , IP['inhibitor_stageX_max'], False)
    inhibitor_antFocus_stageXII = parabola( inhibitor_antFocus_stageXII, IP['anterior_cell'], IP['inhibitor_stageXII_min'], IP['inhibitor_stageXII_max'], False)
    inhibitor_antFocus_stage2 = parabola( inhibitor_antFocus_stage2, IP['anterior_cell'], IP['inhibitor_stage2_min'], IP['inhibitor_stage2_max'], False)

    inhibitor_postFocus_stageX = parabola( inhibitor_postFocus_stageX, IP['posterior_cell'], IP['inhibitor_stageX_min'], IP['inhibitor_stageX_max'], True)
    inhibitor_postFocus_stageXII = parabola( inhibitor_postFocus_stageXII, IP['posterior_cell'], IP['inhibitor_stageXII_min'], IP['inhibitor_stageXII_max'], True)
    inhibitor_postFocus_stage2 = parabola( inhibitor_postFocus_stage2, IP['posterior_cell'], IP['inhibitor_stage2_min'], IP['inhibitor_stage2_max'], True)
    
    IC = {}

    IC["inhibitor_postFocus_stageX"] = inhibitor_postFocus_stageX
    IC["inhibitor_postFocus_stageXII"] = inhibitor_postFocus_stageXII
    IC["inhibitor_postFocus_stage2"] = inhibitor_postFocus_stage2
    IC["inhibitor_antFocus_stageX"] = inhibitor_antFocus_stageX
    IC["inhibitor_antFocus_stageXII"] = inhibitor_antFocus_stageXII
    IC["inhibitor_antFocus_stage2"] = inhibitor_antFocus_stage2

    IC["inducer_postFocus_stageX"] = inducer_postFocus_stageX
    IC["inducer_postFocus_stageXII"] = inducer_postFocus_stageXII
    IC["inducer_postFocus_stage2"] = inducer_postFocus_stage2
    IC["inducer_antFocus_stageX"] = inducer_antFocus_stageX

    initial_concentrations = IC
   
    return initial_concentrations
    

def add_streak(streak_embryo, streak_center, streak_width):

    half_width = int(np.floor(streak_width / 2))
    indices = range(streak_center - half_width, streak_center + half_width + 1)
        
    for idx in indices:
        streak_embryo[idx] = 1
        
    return streak_embryo
   
    
def setup_embryos(list_of_embryos, Model, initial_concentrations):
    
    IC = initial_concentrations
    
    inhibitor_postFocus_stageX = IC["inhibitor_postFocus_stageX"]
    inhibitor_postFocus_stageXII = IC["inhibitor_postFocus_stageXII"]
    inhibitor_postFocus_stage2 = IC["inhibitor_postFocus_stage2"]
    inhibitor_antFocus_stageX = IC["inhibitor_antFocus_stageX"]
    inhibitor_antFocus_stageXII = IC["inhibitor_antFocus_stageXII"]
    inhibitor_antFocus_stage2 = IC["inhibitor_antFocus_stage2"]
    inducer_postFocus_stageX = IC["inducer_postFocus_stageX"]
    inducer_postFocus_stageXII = IC["inducer_postFocus_stageXII"]
    inducer_postFocus_stage2 = IC["inducer_postFocus_stage2"]
    inducer_antFocus_stageX = IC["inducer_antFocus_stageX"]
    
    cell_pellet_width = Model.bead_params['cell_pellet_width']
    cell_pellet_spread = Model.bead_params['cell_pellet_spread']
    vg1_cell_conc = Model.bead_params['vg1_cell_conc']
    bmp4_cell_conc = Model.bead_params['bmp4_cell_conc']
    
    afigel_width = Model.bead_params['afigel_width']
    afigel_50_spread = Model.bead_params['afigel_50_spread']
    afigel_25_spread = Model.bead_params['afigel_25_spread']
    afigel_12_spread = Model.bead_params['afigel_12_spread']
    afigel_6_spread = Model.bead_params['afigel_6_spread']
    bmp4_50_conc = Model.bead_params['bmp4_50_conc']
    bmp4_25_conc = Model.bead_params['bmp4_25_conc']
    bmp4_12_conc = Model.bead_params['bmp4_12_conc']
    bmp4_6_conc = Model.bead_params['bmp4_6_conc']
    
    heparin_width = Model.bead_params['heparin_width']
    heparin_10_spread = Model.bead_params['heparin_10_spread']
    heparin_2_spread = Model.bead_params['heparin_2_spread']
    activin_10_conc = Model.bead_params['activin_10_conc']
    activin_2_conc = Model.bead_params['activin_2_conc']
    
    AG1X2_spread = Model.bead_params['AG1X2_spread']
    AG1X2_width = Model.bead_params['AG1X2_width']
    DM_conc = Model.bead_params['DM_conc']
    SB_conc = Model.bead_params['SB_conc']
    
    noc = list_of_embryos[0].number_of_cells
    streak_width = 13 # probably should be odd
    no_streak = [-1 for cell in range(noc)]
    
    
    list_of_embryos[0].index = '00_'
    list_of_embryos[0].name = 'Stage EG&K XII, antfocus'
    list_of_embryos[0].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[0].target = copy.deepcopy(no_streak)
    list_of_embryos[0].target = add_streak(list_of_embryos[0].target, 0, streak_width)  
    
    
    list_of_embryos[1].index = '01_'
    list_of_embryos[1].name = '2 Vg1 pellets, ant'
    list_of_embryos[1].fig_title = '2 cVG1 pellets'
    list_of_embryos[1].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[1].inducer.add_bead(-15, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[1].inducer.add_bead(15, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[1].target = copy.deepcopy(no_streak)
    list_of_embryos[1].target = add_streak(list_of_embryos[1].target, 0, streak_width)
    list_of_embryos[1].target = add_streak(list_of_embryos[1].target, 300, streak_width)

    list_of_embryos[2].index = '02_'
    list_of_embryos[2].name = '4 Vg1 pellets, ant'
    list_of_embryos[2].fig_title = '4 cVG1 pellets'
    list_of_embryos[2].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[2].inducer.add_bead(-45, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[2].inducer.add_bead(-15, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[2].inducer.add_bead(15, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[2].inducer.add_bead(45, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[2].target = copy.deepcopy(no_streak)
    list_of_embryos[2].target = add_streak(list_of_embryos[2].target, 0, streak_width)

    list_of_embryos[3].index = '03_'
    list_of_embryos[3].name = '4 Vg1 pellets, ant, w ctrl'
    list_of_embryos[3].fig_title = '4 cVG1 pellets + space'
    list_of_embryos[3].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[3].inducer.add_bead(-60, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[3].inducer.add_bead(-30, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[3].inducer.add_bead(30, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[3].inducer.add_bead(60, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[3].target = copy.deepcopy(no_streak)
    list_of_embryos[3].target = add_streak(list_of_embryos[3].target, 0, streak_width)
    list_of_embryos[3].target = add_streak(list_of_embryos[3].target, 255, streak_width)
    list_of_embryos[3].target = add_streak(list_of_embryos[3].target, 345, streak_width)

    list_of_embryos[4].index = '04_'
    list_of_embryos[4].name = '4 Vg1 pellets, ant, w BMP4'
    list_of_embryos[4].fig_title = '4 cVG1 pellets + BMP4'
    list_of_embryos[4].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[4].inducer.add_bead(-60, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[4].inducer.add_bead(-30, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[4].inhibitor.add_bead(0, cell_pellet_spread, bmp4_cell_conc, cell_pellet_width)
    list_of_embryos[4].inducer.add_bead(30, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[4].inducer.add_bead(60, cell_pellet_spread, vg1_cell_conc, cell_pellet_width)
    list_of_embryos[4].target = copy.deepcopy(no_streak)
    list_of_embryos[4].target = add_streak(list_of_embryos[4].target, 0, streak_width)
    list_of_embryos[4].target = add_streak(list_of_embryos[4].target, 255, streak_width)
    list_of_embryos[4].target = add_streak(list_of_embryos[4].target, 345, streak_width)
    
    
    bead_disp = 50 # microbead displacement from anterior centre (positive and negative)
    list_of_embryos[5].index = '05_'
    list_of_embryos[5].name = 'C-A-C, 0-10-0'
    list_of_embryos[5].fig_title = 'C-A-C'
    list_of_embryos[5].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[5].inducer.add_bead(0, heparin_10_spread, activin_10_conc, heparin_width)
    list_of_embryos[5].target = copy.deepcopy(no_streak)
    list_of_embryos[5].target = add_streak(list_of_embryos[5].target, 300, streak_width)
    list_of_embryos[5].target = add_streak(list_of_embryos[5].target, 0, streak_width)
    
    list_of_embryos[6].index = '06_'
    list_of_embryos[6].name = 'A-A-A, 10-10-10'
    list_of_embryos[6].fig_title = 'A-A-A'
    list_of_embryos[6].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[6].inducer.add_bead(-bead_disp, heparin_10_spread, activin_10_conc, heparin_width)
    list_of_embryos[6].inducer.add_bead(0, heparin_10_spread, activin_10_conc, heparin_width)
    list_of_embryos[6].inducer.add_bead(bead_disp, heparin_10_spread, activin_10_conc, heparin_width)
    list_of_embryos[6].target = copy.deepcopy(no_streak)
    list_of_embryos[6].target = add_streak(list_of_embryos[6].target, 0, streak_width)
    
    list_of_embryos[7].index = '07_'
    list_of_embryos[7].name = 'A-C-A, 10-0-10'
    list_of_embryos[7].fig_title = 'A-C-A'
    list_of_embryos[7].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[7].inducer.add_bead(-bead_disp, heparin_10_spread, activin_10_conc, heparin_width)
    list_of_embryos[7].inducer.add_bead(bead_disp, heparin_10_spread, activin_10_conc, heparin_width)
    list_of_embryos[7].target = copy.deepcopy(no_streak)
    list_of_embryos[7].target = add_streak(list_of_embryos[7].target, 300 - bead_disp , streak_width)
    list_of_embryos[7].target = add_streak(list_of_embryos[7].target, 300 + bead_disp , streak_width)
    list_of_embryos[7].target = add_streak(list_of_embryos[7].target, 0, streak_width)

    list_of_embryos[8].index = '08_'
    list_of_embryos[8].name = 'A-B-A, 10-12-10'
    list_of_embryos[8].fig_title = 'A-B-A'
    list_of_embryos[8].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[8].inducer.add_bead(-bead_disp, heparin_10_spread, activin_10_conc, heparin_width)
    list_of_embryos[8].inhibitor.add_bead(0, afigel_12_spread, bmp4_12_conc, afigel_width)
    list_of_embryos[8].inducer.add_bead(bead_disp, heparin_10_spread, activin_10_conc, heparin_width)
    list_of_embryos[8].target = copy.deepcopy(no_streak)
    list_of_embryos[8].target = add_streak(list_of_embryos[8].target, 300 - bead_disp, streak_width)
    list_of_embryos[8].target = add_streak(list_of_embryos[8].target, 300 + bead_disp, streak_width)
    list_of_embryos[8].target = add_streak(list_of_embryos[8].target, 0, streak_width)
    
    list_of_embryos[9].index = '09_'
    list_of_embryos[9].name = 'A-B-A, 10-25-10'
    list_of_embryos[9].fig_title = 'A-B*-A'
    list_of_embryos[9].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[9].inducer.add_bead(-bead_disp, heparin_10_spread, activin_10_conc, heparin_width)
    list_of_embryos[9].inhibitor.add_bead(0, afigel_25_spread, bmp4_25_conc, afigel_width)
    list_of_embryos[9].inducer.add_bead(bead_disp, heparin_10_spread, activin_10_conc, heparin_width)
    list_of_embryos[9].target = copy.deepcopy(no_streak)
    list_of_embryos[9].target = add_streak(list_of_embryos[9].target, 0, streak_width)
    
    
    list_of_embryos[10].index = '10_'
    list_of_embryos[10].name = 'DM bead, ant'
    list_of_embryos[10].fig_title = 'DM bead'
    list_of_embryos[10].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[10].inhibitor.add_bead(40, AG1X2_spread, DM_conc, AG1X2_width)
    list_of_embryos[10].target = copy.deepcopy(no_streak)
    list_of_embryos[10].target = add_streak(list_of_embryos[10].target, 340, streak_width)
    list_of_embryos[10].target = add_streak(list_of_embryos[10].target, 0, streak_width)
    
    list_of_embryos[11].index = '11_'
    list_of_embryos[11].name = 'BMP4 bead, ant'
    list_of_embryos[11].fig_title = 'BMP4 bead'
    list_of_embryos[11].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[11].inhibitor.add_bead(0, afigel_12_spread, bmp4_12_conc, afigel_width)
    list_of_embryos[11].target = copy.deepcopy(no_streak)
    list_of_embryos[11].target = add_streak(list_of_embryos[11].target, 0, streak_width)
    
    
    list_of_embryos[12].index = '12_'
    list_of_embryos[12].name = 'B-C-B, 6-0-6'
    list_of_embryos[12].fig_title = 'B-C-B'
    list_of_embryos[12].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[12].inhibitor.add_bead(-bead_disp, afigel_6_spread, bmp4_6_conc, afigel_width)
    list_of_embryos[12].inhibitor.add_bead(bead_disp, afigel_6_spread, bmp4_6_conc, afigel_width)
    list_of_embryos[12].target = copy.deepcopy(no_streak)
    list_of_embryos[12].target = add_streak(list_of_embryos[12].target, 0, streak_width)

    list_of_embryos[13].index = '13_'
    list_of_embryos[13].name = 'C-A-C, 0-2-0'
    list_of_embryos[13].fig_title = 'C-A-C'
    list_of_embryos[13].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[13].inducer.add_bead(0, heparin_2_spread, activin_2_conc, heparin_width)
    list_of_embryos[13].target = copy.deepcopy(no_streak)
    list_of_embryos[13].target = add_streak(list_of_embryos[13].target, 0, streak_width)
    
    list_of_embryos[14].index = '14_'
    list_of_embryos[14].name = 'B-A-B, 6-2-6'
    list_of_embryos[14].fig_title = 'B-A-B'
    list_of_embryos[14].set_starting_conc(inducer_postFocus_stageXII, inhibitor_antFocus_stageXII)
    list_of_embryos[14].inhibitor.add_bead(-bead_disp, afigel_6_spread, bmp4_6_conc, afigel_width)
    list_of_embryos[14].inducer.add_bead(0, heparin_2_spread, activin_2_conc, heparin_width)
    list_of_embryos[14].inhibitor.add_bead(bead_disp, afigel_6_spread, bmp4_6_conc, afigel_width)
    list_of_embryos[14].target = copy.deepcopy(no_streak)
    list_of_embryos[14].target = add_streak(list_of_embryos[14].target, 300, streak_width)
    list_of_embryos[14].target = add_streak(list_of_embryos[14].target, 0, streak_width)
    
    
    # check none of the protein concentrations drop below 0
    # if they do, fix the concentration at 0.02, to avoid inf problems
    for embryo in list_of_embryos:
        embryo.inhibitor.check_below_zero()
        embryo.inducer.check_below_zero()

    return list_of_embryos
    

def run_model(Embryo, Model):
    
    # unload parameter values
    noc = Embryo.number_of_cells
    
    # calculate values of F_i and G_i (inducer- and inhibitor-bound fraction of SMAD4)
    Embryo.inducer_fraction = [(Model.inducer_scaling * value) / (1 + Model.inducer_scaling * value + Model.inhibitor_scaling * Embryo.inhibitor.conc[idx]) for idx, value in enumerate(Embryo.inducer.conc)]
    Embryo.inhibitor_fraction = [(Model.inhibitor_scaling * value) / (1 + Model.inducer_scaling * Embryo.inducer.conc[idx] + Model.inhibitor_scaling * value) for idx, value in enumerate(Embryo.inhibitor.conc)]
    
    model_string_options = ['inducer_SMAD_nbhd', 'inducer_SMAD']
    #['inhibitor_gradient', 'inducer_bound_fraction', 'inducer_bound_fraction_peak']
    if Model.name == 'inducer_SMAD_nbhd':
        
        half_nbhd = int(np.floor(Model.nbhd_size / 2))
    
        inducer_fraction_nbhd_avg = copy.deepcopy(Embryo.inducer_fraction)
        for center in range(noc):
            
            nbhd_indices = list(range(center - half_nbhd, center + half_nbhd + 1))
            nbhd_indices.remove(center)
            nbhd_indices = [index % noc for index in nbhd_indices]
        
            inducer_fraction_nbhd = [Embryo.inducer_fraction[index] for index in nbhd_indices]
            inducer_fraction_nbhd_avg[center] = np.sum(inducer_fraction_nbhd) / (Model.nbhd_size - 1)
    
        Embryo.model_value = [(value - inducer_fraction_nbhd_avg[idx]) / value for idx, value in enumerate(Embryo.inducer_fraction)]
        Embryo.brachyury = [int(value > Model.threshold) for value in Embryo.model_value]
        
    elif Model.name == 'inducer_SMAD':
        
        Embryo.model_value = Embryo.inducer_fraction
        Embryo.brachyury = [int(value > Model.threshold) for value in Embryo.model_value]
    
    else:
        message_string = ['The input argument model_string must be one of:\n' + i + '\n' for i in model_string_options]
        print(message_string)
        
    max_model_value = np.max([np.max(Embryo.model_value), Model.threshold])
    model_yheight = max_model_value - np.min(Embryo.model_value)
    Embryo.auto_model_ylim = [np.min(Embryo.model_value) - 0.05*model_yheight, max_model_value + 0.05*model_yheight]
    
def check_embryos_success(list_of_embryos):
    
    successN = 0
    failureN = 0
    for embryo in list_of_embryos:
        embryo.success = False
        
        
    if list_of_embryos[0].streakN == 1 and 'posterior' in list_of_embryos[0].streak_positions:
        list_of_embryos[0].success = True
        successN += 1
    else:
        failureN += 1
    
    
    if 'anterior' in list_of_embryos[1].streak_positions:
        list_of_embryos[1].success = True
        successN += 1
    else:
        failureN += 1
    
    if list_of_embryos[2].streakN == 0:
        list_of_embryos[2].success = True
        successN += 1
    else:
        failureN += 1
        
    if list_of_embryos[3].streakN == 3 and 'posterior' in list_of_embryos[3].streak_positions:
        list_of_embryos[3].success = True
        successN += 1
    elif list_of_embryos[3].streakN == 2 and 'posterior' not in list_of_embryos[3].streak_positions:
        list_of_embryos[3].success = True
        successN += 1
    else:
        failureN += 1
        
    if list_of_embryos[4].streakN == 3 and 'posterior' in list_of_embryos[4].streak_positions:
        list_of_embryos[4].success = True
        successN += 1
    elif list_of_embryos[4].streakN == 2 and 'posterior' not in list_of_embryos[4].streak_positions:
        list_of_embryos[4].success = True
        successN += 1
    else:
        failureN += 1
        
        
        
    if list_of_embryos[5].streakN == 1 and 'anterior' in list_of_embryos[5].streak_positions:
        list_of_embryos[5].success = True
        successN += 1
    elif list_of_embryos[5].streakN == 2 and 'posterior' in list_of_embryos[5].streak_positions:
        list_of_embryos[5].success = True
        successN += 1
    else:
        failureN += 1
    
    if list_of_embryos[6].streakN == 0:
        list_of_embryos[6].success = True
        successN += 1
    elif list_of_embryos[6].streakN == 1 and 'posterior' in list_of_embryos[6].streak_positions:
        list_of_embryos[6].success = True
        successN += 1
    else:
        failureN += 1
        
    if list_of_embryos[7].streakN == 3 and 'posterior' in list_of_embryos[7].streak_positions:
        list_of_embryos[7].success = True
        successN += 1
    elif list_of_embryos[7].streakN == 2 and 'posterior' not in list_of_embryos[7].streak_positions:
        list_of_embryos[7].success = True
        successN += 1
    else:
        failureN += 1
        
    if list_of_embryos[8].streakN == 3 and 'posterior' in list_of_embryos[8].streak_positions:
        list_of_embryos[8].success = True
        successN += 1
    elif list_of_embryos[8].streakN == 2 and 'posterior' not in list_of_embryos[8].streak_positions:
        list_of_embryos[8].success = True
        successN += 1
    else:
        failureN += 1
    
    if list_of_embryos[9].streakN == 0:
        list_of_embryos[9].success = True
        successN += 1    
    elif list_of_embryos[9].streakN == 1 and 'posterior' in list_of_embryos[9].streak_positions:
        list_of_embryos[9].success = True
        successN += 1
    else:
        failureN += 1
        
    
    if list_of_embryos[10].streakN == 1 and 'posterior' not in list_of_embryos[10].streak_positions:
        list_of_embryos[10].success = True
        successN += 1
    elif list_of_embryos[10].streakN == 2 and 'posterior' in list_of_embryos[10].streak_positions:
        list_of_embryos[10].success = True
        successN += 1
    else:
        failureN += 1
    
    if list_of_embryos[11].streakN == 0:
        list_of_embryos[11].success = True
        successN += 1
    elif list_of_embryos[11].streakN == 1 and 'posterior' in list_of_embryos[11].streak_positions:
        list_of_embryos[11].success = True
        successN += 1
    else:
        failureN += 1
        
    

    if list_of_embryos[12].streakN == 0:
        list_of_embryos[12].success = True
        successN += 1
    elif list_of_embryos[12].streakN == 1 and 'posterior' in list_of_embryos[12].streak_positions:
        list_of_embryos[12].success = True
        successN += 1
    else:
        failureN += 1
        
    if list_of_embryos[13].streakN == 0:
        list_of_embryos[13].success = True
        successN += 1
    elif list_of_embryos[13].streakN == 1 and 'posterior' in list_of_embryos[13].streak_positions:
        list_of_embryos[13].success = True
        successN += 1
    else:
        failureN += 1
        
    if list_of_embryos[14].streakN == 2 and 'posterior' in list_of_embryos[14].streak_positions:
        list_of_embryos[14].success = True
        successN += 1
    elif list_of_embryos[14].streakN == 1 and 'posterior' not in list_of_embryos[14].streak_positions:
        list_of_embryos[14].success = True
        successN += 1
    else:
        failureN += 1
        


    if successN + failureN != len(list_of_embryos):
        print('check_embryos_success something went wrong')
    
    return successN, failureN


def define_experiment_groups(list_of_embryos):
    
    expN = 5
    experiments = [Experiment('title') for i in range(expN)]
    
    experiments[0].name = 'normal embryo'
    experiments[0].embryos = [list_of_embryos[0]]

    experiments[1].name = 'cVg1 cell pellets in ant.'
    experiments[1].embryos = list_of_embryos[1:5]
    
    experiments[2].name = 'activin repeat of cell pellets in ant.'
    experiments[2].embryos = list_of_embryos[5:10]
    
    experiments[3].name = 'BMP bead in anterior'
    experiments[3].embryos = list_of_embryos[10:12]
    
    experiments[4].name = 'activin and BMP in ant., threshold'
    experiments[4].embryos = list_of_embryos[12:15]
    
    return experiments
    
def set_params_from_df(df, current_model):
    
    if '$b_B$' in df.columns:
        current_model.inhibitor_scaling = 10.0 ** df.iloc[0]['$b_B$']
    if '$b_V$' in df.columns:
        current_model.inducer_scaling = 10.0 ** df.iloc[0]['$b_V$']
    if 'threshold' in df.columns:
        current_model.threshold = df.iloc[0]['threshold']
    if 'n' in df.columns:
        current_model.nbhd_size = 2*np.floor(df.iloc[0]['n']) + 1
    if 'activin_conc' in df.columns:
        current_model.bead_params['activin_2_conc'] = df.iloc[0]['activin_conc']
        current_model.bead_params['activin_10_conc'] = df.iloc[0]['activin_conc']
    if 'activin_2_conc' in df.columns:
        current_model.bead_params['activin_2_conc'] = df.iloc[0]['activin_2_conc']
    if 'activin_10_conc' in df.columns:
        current_model.bead_params['activin_10_conc'] = df.iloc[0]['activin_10_conc']
    if 'activin_spread' in df.columns:
        current_model.bead_params['heparin_2_spread'] = df.iloc[0]['activin_spread']
        current_model.bead_params['heparin_10_spread'] = df.iloc[0]['activin_spread']
    if 'activin_2_spread' in df.columns:
        current_model.bead_params['heparin_2_spread'] = df.iloc[0]['activin_2_spread']
    if 'activin_10_spread' in df.columns:
        current_model.bead_params['heparin_10_spread'] = df.iloc[0]['activin_10_spread']
    if 'bmp4_conc' in df.columns:
        current_model.bead_params['bmp4_50_conc'] = df.iloc[0]['bmp4_conc']
        current_model.bead_params['bmp4_25_conc'] = df.iloc[0]['bmp4_conc']
        current_model.bead_params['bmp4_12_conc'] = df.iloc[0]['bmp4_conc']
        current_model.bead_params['bmp4_6_conc'] = df.iloc[0]['bmp4_conc']
    if 'bmp4_50_conc' in df.columns:
        current_model.bead_params['bmp4_50_conc'] = df.iloc[0]['bmp4_50_conc']
    if 'bmp4_25_conc' in df.columns:
        current_model.bead_params['bmp4_25_conc'] = df.iloc[0]['bmp4_25_conc']
    if 'bmp4_12_conc' in df.columns:
        current_model.bead_params['bmp4_12_conc'] = df.iloc[0]['bmp4_12_conc']
    if 'bmp4_6_conc' in df.columns:
        current_model.bead_params['bmp4_6_conc'] = df.iloc[0]['bmp4_6_conc']
    if 'bmp4_spread' in df.columns:
        current_model.bead_params['afigel_50_spread'] = df.iloc[0]['bmp4_spread']
        current_model.bead_params['afigel_25_spread'] = df.iloc[0]['bmp4_spread']
        current_model.bead_params['afigel_12_spread'] = df.iloc[0]['bmp4_spread']
        current_model.bead_params['afigel_6_spread'] = df.iloc[0]['bmp4_spread']
    if 'bmp4_50_spread' in df.columns:
        current_model.bead_params['afigel_50_spread'] = df.iloc[0]['bmp4_50_spread']
    if 'bmp4_25_spread' in df.columns:
        current_model.bead_params['afigel_25_spread'] = df.iloc[0]['bmp4_25_spread']
    if 'bmp4_12_spread' in df.columns:
        current_model.bead_params['afigel_12_spread'] = df.iloc[0]['bmp4_12_spread']
    if 'bmp4_6_spread' in df.columns:
        current_model.bead_params['afigel_6_spread'] = df.iloc[0]['bmp4_6_spread']
    if 'DM_conc' in df.columns:
        current_model.bead_params['DM_conc'] = df.iloc[0]['DM_conc']
    if 'AG1X2_spread' in df.columns:
        current_model.bead_params['AG1X2_spread'] = df.iloc[0]['AG1X2_spread']
    if 'cell_pellet_spread' in df.columns:
        current_model.bead_params['cell_pellet_spread'] = df.iloc[0]['cell_pellet_spread']
    if 'vg1_cell_conc' in df.columns:
        current_model.bead_params['vg1_cell_conc'] = df.iloc[0]['vg1_cell_conc']
    if 'bmp4_cell_conc' in df.columns:
        current_model.bead_params['bmp4_cell_conc'] = df.iloc[0]['bmp4_cell_conc']
        
    return current_model
    
        
    