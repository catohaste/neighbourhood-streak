import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats
import copy
import os.path

class Model:
    def __init__(self):
        self.plot_color = 'C0'
        self.index_string = 'X_'
        self.inducer_scaling = 1
        self.inhibitor_scaling = 1
        self.threshold = 0.1
        self.nbhd_size = 31 # odd number !!!!
        
class Embryo:
    
    class Protein:
        def __init__(self, start_conc):
            self.conc = start_conc.copy()
            self.number_of_cells = len(self.conc)
        
        def add_bead(self, cells_from_anterior, spread_cells, concentration, bead_width):
            
            anterior_cell = int(self.number_of_cells / 2)
            bead_centre = anterior_cell + cells_from_anterior
            bead_half_width = int(np.floor(bead_width / 2))
            if bead_width % 2 == 0:
                bead_cells = range(bead_centre - bead_half_width, bead_centre + bead_half_width)
            else:
                bead_cells = range(bead_centre - bead_half_width, bead_centre + bead_half_width + 1)
            bead_right = max(bead_cells)
            bead_left = min(bead_cells)

            dummy_cells = range(0,int(self.number_of_cells/2))
            bead_distribution_right = [concentration * np.exp((-1/spread_cells) * cell) for cell in dummy_cells]
            # bead_effect_limit = int(np.ceil(- spread_cells * np.log(0.005)))
            bead_effect_limit = int(np.ceil((self.number_of_cells - bead_width) / 2))
            
            # print(bead_centre, bead_effect_limit)
                       
            # set concentration for bead cells
            for bead_cell in bead_cells:
                bead_cell_idx = bead_cell % self.number_of_cells
                self.conc[bead_cell_idx] = self.conc[bead_cell_idx] + concentration
            for idx in range(1, bead_effect_limit):
        
                pos_cell = (bead_right + idx) % self.number_of_cells
                neg_cell = (bead_left - idx) % self.number_of_cells
        
                self.conc[pos_cell] =  self.conc[pos_cell] + bead_distribution_right[idx]
                self.conc[neg_cell] =  self.conc[neg_cell] + bead_distribution_right[idx]
                
            # this only really applies if bead_effect_limit is set to whole embryo - bead_width
            # however it won't make a massive difference otherwise
            if bead_width % 2 != 0:
                pos_cell = (bead_right + bead_effect_limit) % self.number_of_cells
                self.conc[pos_cell] =  self.conc[pos_cell] + bead_distribution_right[bead_effect_limit]
        
        def check_below_zero(self):
            self.conc = (self.conc <= 0) * 0.02 + np.multiply((self.conc > 0), self.conc)
        
    def __init__(self, name, number_of_cells):
        self.name = name
        self.index = ''
        self.number_of_cells = number_of_cells
        self.auto_model_ylim = [-0.2,0.15]
        self.plot_model_ylim = [-0.2,0.15]
        self.fig_title = name
        self.success = False
        
        self.streak_positions = []
        self.streakN = int
        
    def set_starting_conc(self, inducer_start_conc, inhibitor_start_conc):
        self.inducer = self.Protein(inducer_start_conc)
        self.inhibitor = self.Protein(inhibitor_start_conc)
    
    def find_streaks(self):
        # clear/initiate streak data
        self.streakN = int
        self.streak_start_indices = []
        self.streak_lengths = []
        self.streak_centers = []
        self.streak_positions = []
        
        streak_cell_number = np.count_nonzero(self.brachyury)
        # indices where value changes
        change_indices = np.where(np.roll(self.brachyury,1)!=self.brachyury)[0]
        if np.count_nonzero(self.brachyury) is 0:
            # print('no streaks')
            self.streakN = 0
        elif len(change_indices) is 0:
            # print('ring streak')
            self.streakN = 1
            self.streak_positions = ['ring']
        else:
            streak_start_indices = [idx for idx in change_indices if self.brachyury[idx] == 1]
            streak_lengths = copy.deepcopy(streak_start_indices)
            streak_centers = copy.deepcopy(streak_start_indices)
            for idx, start_index in enumerate(streak_start_indices):
                counter = 1
                while self.brachyury[(start_index + counter) % self.number_of_cells] == 1:
                    counter += 1
                streak_lengths[idx] = counter
                streak_centers[idx] = (start_index + int(np.floor(counter / 2))) % self.number_of_cells
            streak_positions = copy.deepcopy(streak_centers)
            for idx, center in enumerate(streak_centers):
                if center == 0:
                    streak_positions[idx] = 'posterior'
                elif 0 < center < self.number_of_cells / 4:
                    streak_positions[idx] = 'post-left'
                elif center == self.number_of_cells / 4 :
                    streak_positions[idx] = 'left'
                elif self.number_of_cells / 4 < center < self.number_of_cells / 2:
                    streak_positions[idx] = 'ant-left'
                elif center == self.number_of_cells / 2:
                    streak_positions[idx] = 'anterior'
                elif self.number_of_cells / 2 < center < (3*self.number_of_cells) / 4:
                    streak_positions[idx] = 'ant-right'
                elif center == (3*self.number_of_cells) / 4:
                    streak_positions[idx] = 'right'
                elif (3*self.number_of_cells) / 4 < center < self.number_of_cells:
                    streak_positions[idx] = 'post-right'
                else:
                    print('find_streaks Something has gone wrong', self.index)
                    
            self.streakN = len(streak_start_indices)
            self.streak_start_indices = streak_start_indices
            self.streak_lengths = streak_lengths
            self.streak_centers = streak_centers
            self.streak_positions = streak_positions
        
class Experiment:
    def __init__(self, name):
        self.name = name
        self.embryos = []

    def find_plot_model_ylim(self):

        plot_model_ymin = np.min([embryo.auto_model_ylim[0] for embryo in self.embryos])
        plot_model_ymax = np.max([embryo.auto_model_ylim[1] for embryo in self.embryos])

        for embryo in self.embryos:
            embryo.plot_model_ylim = [plot_model_ymin,plot_model_ymax]
 