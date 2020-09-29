
initial_params = {}

IP = initial_params

IP['number_of_cells'] = 600 # should be divisible by 4

IP['posterior_cell'] = 0
IP['anterior_cell'] = int(IP['number_of_cells'] / 2)
IP['left_lateral'] = int(IP['number_of_cells'] / 4)
IP['right_lateral'] = int( (3 * IP['number_of_cells']) / 4)


""" Set concentrations for different stages """

'''stage X '''
IP['inhibitor_stageX_max'] = 0.4
IP['inhibitor_stageX_min'] = 0.1

IP['inducer_stageX_max'] = 0.7
IP['inducer_stageX_min'] = 0.1

'''stage XII '''
IP['inhibitor_stageXII_max'] = 0.65
IP['inhibitor_stageXII_min'] = 0.1

IP['inducer_stageXII_baseline'] = 0.05
IP['inducer_stageXII_conc_coeff'] = 70  
IP['inducer_stageXII_std'] = 35

'''stage 2 '''
IP['inhibitor_stage2_max'] = 0.9
IP['inhibitor_stage2_min'] = 0.1

IP['inducer_stage2_baseline'] = 0.02
IP['inducer_stage2_conc_coeff'] = 25  
IP['inducer_stage2_std'] = 30

initial_params = IP