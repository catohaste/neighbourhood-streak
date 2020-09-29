from pydream.parameters import SampledParam
import scipy.stats as stats


param_priors_dict = {
    '$b_B$': SampledParam(stats.uniform, loc=-3.0, scale=6.0),
    '$b_V$': SampledParam(stats.uniform, loc=-3.0, scale=6.0),
    'threshold': SampledParam(stats.uniform, loc=0.0, scale=0.5),
    'n': SampledParam(stats.uniform, loc=2.0, scale=73.0),
    
    '$b_B^A$': SampledParam(stats.uniform, loc=-3.0, scale=6.0),
    '$b_B^B$': SampledParam(stats.uniform, loc=-3.0, scale=6.0),
    '$b_V^A$': SampledParam(stats.uniform, loc=-3.0, scale=6.0),
    '$b_V^B$': SampledParam(stats.uniform, loc=-3.0, scale=6.0),
    'threshold$^A$': SampledParam(stats.uniform, loc=0.0, scale=0.5),
    'threshold$^B$': SampledParam(stats.uniform, loc=0.0, scale=0.5),
    
    'activin_2_lambda': SampledParam(stats.uniform, loc=0.025, scale=1.0),
    'activin_10_lambda': SampledParam(stats.uniform, loc=0.025, scale=1.0),
    'activin_2_conc': SampledParam(stats.uniform, loc=0.0, scale=1.0),
    'activin_10_conc': SampledParam(stats.uniform, loc=0.0, scale=1.0),
    'activin_spread': SampledParam(stats.uniform, loc=1.0, scale=42.0),
    'activin_2_spread': SampledParam(stats.uniform, loc=1.0, scale=42.0),
    'activin_10_spread': SampledParam(stats.uniform, loc=1.0, scale=42.0),
    'bmp4_50_lambda': SampledParam(stats.uniform, loc=0.025, scale=1.0),
    'bmp4_25_lambda': SampledParam(stats.uniform, loc=0.025, scale=1.0),
    'bmp4_12_lambda': SampledParam(stats.uniform, loc=0.025, scale=1.0),
    'bmp4_50_conc': SampledParam(stats.uniform, loc=0.0, scale=1.0),
    'bmp4_25_conc': SampledParam(stats.uniform, loc=0.0, scale=1.0),
    'bmp4_12_conc': SampledParam(stats.uniform, loc=0.0, scale=1.0),
    'bmp4_6_conc': SampledParam(stats.uniform, loc=0.0, scale=1.0),
    'bmp4_spread': SampledParam(stats.uniform, loc=1.0, scale=42.0),
    'bmp4_50_spread': SampledParam(stats.uniform, loc=1.0, scale=42.0),
    'bmp4_25_spread': SampledParam(stats.uniform, loc=1.0, scale=42.0),
    'bmp4_12_spread': SampledParam(stats.uniform, loc=1.0, scale=42.0),
    'bmp4_6_spread': SampledParam(stats.uniform, loc=1.0, scale=42.0)
}

param_lims_dict = {
    '$b_B$': (-3, 3),
    '$b_V$': (-3, 3),
    'threshold': (0,0.5),
    'n': (2, 75),
    
    '$b_B^A$': (-3, 3),
    '$b_B^B$': (-3, 3),
    '$b_V^A$': (-3, 3),
    '$b_V^B$': (-3, 3),
    'threshold$^A$': (0,0.5),
    'threshold$^B$': (0,0.5),
    
    'activin_2_lambda': (0,1),
    'activin_10_lambda': (0,1),
    'activin_2_conc': (0,1),
    'activin_10_conc': (0,1),
    'activin_spread': (1,43),
    'activin_2_spread': (1,43),
    'activin_10_spread': (1,43),
    'bmp4_50_lambda': (0,1),
    'bmp4_25_lambda': (0,1),
    'bmp4_12_lambda': (0,1),
    'bmp4_50_conc': (0,1),
    'bmp4_25_conc': (0,1),
    'bmp4_12_conc': (0,1),
    'bmp4_6_conc': (0,1),
    'bmp4_spread': (1,43),
    'bmp4_50_spread': (1,43),
    'bmp4_25_spread': (1,43),
    'bmp4_12_spread': (1,43),
    'bmp4_6_spread': (1,43)
}

axes_labels_dict = {
    '$b_B$': '$b_B$',
    '$b_V$': '$b_V$',
    'threshold': 'threshold',
    'n': 'n',
    
    '$b_B^A$': '$b_B^A$',
    '$b_B^B$': '$b_B^B$',
    '$b_V^A$': '$b_V^A$',
    '$b_V^B$': '$b_V^B$',
    'threshold$^A$': 'threshold$^A$',
    'threshold$^B$': 'threshold$^B$',

    'activin_2_lambda': 'activin\n2\nlambda',
    'activin_10_lambda': 'activin\n10\nlambda',
    'activin_2_conc': 'activin\n2\nconc',
    'activin_10_conc': 'activin\n10\nconc',
    'activin_spread': 'activin\nspread',
    'activin_2_spread': 'activin\n2\nspread',
    'activin_10_spread': 'activin\n10\nspread',
    'bmp4_50_lambda': 'bmp4\n50\nlambda',
    'bmp4_25_lambda': 'bmp4\n25\nlambda',
    'bmp4_12_lambda': 'bmp4\n12\nlambda',
    'bmp4_50_conc': 'bmp4\n50\nconc',
    'bmp4_25_conc': 'bmp4\n25\nconc',
    'bmp4_12_conc': 'bmp4\n12\nconc',
    'bmp4_6_conc': 'bmp4\n6\nconc',
    'bmp4_spread': 'bmp4\nspread',
    'bmp4_50_spread': 'bmp4\n50\nspread',
    'bmp4_25_spread': 'bmp4\n25\nspread',
    'bmp4_12_spread': 'bmp4\n12\nspread',
    'bmp4_6_spread': 'bmp4\n6\nspread'
}
