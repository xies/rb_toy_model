#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:16:39 2020

@author: xies
"""


import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pylab as plt
from numpy import random
from scipy import ndimage, stats
from scipy.interpolate import UnivariateSpline
import statsmodels.api as sm
import simulation

# 1) Initial conditions / parameters
Ncells = 1000
Niter = 100

emparams = pd.read_pickle('/Users/xies/Box/HMECs/RB toy model/emp_parameters.pkl')
trans_params = pd.read_pickle('/Users/xies/Box/HMECs/RB toy model/trans_params.pkl')

emparams = emparams.join(trans_params)
emparams.at['RB conc','G1S trans m'] = trans_params.loc['RB conc','G1S trans m']
emparams.at['RB conc','G1S trans b'] = trans_params.loc['RB conc','G1S trans b']
emparams.at['Time','Mean G2 duration'] = trans_params.loc['Time','Mean G2 duration']
emparams.at['Time','Std G2 duration'] = trans_params.loc['Time','Std G2 duration']

#%%
# @todo: initialize ag G1/S
# Initialize each cell as a DataFrame
# population is handled currently as a list of DFs
population = []
for i in range(Ncells):

    # @todo: fix zero size cells
    init_size = random.lognormal(mean = emparams.loc['Size','Mean G1S'],sigma = emparams.loc['Size','Std G1S'])
    init_rb = random.lognormal(mean = emparams.loc['RB','Mean G1S'],sigma = emparams.loc['RB','Std G1S'])
    rb_conc = init_rb / init_size
    cell = pd.DataFrame( columns = ['CellID','Time','Size','RB','RB conc','Phase',
                                    'Birth size','Div size','G1S size','G1S time','Div time','G1S RB conc'],
                        index=pd.RangeIndex(Niter))
    cell = cell.fillna(np.nan)
    
    init_cell = pd.Series({'CellID':i,'Time':0,'Size':init_size,'RB':init_rb,'RB conc':rb_conc,
                           'Phase':'G2','G1S size':init_size})
    cell.at[0] = init_cell
    
    population.append(cell)
    
    
#%% 

# Start stimulation
for t in np.arange(Niter)+1:
    for i in range(Ncells):
        
        # Go to previous time point and grab cell to work on
        this_cell_now = population[i].loc[t-1].copy()
        
        # Check if cell divided
        if this_cell_now['Phase'] == 'None':
            this_cell_next = this_cell_now
        else:
            this_cell_next = simulation.advance_dt(this_cell_now,emparams)
        
        population[i].at[t] = this_cell_next
        
        
             
#%% Retrieve each datafield into dataframe
        
size = np.array([ cell['Size'] for cell in population])
rb = np.array([ cell['RB'] for cell in population])
rb_conc = np.array([ cell['RB conc'] for cell in population])

birth_size = np.array([ cell.iloc[-1]['Birth size'] for cell in population ])
birth_size = np.array([ cell.iloc[0]['RB conc'] for cell in population ])
g1s_size = np.array([ cell.iloc[-1]['G1S size'] for cell in population ])
div_size = np.array([ cell.iloc[-1]['Div size'] for cell in population ])
g1s_time = np.array([ cell.iloc[-1]['G1S time'] for cell in population ])
div_time = np.array([ cell.iloc[-1]['Div time'] for cell in population ])
g1s_rb_conc = np.array([ cell.iloc[-1]['G1S RB conc'] for cell in population ]) 

        





