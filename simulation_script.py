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
Ncells = 1
Niter = 100

emparams = pd.read_pickle('/Users/xies/Box/HMECs/RB toy model/emp_parameters.pkl')
trans_params = pd.read_pickle('/Users/xies/Box/HMECs/RB toy model/trans_params.pkl')

emparams = emparams.join(trans_params)
emparams.at['RB conc','G1S trans m'] = trans_params.loc['RB conc','G1S trans m']
emparams.at['RB conc','G1S trans b'] = trans_params.loc['RB conc','G1S trans b']


#%%
# Initialize cells as dictionaries
population = np.empty((Ncells,Niter+1),dtype=dict)
for i in range(Ncells):

    # @todo: fix zero size cells
    init_size = random.normal(loc = emparams.loc['Size','Mean birth'],scale = emparams.loc['Size','Std birth']/3)
    init_rb = random.normal(loc = emparams.loc['RB','Mean birth'],scale = emparams.loc['RB','Std birth'])
    rb_conc = init_rb / init_size
    cell = pd.DataFrame(data = [[i,0,init_size,init_rb,rb_conc,'G1',init_size]],
                         columns = ['CellID','Time','Size','RB','RB conc','Phase','Birth size'])
    population[i,0] = cell
    
#%%
    
# Start stimulation
for t in np.arange(Niter)+1:
    for i in range(Ncells):
        
        print t
        
        # Go to previous time point and grab cell to work on
        this_cell_now = population[i,t-1].copy()
        
        # Check if cell divided
        if this_cell_now['Phase'].values == 'None':
            this_cell_next = this_cell_now
        else:
            this_cell_next = simulation.advance_dt(this_cell_now,emparams)
        
        population[i,t] = this_cell_next

        
        
        
        