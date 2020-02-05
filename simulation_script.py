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


# 1) Initial conditions / parameters
Ncells = 10
Niter = 50

emparams = pd.read_pickle('/Users/xies/Box/HMECs/RB toy model/emp_parameters.pkl')
trans_params = pd.read_pickle('/Users/xies/Box/HMECs/RB toy model/trans_params.pkl')
emparams = emparams.append(trans_params)


# Initialize cells as dictionaries
population = np.empty((Ncells,Niter+1),dtype=dict)
for i in range(Ncells):

    # @todo: fix zero size cells
    init_size = random.normal(loc = emparams.loc['Size','Mean birth'],scale = emparams.loc['Size','Std birth'])
    init_rb = random.normal(loc = emparams.loc['RB','Mean birth'],scale = emparams.loc['RB','Std birth'])
    cell = pd.DataFrame(data = [[i,0,init_size,init_rb,'G1',init_size]],
                         columns = ['CellID','Time','Size','RB','Phase','Birth size'])
    population[i,0] = cell
    
# Start stimulation
for t in np.arange(Niter)+1:
    for i in range(Ncells):

        this_cell_now = population[i,t-1]
        this_cell_next = advance_dt(cell,emparams)
        
        population[i,t] = this_cell_next



