#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:16:39 2020

@author: xies
"""

import pandas as pd
import numpy as np
from numpy import random
import seaborn as sb
import matplotlib.pylab as plt
import statsmodels.api as sm
import simulation
import copy

# 1) Initial conditions / parameters
Ncells = 10
max_iter = 200 # in frames
dt = 1/6. # in hours 

print(max_iter, ' iterations = ', max_iter * dt / 24, 'days')
print('Expected final population size: ', Ncells * 2**( max_iter * dt / 15))

# emparams = pd.read_pickle('/Users/xies/Box/HMECs/RB toy model/emp_parameters.pkl')
# trans_params = pd.read_pickle('/Users/xies/Box/HMECs/RB toy model/trans_params.pkl')

# emparams = emparams.join(trans_params)
# emparams.at['RB conc','G1S trans m'] = trans_params.loc['RB conc','G1S trans m']
# emparams.at['RB conc','G1S trans b'] = trans_params.loc['RB conc','G1S trans b']
# emparams.at['Time','Mean G2 duration'] = trans_params.loc['Time','Mean G2 duration']
# emparams.at['Time','Std G2 duration'] = trans_params.loc['Time','Std G2 duration']

'''
Parameter sets should be given in hours and mean-divided quantities of RB or Size
'''
emparams = pd.read_csv('empirical_params.csv',index_col=0)

sim_clock = {}
sim_clock['Max frame'] = max_iter
sim_clock['Max time'] = max_iter * dt
sim_clock['dt'] = dt

#%% Initialize
# Initialize each cell as a DataFrame at G1/S transition so we can specify Size and RB independently

sim_clock['Current time'] = 0
sim_clock['Current frame'] = 0

next_cellID = 0
population = {}
for i in range(Ncells):
    # Initialize cells de novo
    cell = simulation.Cell(i, sim_clock, emparams)
    population[i] = cell
    next_cellID += 1
    
# Deepcopy original
initial_pop = copy.deepcopy(population)

#%% 2. Simulation steps

next_cellID = len(initial_pop)
population = copy.deepcopy(initial_pop)
sim_clock['Current time'] = 0
sim_clock['Current frame'] = 0

for t in np.arange(sim_clock['Max frame'] - 1):
    
    # Advance time step by one
    sim_clock['Current frame'] += 1
    sim_clock['Current time'] += sim_clock['dt']
    
    newly_borns = {}
    for this_cell in population.values():
        
        # Skip cell if divided already
        if this_cell.divided:
            continue
        else:
            
            this_cell.advance_dt(sim_clock,emparams)
            
            if this_cell.divided:
                # Newly divided cell: make daughter cells
                print(f'CellID #{this_cell.cellID} has divided at frame {t}')
                # Randomly draw an asymmettry
                a = np.abs( random.randn() * 0.05 ) # 5 percent
                daughters = this_cell.divide(next_cellID, sim_clock, asymmetry=a)
                next_cellID += 2
                # Put daughters into the population
                newly_borns[daughters[0].cellID] = daughters[0]
                newly_borns[daughters[1].cellID] = daughters[1]
    
    population.update(newly_borns)

             
#%% Retrieve each datafield into dataframe
        
btime = np.vstack( [ cell.ts['Birth time'].astype(np.float) for cell in population.values() ])
size = np.vstack( [ cell.ts['Size'].astype(np.float) for cell in population.values() ])
rb = np.vstack( [ cell.ts['RB'].astype(np.float) for cell in population.values() ])
rb_conc = np.vstack( [ cell.ts['RB conc'].astype(np.float) for cell in population.values() ])


birth_size = np.array([ cell.birth_size for cell in population.values() ])
birth_rb = np.array([ cell.birth_rb for cell in population.values() ])
g1s_rb_conc = np.array([ cell.g1s_rb_conc for cell in population.values() ])
g1s_size = np.array([ cell.g1s_size for cell in population.values() ])
div_size = np.array([ cell.div_size for cell in population.values() ])
g1s_duration = np.array([ cell.g1s_time - cell.birth_time for cell in population.values() ])
div_duration = np.array([ cell.div_time - cell.birth_time for cell in population.values() ])


birth_time = np.array([ cell.birth_time for cell in population.values() ])


#%% Plotting

plt.figure()
plt.plot(size.T); plt.xlabel('Birth time (h)'); plt.ylabel('Size')
plt.figure()
plt.plot(btime.T,rb_conc.T); plt.xlabel('Birth time (h)'); plt.ylabel('RB conc')




