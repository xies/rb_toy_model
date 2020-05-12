#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 01:28:07 2020

@author: xies
"""

import numpy as np
from numpy import random

#def initialize_cell_empirical_kde(kde):
#    '''
#    Initialize values (size or RB) for the simulated cell based on an empirical KDE
#    
#    '''
#    from numpy.random import randn
#    
    

def advance_dt(cell,params):
    '''
    Advance a single cell forward one 'time point' (set by DFB dataset dt)
    
    '''
    # @todo: 1 grow size 2) synthesize rb
    next_cell = cell.copy()
    next_cell['Time'] += 1/6.
    
    # First calculate growths
    drb = rb_synthesis_lut(cell,params)
    ds = size_synthesis_lut(next_cell,params)
    
    next_cell['Size'] += ds
    next_cell['RB'] += drb
    next_cell['RB conc'] = next_cell['RB'] / next_cell['Size']
    
    if cell['Phase'] == 'G2':
        
        # check division
        if division(cell,params):
            print('Divided at time ',str(next_cell.Time))
            next_cell['RB'] = np.nan
            next_cell['Size'] = np.nan
            next_cell['RB conc'] = np.nan
            next_cell['Phase'] = 'None'
            next_cell['Div time'] = cell['Time']
            next_cell['Div size'] = cell['Size']
            
            return next_cell
            
    elif cell['Phase'] == 'G1':
        # check g1/s
        if g1s_transition(cell,params):
            print('G1S transitioned with RB_conc ',str(next_cell['RB conc']))
            next_cell['Phase'] = 'G2'
            next_cell['G1S time'] = cell['Time']
            next_cell['G1S size'] = cell['Size']
            next_cell['G1S RB conc'] = cell['RB conc']
    
    return next_cell


def rb_synthesis_lut(cell,params):
    from numpy.random import randn
    #@todo: build noise estimator + randn
    if cell['Phase'] == 'G1':
        m = params.loc['RB','G1 m']
        b = params.loc['RB','G1 b']
        sigma = np.sqrt(params.loc['RB','G1 msr'])
    elif cell['Phase'] == 'G2':
        m = params.loc['RB','SG2 m']
        b = params.loc['RB','SG2 b']
        sigma = np.sqrt(params.loc['RB','SG2 msr'])
    return np.polyval([m,b],cell['RB'])  + randn(1) * sigma

def size_synthesis_lut(cell,params):
    from numpy.random import randn
    if cell['Phase'] == 'G1':
        m = params.loc['Size','G1 m']
        b = params.loc['Size','G1 b']
        sigma = np.sqrt(params.loc['Size','G1 msr'])
    elif cell['Phase'] == 'G2':
        m = params.loc['Size','SG2 m']
        b = params.loc['Size','SG2 b']
        sigma = np.sqrt(params.loc['Size','SG2 msr'])
    return np.polyval([m,b],cell['Size']) + randn(1) * sigma
    
    
def g1s_transition(cell,params):
    
    transitioned = False
    
    if cell['Phase'] == 'G1':
        # @todo: build logit lookup table baesed on [RB]
        rb_conc = cell['RB'] / cell['Size']
        m = params.loc['RB conc']['G1S trans m']
        b = params.loc['RB conc']['G1S trans b']
        p = np.polyval([m,b], rb_conc)
        
        q = random.rand() # Metropolis acceptance
        if q < p:
            transitioned = True
        else:
            transitioned = False
    elif cell['Phase'] == 'G2':
        transitioned = True
        
        
    return transitioned

def division(cell,params):
    # @todo: currently don't implement daughters div == death
    # Use timer for G2 duration -- need to extractt timing params
    
    mean_g2 = params.loc['Time','Mean G2 duration']
    std_g2 = params.loc['Time','Std G2 duration']
    
    q = random.normal(loc=mean_g2,scale=std_g2)
    g2_duration = cell.Time - cell['G1S time']
    
    if g2_duration > q:
        return True
    else:
        return False



