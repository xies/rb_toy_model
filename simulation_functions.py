#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 01:28:07 2020

@author: xies
"""

import numpy as np

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
    next_cell['Time'] += 1
    
    # check g1/s
    if g1s_transition(cell,params):
        next_cell['Phase'].values = 'G2'
    
    # check division
    if division(cell,params):
        next_cell['RB'] = np.nan
        next_cell['Size'] = np.nan
        next_cell['Phase'] = 'None'
    else:
        # If no division, then calculate next step
        drb = rb_synthesis_lut(next_cell,params)
        ds = size_synthesis_lut(next_cell,params)
        
        next_cell['Size'] += ds
        next_cell['RB'] += drb
    
    return next_cell


def rb_synthesis_lut(cell,params):
    from numpy.random import randn
    #@todo: build noise estimator + randn
    if cell['Phase'].values == 'G1':
        m = params.loc['RB','G1 m']
        b = params.loc['RB','G1 b']
        sigma = np.sqrt(params.loc['RB','G1 msr'])
    elif cell['Phase'].values == 'G2':
        m = params.loc['RB','SG2 m']
        b = params.loc['RB','SG2 b']
        sigma = np.sqrt(params.loc['RB','SG2 msr'])
    return np.polyval([m,b],cell['RB'])  + randn(1) * sigma

def size_synthesis_lut(cell,params):
    from numpy.random import randn
    if cell['Phase'].values == 'G1':
        m = params.loc['Size','G1 m']
        b = params.loc['Size','G1 b']
        sigma = np.sqrt(params.loc['Size','G1 msr'])
    elif cell['Phase'].values == 'G2':
        m = params.loc['Size','SG2 m']
        b = params.loc['Size','SG2 b']
        sigma = np.sqrt(params.loc['Size','SG2 msr'])
    return np.polyval([m,b],cell['Size']) + randn(1) * sigma
    
    
def g1s_transition(cell,params):
    if cell['Phase'].values == 'G1':
        # @todo: build logit lookup table baesed on [RB]
        if cell['RB'].values/cell['Size'].values:
            transitioned = True
        else:
            transitioned = False
    elif cell['Phase'].values == 'G2':
        transitioned = True
        
    return transitioned

def division(cell,params):
    # @todo: currently don't implement daughters div == death
    
    if cell['Size'] > cell['Birth size']/2:
        return True
    else:
        return False


def mc_accept(prob_lut,trial_value):
    """
    Accept or reject the current trial_value given the probability density LUT
    (assuming near-analytic resolution ala KDE estimates)
    
    """
    from numpy.random import rand # Uniform [0,1]
    
    p = rand(1)[0]
    
    I = find_nearest_idx(prob_lut,trial_value)
    
    return p < I # acceptance rate