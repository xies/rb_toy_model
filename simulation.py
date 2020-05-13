#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 01:28:07 2020

@author: xies
"""

import numpy as np
from numpy import random
import pandas as pd

class Cell():
    
    def __init__(self, cellID, sim_clock, params=None, mother=None, inheritance=0.5):
        '''
        Cell for simulating cell growth/division coupling

        Parameters
        ----------
        cellID : int
            Unique cell ID
        sim_clock : dict
            Dictionary containing the time information of simulation
        params : pd.DataFrame
            empirical parameters extracted from DFB dataset
        mother : Cell, optional
            If inheriting properties from a mother. The default is None.
        inheritance : TYPE, optional
            fraction of inheritance from mother. The default is 0.5.

        Returns
        -------
        Cell object.

        '''
        
        self.cellID = cellID
        # ts -- time series
        self.ts = pd.DataFrame( columns = ['Time','Birth time','Size','RB','RB conc','Phase'],
                            index=pd.RangeIndex( sim_clock['Max frame'] ), dtype=np.float)
        # scalar results of siulation
        self.birth_size = np.nan
        self.div_size = np.nan
        self.g1s_size = np.nan
        self.birth_time = np.nan
        self.div_time = np.nan
        self.g1s_time = np.nan
        self.birth_frame = np.nan
        self.div_frame = np.nan
        self.g1s_frame = np.nan
        self.birth_rb = np.nan
        self.g1s_rb = np.nan
        self.div_rb = np.nan
        self.birth_rb_conc = np.nan
        self.g1s_rb_conc = np.nan
        self.div_rb_conc = np.nan
        # Flag for if cell is has already divided
        self.divided = False
        
        # Check if a mother cell was passed in
        if mother == None:
            # If Mother is not provided, must provide empirical parameters to initialize de novo
            # Will initialize a cell at the G1/S transition
            assert(params is not None)
            # Create cell de novo using empirical paramters
            g1s_size = random.lognormal(mean = params.loc['Size','Mean G1S'],sigma = params.loc['Size','Std G1S'])
            g1s_rb = random.lognormal(mean = params.loc['RB','Mean G1S'],sigma = params.loc['RB','Std G1S'])
            rb_conc = g1s_rb / g1s_size
            
            init_cell = pd.Series({'Time':sim_clock['Current time'],'Birth time':0,
                                   'Size':g1s_size,'RB':g1s_rb,'RB conc':rb_conc,'Phase':'G2'})
            self.parentID = np.nan
            self.g1s_time = g1s_size
            self.g1s_frame = sim_clock['Current frame']
            self.g1s_time = sim_clock['Current time']
            self.g1s_rb = g1s_rb
            self.g1s_rb_conc = rb_conc
            
        else: # Create daughter cell via symmetric dividion (ignores params)
            #NB: the "halving" is taken care of in @divide function
            init_size = mother.div_size * inheritance
            init_rb = mother.div_rb * inheritance
            rb_conc = init_rb / init_size
            
            init_cell = pd.Series({'Time':sim_clock['Current time'],'Birth time':0,
                                   'Size':init_size,'RB':init_rb,'RB conc':rb_conc,'Phase':'G1'})
            self.parentID = mother.cellID
            self.birth_size = init_size
            self.birth_frame = sim_clock['Current frame']
            self.birth_time = sim_clock['Current time']
            self.birth_rb = init_rb
            self.birth_rb_conc = rb_conc
        
        self.ts.at[ sim_clock['Current frame'] ] = init_cell
        
    
    def divide( self, cellID_beginning, sim_clock, asymmetry = 0.0):
        '''
        Divides the current mother cell into two daughters.

        Parameters
        ----------
        cellID_beginning : int
            the cellID of first daughter, +1 will be second daughter
        sim_clock : dict
            Dictionary containing the time information of simulation
        asymmetry : float, optional
            The difference between larger daughter + smaller daughter normalized by mother.
            sym = (D_L - D_s) / M
            The default is 0.0.

        Returns
        -------
        daughter_a : Cell
            Larger daughter
        daughter_b : Cell
            Smaller daughter

        '''
        current_frame = sim_clock['Current frame']
        
        # Onlly allow G2 divisions
        assert(self.ts.iloc[current_frame]['Phase'] == 'None')
        
        # Calculate respective inheritance fractions 
        assert(asymmetry < 1.0)
        inh_a = (asymmetry + 1.) / 2
        inh_b = 1 - (asymmetry + 1.) / 2
        
        daughter_a = Cell(cellID_beginning, sim_clock, mother=self, inheritance=inh_a)
        daughter_b = Cell(cellID_beginning+1, sim_clock, mother=self, inheritance=inh_b)
        
        return (daughter_a, daughter_b)

# --- SIMULATION METHODS ----

    def advance_dt(self,clock,params):
        '''
        
        Parameters
        ----------
        clock : dict
            Dictionary containing current simulation time information
        params : pd.Dataframe
            Empirical parameters extracted from DFB

        '''
        
        assert(self.divided == False)
        
        prev_frame = clock['Current frame'] - 1
        prev_values = self.ts.iloc[prev_frame] # pd.Series
        
        # Initialize current cell's values as pd.Series from prev. time point
        current_values = prev_values.copy()
        current_values['Time'] += clock['dt']
        current_values['Birth time'] += clock['dt']
        
        # First calculate growths based on previous values
        drb = self.rb_synthesis_lut(prev_frame,params)
        ds = self.size_synthesis_lut(prev_frame,params)
        # Add to current values
        current_values['Size'] += ds
        current_values['RB'] += drb
        current_values['RB conc'] = current_values['RB'] / current_values['Size']
        
        # Now check for transitions
        if current_values['Phase'] == 'G2':
            
            # check division
            if self.decide2divide(prev_frame,params):
                current_values['RB'] = np.nan
                current_values['Size'] = np.nan
                current_values['RB conc'] = np.nan
                current_values['Phase'] = 'None'
                self.div_time = prev_values['Time']
                self.div_size = prev_values['Size']
                self.div_frame = prev_frame
                self.div_rb = prev_values['RB']
                self.div_rb = prev_values['RB conc']
                self.divided = True
                
        elif prev_values['Phase'] == 'G1':
            # check g1/s
            if self.g1s_transition(prev_frame,params):
                current_values['Phase'] = 'G2'
                self.g1s_time = prev_values['Time']
                self.g1s_frame = prev_frame
                self.g1s_size = prev_values['Size']
                self.g1s_rb_conc = prev_values['RB conc']
                self.g1s_rb = prev_values['RB']
              
        # Put the current values into the timeseries
        self.ts.at[prev_frame + 1] = current_values

    
    def rb_synthesis_lut(self, frame, params):
        #Look up RB synthesis amount based on a given frame's RB amount
        cell = self.ts.iloc[frame]
        if cell['Phase'] == 'G1':
            m = params.loc['RB','G1 m']
            b = params.loc['RB','G1 b']
            sigma = np.sqrt(params.loc['RB','G1 msr'])
        elif cell['Phase'] == 'G2':
            m = params.loc['RB','SG2 m']
            b = params.loc['RB','SG2 b']
            sigma = np.sqrt(params.loc['RB','SG2 msr'])
            
        return np.polyval([m,b],cell['RB']) + random.randn() * sigma
    
    
    def size_synthesis_lut(self, frame, params):
        cell = self.ts.iloc[frame]
        if cell['Phase'] == 'G1':
            m = params.loc['Size','G1 m']
            b = params.loc['Size','G1 b']
            sigma = np.sqrt(params.loc['Size','G1 msr'])
        elif cell['Phase'] == 'G2':
            m = params.loc['Size','SG2 m']
            b = params.loc['Size','SG2 b']
            sigma = np.sqrt(params.loc['Size','SG2 msr'])
        
        ds = np.polyval([m,b],cell['Size'])
        noise = random.randn() * sigma
        return ds + noise
        
        
    def g1s_transition(self,frame,params):
        
        cell = self.ts.iloc[frame]
        transitioned = False
        
        assert(cell['Phase'] == 'G1')
        
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
            
        return transitioned
    
    def decide2divide(self,frame,params):
        # Use timer for G2 duration
        
        cell = self.ts.iloc[frame]
        mean_g2 = params.loc['Time','Mean G2 duration']
        std_g2 = params.loc['Time','Std G2 duration']
        
        q = random.normal(loc=mean_g2,scale=std_g2)
        g2_duration = cell.Time - self.g1s_time
        
        if g2_duration > q:
            return True
        else:
            return False

# --- DATA METHODS ----
    # def decatenate_fields()
    

    def __repr__(self):
        string = 'Cell ID = ' + str(self.cellID) + '\n'
        string += 'Born at frame ' + str(self.birth_frame) + '\n'
        string += 'Divided : ' + str(self.divided) + '\n\n\n'
        return string
