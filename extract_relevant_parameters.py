#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:43:22 2019

@author: xies
"""

import numpy as np
import scipy as sp
import pandas as pd
import pickle as pkl

# @todo
# Paramters: gamma1- gamma2 beta1 beta2, Pg1->[RB], Pdivision (timer)
# Parameters: Birth size distribution (mu_b,sigma_b)

# 1-- size growth parameters
# Gamma1/2 - d Size / dt = gamma * Size

filename = '/Users/xies/Box/HMECs/DFB/individuals_pass_g1s.xlsx.synthesis_linear_fits.csv'
lin_fits = pd.read_csv(filename,index_col=0)

