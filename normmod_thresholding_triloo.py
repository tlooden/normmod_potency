#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:05:03 2019

@author: triloo
"""

import scipy.stats
import nispat
import numpy as np
import statsmodels.stats.multitest

"""
Function that takes as input Z.txt files (features x observations) from a 
normative modeling analysis.
Returns two binarized matrices. One containing the significant positive 
deviations, and one containing the significant negative deviations - in 
numpy format.


Parameters -----------------

inputfile:  Enter the location of your Z.txt file

threshold:  you can use a predetermined numerical Z-threshold, otherwise  
            it will default to fdr-correction.
            
persub:     Enter 'True' for this parameter for the function to return one 
            value per participant wrt their # deviant edges.
            
variant"    Enter 'pos', to use onnly the positive deviations. 'neg' for the 
            negative deviations, and 'comb' to add them together.

"""



def threshmat(inputfile, threshold='fdr', variant='comb'):
    
    inputfile=np.array(nispat.fileio.load(inputfile))
    
    if threshold=='fdr':
        #preload matrices for new entries
        posthreshmat=np.zeros(np.shape(inputfile))
        negthreshmat=np.zeros(np.shape(inputfile))
        
        #FDR correction is done on a per feature basis
        for index,row in enumerate(inputfile):
            
            #Transform to p-values for positive and negative deviations
            pospvalrow= scipy.stats.norm.sf(row)
            negpvalrow= scipy.stats.norm.sf(-row)
            
            #FDR-correcting the p-values and returning binarized matrix
            posthreshmat[index] = statsmodels.stats.multitest.fdrcorrection(pospvalrow)[0]
            negthreshmat[index] = statsmodels.stats.multitest.fdrcorrection(negpvalrow)[0]
            
    else:
        #Check simply if Z-value exceeds predetermined threshold
        posthreshmat=np.where(inputfile>threshold, 1, 0)
        negthreshmat=np.where(inputfile<(threshold*-1), 1, 0)
    
    if variant=='pos':
        return posthreshmat
    elif variant=='neg':
        return negthreshmat
    elif variant=='comb':
        # Combine both matrices and sum over edges to get one value per subject
        combthreshmat=posthreshmat+negthreshmat
        return combthreshmat
        