#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 13:51:40 2018

@author: marzab
"""
# process_dir: Z.txt
import numpy as np
from scipy import stats
def EVD(process_dir):
    thr=0.01
    Z=np.loadtxt(process_dir)

    m=Z.shape
    l=np.int(m[0]*thr)
    
    T=np.sort(np.abs(Z),axis=0)[m[0]-l:m[0],:]
    E=stats.trim_mean(T, 0.1)
    return E