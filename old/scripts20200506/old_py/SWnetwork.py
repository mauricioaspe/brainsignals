# -*- coding: utf-8 -*-
"""
Created on Sat Jul 28 03:02:44 2018

@author: CNSLABA
"""
import numpy as np
# import matplotlib.pyplot as plt
#import bct

def SW_network(KK, pij, nnodes):

    CM = np.zeros((nnodes,nnodes)) 
    for n in range(nnodes):
        for i in range(1,KK+1):
            if np.random.uniform() < pij:
                ii = np.random.randint(n+1,nnodes+n)
                CM[n,ii%nnodes]=1
            else:
                CM[n,n-i]=1
                
    CM+=CM.T
    CM=np.minimum(CM,np.ones_like(CM))
    return(CM)
                