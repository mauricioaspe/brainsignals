#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 17:51:18 2020

@author: diogenes
"""
import sys
sys.path.append('/home/diogenes/Projects/brainsignals/modules/')

#import pickle
import mouseInfo


mice_list = mouseInfo.getMice(path='/home/diogenes/Projects/brainsignals/DATA/MICE/')

mouseInfo.getInfo(mice_list, save=True)




