#!/usr/bin/env python
# coding: utf-8

import glob
import os
from shutil import copyfile

path = '.'

os.chdir(path)
os.getcwd()

files = sorted(glob.glob('*.openephys'))
print(files)

for thisfile in files:
    copyfile('master.params', thisfile[:4] + '.params')

