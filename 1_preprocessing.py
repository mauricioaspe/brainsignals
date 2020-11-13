#!/usr/bin/env python3
# coding: utf-8

import openephys2brainsig as o2b


#IDs = input("IDs = ['ID1597', 'ID1659'] (y/n)? ")
IDs = ['ID1597']#, 'ID1659']

#option = input("Substract = median (y/n)? ")
substract = 'median'


structure = ['mPFC']
print(structure)
option = input("Structure: mPFC (y/n)? ").lower()
if option == 'n':
    structure = input("Select structure: ")
    print(structure)
        

save_data = False
option = input("Save data = False (y/n)? ").lower()
if option == "y":
    save_data = True


o2b.raw_to_npy(IDs=IDs, structures=structure, only=[], filter_order=9, detrend=False, substract=substract, downsampling=False, load_triggers=True, load_accelerometers=False, save_data=save_data)



