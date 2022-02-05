#!/usr/bin/python3

############ By Tülin Stapke - 2022 ############

import os
import sys
import pdb_atoms as atms
import numpy as np

#calculates the center of masses of each atoms in each ligand pdb file, subtract them from the original coordinates,writes the new ligand pdbs with new coordinates into CMS folder#

path_to_ligands=('/home/user/Desktop/small_test/ligands/')
path_to_CMS=('/home/user/Desktop/small_test/CMS/')


with os.scandir(path_to_ligands) as lig:
    for each_file in lig:       
        new_name = path_to_CMS + each_file.name.replace('ligand', 'lig_CMS')
        with open(new_name, 'w') as w:
            text=open(each_file, "r")
            lines=text.readlines()
            CMS_coord= np.array( atms.CMS(each_file) )  #getting the coordinates of each atoms center of masses
            for line in lines:
                x =[]
                x.append( float(line[30:38]) - CMS_coord[0])
                x.append( float(line[38:46]) - CMS_coord[1])
                x.append( float(line[46:54]) - CMS_coord[2])
                line = line.strip()
                text_to_write=line[:30] + ("%8.3f%8.3f%8.3f" % (x[0],x[1],x[2])) +line[54:]
                w.write(text_to_write + '\n')
