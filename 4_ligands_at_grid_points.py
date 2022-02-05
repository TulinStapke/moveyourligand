#!/usr/bin/python3

############ By Tülin Stapke - 2022 ############

import os
import sys
import read_receptor as rr
import numpy as np

#takes the min, max values of the coordinates, step size delta, and number of bins from the read_receptor.py, reads the new coordinates of the ligand that is reduced by its center of mass coordinates. adds the min values of the coordinates to visit, then keeps adding new coordinates increasing by delta. writes a new ligand pdb file for each increment in the coordinates.#

path_to_CMS=('/home/user/Desktop/small_test/CMS/')
path_to_grid=('/home/user/Desktop/small_test/grid/')

grid_min_x, nr_bins_x,  grid_min_y, nr_bins_y, grid_min_z, nr_bins_z, delta = rr.ReadReceptor()

with os.scandir(path_to_CMS) as cm:
    for ligand in cm:
        new_name = path_to_grid + ligand.name.replace('lig_CMS', 'grid')[:-4]
        with open(ligand, 'r') as text:
            lines=text.readlines()
            pos = [] # all positions of ligand
            for line in lines:                
                X=float(line[30:38])
                Y=float(line[38:46])
                Z=float(line[46:54])
                pos.append( np.array( [X,Y,Z] ) )
            for i in range(0, nr_bins_x):
                X = i * delta + grid_min_x
                for j in range( 0, nr_bins_y):
                    Y = j * delta + grid_min_y
                    for k in range( 0, nr_bins_z):
                        Z = k * delta + grid_min_z
                        file_name = new_name + '_' + str(i)+ '_' + str(j)+ '_' + str(k) + '.pdb'
                        with open(file_name , 'w') as w:
                            for p,l in zip( pos, lines):
                                text_to_write=l[:30] + ("%8.3f%8.3f%8.3f" % (X+p[0],Y+p[1],Z+p[2]) ) +l[54:]
                                w.write( text_to_write )
