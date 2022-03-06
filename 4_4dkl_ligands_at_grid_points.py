#!/usr/bin/python3

import os
import sys
import read_receptor as rr
import numpy as np

path_to_CMS=('/path/CMS/')
path_to_grid=('/path/grid/')



grid_min_x_rec, nr_bins_x_rec,  grid_min_y_rec, nr_bins_y_rec, grid_min_z_rec, nr_bins_z_rec, delta = rr.ReadReceptor()
#print( 'break the limits:', grid_min_x_rec, nr_bins_x_rec,  grid_min_y_rec, nr_bins_y_rec, grid_min_z_rec, nr_bins_z_rec, delta )

with os.scandir(path_to_CMS) as cm:
    for ligand in cm:
        new_name = path_to_grid + ligand.name.replace('lig_CMS', 'grid')[:-4]
#        print( ligand,new_name)
        with open(ligand, 'r') as text:
            lines=text.readlines()
 #           print( 'ligand', len(lines), 'atoms')
            pos = [] # all positions of ligand
            for line in lines:                
                X_lig=float(line[30:38])
                Y_lig=float(line[38:46])
                Z_lig=float(line[46:54])
                pos.append( np.array( [X_lig, Y_lig, Z_lig] ) )
#            print( 'nr pos:', len(pos))
            for i in range(0, nr_bins_x_rec):
                X_lig = i * delta + grid_min_x_rec
 #               print( i, X)
                for j in range( 0, nr_bins_y_rec):
                    Y_lig = j * delta + grid_min_y_rec
                    for k in range( 0, nr_bins_z_rec):
                        Z_lig = k * delta + grid_min_z_rec
                        file_name = new_name + '_' + str(i)+ '_' + str(j)+ '_' + str(k) + '.pdb'
 #                       print( file_name )
                        with open(file_name , 'w') as w:
                            for p,l in zip( pos, lines):
                                text_to_write=l[:30] + ("%8.3f%8.3f%8.3f" % (X_lig+p[0],Y_lig+p[1],Z_lig+p[2]) ) +l[54:]
                                w.write( text_to_write )
