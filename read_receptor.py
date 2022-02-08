#!/usr/bin/python3

############ By Tülin Stapke - 2022 ############

import os
import sys
import numpy as np

#reads the coordinates of atoms in the receptor, defines delta:step length, gives the min, max values of the x,y coordinates, defines a min and max value for the z coordinate, and number of bins for each coordinates#

path_to_receptor=('/home/user/Desktop/small_test/receptor/')

def ReadReceptor( ):
    with os.scandir(path_to_receptor) as rp:
        for receptor in rp:            
            with open(receptor, 'r') as text:
                lines=text.readlines()
                x = []
                y = []
                z = []
                for line in lines:
                    if line[:4] != "ATOM" and line[:6] != "HETATM": continue
                    x.append(float(line[30:38]))
                    y.append(float(line[38:46]))
                    z.append(float(line[46:54])) 

                grid_min_x=float(min(x))
                grid_max_x=float(max(x))
                grid_min_y=float(min(y))
                grid_max_y=float(max(y))
                grid_min_z=float((min(z)+max(z))*2/3)
                grid_max_z=float((max(z))+10)
                delta=2
                nr_bins_x=int((grid_max_x - grid_min_x)/delta)
                nr_bins_y=int((grid_max_y - grid_min_y)/delta)
                nr_bins_z=int((grid_max_z - grid_min_z)/delta)
                #print(nr_bins_x)
                #print( "min:", min(x),min(y),min(z))
               # print( 'max:', max(x),max(y),max(z))
    return grid_min_x, nr_bins_x,  grid_min_y, nr_bins_y, grid_min_z, nr_bins_z, delta
