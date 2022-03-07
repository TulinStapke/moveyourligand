#!/usr/bin/python3

############ By Tülin Stapke - 2022 ############



import os
import sys
import numpy as np

path_to_receptor=('/path/receptor/')


def ReadReceptor( ):
    with os.scandir(path_to_receptor) as rp:
        for receptor in rp:            
            with open(receptor, 'r') as text:
                lines=text.readlines()
                x_rec = []
                y_rec = []
                z_rec = []
                for line in lines:
                    if line[:4] != "ATOM" and line[:6] != "HETATM": continue
                    x_rec.append(float(line[30:38]))
                    y_rec.append(float(line[38:46]))
                    z_rec.append(float(line[46:54])) 


                grid_min_x_rec=float(min(x_rec))
                grid_max_x_rec=float(max(x_rec))
                # grid_min_y_rec=float(min(y_rec))
                # grid_max_y_rec=float(max(y_rec)) #y changed to z, bcs z and y are shifted in the receptor
                grid_min_z_rec=float(min(z_rec))
                grid_max_z_rec=float(max(z_rec))
                
                # grid_min_z_rec=float((min(z_rec)+max(z_rec))*2/3)
                # grid_max_z_rec=float((max(z_rec))+30)
                grid_min_y_rec=float((min(y_rec)+max(y_rec))*2/3) #bcs z and y are shifted in the receptor
                grid_max_y_rec=float((max(y_rec))+30)
                
                delta=3
                nr_bins_x_rec=int((grid_max_x_rec - grid_min_x_rec)/delta)
                nr_bins_y_rec=int((grid_max_y_rec - grid_min_y_rec)/delta)
                nr_bins_z_rec=int((grid_max_z_rec - grid_min_z_rec)/delta)
                #print(nr_bins_x_rec)
                #print( "min:", min(x_rec),min(y_rec),min(z_rec))
               # print( 'max:', max(x_rec),max(y_rec),max(z_rec))
    return grid_min_x_rec, nr_bins_x_rec,  grid_min_y_rec, nr_bins_y_rec, grid_min_z_rec, nr_bins_z_rec, delta
   # print(nr_bins_z_rec)


