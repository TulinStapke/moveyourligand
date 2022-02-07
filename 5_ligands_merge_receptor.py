#!/usr/bin/python3

############ By Tülin Stapke - 2022 ############

import os
import sys
import numpy as np
import glob

#reads the coordinates from the receptor, and ligands. if the minimum distance between any coordinates of them >=0.3,
#   adds all the lines of the receptor under the lines of the each ligand file and writes it as a new output file. 
# if the ligand and receptor are getting closer than 0.3 at any point skips that ligand file

path_to_grid=('/home/user/Desktop/small_test/grid/')
path_to_receptor=('/home/user/Desktop/small_test/receptor/')
path_to_whole_pdbs=('/home/user/Desktop/small_test/whole_pdbs/')

### Read coordinates from the receptor

#with os.scandir(path_to_receptor) as rp:        # rp is a list of all receptor PDBs

rp = glob.glob(f"{path_to_receptor}*")[0]

print(f"DEBUG: {rp = }")

with open(rp, 'r') as text:
    lines=text.readlines()
    r_x = []
    r_y = []
    r_z = []                            # Open the PDB file and initialize the lists for coordinates
    for line in lines:

        if line[:4] != "ATOM" and line[:6] != "HETATM": # Skip lines we dont want
            continue

        r_x.append(float(line[30:38]))
        r_y.append(float(line[38:46]))
        r_z.append(float(line[46:54])) 
      
receptor_coords = np.empty([3,len(r_x)], dtype=float)
receptor_coords[0,:] = r_x
receptor_coords[1,:] = r_y
receptor_coords[2,:] = r_z              # receptor-coords is a two-dimensional array containing x,y,z for each atom.
                                        # order of dimensions: xyz, atom [:,:]

### Read the data from each ligand file


g = glob.glob(f"{path_to_grid}*")    # g is a list of ligand files in the corresponding path

ligand_atom_count = 36
ligand_count = len(g)

ligand_coords = np.empty([ligand_count, 3, ligand_atom_count], dtype=float)  # order of dimensions: ligand file , xyz, atom [:,:,:]

for ligand_index, ligand in enumerate(g):                # Do this for each ligand file
    print(ligand)      
    # new_name = path_to_whole_pdbs + ligand.name.replace('grid', 'score') # This we need when we write stuff
    # print(new_name)
    # with open(new_name, 'a') as w, open(ligand, 'r') as text:

    with open(ligand, 'r') as text:
        lig_lines = text.readlines()
        #overlap = False
        
        for atom_index, lig_line in enumerate(lig_lines):                

            ligand_coords[ligand_index, :, atom_index] = [float(lig_line[30:38]), float(lig_line[38:46]), float(lig_line[46:54])]
            #   The xyz of the atom are read into the array at the corresponding ligand index and atom index

### LOGIC BLOCK

# We want to subtract the coords of each atom of the ligand from each atom of the receptor
# If the amount of difference (abs()) is lower than 0.3 at ANY of these positions, we DONT want to write it down
# If it is ALWAYS larger, we want to write both ligand and receptor to a new file.

def evaluate_distance(receptor_coords, ligand_coords, cutoff=0.3): 

    '''Input ONE receptor and ONE ligand, the function calculates pairwise distances, as soon as one atom pair is below the cutoff value, returns False
    Otherwise, if the ligand is sufficiently far away, returns True'''

    # data structure: [3, number_of_atoms]
    mydims = receptor_coords.shape # (3, number_of_atoms)
    receptor_atoms = mydims[1] # take the value from the dimensions of the array
    ligand_atoms = 36

    for receptor_atom in range(receptor_atoms):
        rec_atom_coords = receptor_coords[:,receptor_atom]

        for ligand_atom in range(ligand_atoms):
            lig_atom_coords = ligand_coords[:,ligand_atom]

            dx2 = (rec_atom_coords[0]-lig_atom_coords[0])**2
            dy2 = (rec_atom_coords[1]-lig_atom_coords[1])**2
            dz2 = (rec_atom_coords[2]-lig_atom_coords[2])**2

            distance = np.sqrt(dx2 + dy2 + dz2)

            if distance < cutoff:
                return False

    return True

### Iteration and Writing of the new files, if the receptor and ligand are far enough apart

# For each ligand, check if it is far enough apart from the receptor,
# If yes , combine both coords into one file.
# Write ligand first
# Put a TER statement
# Write receptor
# Put a TER
# Put END

with open(rp, "r") as receptor_file:    # Open receptor file and copy manually, eliminating the empty lind and other erroneus statements
    receptor_data = receptor_file.readlines()

for ligand_index in range(ligand_count): # For each ligand, perform check, if thats true, write it down

    check_flag = evaluate_distance(receptor_coords, ligand_coords[ligand_index,:,:], cutoff=0.35)

    if check_flag == True:
        ligand_name = g[ligand_index] # Get the filename from the corresponding position in the filelist
        ligand_filename = ligand_name.split("/")[-1] # Extract only the filename from the whoel path
        new_name = path_to_whole_pdbs + ligand_filename.replace('grid', 'score')

        newfile = open(new_name,'w')

        with open(ligand_name,'r') as ligand_file:  # Open corresponding ligand file and copy everything
            data = ligand_file.read()
            newfile.write(data)

        newfile.write("TER\n")

        for line in receptor_data:                  # Use the previously read receptor data to copy contents
            if line.startswith("ATOM") or line.startswith("HETATM"):
                newfile.write(line)

        newfile.write("TER\nEND\n")
        newfile.close()
