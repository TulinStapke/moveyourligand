#!/usr/bin/python3

############ By Tülin Stapke - 2022 ############

import os
import sys
import shutil
import glob
import pathlib
import pdb_filter_chains_modified as flt
import read_receptor as rr
import pdb_atoms as atms
import numpy as np

#This part is to collect all the cluster.pdb files except for the ones in the summary folders, move them all into a new folder named clusters and rename them#

start_folder_path=pathlib.Path('/home/tulin/Desktop/test/')

end_folder_path=pathlib.Path('/home/tulin/Desktop/test/clusters/')

paths_to_clusters=[path for path in start_folder_path.rglob('cluster*pdb')]

for target_file_path in paths_to_clusters:
	if target_file_path.parent.name == glob.glob('/home/tulin/Desktop/test/**/*_summary'):
		continue
	else:
		new_file_name=((target_file_path.parent.parent.name)+'_'+(target_file_path.name))

		new_file_path=end_folder_path.joinpath(new_file_name)

		shutil.copy(target_file_path, new_file_path)




#This part is to extract ligand chain L from each files in the cluster folder, move them into ligands folder and rename them#

path_to_clusters = ('/home/tulin/Desktop/small_test/here/')
path_to_ligands=('/home/tulin/Desktop/small_test/ligands/')
chains = ["L"]
with os.scandir( path_to_clusters ) as itr:
	for file_name in itr:
		lines = flt.filter_chain( path_to_clusters + file_name.name, chains)
		new_name = path_to_ligands + file_name.name.replace('cluster', 'ligand')
		print(new_name)
		with open( new_name, 'w') as w:
			for l in lines:
				w.write(l + '\n')


							
#This part calculates the center of masses of each atoms in each ligand pdb file, subtract them from the original coordinates,writes the new ligand pdbs with new coordinates into CMS folder#

path_to_CMS=('/home/tulin/Desktop/small_test/CMS/')

with os.scandir( path_to_ligands) as lig:
	for each_file in lig:
		text=open(each_file, "r")
		lines=text.readlines()
		result=[]
		for x in lines:
			result.append( [ float(x[30:38]) ,float(x[39:46]) ,float(x[47:54]) ] )	
		new_coord=[]
		print( "CMS all: ", atms.CMS( each_file) )
		array_1=np.array(result)
		array_2=np.array(atms.CMS( each_file))
		subtracted_array = np.subtract(array_1, array_2)
		subtracted = list(subtracted_array)
		for i in subtracted:
			new_coord=i
			print(i)		



#This part takes the min, max values of the coordinates, step size delta, and number of bins from the read_receptor.py, reads the new coordinates of the ligand that is reduced by its center of mass coordinates. adds the min values of the coordinates to visit, then keeps adding new coordinates increasing by delta. writes a new ligand pdb file for each increment in the coordinates.#
		
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
                                

# This part reads the coordinates from the receptor, and ligands. if the minimum distance between any coordinates of them >=0.3,
#   adds all the lines of the receptor under the lines of the each ligand file and writes it as a new output file. 
# if the ligand and receptor are getting closer than 0.3 at any point skips that ligand file

path_to_grid=('/home/tulin/Desktop/small_test/grid/')
path_to_receptor=('/home/tulin/Desktop/small_test/receptor/')
path_to_whole_pdbs=('/home/tulin/Desktop/small_test/whole_pdbs/')

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
