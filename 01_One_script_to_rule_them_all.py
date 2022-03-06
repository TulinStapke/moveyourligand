#!/usr/bin/python3

import os
import sys
import shutil
import glob
import pathlib
import math
import numpy as np
import pdb_filter_chains_modified as flt
import pdb_atoms as atms
import read_receptor as rr


path_to_clusters = ('/path/clusters/')
path_to_ligands=('/path/ligands/')
path_to_CMS=('/path/CMS/')
path_to_grid=('/path/grid/')
path_to_receptor=('/path/receptor_4dkl/')
path_to_whole_pdbs=('/path/whole_pdbs/')

#This part is to collect all the cluster.pdb files except for the ones in the summary folders, move them all into a new folder named clusters and rename them#

start_folder_path=pathlib.Path('/path/test_4dkl')

end_folder_path=pathlib.Path('/path/clusters/')

paths_to_clusters=[path for path in start_folder_path.rglob('cluster*pdb')]

for target_file_path in paths_to_clusters:
	if target_file_path.parent.name == glob.glob('/home/user/Desktop/small_test/test_4dkl/**/*_summary'):
		continue
	else:
		new_file_name=((target_file_path.parent.parent.name)+'_'+(target_file_path.name))

		new_file_path=end_folder_path.joinpath(new_file_name)

		shutil.copy(target_file_path, new_file_path)

#This part is to extract ligand chain L from each files in the cluster folder, move them into ligands folder and rename them#

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

#This part puts ligands at grid points that are defined around the receptor

grid_min_x, nr_bins_x,  grid_min_y, nr_bins_y, grid_min_z, nr_bins_z, delta = rr.ReadReceptor()
print( 'break the limits:', grid_min_x, nr_bins_x,  grid_min_y, nr_bins_y, grid_min_z, nr_bins_z, delta )

with os.scandir(path_to_CMS) as cm:
    for ligand in cm:
        new_name = path_to_grid + ligand.name.replace('lig_CMS', 'grid')[:-4]
        print( ligand,new_name)
        with open(ligand, 'r') as text:
            lines=text.readlines()
            print( 'ligand', len(lines), 'atoms')
            pos = [] # all positions of ligand
            for line in lines:                
                X=float(line[30:38])
                Y=float(line[38:46])
                Z=float(line[46:54])
                pos.append( np.array( [X,Y,Z] ) )
            print( 'nr pos:', len(pos))
            for i in range(0, nr_bins_x):
                X = i * delta + grid_min_x
                print( i, X)
                for j in range( 0, nr_bins_y):
                    Y = j * delta + grid_min_y
                    for k in range( 0, nr_bins_z):
                        Z = k * delta + grid_min_z
                        file_name = new_name + '_' + str(i)+ '_' + str(j)+ '_' + str(k) + '.pdb'
                        #print( file_name )
                        with open(file_name , 'w') as w:
                            for p,l in zip( pos, lines):
                                text_to_write=l[:30] + ("%8.3f%8.3f%8.3f" % (X+p[0],Y+p[1],Z+p[2]) ) +l[54:]
                                w.write( text_to_write )


#This part connects receptor pdbs under the ligand pdbs, avoiding clashes between them

### Read coordinates from the receptor

#with os.scandir(path_to_receptor) as rp:        # rp is a list of all receptor PDBs

rp = glob.glob(f"{path_to_receptor}*")[0]

#print(f"DEBUG: {rp = }")

with open(rp, 'r') as text:
    lines=text.readlines()
    x_rec = []
    y_rec = []
    z_rec = []                            # Open the PDB file and initialize the lists for coordinates
    for line in lines:

        if line[:4] != "ATOM" and line[:6] != "HETATM": # Skip lines we dont want
            continue
# in my receptor the y and z axes are shifted.
        x_rec.append(float(line[30:38]))
        y_rec.append(float(line[38:46]))
        z_rec.append(float(line[46:54])) 
        receptor_coords = np.empty([3,len(x_rec)], dtype=float)
receptor_coords[0,:] = x_rec
receptor_coords[1,:] = y_rec
receptor_coords[2,:] = z_rec    


y_rec_mean = (max(y_rec)+min(y_rec))/2 #

z_rec_mean = (max(z_rec)+min(z_rec))/2 # define the cylinder center line that will contain the ligands, use z instead of y, because y and z axes are shifted
x_rec_mean = (max(x_rec)+min(x_rec))/2


cutoff_cylinder = int(max(y_rec) - y_rec_mean) #how far should ligands be fm the center line


#z_rec_ceiling = max(z_rec) +10 
y_rec_ceiling = max(y_rec) +15 
y_rec_bottom =  (min(y_rec)+max(y_rec))*2/3   
#print(z_rec_mean, y_rec_mean, cutoff_cylinder, y_rec_ceiling)


### Read the data from each ligand file


g = glob.glob(f"{path_to_grid}*")    # g is a list of ligand files in the corresponding path

ligand_atom_count = 37 #for 4dkl 37, for 5c1m 36
ligand_count = len(g)

ligand_coords = np.empty([ligand_count, 3, ligand_atom_count], dtype=float)  # order of dimensions: ligand file , xyz, atom [:,:,:]

for ligand_index, ligand in enumerate(g):                # Do this for each ligand file

    with open(ligand, 'r') as text:
        lig_lines = text.readlines()
        #overlap = False
        
        for atom_index, lig_line in enumerate(lig_lines):                

            ligand_coords[ligand_index, :, atom_index] = [float(lig_line[30:38]), float(lig_line[38:46]), float(lig_line[46:54])]


### LOGIC BLOCK

# we want the ligands inside the receptor, until the top of the receptor. if they are outside, dont write it in the grids. 
#ligands should be in a cutoff distance from the cylinder shaped receptors center line

#AND

#ligands should be at a max distance from top of the receptor

#AND

# We want to subtract the coords of each atom of the ligand from each atom of the receptor
# If the amount of difference (abs()) is lower than 1.0at ANY of these positions, we DONT want to write it down
# If it is ALWAYS larger

#we want to write both ligand and receptor to a new file.

def evaluate_distance(receptor_coords, ligand_coords, cutoff=1.0): 

    '''Input ONE receptor and ONE ligand, the function calculates pairwise distances, as soon as one atom pair is below the cutoff value, returns False
    Otherwise, if the ligand is sufficiently far away, returns True'''

    # data structure: [3, number_of_atoms]
    mydims = receptor_coords.shape # (3, number_of_atoms)
    receptor_atoms = mydims[1] # take the value from the dimensions of the array

    for receptor_atom in range(receptor_atoms):
        rec_atom_coords = receptor_coords[:, receptor_atom]
        for ligand_atom in range(ligand_atom_count):
            lig_atom_coords = ligand_coords[:, ligand_atom]
            
            dx2 = (rec_atom_coords[0]-lig_atom_coords[0])**2
            dy2 = (rec_atom_coords[1]-lig_atom_coords[1])**2
            dz2 = (rec_atom_coords[2]-lig_atom_coords[2])**2          
            distance = np.sqrt(dx2 + dy2 + dz2)
       
            if distance < cutoff:
                return False

    return True


def evaluate_ceiling(y_rec_ceiling, ligand_coords):
    y_max = max(ligand_coords[1,:])
    if y_max > y_rec_ceiling :
        return False
    else:
        return True
  
def evaluate_cone(x_mean, y_mean, z_mean, ligand_coords, slope=0.39):
    # The mean values define the ABSOLUTE CENTER of the receptor
    # Define a cone originating from the CENTER, which extends upwards in x direction (x value increasing)
    # we can define a slope (in radians) which can be used to vary the cone extent
    #    a slop of pi/4 is again an infinite cylinder, we want to keep the angle around  45° (pi/8 ~ 0.39)

    for ligand_atom in range(ligand_atom_count):
            lig_atom_coords = ligand_coords[:, ligand_atom]          
            dx_c2 = (x_rec_mean - lig_atom_coords[0])**2
            dz_c2 = (z_rec_mean - lig_atom_coords[2])**2 
            distance_xz_plane = np.sqrt(dx_c2 + dz_c2 )
            
            # the cone_boundary is dependent on the height (y-dimension) of the corresponding atom
            print(f"DEBUG: {ligand_atom = }, {ligand_coords.shape = } {lig_atom_coords = } {lig_atom_coords.shape = }")
            print(f"DEBUG {lig_atom_coords[1] = }, {y_mean = }, {type(lig_atom_coords[1]) = }, {type(y_mean) = }")
            dy = lig_atom_coords[1] - y_rec_mean

            # since the cone ONLY extends UPWARDS, the Y coord. of the atom needs to be LARGER than y_mean
            if dy <= 0:
                return False
            cone_boundary = dy * math.sin(slope)
            #print(f"DEBUG: {ligand_atom = }, {ligand_coords.shape = } {lig_atom_coords = } {lig_atom_coords.shape = }")
            if distance_xz_plane > cone_boundary:
                return False

    return True
    
              
### Iteration and Writing of the new files, if the receptor and ligand are far enough apart


with open(rp, "r") as receptor_file:    # Open receptor file and copy manually, eliminating the empty lind and other erroneus statements
    receptor_data = receptor_file.readlines()

for ligand_index in range(ligand_count): # For each ligand, perform check, if thats true, write it down
    print("*"*30, "\n", "CHECKING LIGAND NUMBER:", ligand_index) 

    if evaluate_distance(receptor_coords, ligand_coords[ligand_index,:,:], cutoff=1.0) == False:
        print("[SKIP]: 1st Distance threshold not satisfied.")
        continue
    print("1st Distance satisfied.")

    if evaluate_ceiling(y_rec_ceiling, ligand_coords[ligand_index,:,:]) == False:
        print("[SKIP]: 2nd Ceiling criterion not satisfied.")
        continue
    print("2nd Ceiling satisfied.")    
    
    if evaluate_cone(x_rec_mean, y_rec_mean, z_rec_mean, ligand_coords[ligand_index,:,:] , slope=0.39) == False:
        print(f"[SKIP]:3rd CONE criterion not satisfied.")
        continue   
    print("3rd Cone satisfied. All Criteria are met. Continuing to writing File.")
   
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



