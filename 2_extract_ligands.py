#!/usr/bin/python3

############ By Tülin Stapke - 2022 ############

import pdb_filter_chains_modified as flt
import os
import sys

#This part is to extract ligand chain L from each files in the cluster folder, move them into ligands folder and rename them#

path_to_clusters = ('/home/user/Desktop/small_test/here/')
path_to_ligands=('/home/user/Desktop/small_test/ligands/')
chains = ["L"]
with os.scandir( path_to_clusters ) as itr:
	for file_name in itr:
		lines = flt.filter_chain( path_to_clusters + file_name.name, chains)
		new_name = path_to_ligands + file_name.name.replace('cluster', 'ligand')
		print(new_name)
		with open( new_name, 'w') as w:
			for l in lines:
				w.write(l + '\n')
			
			

