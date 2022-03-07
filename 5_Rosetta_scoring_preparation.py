#!/usr/bin/python3

############ By Tülin Stapke - 2022 ############



import os

path_to_ligands=('/path/ligands/')
# When scoring with the rosetta, ligand chain needs to be named as X and Ligand "ATOM" needs to be replaced with "HETATM"

with os.scandir(path_to_ligands) as lig:
    for each_file in lig:
        with open(each_file, 'r') as file:
            file_data=file.read()    
            file_data =file_data.replace('L', 'X')
            file_data=file_data.replace( 'ATOM', 'HETATM')
        with open(each_file, 'w') as file:
            file.write(file_data)
# When changing from ATOM to HETATM, one more space is generated which caused the error woth rosetta. Fix this error with running the command
#line on the teminal: for f in *.pdb; do sed -i 's/HETATM  /HETATM/g' "$f" ; done          
            
