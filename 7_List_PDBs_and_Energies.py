#!/usr/bin/python3

import os 
import sys

#open and read each file in the path, if "TOTAL_ENERGY:" is written in a line extract 
#the total energy value, which is seperated by space, and in a txt file make a list
#of the file names and the energy value next to it. sort the list by the energy value

filename= os.listdir("/path/to/files/")
# print(file_list)
filename_and_energy=[]
text_file =  open('list_pdbs_energies.txt' , 'w')
for f in filename:
    with open("{}/{}".format("/path/to/files/", f), 'r') as text:
        lines=text.readlines()
        for l in lines:
            if 'TOTAL_ENERGY:' in l:
                cols =l.split()
                energy = float(cols[1])
                filename_and_energy.append((f, energy))
#print(filename_and_energy )
filename_and_energy.sort(key=lambda x: x[1])
#print(filename_and_energy )
for e in filename_and_energy: 
#    text_file.write("{}\n".format(e))
    text_file.write("{} {}\n".format(e[0], e[1]))
#text_file.write(str(filename_and_energy))
text_file.close() 



