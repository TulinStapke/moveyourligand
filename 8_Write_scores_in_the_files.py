#!/usr/bin/python3

#open and read the list of the pdb files and the scores from the text file, names and 
#energy values are separated with space. in the same folder, the pdb
#files exist without energy values in it. append the pdb files with the corresponding energy
#value from the list, write it as "score: energy value"

import os
import pandas as pd

path = "/home/tulin/Desktop/test/"
dict_energy = {}
for file in os.listdir(path):
    with open(path + "file_names_energies.txt","r") as l:
        data = pd.read_csv(path + "file_names_energies.txt", sep=" ", header=None)
        data.columns = ["file", "energy"]
        for i in range(len(data.file)):
            dict_test = {data.file[i]:data.energy[i]}
            dict_energy.update(dict_test)
    if file != "file_names_energies.txt":
        with open(path + file, "a") as file_object:
            file_object.write('\n' + "haddock_score:" + str(dict_energy.get(file)))

