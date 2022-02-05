#!/usr/bin/python3

############ By Tülin Stapke - 2022 ############


import os
import sys
import shutil
import glob
import pathlib


#This part is to collect all the cluster.pdb files except for the ones in the summary folders, (since in that folder there are only copies of some poses from another folder), move them all into a new folder named clusters and rename them#

start_folder_path=pathlib.Path('/home/user/Desktop/small_test/')

end_folder_path=pathlib.Path('/home/user/Desktop/small_test/clusters/')

paths_to_clusters=[path for path in start_folder_path.rglob('cluster*pdb')]

for target_file_path in paths_to_clusters:
	if target_file_path.parent.name == glob.glob('/home/user/Desktop/small_test/**/*_summary'):
		continue
	else:
		new_file_name=((target_file_path.parent.parent.name)+'_'+(target_file_path.name))

		new_file_path=end_folder_path.joinpath(new_file_name)

		shutil.copy(target_file_path, new_file_path)

				
				
