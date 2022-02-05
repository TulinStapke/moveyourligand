#!/usr/bin/python3


#######################################
###  COPYRIGHT: RENE STARITZBICHLER  ##
###             02.02.2020           ##
#######################################


import sys

def filter_chain( file_name, chains):
	#change side chain	
	lines = []
	with open( file_name ) as f: 
		# gehe durch zeilen	
		for l in f:                
        		# filtern von zeilen, die mit "ATOM" oder "HETATM"  anfangen
     		        if ( l[:4] == "ATOM" or  l[:6] == "HETATM" ) and l[21] in chains:
                                lines.append( l.strip())
	return lines
	
if __name__ == "__main__":
	# pruefe laenge der liste der argumente wenn ungleich print
	if len( sys.argv ) < 3:
		print( "USAGE:", sys.argv[0], "PDB CHAIN1 ..")
		exit(1)


	chains = sys.argv[2:]
	file_name = sys.argv[1]
	for line in filter_chain(file_name, chains):
		print( line, end="")
