import numpy as np
import Bio.PDB


import vector_functions as vf


# to make this safe from HETATM section, use
#from Bio import Struct
#s = Struct.read('protein_A.pdb')
# p = s.as_protein()  # strip off hetatms




def RMSD( first_file, second_file, chains = []):
    pdb_parser = Bio.PDB.PDBParser( QUIET=True )

    first_structure = pdb_parser.get_structure( "first", first_file )[0]
    second_structure = pdb_parser.get_structure( "second", second_file )[0]

    first_atoms = []
    second_atoms = []

    for chain in first_structure:
        #print( chain.id )
        if len(chains) > 0 and chain.id not in chains:
            continue
        for residue in chain:
            for atom in residue:
                first_atoms.append( atom )

    for chain in second_structure:
        if len(chains) > 0 and chain.id not in chains:
            continue
        for residue in chain:
            for atom in residue:
                second_atoms.append( atom )

    if len(second_atoms) != len(first_atoms):
        print( "WARNING: number of atoms do not match!", len(first_atoms),len(second_atoms) )

    if len( first_atoms) == 0:
        #sys.stderr.write( "WARNING: empty mol " + first_file + ' ' + second_file + '\n')
        print( "WARNING: no atoms found!")
        return 0.0
    
    aligner = Bio.PDB.Superimposer()
    aligner.set_atoms( second_atoms, first_atoms )
    return aligner.rms

def CMS( pdb_file, chains = [], atoms = [ ] ):
    chain_coo = ReadChainCoordinates( pdb_file , atoms )

    pos = [0,0,0]
    length = 0
    for chain,coo in chain_coo.items():
        if len(chains) > 0 and chain not in chains:
            continue
        for x in coo:
            pos += x
        length += len(coo)
    return np.array(pos) / float(length)

def MinMax( pdb_file, chains = [], atoms = [] ):
    chain_coo = ReadChainCoordinates( pdb_file , atoms )

    maxi = [-999999,-999999,-999999]
    mini = [999999,999999,999999]

    for chain,coo in chain_coo.items():
        if len(chains) > 0 and chain not in chains:
            continue
        for x in coo:
            for i in range(0,3):
                mini[i] = min( mini[i], x[i] )
                maxi[i] = max( maxi[i], x[i] )

    return mini, maxi


def ReadSequencesFromCoordinates( filename, chains = [] ):
    import sequence_functions as sf
    seq = ""
    prev = -99999
    prev_chain = "y908"
    sequences = {}
    #print len(chains), chains
    with open( filename) as f:
        for l in f:
            if l[0:4] != "ATOM": continue
            atom = PDBAtom( l)
            #print atom.atom_name, atom.chain
            if atom.atom_name != "CA" or (len(chains) > 0 and atom.chain not in chains):
                continue
            #print l,
            if atom.chain != prev_chain:
                if len(seq) > 0:
                    sequences[prev_chain] = seq
                    seq = ""
                prev_chain = atom.chain
                prev = -9999
            if atom.residue_id != prev:
                seq += sf.aa321[ atom.residue_name]
            prev_id = atom.residue_id
    if len(seq) > 0:
        sequences[prev_chain] = seq                
    return sequences


def MatchingPositions( first_file, second_file, chains = []):
    from Bio import pairwise2 as pw
    from Bio.pairwise2 import format_alignment as alig
    #print "read!",chains
    seqs1 = ReadSequencesFromCoordinates( first_file, chains)
    seqs2 = ReadSequencesFromCoordinates( second_file, chains)
    #print first_file, seqs1
    #print second_file, seqs2
    #print "align"
    positions = {}
    for c,seq1 in seqs1.items():
        #print "align", c
        seq2 = seqs2[c]
        alignment = pw.align.globalxx( seq1, seq2 )
        alig1 = alignment[0][0]
        alig2 = alignment[0][1]
        #print alig1
        #print
        #print alig2
        c1 = 0
        c2 = 0
        pos = []
        for i in range( 0, len( alig1 )):
            if alig1[i] == '-' or alig2[i] == '-':
                if alig1[i] != '-': c1 += 1
                if alig2[i] != '-': c2 += 1
                continue
            pos.append( [c1,c2] )
            c1 += 1
            c2 += 1
        positions[c] = pos
    return positions

def ReadChainCoordinates( filename, atoms = [] ):
    chains = {}
    positions = []
    prev_chain = 'sdaadf'
    
    with open( filename) as f:
        for l in f:
            if l[0:4] != "ATOM" and not l[0:6] == "HETATM": continue
            atom = PDBAtom( l)
            #print atom.atom_name, atom.chain
            if len(atoms) > 0 and atom.atom_name not in atoms:
                continue
            #print l,
            if atom.chain != prev_chain:
                if len(positions) > 0:
                    chains[prev_chain] = positions
                    positions = []
                prev_chain = atom.chain
            positions.append( atom.position )
    if len(positions) > 0:
        chains[prev_chain] = positions
    return chains

            
            

def MatchingPositionsRMSD( first_file, second_file, pos):
    from collections import defaultdict
    pos1 = ReadChainCoordinates( first_file, ['CA'] )
    pos2 = ReadChainCoordinates( second_file, ['CA'] )

    total = 0
    c = 0
    for chain, positions in pos.items():
        rmsd = 0.0
        #print chain, len(positions)
        for p in positions:
            rmsd += np.power( pos1[chain][p[0]] - pos2[chain][p[1]] , 2 ).sum()
            c += 1
        #print chain,  np.sqrt( rmsd / float(len(positions)))
        total += rmsd
    return np.sqrt(total / float(c) )

    
#def BasicRMSD( first_file, second_file, chains = []):
#    pdb_parser = Bio.PDB.PDBParser( QUIET=True )
#
#    first_structure = pdb_parser.get_structure( "first", first_file )[0]
#    second_structure = pdb_parser.get_structure( "second", second_file )[0]
#
#    first_atoms = []
#    second_atoms = []
#
#    for chain in first_structure:
#        if len(chains) > 0 and chain not in chains:
#            continue
#        for residue in chain:
#            first_atoms.append( residue['CA'] )
#
#    for chain in second_structure:
#        if len(chains) > 0 and chain not in chains:
#            continue
#        for residue in chain:
#            second_atoms.append( residue['CA'] )
#
#    if len(second_atoms) != len(first_atoms):
#        print "WARNING: number of atoms do not match!", len(first_atoms),len(second_atoms)
#
#        
#    #chain by chain!!!
#
#        
#    #alignment = \
#
#
#    length = len( alignment[0] )
#    for i in range(0,length):
#        if alignment[0][i] != '-' and alignment[1][i] != '-':
#            ids.append( i )
#        
#    for record in alignment:
#        print record.seq
#        for l in 
#
#    
#    for a1, a2 in zip( first_atoms, second_atoms):
#        if a1.get_name() != a2.get_name():
#            print "WARNING: mismatch:", a1,a2
#        rms += np.power( a1 - a2 , 2 )
#    return rms / float(len(first_atoms))


class PDBAtom:
    atom_name = ""
    residue_name = ""
    chain = ''
    residue_id = -1
    position = [0.0] * 3

    def __init__( self, line = None):
        if line == None:
            self.atom_name = ""
            self.residue_name = ""
            self.chain = ""
            self.residue_id = -9999
            self.position = np.array( [ -9999,-9999,-9999 ] )
        else:
            self.atom_name = line[12:16].strip()
            self.residue_name = line[17:20].strip()
            self.chain = line[21]
            self.residue_id = int( line[22:26] )
            self.position = np.array( [  float( line[30:38] ), float( line[38:46] ), float( line[46:54] ) ] )

def Distance( a1, a2):
    return vf.Distance( a1.position, a2.position )

def Angle( a1, a2, a3):
    return vf.Angle( a1.position - a2.position, a3.position - a2.position )

def Dihedral( a1, a2, a3, a4):
    return vf.Dihedral( a1.position, a2.position, a3.position, a4.position)



