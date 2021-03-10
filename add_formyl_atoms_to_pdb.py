"""
Goal:
    Insert extra atoms at Cytosines in DNA - need to add C, H and O for a formyl group
Input:
    (1) A pdb file that describes our nucleosome structure
Actions:
    (1) Add extra lines to the pdb - need C, H, O atoms for each formyl group on each cytosine
    (BTW For the line to be nice and complete, need to extract/calculate lots of thingies)
Output:
    (1) A pdb file that contains lines for C, H, O atoms in formyl groups
Notes:
    
"""


#############################################
# part 1 = common to all
# set directory, import tools, read in your pdb file

# set working directory to where your input pdb file is
import sys, os
os.chdir('C:\Kotryna\Python\Lab projects\pdb file handling\Input_Output_files')
cwd = os.getcwd()
print("cwd is: ", cwd)

# import tools you'll need
import pandas as pd
from biopandas.pdb import PandasPdb as bpd
import numpy as np



# read in a pdb file
file1 = bpd().read_pdb('nuc_with_hydrogens.pdb')

# create a dataframe that contains all atoms
atoms = file1.df['ATOM']
atoms
# write atoms out to csv for reference
#atoms.to_csv("../Output_files/original_atoms.csv", encoding='utf-8', index=False)

# what are the column names
atoms.columns


#############################################
# part 2 = specific DNA chain
# get the cytosines & their C5 atoms

# loop to go through DNA chains

DNA_chains = ['I', 'J']
for one_chain in DNA_chains:
    
    chain_name = one_chain
    
    chain = atoms[atoms.chain_id == chain_name]
    chain.head()
    
    # get cytosines
    cytosines = chain[chain.residue_name == 'DC']
    cytosines.head(20)
    
    # get C5 atoms from cytosines
    c5_atoms = cytosines[cytosines.atom_name == 'C5']
    c5_atoms.head()
    
    c5_atom_numbers = c5_atoms.atom_number.tolist()
    c5_atom_numbers
    type(c5_atom_numbers)    # just to check
    
    # check: how many C5 atoms there are (ie how many cytosines we have)
    c5_atoms.residue_name.count()
    # alternative, less neat counting of everything
    c5_atoms.count()
    
    
    #############################################
    # part 3 = for each C5 atom
    
    # functions to extract xyz coordinates given atom data-row
    def get_coords(data_row):
    
        # extract coordinates
        x_coord = float(data_row.x_coord)
        y_coord = float(data_row.y_coord)
        z_coord = float(data_row.z_coord)
           
        # output a dictionary
        coords = {"x":x_coord, "y":y_coord, "z":z_coord}
        return coords
    
    
    # the actual loop through all C5 atoms
    for c5_no in c5_atom_numbers:
        # NB the atom number of the C5 atom is c5_no
    
        c5_atom = c5_atoms[c5_atoms.atom_number == c5_no]
        #c5_atom
      
        # create a relevant range of atom numbers to consider around this C5 atom
        range_up = 12
        range_down = 7
        atom_range = range(c5_no - range_up, c5_no + range_down, 1)
        #atom_range
    
        # extract those relevant atoms
        within_range = cytosines[cytosines.atom_number.isin(atom_range)]
        #within_range
    
    
        ### extract relevant info from atoms
        ### for C5 atom
    
        # extract the index of the C5 atom & the line number
        c5_index = c5_atom.index.item()
        #c5_index
        
        c5_line_idx = c5_atom.line_idx.item()
        #c5_line_idx
    
        # extract relevant "general" info
        c5_residue_number = c5_atom.residue_number.item()
        #c5_residue_number
        
        c5_chain_id = c5_atom.chain_id.item()
        #c5_chain_id
        
        c5_segment_id = c5_atom.segment_id.item()
        #c5_segment_id
    
        c5_charge = c5_atom.charge.item()
        #c5_charge
        
        # extract the coordinates of the C5 atom (output = dictionary)
        c5_coords = get_coords(c5_atom)
        #c5_coords
    
    
        ### for N1 atom  
        # find the N1 atom in the range
        n1_atom = within_range[within_range.atom_name == 'N1']
        #n1_atom    
        # extract the coordinates of the N1 atom (output = dictionary)
        n1_coords = get_coords(n1_atom)
        #n1_coords
    
    
        ### for C6 atom    
        # find the C6 atom in the range
        c6_atom = within_range[within_range.atom_name == 'C6']
        #c6_atom    
        # extract the coordinates of the C6 atom (output = dictionary)
        c6_coords = get_coords(c6_atom)
        #c6_coords


        ### for C4 atom    
        # find the C4 atom in the range
        c4_atom = within_range[within_range.atom_name == 'C4']
        #c6_atom    
        # extract the coordinates of the C4 atom (output = dictionary)
        c4_coords = get_coords(c6_atom)
        #c4_coords

    
    
    #    ### calculate distance vectors
    #    ### vector from N1 to C6
    #    n1_c6_x = c6_coords["x"] - n1_coords["x"]
    #    n1_c6_y = c6_coords["y"] - n1_coords["y"]
    #    n1_c6_z = c6_coords["z"] - n1_coords["z"]
    #
    #    # neat output dictionary
    #    n1_c6 = {"dist_x":n1_c6_x, "dist_y":n1_c6_y, "dist_z":n1_c6_z}    
    #    
    #    ### vector from C6 to C5
    #    c6_c5_x = c5_coords["x"] - c6_coords["x"]
    #    c6_c5_y = c5_coords["y"] - c6_coords["y"]
    #    c6_c5_z = c5_coords["z"] - c6_coords["z"]
    #
    #    # neat output dictionary
    #    c6_c5 = {"dist_x":c6_c5_x, "dist_y":c6_c5_y, "dist_z":c6_c5_z}   
    
    
        ### create xyz coords for C5A
        c5a_coords = {}
        c5a_coords["x"] = c5_coords["x"] + c6_coords["x"] - n1_coords["x"]
        c5a_coords["y"] = c5_coords["y"] + c6_coords["y"] - n1_coords["y"]
        c5a_coords["z"] = c5_coords["z"] + c6_coords["z"] - n1_coords["z"]
        c5a_coords
        
        ### create xyz coords for O5A
        o5a_coords = {}
        o5a_coords["x"] = c5a_coords["x"] + c5_coords["x"] - c6_coords["x"]
        o5a_coords["y"] = c5a_coords["y"] + c5_coords["y"] - c6_coords["y"]
        o5a_coords["z"] = c5a_coords["z"] + c5_coords["z"] - c6_coords["z"]
        o5a_coords

        ### create xyz coords for C5A
        h5a_coords = {}
        h5a_coords["x"] = c5a_coords["x"] + c5_coords["x"] - c4_coords["x"]
        h5a_coords["y"] = c5a_coords["y"] + c5_coords["y"] - c4_coords["y"]
        h5a_coords["z"] = c5a_coords["z"] + c5_coords["z"] - c4_coords["z"]
        h5a_coords

    
        # what index does C5 atom have in atoms dataframe
        # insert C5A and O5A right after C5 in atoms
        # renumber things in atoms
      
        # insert C5A & O5A
        # create a dummy b-factor for the new atoms
        dummy_b = c5_atom.b_factor.item()
        dummy_b
        type(dummy_b)   # just to check
        
        # create a data-row for C5A
        c5a_atom = pd.DataFrame({'record_name': 'ATOM', 'atom_number': '', 'blank_1': '', 'atom_name': 'C5A', 'alt_loc': '', 'residue_name': 'DC', 'blank_2': '', 'chain_id': chain_name, 'residue_number': c5_residue_number, 'insertion': '', 'blank_3': '', 'x_coord': c5a_coords["x"], 'y_coord': c5a_coords["y"], 'z_coord': c5a_coords["z"], 'occupancy': 1, 'b_factor': dummy_b, 'blank_4': '', 'segment_id': c5_segment_id, 'element_symbol': 'C', 'charge': c5_charge, 'line_idx': ''}, index = [c5_index + 5])
        c5a_atom
        # create a data-row for O5A
        o5a_atom = pd.DataFrame({'record_name': 'ATOM', 'atom_number': '', 'blank_1': '', 'atom_name': 'O5A', 'alt_loc': '', 'residue_name': 'DC', 'blank_2': '', 'chain_id': chain_name, 'residue_number': c5_residue_number, 'insertion': '', 'blank_3': '', 'x_coord': o5a_coords["x"], 'y_coord': o5a_coords["y"], 'z_coord': o5a_coords["z"], 'occupancy': 1, 'b_factor': dummy_b, 'blank_4': '', 'segment_id': c5_segment_id, 'element_symbol': 'O', 'charge': c5_charge, 'line_idx': ''}, index = [c5_index + 6])
        o5a_atom          
        # create a data-row for H5A
        h5a_atom = pd.DataFrame({'record_name': 'ATOM', 'atom_number': '', 'blank_1': '', 'atom_name': 'H5A', 'alt_loc': '', 'residue_name': 'DC', 'blank_2': '', 'chain_id': chain_name, 'residue_number': c5_residue_number, 'insertion': '', 'blank_3': '', 'x_coord': h5a_coords["x"], 'y_coord': h5a_coords["y"], 'z_coord': h5a_coords["z"], 'occupancy': 1, 'b_factor': dummy_b, 'blank_4': '', 'segment_id': c5_segment_id, 'element_symbol': 'H', 'charge': c5_charge, 'line_idx': ''}, index = [c5_index + 7])
        h5a_atom       
        
        # insert C5A and O5A into the atoms dataframe    
        atoms = pd.concat([atoms.loc[:c5_index+5], c5a_atom, o5a_atom, h5a_atom, atoms.loc[(c5_index + 6):]])
    
        # rearrange column positions to match pdb format
        atoms = atoms [['record_name', 'atom_number', 'blank_1', 'atom_name', 'alt_loc', 'residue_name', 'blank_2', 'chain_id', 'residue_number', 'insertion', 'blank_3', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'b_factor', 'blank_4', 'segment_id', 'element_symbol', 'charge', 'line_idx']]
        atoms
        
        
    #############################################
    # part 4 = common to all
    # renumber things
    # write out a csv file of atoms, just for fun
    # save into a new pdb file
    
    # redo indexing
    atoms.index = range(1, 1+len(atoms.index), 1)
    atoms.index
    
    # redo atom numbers
    atoms.atom_number = range(1, 1+len(atoms.index), 1)
    atoms.atom_number
    
    # redo line idx
    atoms.line_idx = range(1, 1+len(atoms.index), 1)
    atoms.line_idx
    
    # final product
    atoms
    atoms[atoms.residue_name == 'DC']

# rename pdb identifiers for cytosines: change "record_name" from ATOM into HETATM
# and change "residue_name" from DC into 5FC
## df.loc[<row selection>, <column selection>]

fiddle_atoms = atoms
fiddle_atoms

fiddle_atoms.loc[fiddle_atoms.residue_name == 'DC', 'residue_name'] = '5FC'
fiddle_atoms[fiddle_atoms.residue_name == '5FC']

fiddle_atoms.loc[fiddle_atoms.residue_name == '5FC', 'record_name'] = 'HETATM'
fiddle_atoms[fiddle_atoms.residue_name == '5FC']
        
atoms = fiddle_atoms         # assuming stuff worked, let's modify atoms dataframe


# random check for other purpose
# atoms[(atoms.atom_name == 'OP2') & (atoms.residue_name == '5FC')]



# write out all atoms into a csv file, just for fun
#fiddle_atoms.to_csv("../Output_files/updated_atoms.csv", encoding='utf-8', index=False)

### create a new pdb file

# check if pdb file creation works
#file1.to_pdb(path='./bur_pdb.pdb', records=None, gz=False, append_newline=True)

# write out our updated atoms into a new pdb file
new_file = file1
new_file.df['ATOM'] = atoms
new_file.df['ATOM']

new_file.to_pdb(path='../Input_Output_files/nuc_with_5FC.pdb', records=None, gz=False, append_newline=True)