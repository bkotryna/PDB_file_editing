### get started: set directory, import tools, read in your pdb file

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
file = bpd().read_pdb('pdb_with_5FC.pdb')

# create a dataframe that contains all atoms
atoms = file.df['ATOM']
atoms

hetatm = file.df['HETATM']
hetatm

all_atoms = atoms.append(hetatm)
all_atoms

all_atoms.tail(50)


## write atoms out to csv for reference
#atoms.to_csv("../Output_files/original_atoms.csv", encoding='utf-8', index=False)

## what are the column names
#atoms.columns

### extract info from the pdb file for use in creating CONECT lines

## create a file to store CONECT info
#with open('../Output_files/conect_info.txt', "w") as f:
#    f.write("atom" + ','+ "number")













c5_no = 19696



# NB the atom number of the C5 atom is c5_no

# get the pdb line for this particular C5 atom
c5_atom = all_atoms[all_atoms.atom_number == c5_no]
c5_atom
  
# create a relevant range of atom numbers to consider around this C5 atom
# this "relevant range" covers all 5FC atoms (and no other atoms??? Ig)
range_up = 12
range_down = 21
atom_range = range(c5_no - range_up, c5_no + range_down, 1)
atom_range

# extract atoms from this "relevant range"
within_range = all_atoms[all_atoms.atom_number.isin(atom_range)]
within_range

# extract atom_name and atom_number from the "relevant range"
# write it out in a dictionary
# will need to use this to write out CONECT lines
range_atom_nos = within_range.atom_number.tolist()
range_atom_nos

# make a dictionary where keys =  atom names and vbalues =  atom positions, within this "relevant range"
positions_dict = {}
for atom_no in range_atom_nos:
    atom_nm = within_range[within_range.atom_number == atom_no].atom_name.item()
    positions_dict.update({atom_nm : atom_no})
positions_dict

# add atoms before & after the 5FC group, also needed as input for making CONECT lines
# for this, first generate dictionary entries and then append them to our name-position dictionary
# for each entry, we give a manual unique name and look up the relrevant atom number based on empirical observations

# generate entry for "O3' before"
O3_before_atom_no = c5_no - 25  # the 25 is the empirically derived number :)
O3_before_atom_no
O3_before_atom_nm_check = all_atoms[all_atoms.atom_number == O3_before_atom_no].atom_name.item()
O3_before_atom_nm_check
O3_before_atom_nm = "O3' before"

# generate entry for "P after"
P_after_atom_no = int(c5_no + 9) # the 9 is the empirically derived number :)
P_after_atom_no
P_after_atom_nm_check = all_atoms[all_atoms.atom_number == P_after_atom_no].atom_name.item()
P_after_atom_nm_check
P_after_atom_nm = "P after"

# add both entries to the dictionary
positions_dict.update({O3_before_atom_nm : O3_before_atom_no, P_after_atom_nm : P_after_atom_no})
positions_dict


   
#        # dump dictionary items via json for potential use outside this script
#        import json 
#
#        dict_filename = '../Output_files/chain_' + chain_name + '_' + 'dict_c5_' + str(c5_no) + '_conect_info.csv'
#        with open(dict_filename, 'w') as file:
#            file.write(json.dumps(positions_dict)) # use `json.loads` to do the reverse
#        
#              
#        # append info to a text file      
#        with open('../Output_files/connect_info.csv', "a") as f:
#            for k, v in positions_dict.items():
#                f.write('\n' + str(k) + ','+ str(v))
#                
#        # write out individual text files for each 5FC group
#        filename = '../Output_files/c5_' + str(c5_no) + '_conect_info.csv'
#        with open(filename, "a") as file:
#            for k, v in positions_dict.items():
#                file.write('\n' + str(k) + ','+ str(v))
#
#       # check if we can read in the new text file as a pandas dataframe
#       check_conect_info = pd.read_csv('../Output_files/conect_info.txt', encoding='utf-8')
#       check_conect_info


### use CONECT line template (in the form of a nested list)
# as a template to create CONECT lines for particular 5FC group in hadn

# acquire atom name dictionary for CONECT lines
import pickle
with open('C:/Kotryna/Python/Lab projects/pdb file handling/Input_Output_files/atom_name_list_for_conect', 'rb') as file:
    conect_lines = pickle.load(file)
conect_lines    # just print it out to inspect

# generate CONECT lines for this one 5FC group

### obtain CONECT position coordinates
# read out stuff from the parent list element-by-element (ie child-list-by-child-list) from the CONECT line template
# within each element (ie child list), take each element (ie atom name) and convert it to atom number
# (via the dictionary that we've just created)



number_lines = []

# go through parent-list element-by-element
for list in conect_lines:
    print(list)
    # go through the child-list element-by-element
    positions = []
    for name in list:
        print(name)
        # retrieve the atom position that corresponds to each atom name
        print(int(positions_dict[name]))
        # store atom position in a list
        positions.append(int(positions_dict[name]))
        print(positions)
    # add the list containing atom names (aka one CONECT line) to the list of CONECT lines
    number_lines.append(positions)

#list of lines, expressed in atom numbers rather than positions (each line is a list)
number_lines

#        # pickle number_lines for use in other scripts
#        
#        import pickle
#        
#        import sys, os
#        os.chdir('C:\Kotryna\Python\Lab projects\pdb file handling\Output_files')
#        cwd = os.getcwd()
#        print("cwd is: ", cwd)
#        
#        with open('atom_positions_list_for_conect', 'wb') as file:
#            pickle.dump(number_lines, file)        

### make nice txt output file for CONECT info for each 5FC group
# change directory to Output_files
import sys, os
os.chdir('C:\Kotryna\Python\Lab projects\pdb file handling\Input_Output_files')
cwd = os.getcwd()
print("cwd is: ", cwd)
       
conect_filename = 'neat_conect_lines_c5_' + str(c5_no) + '.txt'
with open(conect_filename, 'w') as file:    # this way of opening it deletes file contents first, before writing stuff
# go through parent-list element-by-element
    for list in number_lines:
        text = "CONECT "
        for number in list:
            text += str(number) + "  "
        text += "\n"
        print(text)
        file.write(text)
    
### append CONECT lines to one joint txt file for the whole pdb
joint_conect_filename = 'joint_conect_lines.txt'
with open(joint_conect_filename, 'a') as file:    # this way of opening it just keeps adding stuff to the file without deleting anything
# go through parent-list element-by-element
    for list in number_lines:
        text = "CONECT "
        for number in list:
            text += str(number) + "  "
        text += "\n"
        print(text)
        file.write(text)     





### check the output: express it in atom names again, see if they match the CONECT line template

# extract atom_name and atom_number from the "relevant range"
range_atom_nms = within_range.atom_name.tolist()
range_atom_nms
        
# make a reverse dictionary where keys =  atom positions and values =  atom names, within your "relevant range"
names_dict = {}
for atom_nm in range_atom_nms:
    atom_no = within_range[within_range.atom_name == atom_nm].atom_number.item()
    names_dict.update({atom_no : atom_nm})
names_dict

# add atoms before & after the 5FC group, also needed as input for making CONECT lines
# for this, first generate dictionary entries and then append them to our name-position dictionary
# for each entry, we give a manual unique name and look up the relrevant atom number based on empirical observations

# add atoms before & after the 5FC group, also needed as input for making CONECT lines
# for this, first generate dictionary entries and then append them to our name-position dictionary
# for each entry, we give a manual unique name and look up the relrevant atom number
# this is found by scrolling through atoms above 5FC group for "O3' before"
# and based on empirical observations for "P after"

start = c5_no - range_up - 1    # this gives the first atom above 5FC group
start
start_check = all_atoms[all_atoms.atom_number == start].atom_name.item()
start_check

check_atom_name = None
i = start
while check_atom_name != "O3'":
    check_atom_name = all_atoms[all_atoms.atom_number == i].atom_name.item()
    i -= 1
O3_before_atom_no = i + 1 
O3_before_atom_no

O3_before_atom_nm_check = all_atoms[all_atoms.atom_number == O3_before_atom_no].atom_name.item()   
O3_before_atom_nm_check
if O3_before_atom_nm_check != "O3'":
    print("Stuff went wrong!! \"O3' before\" atom position not identified correctly")
O3_before_atom_nm = "O3' before"

# generate entry for "P after"
P_after_atom_no = int(c5_no + range_down) # the range_down is the empirically derived number :)
P_after_atom_no
P_after_atom_nm_check = all_atoms[all_atoms.atom_number == P_after_atom_no].atom_name.item()
P_after_atom_nm_check
if P_after_atom_nm_check != "P":
    print("Stuff went wrong!! \"P after\" atom position not identified correctly")
P_after_atom_nm = "P after"

# add both entries to the dictionary
names_dict.update({O3_before_atom_no : O3_before_atom_nm, P_after_atom_no : P_after_atom_nm})
names_dict



# now, time for the real check: translate CONECT lines from atom positions into atom names
name_lines = []

# go through parent-list element-by-element
for list in number_lines:
    print(list)
    # go through the child-list element-by-element
    names = []
    for position in list:
        print(position)
        # retrieve the atom position that corresponds to each atom name
        print(names_dict[position])
        # store atom position in a list
        names.append(names_dict[position])
        print(names)
    # add the list containing atom names (aka one CONECT line) to the list of CONECT lines
    name_lines.append(names)

#list of lines, expressed in atom numbers rather than positions (each line is a list)
name_lines

mismatch1 = [x for x in name_lines if x not in conect_lines]
mismatch1

mismatch2 = [x for x in conect_lines if x not in name_lines]
mismatch2

if name_lines != conect_lines:
    print("Stuff went wrong!! Atom positions and names don't match in CONECT lines.")
else: print("Good stuff, all okay. CONECT line atom positions match the expected atom names.")



**********************************************************************************************************
**********************************************************************************************************
**********************************************************************************************************



### create a good range around c5 atom


### get started: set directory, import tools, read in your pdb file

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
file = bpd().read_pdb('pdb_with_5FC.pdb')

# create a dataframe that contains all atoms
atoms = file.df['ATOM']
atoms

hetatm = file.df['HETATM']
hetatm

all_atoms = atoms.append(hetatm)
all_atoms

all_atoms.tail(50)










c5_no = 19696



# NB the atom number of the C5 atom is c5_no

# get the pdb line for this particular C5 atom
c5_atom = all_atoms[all_atoms.atom_number == c5_no]
c5_atom


### establish a range where the 5FC residue is

# find the first P above c5 (ie beginning of the 5FC nucleotide)
first_P_name = None
count_P = c5_no
while first_P_name != "P":
    first_P_name = all_atoms[all_atoms.atom_number == count_P].atom_name.item()
    count_P -= 1
first_P_no = count_P + 1 
first_P_no

first_P_name_check = all_atoms[all_atoms.atom_number == first_P_no].atom_name.item()
first_P_name_check
if first_P_name_check != "P":
    print("Stuff went wrong!! \"first P\" atom position not identified correctly")

range_start = first_P_no
range_start

# find the first P below c5 (ie beginning of a new nucleotide)
next_P_name = None
count_Pn = c5_no
while next_P_name != "P":
    next_P_name = all_atoms[all_atoms.atom_number == count_Pn].atom_name.item()
    count_Pn += 1
next_P_no = count_Pn - 1 
next_P_no

next_P_name_check = all_atoms[all_atoms.atom_number == next_P_no].atom_name.item()
next_P_name_check
if next_P_name_check != "P":
    print("Stuff went wrong!! \"first P\" atom position not identified correctly")

range_finish = next_P_no    # this is the first entry outside the 5FC nucleotide
range_finish

  
# create a relevant range of atom numbers to consider around this C5 atom
# this "relevant range" covers all 5FC atoms (and no other atoms??? Ig)
atom_range = range(range_start, range_finish, 1)
atom_range

# extract atoms from this "relevant range"
within_range = all_atoms[all_atoms.atom_number.isin(atom_range)]
within_range



##############################################################
### obtain precise range of c5 atoms that need formyl groups to be added onto them

### get started: set directory, import tools, read in your pdb file

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
file = bpd().read_pdb('nuc.pdb')

# create a dataframe that contains all atoms
atoms = file.df['ATOM']
atoms

hetatm = file.df['HETATM']
hetatm

all_atoms = atoms.append(hetatm)
all_atoms

all_atoms.tail(50)

chain_name = "I"

chain = atoms[atoms.chain_id == chain_name]
chain.head()

# get cytosines
cytosines = chain[chain.residue_name == 'DC']
cytosines.head(20)

# get C5 atoms from cytosines
c5_atoms = cytosines[cytosines.atom_name == 'C5']
c5_atoms.head(10)
c5_atoms.tail(10)


# narrow down to the c5 atoms that actually are supposed to have formyl groups
# ie skip some c5 atoms in the beginning and end
# BTW "some" = look up which ones in nucleosome info word.doc file
# ie first six and last six c5 atoms don't need a formyl group


# check: how many C5 atoms there are (ie how many cytosines we have)
starting_no_of_c5 = c5_atoms.residue_name.count()

# the narrowed down array of c5 atoms
c5_atoms = c5_atoms.iloc[6:(starting_no_of_c5 - 6)] 

# inspect what's happened
c5_atoms.head(10)
c5_atoms.tail(10)

# count number of C5 atoms in 5FC groups (aka number of 5FC groups) - just to check
c5_in_5fc = atoms.loc[(atoms['residue name'] == '5FC') & (atoms['atom name'] == 'C5')]
c5_in_5fc_count = c5_in_5fc.count()
print(c5_in_5fc_count)



###################################################

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

###################################################


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

    
chain_name = "I"

chain = atoms[atoms.chain_id == chain_name]
chain.head()

# get cytosines
cytosines = chain[chain.residue_name == 'DC']
cytosines.head(20)

# get C5 atoms from cytosines
c5_atoms = cytosines[cytosines.atom_name == 'C5']
c5_atoms.head(10)
c5_atoms.tail(10)

# narrow down to the c5 atoms that actually are supposed to have formyl groups
# ie skip some c5 atoms in the beginning and end
# BTW "some" = look up which ones in nucleosome info word.doc file
# ie first six and last six c5 atoms don't need a formyl group


# check: how many C5 atoms there are (ie how many cytosines we have)
starting_no_of_c5 = c5_atoms.residue_name.count()

# the narrowed down array of c5 atoms (lacks six atoms at the beginning & six atoms at the end)
c5_atoms = c5_atoms.iloc[6:(starting_no_of_c5 - 6)] 

# inspect what's happened
c5_atoms.head(10)
c5_atoms.tail(10)

# proceed
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








    # rename pdb identifiers for cytosines: change "record_name" from ATOM into HETATM
    # and change "residue_name" from DC into 5FC
    atoms.loc[(atoms.atom_number.isin(atom_range)), "residue_name"] = "5FC"
    atoms.loc[(atoms.atom_number.isin(atom_range)), "record_name"] = "HETATM"
#        # rename residue within range - go from "DC" to "5FC"
#        atoms.loc[atoms.residue_name == 'DC', 'residue_name'] = '5FC'
#        # rename record name from "ATOM" to "HETATM"
#        atoms.loc[atoms.residue_name == '5FC', 'record_name'] = 'HETATM'






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
    c5a_atom = pd.DataFrame({'record_name': 'HETATM', 'atom_number': '', 'blank_1': '', 'atom_name': 'C5A', 'alt_loc': '', 'residue_name': '5FC', 'blank_2': '', 'chain_id': chain_name, 'residue_number': c5_residue_number, 'insertion': '', 'blank_3': '', 'x_coord': c5a_coords["x"], 'y_coord': c5a_coords["y"], 'z_coord': c5a_coords["z"], 'occupancy': 1, 'b_factor': dummy_b, 'blank_4': '', 'segment_id': c5_segment_id, 'element_symbol': 'C', 'charge': c5_charge, 'line_idx': ''}, index = [c5_index + 5])
    c5a_atom
    # create a data-row for O5A
    o5a_atom = pd.DataFrame({'record_name': 'HETATM', 'atom_number': '', 'blank_1': '', 'atom_name': 'O5A', 'alt_loc': '', 'residue_name': '5FC', 'blank_2': '', 'chain_id': chain_name, 'residue_number': c5_residue_number, 'insertion': '', 'blank_3': '', 'x_coord': o5a_coords["x"], 'y_coord': o5a_coords["y"], 'z_coord': o5a_coords["z"], 'occupancy': 1, 'b_factor': dummy_b, 'blank_4': '', 'segment_id': c5_segment_id, 'element_symbol': 'O', 'charge': c5_charge, 'line_idx': ''}, index = [c5_index + 6])
    o5a_atom          
    # create a data-row for H5A
    h5a_atom = pd.DataFrame({'record_name': 'HETATM', 'atom_number': '', 'blank_1': '', 'atom_name': 'H5A', 'alt_loc': '', 'residue_name': '5FC', 'blank_2': '', 'chain_id': chain_name, 'residue_number': c5_residue_number, 'insertion': '', 'blank_3': '', 'x_coord': h5a_coords["x"], 'y_coord': h5a_coords["y"], 'z_coord': h5a_coords["z"], 'occupancy': 1, 'b_factor': dummy_b, 'blank_4': '', 'segment_id': c5_segment_id, 'element_symbol': 'H', 'charge': c5_charge, 'line_idx': ''}, index = [c5_index + 7])
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
atoms[atoms.residue_name == '5FC']

atoms.residue_name.loc[(atoms['residue_name'] == '5FC') & (atoms['atom_name'] == 'C5')].count()

# count number of C5 atoms in 5FC groups (aka number of 5FC groups) - just to check
c5_in_5fc = atoms.loc[(atoms['residue_name'] == '5FC') & (atoms['atom_name'] == 'C5')]
c5_in_5fc_count = c5_in_5fc.count()
print(c5_in_5fc_count)


atoms.loc[(atoms['residue_name'] == 'DC') & (atoms['atom_name'] == 'C5')].count()
atoms.loc[(atoms['residue_name'] == 'DG') & (atoms['atom_name'] == 'C5')].count()

