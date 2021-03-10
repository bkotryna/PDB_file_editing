"""
Goal:
    Generate CONECT lines that need to be added to the end of a pdb file (to describe formyl group atoms)
Input:
    (1) A pdb file with 5FC atoms already inserted
    (2) A template for CONECT lines
    (in the form of a pickled nested list with atom names
    that will need to be translated into atom positions)
Actions:
    (1) Extract individual 5FC groups
    (2) Look up the positions that correspond to the atom names and "translate" the lines into atom positions
    using the pickeld nested list as the template
    (3) Assemble CONECT line text, write it out into a file.
    Could be a separate file for each formyl group and then also one joint file for the whole pdb file
    
Output:
    (1) A txt file for each 5FC group where each line is a CONECT line with atom positions
    (2) A txt file where each line is a CONECT line with atom positions,
    covering the whole pdb file (ie ~140 5FC groups)
Notes:
    (1) I've attempted to tidy things up a bit by using functions (instead of one massive roll of code)
"""


def clear_files():
    # clean any previously generated files
    import sys, os
    os.chdir('C:\Kotryna\Python\Lab projects\pdb file handling\Input_Output_files')
    cwd = os.getcwd()
    print("cwd is: ", cwd)
    
    with open('run_log.txt', 'w') as file:
        file.close()
    with open('joint_conect_lines.txt', 'w') as file:
        file.close()    



def get_conect_lines():
    # acquire atom name dictionary for CONECT lines
    import pickle
    with open('C:/Kotryna/Python/Lab projects/pdb file handling/Input_Output_files/atom_name_list_for_conect', 'rb') as file:
        conect_lines = pickle.load(file)
    conect_lines    # just print it out to inspect    
    return(conect_lines)



def get_all_atoms():
    # read in a pdb file
    file = bpd().read_pdb('nuc_with_5FC.pdb')
    
    # create a dataframe that contains all atoms
    atoms = file.df['ATOM']
    atoms
    
    hetatm = file.df['HETATM']
    hetatm
    
    all_atoms = atoms.append(hetatm)
    all_atoms
    all_atoms.head(50)
    all_atoms.tail(50)
    
    # clean up the order - make sure atom numbers are nicely ascending
    all_atoms = all_atoms.sort_values(['atom_number'], ascending=[True])
    all_atoms
    all_atoms.head(50)
    all_atoms.tail(50)
    
    return(all_atoms)



def find_range_start(c5_no, all_atoms):
 
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
    
    return(range_start)



def find_range_finish(c5_no, all_atoms):

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
    
    return(range_finish)



def find_O3_before(range_start, all_atoms):

    check_O3_atom_name = None
    i = range_start
    while check_O3_atom_name != "O3'":
        check_O3_atom_name = all_atoms[all_atoms.atom_number == i].atom_name.item()
        i -= 1
    O3_before_atom_no = i + 1 
    O3_before_atom_no
    
    O3_before_atom_nm_check = all_atoms[all_atoms.atom_number == O3_before_atom_no].atom_name.item()   
    O3_before_atom_nm_check
    if O3_before_atom_nm_check != "O3'":
        print("Stuff went wrong!! \"O3' before\" atom position not identified correctly")
    
    return(O3_before_atom_no)



def find_positions_dict(c5_no, all_atoms, range_start, range_finish):
                  
    # create a relevant range of atom numbers to consider around this C5 atom
    # this "relevant range" covers all 5FC atoms (and no other atoms??? Ig)
    atom_range = range(range_start, range_finish, 1)
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
    
    # add atoms before & after the 5FC group, also needed as input for making CONECT lines
    # for this, first generate dictionary entries and then append them to our name-position dictionary
    # for each entry, we give a manual unique name and look up the relrevant atom number based on empirical observations
    
    # add atoms before & after the 5FC group, also needed as input for making CONECT lines
    # for this, first generate dictionary entries and then append them to our name-position dictionary
    # for each entry, we give a manual unique name and look up the relrevant atom number
    # this is found by scrolling through atoms above 5FC group for "O3' before"
    # and based on empirical observations for "P after"
    
    O3_before_atom_no = find_O3_before(range_start, all_atoms)
    O3_before_atom_nm = "O3' before"
    
    # generate entry for "P after"
    P_after_atom_no = range_finish
    P_after_atom_nm = "P after"
    
    # add both entries to the dictionary
    positions_dict.update({O3_before_atom_nm : O3_before_atom_no, P_after_atom_nm : P_after_atom_no})
    positions_dict

    return(positions_dict) 



def find_number_lines(c5_no, all_atoms, positions_dict, conect_lines):
    
    ### use CONECT line template (in the form of a nested list)
    # as a template to create CONECT lines for particular 5FC group in hadn
       
    # generate CONECT lines for this one 5FC group
    
    ### obtain CONECT position coordinates
    # read out stuff from the parent list element-by-element (ie child-list-by-child-list) from the CONECT line template
    # within each element (ie child list), take each element (ie atom name) and convert it to atom number
    # (via the dictionary that we've just created)
    
    number_lines = []
    # go through parent-list element-by-element
    for list in conect_lines:
        #print(list)
        # go through the child-list element-by-element
        positions = []
        for name in list:
            #print(name)
            # retrieve the atom position that corresponds to each atom name
            #print(int(positions_dict[name]))
            # store atom position in a list
            positions.append(int(positions_dict[name]))
            #print(positions)
        # add the list containing atom names (aka one CONECT line) to the list of CONECT lines
        number_lines.append(positions)
    
    return(number_lines)



def write_conect_txt(c5_no, number_lines):

    ### make nice txt output file for CONECT info for each 5FC group
    # change directory to Output_files
    import sys, os
    os.chdir('C:\Kotryna\Python\Lab projects\pdb file handling\Input_Output_files')
    cwd = os.getcwd()
    #print("cwd is: ", cwd)
           
#    conect_filename = 'neat_conect_lines_c5_' + str(c5_no) + '.txt'
#    with open(conect_filename, 'w') as file:    # this way of opening it deletes file contents first, before writing stuff
#    # go through parent-list element-by-element
#        for list in number_lines:
#            text = "CONECT "
#            for number in list:
#                text += str(number) + "  "
#            text += "\n"
#            #print(text)
#            file.write(text)
        
    ### append CONECT lines to one joint txt file for the whole pdb
    joint_conect_filename = 'joint_conect_lines.txt'
    with open(joint_conect_filename, 'a') as file:    # this way of opening it just keeps adding stuff to the file without deleting anything
    # go through parent-list element-by-element
        for list in number_lines:
            text = "CONECT "
            for number in list:
                text += str(number) + "  "
            text += "\n"
            #print(text)
            file.write(text) 



def find_names_dict(c5_no, all_atoms, range_start, range_finish):

    # create a relevant range of atom numbers to consider around this C5 atom
    # this "relevant range" covers all 5FC atoms (and no other atoms??? Ig)
    atom_range = range(range_start, range_finish, 1)
    atom_range
    
    # extract atoms from this "relevant range"
    within_range = all_atoms[all_atoms.atom_number.isin(atom_range)]
    within_range
    
    
    
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
    
    O3_before_atom_no = find_O3_before(range_start, all_atoms)
    O3_before_atom_nm = "O3' before"
    
    # generate entry for "P after"
    P_after_atom_no = range_finish
    P_after_atom_nm = "P after"
    
    # add both entries to the dictionary
    names_dict.update({O3_before_atom_no : O3_before_atom_nm, P_after_atom_no : P_after_atom_nm})
    names_dict

    return(names_dict)



def find_name_lines(c5_no, all_atoms, names_dict, number_lines, conect_lines): 
    
    # now, time for the real check: translate CONECT lines from atom positions into atom names
    name_lines = []
    
    # go through parent-list element-by-element
    for list in number_lines:
        #print(list)
        # go through the child-list element-by-element
        names = []
        for position in list:
            #print(position)
            # retrieve the atom position that corresponds to each atom name
            #print(names_dict[position])
            # store atom position in a list
            names.append(names_dict[position])
            #print(names)
        # add the list containing atom names (aka one CONECT line) to the list of CONECT lines
        name_lines.append(names)
    
    #list of lines, expressed in atom numbers rather than positions (each line is a list)
    name_lines

    return(name_lines)



def write_run_log_txt(c5_no, name_lines, conect_lines):

    # change directory to Output_files
    import sys, os
    os.chdir('C:\Kotryna\Python\Lab projects\pdb file handling\Input_Output_files')
    cwd = os.getcwd()
    #print("cwd is: ", cwd)
    
    mismatch1 = [x for x in name_lines if x not in conect_lines]
    mismatch1
    
    mismatch2 = [x for x in conect_lines if x not in name_lines]
    mismatch2
    
    if name_lines != conect_lines:
        print(str(c5_no), ": Stuff went wrong!! Atom positions and names don't match in CONECT lines.")
        with open('run_log.txt', 'a') as file:
            text = "C5 no is: " + str(c5_no) + ".\nAnd the final back-check of CONECT line info assignments is:\nsth is wrong.\n\n"
            file.write(text)
    else:
        #print(str(c5_no), ": Good stuff, all okay. CONECT line atom positions match the expected atom names.")
        with open('run_log.txt', 'a') as file:
            text = "C5 no is: " + str(c5_no) + ".\nAnd the final back-check of CONECT line info assignments is:\nall good.\n\n"
            file.write(text)


#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************


### get started

clear_files()

# set working directory to where your input pdb file is
import sys, os
os.chdir('C:\Kotryna\Python\Lab projects\pdb file handling\Input_Output_files')
cwd = os.getcwd()
print("cwd is: ", cwd)

# import tools you'll need
import pandas as pd
from biopandas.pdb import PandasPdb as bpd
import numpy as np


conect_lines = get_conect_lines()
all_atoms = get_all_atoms()


# go through the whole big pdb file
# first, loop through DNA chains
DNA_chains = ['I', 'J']
for one_chain in DNA_chains:
    
    chain_name = one_chain
    print("Chain name is: ", chain_name)
    
    # get all atoms that belong to this particular DNA chain
    chain = all_atoms[all_atoms.chain_id == chain_name]
    chain.head()
    
    # get all atoms that belong to 5FC cytosines
    cytosines = chain[chain.residue_name == '5FC']
    cytosines.head(20)
    
    # get all C5 atoms that belong to 5FC cytosines
    c5_atoms = cytosines[cytosines.atom_name == 'C5']
    c5_atoms.head()

    # extract atom numbers for the 5FC-C5 atoms    
    c5_atom_numbers = c5_atoms.atom_number.tolist()
    c5_atom_numbers
    type(c5_atom_numbers)    # just to check
    
    # check: how many C5 atoms there are (ie how many cytosines we have)
    c5_atoms.residue_name.count()
    # alternative, less neat counting of everything
    c5_atoms.count()
    
    # the actual loop through all C5 atoms
    for c5_no in c5_atom_numbers:
        print("c5_no is: ", str(c5_no))
        
        # NB the atom number of the C5 atom is c5_no

        range_start = find_range_start(c5_no, all_atoms)
        range_finish = find_range_finish(c5_no, all_atoms)
        
        positions_dict = find_positions_dict(c5_no, all_atoms, range_start, range_finish)
        number_lines = find_number_lines(c5_no, all_atoms, positions_dict, conect_lines)
        write_conect_txt(c5_no, number_lines)
        
        names_dict = find_names_dict(c5_no, all_atoms, range_start, range_finish)
        name_lines = find_name_lines(c5_no, all_atoms, names_dict, number_lines, conect_lines)
        write_run_log_txt(c5_no, name_lines, conect_lines)

