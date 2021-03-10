"""
Goal:
    Generate a template for CONECT lines
    (basically, CONECT info expressed in atom names, which can then
    be translated into atom positions for each individual formyl group)
Input:
    (1) A snippet of a pdb that contains all atoms from one 5FC group
    (2) An example txt file with all relevant CONECT lines for one 5FC group
    (this will act as a template for CONECT line info)
Actions:
    (1) Take the pdb snippet: Make a dictionary where atom positions can be looked up using their names
    (2) Take example txt file: Convert the numbers in the CONECT info into atom names
Output:
    (1) A nested list that contains CONECT line info, expressed in atom names
Notes:
    
"""


# set working directory to where your input pdb file is
import sys, os
os.chdir('C:\Kotryna\Python\Lab projects\pdb file handling\Input_Output_files')
cwd = os.getcwd()
print("cwd is: ", cwd)

# import tools you'll need
import pandas as pd
from biopandas.pdb import PandasPdb as bpd
import numpy as np

### make a dictionary to match atom numbers to atom names
# open pdb snippet that contains the relevant atom names & numbers
file = bpd().read_pdb('pdb_snippet_for_conect.txt')

# create a dataframe that contains all atoms
atoms = file.df['ATOM']
atoms

hetatm = file.df['HETATM']
hetatm

all_atoms = atoms.append(hetatm)
all_atoms

# extract atom numbers & names into lists
atom_numbers = all_atoms.atom_number.tolist()
atom_names = all_atoms.atom_name.tolist()
length = len(atom_numbers)
length

# create dictionary with atom numbers as keys and atom names as values
names_dict = {}
for i in range(0,length):
    k = atom_numbers[i]
    v = atom_names[i]
    names_dict.update({k:v})
names_dict
len(names_dict)
names_dict[243] # just checking

# add to the dictionary two more atoms that are needed: "O3' before" and "P after"
names_dict.update({198:"O3' before"})
names_dict.update({255:"P after"})

names_dict
len(names_dict)
names_dict[255] # just checking

### obtain CONECT position coordinates
# create a list for CONECT lines. Each item in the list will be info for one CONECT line
conect_lines = []
filename = "template_for_conect.txt"
with open(filename, 'r') as file:
    #count = 1
    # go through file line-by-line
    for line in file:
        # split each line to separate "words"
        words = line.split()
        # delete the first word (which is always "CONECT")
        del words[0]
        print(words)

        # create a list for atom positions from each CONECT line
        names = []
        # add each word from the line into a dictionary: each "word" is atom position
        for word in words:
            print(int(word))
            # retrieve the atom name that corresponds to each atom positions
            print(names_dict[int(word)])
            # store atom name in a list
            names.append(names_dict[int(word)])
        print(names)
        # add the list containing atom names (aka one CONECT line) to the list of CONECT lines
        conect_lines.append(names)

        #count +=1

#list of lines, expressed in atom names rather than numbers (each line is a list)
conect_lines

# pickle conect_lines for use in other scripts

import pickle
with open('C:/Kotryna/Python/Lab projects/pdb file handling/Input_Output_files/atom_name_list_for_conect', 'wb') as file:
    pickle.dump(conect_lines, file)



