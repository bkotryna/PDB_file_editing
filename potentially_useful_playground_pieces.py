"""
Goal:
    Store bits of coding practice that may be useful, but didn't end up being used in the final programmes
Input:
    
Actions:
    
Output:
    
Notes:
    
"""

troubleshooting testing

chain_name = 'I'
    
chain = atoms[atoms.chain_id == chain_name]
chain.head()

# get cytosines
cytosines = chain[chain.residue_name == '5FC']
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







atoms
atoms[atoms.residue_name == '5FC']

c5_no = 6081
c5_atom = atoms[atoms.atom_number == c5_no]
c5_atom
  
# create a relevant range of atom numbers to consider around this C5 atom
range_up = 12
range_down = 9
atom_range = range(c5_no - range_up, c5_no + range_down, 1)
atom_range

# extract those relevant atoms
within_range = atoms[atoms.atom_number.isin(atom_range)]
within_range

# extract atom_name and atom_number for each atom within the 5FC group, plus some before & after
# write it out in a dictionary
# will need to use this to write out CONECT lines
range_atom_nos = within_range.atom_number.tolist()
range_atom_nos

atom_positions = {}
for atom_no in range_atom_nos:
    atom_nm = within_range[within_range.atom_number == atom_no].atom_name.item()
    atom_positions.update({atom_nm : atom_no})
atom_positions

# add atoms before & after the 5FC group, also needed for CONECT lines
O3_before_atom_no = c5_no - 25
O3_before_atom_no
O3_before_atom_nm_check = atoms[atoms.atom_number == O3_before_atom_no].atom_name.item()
O3_before_atom_nm_check
O3_before_atom_nm = "O3' before"

P_after_atom_no = c5_no + 9
P_after_atom_no
P_after_atom_nm_check = atoms[atoms.atom_number == P_after_atom_no].atom_name.item()
P_after_atom_nm_check
P_after_atom_nm = "P after"

atom_positions.update({O3_before_atom_nm : O3_before_atom_no, P_after_atom_nm : P_after_atom_no})
atom_positions

# append info to a text file
with open('../Output_files/conect_info.txt', "a") as f:
    for k, v in atom_positions.items():
        f.write('\n' + str(k) + ','+ str(v))

# check if we can read in the new text file as a pandas dataframe
check_conect_info = pd.read_csv('../Output_files/conect_info.txt', encoding='utf-8')
check_conect_info



        
           
# dump dictionary items via json for potential use outside this script
import json 

dict_filename = '../Output_files/chain_' + chain_name + '_' + 'dict_c5_' + str(c5_no) + '_conect_info.csv'
with open(dict_filename, 'w') as file:
    file.write(json.dumps(positions_dict)) # use `json.loads` to do the reverse

      
# append info to a text file      
with open('../Output_files/connect_info.csv', "a") as f:
    for k, v in positions_dict.items():
        f.write('\n' + str(k) + ','+ str(v))
        
# write out individual text files for each 5FC group
filename = '../Output_files/c5_' + str(c5_no) + '_conect_info.csv'
with open(filename, "a") as file:
    for k, v in positions_dict.items():
        file.write('\n' + str(k) + ','+ str(v))

# check if we can read in the new text file as a pandas dataframe
check_conect_info = pd.read_csv('../Output_files/conect_info.txt', encoding='utf-8')
check_conect_info