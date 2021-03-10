"""
Goal:
    Append CONECT lines to the end of the pdb file
Input:
    (1) pdb file with 5FC
    (2) txt file containing CONECT lines
Actions:
    (1) Delete the "END" line in pdb - write it out as a temp file
    (2) Append CONECT lines to the temp file, creating a new file
    (3) Add "END" line at the end of the new file
    (4) Delete the temp file
Output:
    (1) extended pdb file ("nuc_with_5FC_with_conect.pdb")
Notes:
    
"""

### set working directory to where your input pdb file is
import sys, os
os.chdir('C:\Kotryna\Python\Lab projects\pdb file handling\Input_Output_files')
cwd = os.getcwd()
#print("cwd is: ", cwd)

### add the generated CONECT lines to the bottom of nuc_with_FC file

# remove the "END" line from nuc_with_5FC; save the output in a different file

file1 = open("nuc_with_5FC.pdb")
lines1 = file1.readlines()
file1.close()

file2 = open("nuc_with_5FC_trunc.pdb",'w')
file2.writelines([item for item in lines1[:-1]])
file2.writelines("\n")
file2.close()


filenames = ['nuc_with_5FC_trunc.pdb', 'joint_conect_lines.txt']
with open('nuc_with_5FC_with_conect.pdb', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)


# remove the truncated pdb without "END" line
os.remove("nuc_with_5FC_trunc.pdb")

# add the "END" line at the end of the appended pdb
file4 = open("nuc_with_5FC_with_conect.pdb",'a')
file4.writelines("\nEND")
file4.close()
