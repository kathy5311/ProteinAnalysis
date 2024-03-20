import os
#from Combination_Assign import *
#from protein_sequence import parsing_pdb, amino_name

#requirement
#1.protein_sequence
#2. secondary structure labels

#testing one file
#protein_sequence
amino_dict={
        'ALA':'A',
        'ARG':'R',
        'ASN':'N',
        'ASP':'D',
        'CYS':'C',
        'GLU':'E',
        'GLN':'Q',
        'GLY':'G',
        'HIS':'H',
        'ILE':'I',
        'LEU':'L',
        'LYS':'K',
        'MET':'M',
        'PHE':'F',
        'PRO':'P',
        'SER':'S',
        'THR':'T',
        'TRP':'W',
        'TYR':'Y',
        'VAL':'V'
        
    }

def amino_name(a):
    amino_dict={
        'ALA':'A',
        'ARG':'R',
        'ASN':'N',
        'ASP':'D',
        'CYS':'C',
        'GLU':'E',
        'GLN':'Q',
        'GLY':'G',
        'HIS':'H',
        'ILE':'I',
        'LEU':'L',
        'LYS':'K',
        'MET':'M',
        'PHE':'F',
        'PRO':'P',
        'SER':'S',
        'THR':'T',
        'TRP':'W',
        'TYR':'Y',
        'VAL':'V'
        
    }
    return amino_dict[a]

def parsing_pdb(pdb_file):
    with open(pdb_file, 'r') as file:
        xyz={}
        for line in file:
            if line.startswith("ATOM"):
                atom_info={
                    #'Serial_Num':line[6:11],
                    #'Atom_Name':line[12:16].strip(),
                    #'Alt_location':line[16:17],
                    'Residue_Name':line[17:20].strip(),
                    'Residue_num':int(line[22:26].strip()),
                    #'xyz':line[30:54].strip()
                }

                if atom_info['Residue_num'] not in xyz:
                    xyz[atom_info['Residue_num']]=atom_info['Residue_Name']
                
                
    return xyz


path="/Users/kistintern6/set2.antigen/"
file_list=os.listdir(path)
file_list_pdb=[file for file in file_list if file.endswith('.pdb')]
#print(file_list_pdb)
exception_file=[]
for file in file_list_pdb:
    path_file=path+file
    a=parsing_pdb(path_file)
    
    exception=[]
    for i in a:
        if (a[i] in amino_dict.keys()):
            a[i]=amino_name(a[i])
        else:
            if file not in exception_file:
                exception_file.append(file)
            exception.append(i)
    
    for key in exception:
        del a[key]
    
    if a=={}:
        os.remove(path_file)
    
    print(file)    
    print(a)
    print(exception)
    print()
print(exception_file)



