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
    
'''
def pdb_open(pdb_file):
    amino_list=[]
    with open(pdb_file,'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                pdb_dict={
                    'amino_name':line[17:20].strip()
                }
                amino_list.append(pdb_dict['amino_name'])
    return amino_list

pdb_file='/Users/kistintern6/ProteinAnalysis/Day1/PDBex.pdb'
count=[]
b=pdb_open(pdb_file)
print(b)
print()
for name in b:
    k=amino_name(name)
    if k not in count:
        count.append(k)
    
    print(k, end="")
print(len(count))
'''

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
                    xyz[atom_info['Residue_num']]={
                        'atoms':"A"
                    }
                
                    
                if atom_info['Residue_Name'] not in xyz[atom_info['Residue_num']]['atoms']:
                        xyz[atom_info['Residue_num']]['atoms']=atom_info['Residue_Name']
                
                
    return xyz

pdb_file='/Users/kistintern6/set2.antigen/3ijy_ag.pdb'

b=parsing_pdb(pdb_file)
print(b)
print()

for i in range(1,len(b)+1):
    b[i]['atoms']=amino_name(b[i]['atoms'])
print(b)


    
    
    