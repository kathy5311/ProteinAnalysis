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

pdb_file='/Users/kistintern6/Desktop/PDBex.pdb'
b=pdb_open(pdb_file)
for name in b:
    k=amino_name(name)
    print(k, end="")

            

    
    
    