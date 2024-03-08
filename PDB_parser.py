#pdb 딕셔너리
import numpy as np
def parsing_pdb(pdb_file):
    with open(pdb_file, 'r') as file:
        xyz={}
        for line in file:
            if line.startswith("ATOM"):
                chain=line[21]
                if chain not in xyz:
                    xyz[chain]={}
                
                residue_num=line[22:26].strip()
                if residue_num not in xyz[chain]:
                    xyz[chain][residue_num]={}
                
                atom_name=line[12:16].strip()
                if atom_name not in xyz[chain][residue_num]:
                    xyz[chain][residue_num][atom_name]={}
                
                xyz[chain][residue_num][atom_name]=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
                
    return xyz

pdb_file='/Users/kistintern6/ProteinAnalysis/Day1/PDBex.pdb'

#print(parsing_pdb(pdb_file))
parsing_dict=parsing_pdb(pdb_file)
print(parsing_dict)

    