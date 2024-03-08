#day1. atom 좌표 저장

import numpy as np

data = open("/Users/kistintern6/ProteinAnalysis/Day1/PDBex.pdb") #pdb파일을 열어서 data에 저장
atom_dict = {} #3중 딕셔너리임. chain에 해당하는 residue에 있는 atom에 대한 좌표값

for line in data: #data 한 줄씩 가져오기
    if line[:4] == 'ATOM': #ATOM으로 시작하는 line 찾기
        chain = line[21] #chain identifer
        if chain not in atom_dict:
            atom_dict[chain] = {}
        atom_name = line[12:16] #원자 이름
        atom_name_split = ('').join(atom_name.split()) #원자이름 공백 제거
        num = int(line[22:26])
        if num not in atom_dict[chain]: #atom_dict[chain]에 num 있는지 확인
            atom_dict[chain][num] = {}
        atom_dict[chain][num][atom_name_split] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            #원자 이름에 해당하는 좌표값 저장
data.close()
