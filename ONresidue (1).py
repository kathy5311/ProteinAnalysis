#day4. O,N residue
import numpy as np
from atomxyz import atom_dict

def cal_len(A, B): #길이 계산
	P = A - B
	len = np.sqrt(np.sum(P*P))
	return (len)

residue_dict = {}
for ch in atom_dict: #residue_dict에 [O,N] residue 저장
	if ch not in residue_dict:
		residue_dict[ch] = {}
		residue = []
	for num in atom_dict[ch]:
		for i in atom_dict[ch]:
			if num != i and 'O' in atom_dict[ch][num] and 'N' in atom_dict[ch][i]:
				len = cal_len(atom_dict[ch][num]['O'], atom_dict[ch][i]['N'])
				if 2.5 < len < 3.5 :
					residue.append([num, i])
	residue_dict[ch] = residue

