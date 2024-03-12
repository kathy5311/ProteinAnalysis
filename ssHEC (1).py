#day5. SS에 H, E, C 채우기
 
from atomxyz import atom_dict
from ONresidue import residue_dict
from phipsi import psi_dict
from phipsi import phi_dict

#H, E 둘다 만약 양 끝일 때는 각도 판별 안하는거 이거 안함....


def fill_H(ss_dict, residue, phi_dict, psi_dict):
	for A in residue: #H 넣기 오른손감김 알파 헬릭스
		if A[1] == A[0] + 4:
			for i in range(A[0], A[1] + 1):
				if not ss_dict[i]:
					if i not in phi_dict and (i in psi_dict and -90 < psi_dict[i] < 30):
						ss_dict[i] = 'H'
					elif i not in psi_dict and (i in phi_dict and -150 < phi_dict[i] < -30):
						ss_dict[i] = 'H'
					elif (i in phi_dict and -150 < phi_dict[i] < -30) and (i in psi_dict and -90 < psi_dict[i] < 30):
						ss_dict[i] = 'H'
	return (ss_dict)

def check_phipsi_E(ss_dict, atom_dict, phi_dict, psi_dict):
	for i in atom_dict: # E 넣기 전에 각도로 먼저 판단하여 각도 만족 안하면 C 넣기
		if (i in phi_dict and phi_dict[i] > -20) or (i in psi_dict and psi_dict[i] < 45):
			print(i, '_________')
			if not ss_dict[i]:
				ss_dict[i] = 'C'
				print(i,'C')
	return (ss_dict)

def fill_E(ss_dict, residue):
	for list in residue: #E 넣기
		i = list[0] #i가 O, j가 N일때
		j = list[1]
		if [i - 2, j - 2] in residue: #i - 2, j - 2 있을때
			if not ((i < j and i >= j - 2) or (i > j and i - 2 <= j)):
				for x in range(i, i - 3, -1):
					if x in ss_dict and not ss_dict[x]:
						ss_dict[x] = 'E'
						print(i,'E', 1)
		if [i + 2, j - 2] in residue: #i + 2, j - 2 있을 때
			if not (i < j and i + 2 >= j - 2):
				for x in range(i, i + 3):
					if x in ss_dict and not ss_dict[x]:
						ss_dict[x] = 'E'
						print(i, 'E', 3)
	for list in residue:
		i = list[1] #i가 N, j가 O일때
		j = list[0]
		if [j - 2, i - 2] in residue:
			if not ((i < j and i >= j - 2) or (i > j and i - 2 <= j)):
				for x in range(i, i - 3, -1):
					if x in ss_dict and not ss_dict[x]:
						ss_dict[x] = 'E'
						print(i, 'E', 5)
		if [j - 2, i + 2] in residue:
			if not (i < j and i + 2 >= j - 2):
				for x in range(i, i + 3):
					if x in ss_dict and not ss_dict[x]:
						ss_dict[x] = 'E'
						print(i, 'E', 7)
	return (ss_dict)

def fill_C(ss_dict, atom_dict):
	for i in atom_dict: #안채워진 부분은 C로 채우기 + E가 2개 미만일 때 C로 바꿈
		if not ss_dict[i]:
			ss_dict[i] = 'C'
		if ss_dict[i] == 'E': #E가 1개일 때 안되므로 C로 바꿔줌
			if i - 1 not in ss_dict and i + 1 in ss_dict and ss_dict[i + 1] != 'E': #i - 1에 residue가 없고, i + 1이 E가 아닐 때
				ss_dict[i] = 'C'
			elif i + 1 not in ss_dict and i - 1 in ss_dict and ss_dict[i - 1] != 'E': #i + 1에 residue가 없고, i - 1이 E가 아닐 때
				ss_dict[i] = 'C'
			elif (i + 1 in ss_dict and ss_dict[i + 1] != 'E') and (i - 1 in ss_dict and ss_dict[i - 1] != 'E'): #양 옆이 'E'가 아닐 때
				ss_dict[i] = 'C'
	return (ss_dict)

ss_dict = {}
for ch in atom_dict: # residue 수 만큼 칸 할당하고
	ss_dict[ch] = {}
	for i in atom_dict[ch]:
		ss_dict[ch][i] = {}
	ss_dict[ch] = fill_H(ss_dict[ch], residue_dict[ch], phi_dict[ch], psi_dict[ch])
	ss_dict[ch] = check_phipsi_E(ss_dict[ch], atom_dict[ch], phi_dict[ch], psi_dict[ch])
	ss_dict[ch] = fill_E(ss_dict[ch], residue_dict[ch])
	ss_dict[ch] = dict(sorted(ss_dict[ch].items())) #ss_dict[ch]을 key 값인 residue 순으로 정렬
	ss_dict[ch] = fill_C(ss_dict[ch], atom_dict[ch])

print(ss_dict)

