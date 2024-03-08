#day2-3. phi/psi 계산

import numpy as np
import math

from atomxyz import atom_dict

def return_angle(A, B, C, D):
    AB = B - A
    AC = C - A
    BC = C - B
    BD = D - B
    v1 = np.cross(AB, AC) #np.cross : 벡터의 외적
    v2 = np.cross(BC, BD)
    v1v2 = np.dot(v1, v2) #np.dot : 벡터 점의 곱
    len_v1 = math.sqrt(np.sum(v1*v1)) #np.sqrt : 제곱근, np.sum : 합
    len_v2 = math.sqrt(np.sum(v2*v2))
    radi = math.acos(v1v2 / (len_v1*len_v2)) #math.acos : radian구하기
    x = math.degrees(radi) #math.degrees : 라디안을 각도로 변환
    k = np.dot(np.cross(v1,v2), BC)
    sign = 1
    if k < 0:
        sign = -1
    return (sign*x)

def phi_cal_angle(atom_dict, num): #좌표찾고 phi 계산
    if 'C' in atom_dict[num - 1] and 'N' in atom_dict[num] and 'CA' in atom_dict[num] and 'C' in atom_dict[num]:
        phi = return_angle(atom_dict[num - 1]['C'], atom_dict[num]['N'], atom_dict[num]['CA'], atom_dict[num]['C'])
    else:
        phi = None
    return (phi)

def psi_cal_angle(atom_dict, num): #각각의 좌표 찾아서 psi계산
    if 'N' in atom_dict[num] and 'CA' in atom_dict[num] and 'C' in atom_dict[num] and 'N' in atom_dict[num + 1]:
        psi = return_angle(atom_dict[num]['N'], atom_dict[num]['CA'], atom_dict[num]['C'], atom_dict[num + 1]['N'])
    else:
        psi = None
    return(psi)

phi_dict = {}
psi_dict = {}
for ch in atom_dict:
    if ch not in phi_dict:
        phi_dict[ch] = {}
    if ch not in psi_dict:
        psi_dict[ch] = {}
    for num in atom_dict[ch]:
        if num - 1 in atom_dict[ch]:
            phi = phi_cal_angle(atom_dict[ch], num)
            if phi is not None:
                phi_dict[ch][num] = phi
        if num + 1 in atom_dict[ch]:
            psi = psi_cal_angle(atom_dict[ch], num)
            if psi is not None:
                psi_dict[ch][num] = psi_cal_angle(atom_dict[ch], num)

'''
print("phi---------------------------------------------")
print(phi_dict)
print("psi---------------------------------------------")
print(psi_dict)
print("------------------------------------------------")
'''
