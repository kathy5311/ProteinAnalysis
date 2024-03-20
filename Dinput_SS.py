import numpy as np
import math
import os
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

def cal_len(A, B): #길이 계산
	P = A - B
	len = np.sqrt(np.sum(P*P))
	return (len)

def H_func(ss_dict,residue,phi_dict,psi_dict):
    for L in residue:
        if L[0]+4==L[1]:
            for i in range(L[0],L[1]+1):
                if ((i not in phi_dict['A']) and (-90<psi_dict['A'][i]<30)):
                    ss_dict[i]="H"
                elif ((i not in psi_dict['A']) and (-150<phi_dict['A'][i]<-30)):
                    ss_dict[i]="H"
                elif ((-150<phi_dict['A'][i]<-30) and (-90<psi_dict['A'][i]<30)):
                    ss_dict[i]="H"
    return ss_dict

#beta-helix
#angle
def E_angle(ss_dict,atom_dict,phi_dict,psi_dict):
    for i in atom_dict['A']:
        if (i in phi_dict['A'] and phi_dict['A'][i]> -20) or (i in psi_dict['A'] and psi_dict['A'][i] < 45):
            if not ss_dict[i]:
                ss_dict[i]="C"
    return (ss_dict)
                
def E_func(ss_dict,residue):
    for L in residue:
        i=L[0]
        j=L[1]    
        if [i-2, j-2] in residue:
            for x in range(i-2,i+1):
                if (ss_dict[x]=="H") or (ss_dict[x]=="C"):
                    pass
                elif x in ss_dict:
                    ss_dict[x]="E"
                
                
        elif [i+2, j-2] in residue:
            for x in range(i,i+3):
                if (ss_dict[x]=="H") or (ss_dict[x]=="C"):
                    pass
                elif x in ss_dict:
                    ss_dict[x]="E"
    
    for L in residue:
        i=L[1]
        j=L[0]    
        if [j-2, i-2] in residue:
            for x in range(i-2,i+1):
                if (ss_dict[x]=="H") or (ss_dict[x]=="C"):
                    pass
                elif x in ss_dict:
                    ss_dict[x]="E"
                
                
        elif [j+2, i-2] in residue:
            for x in range(i-2,i+1):
                if (ss_dict[x]=="H") or (ss_dict[x]=="C"):
                    pass
                elif x in ss_dict:
                    ss_dict[x]="E"
                    
    return ss_dict

#coil
def C_func(ss_dict, atom_dict):
    for i in range(1,count+1):
        if (ss_dict[i]=="H") or (ss_dict[i]=="E"):
            pass
        else:
            ss_dict[i]="C"
    
    return ss_dict

def check_C(ss_dict, atom_dict):
    for i in range(1, (count+1)):
         if ss_dict[i]=="E":
             if (i-1 in ss_dict and ss_dict[i-1]!="E") or (i+1 in ss_dict and ss_dict[i+1]!="E"):
                 if (i+3 in ss_dict and ss_dict[i+2]!="E") and (i-3 in ss_dict and ss_dict[i-2]!="E"):
                     ss_dict[i]="C"
    return ss_dict

path="/Users/kistintern6/set2.antigen/"
file_list=os.listdir(path)
file_list_pdb=[file for file in file_list if file.endswith('.pdb')]

for file in file_list_pdb:
    data = open(path+file) 
    
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
    print(file)
    #print(atom_dict)
    print()

    count=0
    for k in atom_dict[chain]:
        count+=1



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
    residue_dict = {}
    for ch in atom_dict:
        if ch not in residue_dict:
            residue_dict[ch]={}
            residue=[]
        for num in atom_dict[ch]:
            for i in atom_dict[ch]:
                if num != i and '0' in atom_dict[ch][num] and 'N' in atom_dict[ch][i]:
                    len = cal_len(atom_dict[ch][num]['O'], atom_dict[ch][i]['N'])
                    if 2.5<len<3.5:
                        residue.append([num,i])
        residue_dict[ch]=residue

#print(len(residue_dict))

    ss_dict={}
    if 'A' not in ss_dict:
        ss_dict['A']={}



    for i in range(1, (count+1)):
            if i not in ss_dict['A']:
                ss_dict['A'][i]={}


    b=H_func(ss_dict[chain],residue,phi_dict,psi_dict)
    ss_dict[chain]=b
#print(ss_dict['A'])
#print()
    c=E_angle(ss_dict[chain],atom_dict,phi_dict,psi_dict)
    ss_dict[chain]=c

    d=E_func(ss_dict[chain],residue)
    ss_dict[chain]=d

    e=check_C(ss_dict[chain], atom_dict['A'])
    ss_dict[chain]=e
    a=C_func(ss_dict[chain],atom_dict[chain])
    ss_dict[chain]=a

    print(ss_dict[chain])
