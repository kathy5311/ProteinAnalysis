import math
import numpy as np
import os
from collections import OrderedDict

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

def H_func(ss_dict,Hdict,phi_dict,psi_dict):

    for L in Hdict:
        if L[0]+4==L[1]:
            for i in range(L[0],L[1]+1):
                if ((i not in phi_dict) and (-90<psi_dict[i]<30)):
                    ss_dict[i]="H"
                elif ((i not in psi_dict) and (-150<phi_dict[i]<-30)):
                    ss_dict[i]="H"
                elif ((-150<=phi_dict[i]<=-30) and (-90<=psi_dict[i]<=30)):
                    ss_dict[i]="H"
                else:
                    pass
    return ss_dict

def E_angle(ss_dict,atom_dict,phi_dict,psi_dict):
    
    for i in atom_dict:
        if (i in phi_dict and phi_dict[i]> -20) or (i in psi_dict and psi_dict[i] < 45):
            if not ss_dict[i]:
                ss_dict[i]="C"
    return (ss_dict)
                
def E_func(ss_dict,Hdict):
    for L in Hdict:
        i=L[0]
        j=L[1]    
        if [i-2, j-2] in Hdict:
            for x in range(i-2,i+1):
                if (ss_dict[x]=="H") or (ss_dict[x]=="C"):
                    pass
                elif x in ss_dict:
                    ss_dict[x]="E"
                
                
        elif [i+2, j-2] in Hdict:
            for x in range(i,i+3):
                if (ss_dict[x]=="H") or (ss_dict[x]=="C"):
                    pass
                elif x in ss_dict:
                    ss_dict[x]="E"
    
    for L in Hdict:
        i=L[1]
        j=L[0]    
        if [j-2, i-2] in Hdict:
            for x in range(i-2,i+1):
                if (ss_dict[x]=="H") or (ss_dict[x]=="C"):
                    pass
                elif x in ss_dict:
                    ss_dict[x]="E"
                
                
        elif [j+2, i-2] in Hdict:
            for x in range(i-2,i+1):
                if (ss_dict[x]=="H") or (ss_dict[x]=="C"):
                    pass
                elif x in ss_dict:
                    ss_dict[x]="E"
                    
    return ss_dict

def C_func(ss_dict, atom_dict):
    for i in atom_dict:
        if (ss_dict[i]=="H") or (ss_dict[i]=="E"):
            pass
        else:
            ss_dict[i]="C"
    
    return ss_dict

def check_C(ss_dict, atom_dict):
    for i in atom_dict:
         if ss_dict[i]=="E":
             if (i-1 in ss_dict and ss_dict[i-1]!="E") or (i+1 in ss_dict and ss_dict[i+1]!="E"):
                 if (i+3 in ss_dict and ss_dict[i+2]!="E") and (i-3 in ss_dict and ss_dict[i-2]!="E"):
                     ss_dict[i]="C"
    return ss_dict

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

def process_pdb(file_path):
    a=parsing_pdb(file_path)
    
    exception=[]
    for i in a:
        try:
            a[i]=amino_name(a[i])
        
        except Exception as e:
            print(e)
            exception.append(i)
        
    
    for key in exception:
        a[key]="-"
    
    if a=={}:
            os.remove(file_path)
            return 1
        
    data=open(file_path)
    atom_dict = {} #3중 딕셔너리임. chain에 해당하는 residue에 있는 atom에 대한 좌표값
    count=0
    error_dict={}
    error_list=[]
    new_atom={}
    seq_dict={}
    for line in data:#data 한 줄씩 가져오기
            
        if line[:4] == 'ATOM': #ATOM으로 시작하는 line 찾기
        
            chain = line[21] #chain identifer
            if chain not in atom_dict:
                count=0
                error_list.clear()
                atom_dict[chain] = {}
            
            if chain not in new_atom:
                new_atom[chain]={}
            
            if chain not in error_dict:
                error_dict[chain]={}
            
            if chain not in seq_dict:
                seq_dict[chain]={}
        
            resn=line[17:20].strip()
            atom_name = line[12:16] #원자 이름
            atom_name_split = ('').join(atom_name.split()) #원자이름 공백 제거
            resi = line[22:27].strip()
            real_index=line[22:27].strip()
                   
        
            if real_index not in new_atom[chain]:
                new_atom[chain][real_index]=0
            
            try:
                resi=int(resi)+count
            except:
                if int(resi[:-1])==1:
                    if resi not in error_list:
                        error_list.append(resi)
                        count+=1
                    resi=int(resi[:-1])
                else:
                    if resi not in error_list:
                        error_list.append(resi)
                        count+=1
                    resi=int(resi[:-1])+count
            
            error_dict[chain]=error_list
        
            if resi not in atom_dict[chain]:
                #atom_dict[chain]에 num 있는지 확인
                atom_dict[chain][resi] = {}
            atom_dict[chain][resi][atom_name_split] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        
            if resi not in seq_dict[chain]:
                seq_dict[chain][resi]=amino_name(resn)
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
    
    Hdict = {}
    for ch in atom_dict: #residue_dict에 [O,N] residue 저장
        if ch not in Hdict:
            Hdict[ch] = {}
            residue = []
        for num in atom_dict[ch]:
            for i in atom_dict[ch]:
                if num != i and 'O' in atom_dict[ch][num] and 'N' in atom_dict[ch][i]:
                    len = cal_len(atom_dict[ch][num]['O'], atom_dict[ch][i]['N'])
                    if 2.5 < len < 3.5 :
                        residue.append([num, i])
        Hdict[ch] = residue

    ss_dict={}
    for ch in atom_dict:
        
        ss_dict[ch]={}
    
        for i in atom_dict[ch]:
            ss_dict[ch][i]={}

        b=H_func(ss_dict[ch],Hdict[ch],phi_dict[ch],psi_dict[ch])
        ss_dict[ch]=b
        c=E_angle(ss_dict[ch],atom_dict[ch],phi_dict[ch],psi_dict[ch])
        ss_dict[ch]=c
        d=E_func(ss_dict[ch],Hdict[ch])
        ss_dict[ch]=d
        a=C_func(ss_dict[ch],atom_dict[ch])
        ss_dict[ch]=a
        
    for ch in new_atom:
        ss_list_index=0
        count_zero=0
        seq_list=[]
        ss_list=[]
        for i in ss_dict[ch]:
            ss_list.append(ss_dict[ch][i])
            count_zero+=1
        #print(ss_list)
        for k in seq_dict[ch]:
            seq_list.append(seq_dict[ch][k])
        #print(seq_list)
        #print()    
        for j in new_atom[ch]:
            if ss_list_index<count_zero:
    
                new_atom[ch][j]=[ss_list[ss_list_index]]
                new_atom[ch][j].append(seq_list[ss_list_index])
                ss_list_index+=1
            else:
                break

    data.close()
    return new_atom

def sequenc_func(new_atom):
    for ch in new_atom:
        new_string=""
        ordered_dict=OrderedDict(new_atom[ch])
        for key,value in ordered_dict.items():
            new_string+=value[1]
    
        print(ch)
        print(new_string)
        print()

def save_result(result, filename):
    npz_filename = filename[:-4] + ".npz"  # .pdb 확장자를 .npz 확장자로 바꿉니다.
    np.savez(npz_filename, result=result)
    print(f"Result saved as {npz_filename}")
       
#Execution part
exception_file=[]
error_list=[]
data_directory = "/Users/kistintern6/set2.antigen/" #pdb파일을 열어서 data에 저장
for filename in os.listdir(data_directory):
    if filename.endswith(".pdb"):
        file_path = os.path.join(data_directory, filename)
        try:
            
            final=process_pdb(file_path)
            save_result(final,filename)
            print(filename)
            print(final)
            sequenc_func(final)
            
        except Exception as e:
            exception_file.append(filename)
            if e not in error_list:
                error_list.append([e,filename])


