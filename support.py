
pairs=[['C', 'C2'], ['C', 'C14'], ['C2', 'C3'], ['C3', 'C4'], ['C4', 'C13'], ['C7', 'C8'], ['C7', 'C12'], ['C8', 'C9'], ['C9', 'C10'], ['C10', 'C11'], ['C11', 'C12'], ['C13', 'C14']]

#save=[]
for i in range(len(pairs)):
    print(i)
    save=[]
    save.append(pairs[i][0])
    save.append(pairs[i][1])
    print(save)
    for j in pairs:
        if i!=j:
            atom1=i[0]
            atom2=i[1]
            
            if (atom1 in j) or (atom2 in j):
                if j[0] not in save:
                    save.append(j[0])
                if j[1] not in save:
                    save.append(j[1])
#print(save)
print()
            
                     
    
'''
#print(pairs[0])
for atoms in pairs[0]:
    save.append(atoms)
print(save)

for pair in pairs:
    atom1=pair[0]  
    atom2=pair[1]
    if atom1 in save or atom2 in save:
        save.append(atom1)
        save.append(atom2)
        save=list(set(save))
    elif atom1 not in save and atom2 not in save:
        new1=[]
        new1.append(atom1)
        new1.append(atom2)
        new1=list(set(new1))

#
print(save)
print(new1)
'''