import numpy as np
import os
import random

path="/Users/kistintern6/ProteinAnalysis/Day4/npz/"
file_list=os.listdir(path)

random.shuffle(file_list)

train_size=int(0.8*len(file_list))
train=np.array(file_list[:train_size])
test = np.array(file_list[train_size:])

#print(train)
#np.save('train.npy',train)
#np.save('valid.npy',test)
#print()
#print(test)
data=np.load('train.npy')
print(data)

        