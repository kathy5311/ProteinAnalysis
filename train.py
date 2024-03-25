from sklearn.metrics import confusion_matrix
import torch
device = torch.device("cuda:0" if (torch.cuda.is_available()) else "cpu")
#import models
import numpy as np

#model
import torch.nn as nn

MAXEPOCH=100
BATCH=1

class CNN(nn.Module):
    def __init__(self, nlayer=4, dropout=0.1):
        super().__init__()
        layers = []

        drop = torch.nn.Dropout(p=dropout)
        conv1 = torch.nn.Conv1d(21,32,3,padding=1) # aa1hot,channel,
        layers = [drop,conv1]
        
        for k in range(nlayer):
            conv2 = torch.nn.Conv1d(32,32,3,padding=1) # aa1hot,channel,
            layers.append(conv2)
            layers.append(nn.BatchNorm1d(32))
            layers.append(nn.ReLU(inplace=True))

        self.layers = nn.ModuleList(layers)

        # 1 x 32 x nres
        self.outlayer = nn.Linear(32,3)
        
    def forward(self, seq):
        #pred = seq # should B x 20 x nres
        for layer in self.layers:
            seq = layer(seq)

        seq = torch.transpose(seq,1,2) # put channel at the last
        
        pred = self.outlayer(seq)
        pred = torch.transpose(pred,2,1)
        return pred

class DataSet(torch.utils.data.Dataset):
    def __init__(self, datalist): #idx제거
        self.tags = [tag for tag in datalist if not tag.endswith('.DS_Store')]

    def __len__(self):
        return len(self.tags)

    def __getitem__(self,index):
        npz = self.tags[index]
        #print(npz)
          
        data = np.load(npz,allow_pickle=True)

        aas = 'ACDEFGHIKLMNPQRSTVWYX'
        SS3 = 'HEC'
    
        seqs = [aas.index(a) for a in  data['sequence']]#변수 명 바꿈
        #print(seqs)
        SSs  = [SS3.index(a) for a in data['SS']]
        
        seq1hot = np.transpose(np.eye(21)[seqs],(1,0)) # 20xnres #tensor size 20으로 바꿈.
        #np.eye(21)[seqs]=> 21x21항등행렬을 만드는데 seqs에 맞게 설정된 항등행렬을 만들어라.
        
        SS1hot = np.transpose(np.eye(3)[SSs],(1,0)) #np.eye는 항등행렬을 생성해줌
        #print(SS1hot)
        return seq1hot, SSs, seq1hot.shape[1] #SSs[index]??->SS1hot
        #seq1hot.shape[1]=시퀀스의 길이다.
    
def collate(samples): # 같은 배치 안에 길이가 가장 긴 input에 맞춰 다른 input들에 임의로 zero-padding.
    seq,SS,nres = map(list, zip(*samples))
    valid = [i for i,n in enumerate(nres) if n > 50]
    #print(valid)
    if len(valid) == 0: return [],[]
    
    seq = [seq[i] for i in valid]
    SS = [SS[i] for i in valid]
    
    nres = max(nres)
    B = len(seq)

    # map into maxres
    seqs = torch.zeros(B,21,nres)
    #print(seqs)
    SSs  = torch.zeros(B,nres,dtype=torch.long)
    for i,s in enumerate(seq): #seqs는 0한 개로만 이루어진 리스트다. 때문에 i는 0만 출력된다.
        #print(i)
        seqs[i][:len(s[1])] = torch.tensor(s) 
    for i,s in enumerate(SS): 
        #print(i)
        SSs[i][:len(s)] = torch.tensor(s)

    return seqs, SSs

model = CNN()
model.to(device)

## load dataset
trainlist = np.load('train.npy')
validlist = np.load('valid.npy')

trainset = DataSet(trainlist)
validset = DataSet(validlist)

generator_params = {
    'shuffle': True,
    'num_workers': 0, #num_worker 오류가 나서 0으로 바꿔줌.
    'pin_memory': True,
    'collate_fn': collate,
    'batch_size': BATCH,
    'worker_init_fn' : np.random.seed()
}
train_generator = torch.utils.data.DataLoader(trainset, **generator_params)
valid_generator = torch.utils.data.DataLoader(validset, **generator_params)
 
optimizer = torch.optim.Adam(model.parameters(), lr=1.0e-4)

lossfunc = torch.nn.CrossEntropyLoss()
for epoch in range(MAXEPOCH):  
    loss_t = []
    for i,(seq,SS) in enumerate(train_generator):
        if len(seq) == 0: continue
        
        # get prediction
        SSpred = model(seq.to(device))
        
        # calculate loss
        SS = SS.to(device)
        loss = lossfunc(SSpred,SS)
        loss.backward(retain_graph=True)
        optimizer.step()

        loss_t.append(loss.cpu().detach().numpy())
    #print("TRAIN:", epoch, float(np.mean(loss_)))
        
    loss_v = []
    for i,(seq,SS) in enumerate(valid_generator):
        if len(seq) == 0: continue        
        # get prediction
        SSpred = model(seq.to(device))
        SS_valid=SSpred
        
        # calculate loss
        SS = SS.to(device)
        loss = lossfunc(SSpred,SS)
        loss_v.append(loss.cpu().detach().numpy())
        
    print("Train/Valid: %3d %8.4f %8.4f"%(epoch, float(np.mean(loss_t)), float(np.mean(loss_v))))

'''
correct_train=0
total_train=0

for i, (seq, SS) in enumerate(valid_generator):
    print(i)
    print(seq,SS)
    
    if len(seq)==0:
        continue
    
    SSpred=model(seq.to(device))
    
    _, predicted = torch.max(SSpred.data,1)
    
    correct_train += (predicted == SS.to(device)).sum().item()
    total_train+=seq.size(0)

accuracy_train=correct_train/total_train

correct_valid=0
total_valid=0

for i, (seq, SS) in enumerate(valid_generator):
    if len(seq)==0:
        continue
    
    SSpred=model(seq.to(device))
    
    _, predicted = torch.max(SSpred.data,1)
    
    correct_valid += (predicted == SS.to(device)).sum().item()
    total_valid+=seq.size(0)

accuracy_valid=correct_valid/total_valid

# 정확도 출력
print("Train Accuracy: {:.2f}%".format(accuracy_train* 100))
print("Valid Accuracy: {:.2f}%".format(accuracy_valid * 100))
'''

for i, (seq, SS) in enumerate(train_generator):
    if len(seq)==0:
        continue
    
    SSpred=model(seq.to(device))
    
    _, predicted = torch.max(SSpred.data,1)
    
    cm=confusion_matrix(SS.to(device).detach().numpy(), SSpred.detach().numpy())
    print(cm)
