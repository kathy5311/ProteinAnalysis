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
        conv1 = torch.nn.Conv1d(20,32,3,padding=1) # aa1hot,channel,
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
    def __init__(self, datalist, idx):
        self.tags = datalist 

    def __len__(self):
        return len(self.tags)

    def __getitem__(self,index):
        npz = 'data/'+self.tags[index]+'.npz'
        data = np.load(npz,allow_pickle=True)

        aas = 'ACDEFGHIKLMNPQRSTVWYX'
        SS3 = 'HEC'
        
        seqs = [AAs.index(a) for a in  data['seqs']]
        SSs  = [SS3.index(a) for a in data['SSs']]
        
        seq1hot = np.transpose(np.eye(21)[seqs],(1,0)) # 20xnres
        SS1hot = np.transpose(np.eye(3)[SSs],(1,0))
        
        return seq1hot, self.SSs[index], seq1hot.shape[1]
    
def collate(samples):
    seq,SS,nres = map(list, zip(*samples))
    nres = max(nres)
    B = len(seq)

    # map into maxres
    seqs = torch.zeros(B,20,nres)
    SSs  = torch.zeros(B,nres,dtype=torch.long)
    for i,s in enumerate(seq): seqs[i][:len(s[1])] = torch.tensor(s) 
    for i,s in enumerate(SS):  SSs[i][:len(s)] = torch.tensor(s)
    
    return seqs, SSs

model = CNN()
model.to(device)

## load dataset
trainlist = np.load('data/train.npy')
validlist = np.load('data/valid.npy')

trainset = DataSet(trainlist)
validset = DataSet(validlist)

generator_params = {
    'shuffle': True,
    'num_workers': 4,
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
        # get prediction
        SSpred = model(seq.to(device))
        # calculate loss
        SS = SS.to(device)
        loss = lossfunc(SSpred,SS)
        loss_v.append(loss.cpu().detach().numpy())
        
    print("Train/Valid: %3d %8.4f %8.4f"%(epoch, float(np.mean(loss_t)), float(np.mean(loss_v))))