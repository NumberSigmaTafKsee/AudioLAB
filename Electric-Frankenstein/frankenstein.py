from torch.utils.data import Dataset,DataLoader
import numpy as np
import torch
import torch.nn as nn
from torch.nn import Module, Parameter
from torch import FloatTensor
from scipy import signal
import math
import matplotlib as mpl
import matplotlib.pyplot as plt

class DIIRDataSet(Dataset):
    def __init__(self, input, target, sequence_length):
        self.input = input
        self.target = target
        self._sequence_length = sequence_length
        self.input_sequence = self.wrap_to_sequences(self.input, self._sequence_length)
        self.target_sequence = self.wrap_to_sequences(self.target, self._sequence_length)
        self._len = self.input_sequence.shape[0]

    def __len__(self):
        return self._len

    def __getitem__(self, index):
        return {'input': self.input_sequence[index, :, :]
               ,'target': self.target_sequence[index, :, :]}

    def wrap_to_sequences(self, data, sequence_length):
        num_sequences = int(np.floor(data.shape[0] / sequence_length))
        truncated_data = data[0:(num_sequences * sequence_length)]
        wrapped_data = truncated_data.reshape((num_sequences, sequence_length, 1))
        return np.float32(wrapped_data)
        
class DOnePoleCell(Module):
    def __init__(self, a1=0.05, b0=1.0, b1=0.0):
        super(DOnePoleCell, self).__init__()
        self.b0 = Parameter(FloatTensor([b0]))
        self.b1 = Parameter(FloatTensor([b1]))
        self.a1 = Parameter(FloatTensor([a1]))

    def init_states(self, size):
        state = torch.zeros(size).to(self.a1.device)
        return state

    def forward(self, input, state):
        self.a1.data = self.a1.clamp(-1, 1)
        output = self.b0 * input + state
        state = self.b1 * input + self.a1 * output
        return output, state

class DOnePole(Module):
    def __init__(self):
        super(DOnePole, self).__init__()
        self.cell = DOnePoleCell()

    def forward(self, input, initial_states=None):
        batch_size = input.shape[0]
        sequence_length = input.shape[1]

        if initial_states is None:
            states = self.cell.init_states(batch_size)
        else:
            states = initial_states

        out_sequence = torch.zeros(input.shape[:-1]).to(input.device)
        for s_idx in range(sequence_length):
            out_sequence[:, s_idx], states = self.cell(input[:, s_idx].view(-1), states)
        out_sequence = out_sequence.unsqueeze(-1)

        if initial_states is None:
            return out_sequence
        else:
            return out_sequence, states        
            
            
fs = 48e3
f0 = 20
f1 = 20e3
t = np.linspace(0, 60, math.ceil(60*fs))

train_input = signal.chirp(t=t, f0=f0, t1=60, f1=f1, method='logarithmic') + np.random.normal(scale=5e-2, size=len(t))

fc = 2e3
sos = signal.butter(N=2, Wn=fc/fs, output='sos')
train_target = signal.sosfilt(sos, train_input)

device = torch.device('cuda')
model = DOnePole().to(device)
batch_size = 1024
sequence_length = 512

loader = DataLoader(dataset=DIIRDataSet(train_input, train_target, sequence_length), batch_size=batch_size, shuffle = True)
n_epochs = 10
lr = 1e-3

optimizer = torch.optim.Adam(model.parameters(), lr=lr, betas=(0.9, 0.999), eps=1e-08, weight_decay=0, amsgrad=False)

criterion = nn.MSELoss()


def train(criterion, model, loader, optimizer):
    model.train()
    device = next(model.parameters()).device
    total_loss = 0
    for batch in loader:
        input_seq_batch = batch['input'].to(device)
        target_seq_batch = batch['target'].to(device)        
        optimizer.zero_grad()
        predicted_output = model(input_seq_batch)
        loss = criterion(target_seq_batch, predicted_output)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()

    total_loss /= len(loader)
    return total_loss
    
for epoch in range(n_epochs):
    loss = train(criterion, model, loader, optimizer)
    print("Epoch {} -- Loss {:3E}".format(epoch, loss))


plt.plot(train_input)
plt.plot(train_target)
train_target = signal.sosfilt(sos, train_input)
print(model.cell.a1,model.cell.b0,model.cell.b1)
print(model.cell.a1[0].item())

