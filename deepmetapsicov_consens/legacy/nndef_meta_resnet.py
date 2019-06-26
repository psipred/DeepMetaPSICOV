import torch.nn as nn
import torch.nn.functional as F
from math import sqrt

NUM_CHANNELS = 441+60

class Maxout2d(nn.Module):

    def __init__(self, in_channels, out_channels, pool_size):
        super(Maxout2d, self).__init__()
        self.in_channels, self.out_channels, self.pool_size = in_channels, out_channels, pool_size
        self.lin = nn.Conv2d(in_channels=in_channels, out_channels=out_channels * pool_size, kernel_size=1)

    def forward(self, inputs):
        shape = list(inputs.size())
        out = self.lin(inputs)
        m = out.view(shape[0], self.out_channels, self.pool_size, shape[2], shape[3]).max(dim=2)[0]
        return m

# ResNet Module
class ResNet(nn.Module):
    def __init__(self,width):
        super(ResNet, self).__init__()
        self.redim = Maxout2d(in_channels=NUM_CHANNELS, out_channels=width, pool_size=3)
        self.firstnorm = nn.InstanceNorm2d(width, affine=True)
        self.resblocks = nn.ModuleList()

        for fsize,dilv in [(5,1), (5,2), (5,1), (5,4), (5,1), (5,8), (5,1), (5,16), (5,1), (5,32), (5,1), (5,64), (5,1), (5,1), (5,1), (5,1), (5,1), (5,1)]:
             if fsize > 0:
                layer = nn.Conv2d(in_channels=width, out_channels=width, kernel_size=5, dilation=dilv, padding=int(dilv*(fsize-1)/2))
                nn.init.xavier_uniform(layer.weight, gain=sqrt(2.0))
                self.resblocks.append(layer)
                self.resblocks.append(nn.InstanceNorm2d(width, affine=True))
                layer = nn.Conv2d(in_channels=width, out_channels=width, kernel_size=5, dilation=dilv, padding=int(dilv*(fsize-1)/2))
                nn.init.xavier_uniform(layer.weight, gain=sqrt(2.0))
                self.resblocks.append(layer)
                self.resblocks.append(nn.InstanceNorm2d(width, affine=True))

        self.lastnorm = nn.InstanceNorm2d(1, affine=True)
        self.outlayer = nn.Conv2d(in_channels=width, out_channels=1, kernel_size=1)
        nn.init.xavier_uniform(self.outlayer.weight)

    def forward(self, x):
        out = self.redim(x)
        out = self.firstnorm(out)
        out = F.dropout(out, p=0.2)
        out = F.dropout2d(out, p=0.2)
        for i in range(int(len(self.resblocks)/4)):
            residual = out
            out = self.resblocks[i*4](out)
            out = F.relu(self.resblocks[i*4+1](out))
            out = self.resblocks[i*4+2](out)
            out = self.resblocks[i*4+3](out)
            out += residual
            out = F.relu(out)
        out = self.outlayer(out)
        out = self.lastnorm(out)
        return out
