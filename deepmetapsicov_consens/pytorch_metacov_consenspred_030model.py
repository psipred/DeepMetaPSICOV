#!/usr/bin/env python

from __future__ import print_function

import sys
import os
#os.environ['CUDA_VISIBLE_DEVICES'] = "0"
import time
import random

from math import sqrt, log

import numpy as np

import torch
#from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F

from nndef_meta_resnet import ResNet

MAP_CHANNELS = 60

# ############################## Main program ################################
# Everything else will be handled in our main program now. We could pull out
# more functions to better separate the code, but it wouldn't make it any
# easier to read.

def main():

    device = torch.device("cpu")
    scriptdir = os.path.dirname(os.path.realpath(__file__))

    # Create neural network model (depending on first command line parameter)
    network = ResNet(64).to(device).eval()

    mapdata = np.fromfile(sys.argv[1], dtype=np.float32)
    length = int(sqrt(mapdata.shape[0]/441))
    inputs = mapdata.reshape(1,441,length,length)

    mapdata = np.fromfile(sys.argv[2], dtype=np.float32)

    inputs = torch.from_numpy(np.concatenate((mapdata.reshape(1,MAP_CHANNELS,length,length), inputs), axis=1)).type(torch.FloatTensor).to(device)

    pretrained_dict = torch.load(scriptdir + "/DMP_model1.pt", map_location=lambda storage, loc: storage)
    model_dict = network.state_dict()
    pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}

    network.load_state_dict(pretrained_dict)

    with torch.no_grad():
        result = torch.sigmoid(network(inputs).data)

        extra_models = ["DMP_model2.pt", "DMP_model3.pt", "DMP_model4.pt", "DMP_model5.pt"]

        for model in extra_models:
            pretrained_dict = torch.load(scriptdir + '/' + model, map_location=lambda storage, loc: storage)
            pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}
            network.load_state_dict(pretrained_dict)
            result += torch.sigmoid(network(inputs).data)

    n = len(extra_models)+1

    result /= n

    for wi in range(0, length-1):
        for wj in range(wi+2, length):
            print("{} {} 0 8 ".format(wi+1,wj+1), end='')
            print(" {}".format(0.5 * (result[0,0,wi,wj] + result[0,0,wj,wi])))

if __name__=="__main__":
    main()
