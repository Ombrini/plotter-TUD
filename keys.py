import os.path as osp
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import os

import mpet.mod_cell as mod_cell
import mpet.geometry as geom
import mpet.props_am as props_am
import mpet.utils as utils

from mpet.config import Config, constants

def keys(resultDir):
    config = Config.from_dicts(resultDir)
    matfile = osp.join(resultDir, 'output_data.mat')
    sim_output = sio.loadmat(matfile)   
    for k in sim_output.keys():
         print(k)
    return