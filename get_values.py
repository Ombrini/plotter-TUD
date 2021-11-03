
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



def get_Npart(resultDir,electrode):
    config = Config.from_dicts(resultDir)   
    Npart =  config["Npart"][electrode]
    return Npart

def get_Nvol(resultDir,electrode):
    config = Config.from_dicts(resultDir) 
    Nvol =  config["Nvol"][electrode]
    return Nvol


def get_C_rate(resultDir):
    #return single value for CC and an array for CCsegments
    config = Config.from_dicts(resultDir)

    if config["profileType"] == "CC":
        Crate = config["Crate"]
        return Crate
    if config["profileType"] == "CCsegments":
        Crate_vec = np.empty([0])
        segments = config["segments"]
        for i in range(len(segments)):
            Crate_i = np.round(segments[i][0]/0.00216464,decimals=1)
            # I don't know why using config on "segments" gives me a decimal number
            # that must be devided by 0.0021 to make it equal to the Crate
            Crate_vec = np.append(Crate_vec,Crate_i)
        return Crate_vec

# def get_times(resultDir):
    # the idea is the have the possibility to zoom in the 
    # times or in the ffrac and plot only those parts

    # for example we set different segments and we want to plot 
    # c_barLine only in during one or more of these segments 
    # using different rates could help us to see when the solid 
    # electrolyte creates diffusion limitations and when instead 
    # the cell is reaction limited



