
import os.path as osp
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

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
        limtrode = config['limtrode']
        theoretical_1C_current = config[limtrode, 'cap'] / 3600
        for i in range(len(segments)):
            Crate_i = np.round(segments[i][0]*config['t_ref']/60)
            Crate_i = np.round(config["segments"][i][0] *config["1C_current_density"] /theoretical_1C_current /config['curr_ref']/1.8e-5, 1)
            # I don't know why using config on "segments" gives me a decimal number
            # that must be devided by 0.0021 to make it equal to the Crate
            Crate_vec = np.append(Crate_vec,Crate_i)
        return Crate_vec


def keys(resultDir):
    config = Config.from_dicts(resultDir)
    matfile = osp.join(resultDir, 'output_data.mat')
    sim_output = sio.loadmat(matfile)   
    for k in sim_output.keys():
         print(k)
    return

def get_last_V(resultDir):
    k = constants.k                      # Boltzmann constant, J/(K Li)
    Tref = constants.T_ref               # Temp, K
    e = constants.e   

    config = Config.from_dicts(resultDir)
    matfile = osp.join(resultDir, 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    trodes = config["trodes"]
    td = config["t_ref"]

    Etheta = {"a": 0.}
    for trode in trodes:
        Etheta[trode] = -(k*Tref/e) * config[trode, "phiRef"]
    data_volt = Etheta["c"] - Etheta["a"] - (k*Tref/e)*sim_output['phi_cell'][0]
    lastV = data_volt[-1]

    return lastV

def get_bulkCon(resultDir,electrode):
    config = Config.from_dicts(resultDir)
    matfile = osp.join(resultDir, 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    bulkCon = config['sigma_s'][electrode]*config['sigma_s_ref']

    return bulkCon

def get_thickness(resultDir,electrode):
    config = Config.from_dicts(resultDir)
    matfile = osp.join(resultDir, 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    L_c = config["L"][electrode]*config['L_ref']*1e6

    return L_c




