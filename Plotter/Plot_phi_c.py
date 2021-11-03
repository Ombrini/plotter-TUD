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

# Why are those always 0 ? bulk should be the e- potential wich can be always 0, what is part ? 
def phi_part_c(resultDir_dic):
    k = constants.k                      # Boltzmann constant, J/(K Li)
    Tref = constants.T_ref               # Temp, K
    e = constants.e                      # Charge of proton, C
    # taking the mat file

    #once a dictionary is created from the folder cotaing a certain set of simulations
    #the function loop inside of it plotting all the graphs, one over the other 
    fig, ax = plt.subplots(1,2, sharey=True, figsize=(12, 6))
    for i in resultDir_dic.values():
        config = Config.from_dicts(i)
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        trodes = config["trodes"]
        # useful constants
        
        Etheta = {"a": 0.}
        for trode in trodes:
            Etheta[trode] = -(k*Tref/e) * config[trode, "phiRef"]

        data_phi_part_c = (k*Tref/e)*sim_output['phi_part_c']
        ffvec = sim_output['ffrac_c'][0]
        times = sim_output['phi_applied_times'][0]

        print('data_phi_part_c: ', data_phi_part_c)
        print('ffvec: ', ffvec)

        ax[0].plot(ffvec, data_phi_part_c, label='elyte')
        ax[0].set_ylabel('phi_part_C (V)')
        ax[0].set_xlabel('Filling fraction')
        ax[0].set_title('V vs SOC')

        ax[1].plot(times, data_phi_part_c, label='elyte')
        ax[1].set_ylabel('phi_part_C (V)')
        ax[1].set_xlabel('Time (s)')
        ax[1].set_title('V vs t')

    return sim_output

def phi_bulk_c(resultDir_dic):
    k = constants.k                      # Boltzmann constant, J/(K Li)
    Tref = constants.T_ref               # Temp, K
    e = constants.e                      # Charge of proton, C
    # taking the mat file

    #once a dictionary is created from the folder cotaing a certain set of simulations
    #the function loop inside of it plotting all the graphs, one over the other 
    fig, ax = plt.subplots(1,2, sharey=True, figsize=(12, 6))
    for i in resultDir_dic.values():
        config = Config.from_dicts(i)
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        trodes = config["trodes"]
        # useful constants
        
        Etheta = {"a": 0.}
        for trode in trodes:
            Etheta[trode] = -(k*Tref/e) * config[trode, "phiRef"]

        data_phi_bulk_c = (k*Tref/e)*sim_output['phi_bulk_c']
        ffvec = sim_output['ffrac_c'][0]
        times = sim_output['phi_applied_times'][0]

        print('data_phi_bulk_c: ', data_phi_bulk_c)
        print('ffvec: ', ffvec)

        ax[0].plot(ffvec, data_phi_bulk_c, label='elyte')
        ax[0].set_ylabel('phi_part_C (V)')
        ax[0].set_xlabel('Filling fraction')
        ax[0].set_title('V vs SOC')

        ax[1].plot(times, data_phi_bulk_c, label='elyte')
        ax[1].set_ylabel('phi_part_C (V)')
        ax[1].set_xlabel('Time (s)')
        ax[1].set_title('V vs t')

    return sim_output


