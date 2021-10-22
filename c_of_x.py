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

def plot_c(resultDir_dic, num_pictures):

# similarly to c_barLine, but the idea is to plot the c(x) for a selected reange of particles
# it is already possible to select the number of frames to plot
# the idea is to have also the possibility to plot in certain range of time
    fig, axes = plt.subplots(Npart_tot,Nvol_tot, sharey=True, figsize=(35, 35))
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        
        #times of the 0,0 particles are the same for all the particles
        times = sim_output['partTrodecvol0part0_cbar_times'][0] #in mpet1.7 exist only phi_times 
        ffvec = sim_output['ffrac_c'][0]
        config = Config.from_dicts(i)

        Nvol = config["Nvol"]
        trodes = config["trodes"]
        Npart = config["Npart"]
        psd_len = config["psd_len"]
        segments = config["segments"]
        dD_e = {}
        ndD_e = {}
        
        Nvol_c = Nvol['c']
        Npart_c = Npart['c']
        
        tottimesteps = len(times)-1
    #     print(tottimesteps)
            
        # Pick out some useful calculated values
        k = constants.k                      # Boltzmann constant, J/(K Li)
        Tref = constants.T_ref               # Temp, K
        e = constants.e                      # Charge of proton, C
        F = constants.F                      # C/mol
        c_ref = constants.c_ref
        
        partStr = 'partTrodecvol0part0_c'
        c = sim_output[partStr]
        num_of_pictures = 50
        # c is an array times x meters 
    #     timestep_vec = np.arange(tottimesteps)
    #     timestep_vec = np.arange(50)
        timestep_vec = np.around(np.linspace(0,tottimesteps,num=num_of_pictures),0)
        timestep_vec = timestep_vec.astype(np.int32)
    #     timestep_vec = np.arange(tottimesteps)[700:900]

        xaxis = np.arange(len(c[0]))
        fig, axes = plt.subplots(1,num_pictures, sharey=True, figsize=(50, 4))
        j = 0
        for i in timestep_vec:
    
            c_of_t = c[i]
            

            ax = axes[j]
            ax.plot(xaxis, c_of_t, label='c')
            ax.set_ylabel('c(x)')
            ax.set_xlabel('x')
            ax.set_title('c during times')
            j = j + 1

    return sim_output



    