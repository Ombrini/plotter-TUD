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

def plot_cbarLine(resultDir_dic, Nvol_iniz, Nvol_fin, Npart_iniz, Npart_fin):
    #the idea is to plot c_barLine of multiples results and also be able to zoom 
    #to specific volumes or specifi times

    #in this scritp there is also the possibility to select a time range and do the 
    #average of the various c_barLine of the particles distinguishing them between 
    # "activate" particles (which are charging) and "not active" ones

    #it's certainly better to create different functions and put them togather
    Npart_tot = Npart_fin-Npart_iniz
    Nvol_tot = Nvol_fin- Npart_iniz
    fig, axes = plt.subplots(Npart_tot,Nvol_tot, sharey=True, figsize=(35, 35))

    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        
        #times of the 0,0 particles are the same for all the particles
        times = sim_output['partTrodecvol0part0_cbar_times'][0]
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
        
        Crate0 = round(segments[0][0]/0.00216464,1)
        Crate1 = round(segments[2][0]/0.00216464,1)

        tottimesteps = len(times)
        timestep_end = round((tottimesteps/(1/Crate0 + 2 + 1/Crate1))*(1/Crate0 + 2))

        # Pick out some useful calculated values
        k = constants.k                      # Boltzmann constant, J/(K Li)
        Tref = constants.T_ref               # Temp, K
        e = constants.e                      # Charge of proton, C
        F = constants.F                      # C/mol
        c_ref = constants.c_ref
        
        # for loop for plotting the desired set of particles
        Nvol_c_vec = np.arange(Nvol_iniz, Nvol_fin)
        Npart_c_vec = np.arange(Npart_iniz, Npart_fin)
        delith = 0
        almost0 = 0
        almost1 = 0
        
        delith_vec = np.empty([0])
        avg_conc_delith = np.empty([0])
        lith_vec = np.empty([0])
        avg_conc_lith = np.empty([0])
        
    #     print('Crate writing phase = ', Crate0)
        for i in Nvol_c_vec:
                for j in Npart_c_vec:
        
                    partStr = 'partTrodecvol'+str(i)+'part'+str(j)+'_cbar'
                    cbar = sim_output[partStr][0]
                    
                    
                    if cbar[timestep_end-2] < 0.4:
                        delith += 1
                    if cbar[timestep_end-2] < 0.4:
                        delith_vec = np.append(delith_vec, cbar[timestep_end-2])
                    elif cbar[timestep_end-2] > 0.6:
                        lith_vec = np.append(lith_vec, cbar[timestep_end-2])

                    #times or ffcev
                    xaxis = times
    #                 xaxis = ffvec

                    # d[Li]/dt vs time
                    ax = axes[j,i]
                    ax.plot(xaxis, cbar, label='cbarLine')
                    ax.set_ylabel('Li conc')
                    ax.set_xlabel('times')

                    ax.set_title('volume: ' +str(i+1)+' particle: '+str(j+1))
                    
        avg_conc_delith = np.average(delith_vec)
        avg_conc_lith = np.average(lith_vec)
    
    return sim_output, Crate0, avg_conc_delith, avg_conc_lith
    
    