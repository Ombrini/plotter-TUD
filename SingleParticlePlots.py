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
from muFunc import *

def plot_cbar(resultDir_dic, xaxis):
    matfile = osp.join(resultDir_dic["sim1"], 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    config = Config.from_dicts(resultDir_dic["sim1"])

    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    print(Npart_c)
    H = Npart_c * 4 
    if H > 25:
        H = 25
    L = Nvol_c * 3
    if L > 20:
        L = 20


    fig, ax = plt.subplots(Npart_c,Nvol_c, sharey=True, figsize=(L, H))
    
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(i)

        times = sim_output['phi_applied_times'][0]
        ffvec = sim_output['ffrac_c'][0]
        if xaxis == 'time':
            xaxis = times
            xlabel = 'time'
        elif xaxis == 'ffrac':
            xaxis = ffvec
            xlabel = 'ffrac'
        else:
            print("!!! Third parameter must be 'time' or 'ffrac' ")

        
        for k in range(Nvol_c):
            for j in range(Npart_c):
                partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_cbar'
                cbar = sim_output[partStr][0]
            
                if Npart_c == 1:
                    ax[k].plot(xaxis, cbar, label='dcbardt')
                    ax[k].set_ylabel('Li')
                    ax[k].set_xlabel(xlabel)
                    ax[k].set_title('volume: ' +str(k+1)+' particle: '+str(j+1))
                elif Nvol_c == 1:
                    ax[j].plot(xaxis, cbar, label='dcbardt')
                    ax[j].set_ylabel('Li')
                    ax[j].set_xlabel(xlabel)
                    ax[j].set_title('volume: ' +str(k+1)+' particle: '+str(j+1))
                else: 
                    ax[j,k].plot(xaxis, cbar, label='dcbardt')
                    ax[j,k].set_ylabel('Li')
                    ax[j,k].set_xlabel(xlabel)
                    ax[j,k].set_title('volume: ' +str(k+1)+' particle: '+str(j+1))

    return sim_output
    #the idea is to plot c_barLine of multiples results and also be able to zoom 
    #to specific volumes or specifi times

    #in this scritp there is also the possibility to select a time range and do the 
    #average of the various c_barLine of the particles distinguishing them between 
    # "activate" particles (which are charging) and "not active" ones

    #it's certainly better to create different functions and put them togather
    

def plot_dcbardt(resultDir_dic):
    matfile = osp.join(resultDir_dic["sim1"], 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    config = Config.from_dicts(resultDir_dic["sim1"])

    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    print(Npart_c)
    H = Npart_c * 3 
    if H > 20:
        H = 20
    L = Nvol_c * 3
    if L > 20:
        L = 20


    fig, ax = plt.subplots(Npart_c,Nvol_c, sharey=True, figsize=(L, H))
    
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(i)

        times = sim_output['phi_applied_times'][0]
        ffvec = sim_output['ffrac_c'][0]
        xaxis = times

        
        for k in range(Nvol_c):
            for j in range(Npart_c):
                partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_dcbardt'
                dcbardt = sim_output[partStr][0]
            
                if Npart_c == 1:
                    ax[k].plot(xaxis, dcbardt, label='dcbardt')
                    ax[k].set_ylabel('dLi/dt')
                    ax[k].set_xlabel('times')
                    ax[k].set_title('volume: ' +str(k+1)+' particle: '+str(j+1))
                elif Nvol_c == 1:
                    ax[j].plot(xaxis, dcbardt, label='dcbardt')
                    ax[j].set_ylabel('dLi/dt')
                    ax[j].set_xlabel('times')
                    ax[j].set_title('volume: ' +str(k+1)+' particle: '+str(j+1))
                else: 
                    ax[j,k].plot(xaxis, dcbardt, label='dcbardt')
                    ax[j,k].set_ylabel('dLi/dt')
                    ax[j,k].set_xlabel('times')
                    ax[j,k].set_title('volume: ' +str(k+1)+' particle: '+str(j+1))

    return sim_output

def plot_mubar(resultDir_dic):

    config = Config.from_dicts(resultDir_dic['sim1'])
    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    H = Npart_c * 3 
    if H > 20:
        H = 20
    L = 1 * 3
    if L > 20:
        L = 20

    fig, axes = plt.subplots(Npart_c,Nvol_c, sharey=True, figsize=(10, 4)) 
    for i in resultDir_dic.values():
        resultDir = i
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        
        #times of the 0,0 particles are the same for all the particles
        times = sim_output['phi_applied_times'][0] #in mpet1.7 exist only phi_times 
        ffvec = sim_output['ffrac_c'][0]
        config = Config.from_dicts(i)

        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        
        tottimesteps = len(times)-1

        for k in range(Nvol_c):
            for j in range(Npart_c):

                partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_c'
                c = sim_output[partStr]

                partStr_bar = 'partTrodecvol'+str(k)+'part'+str(j)+'_cbar'
                cbar = sim_output[partStr_bar]
        
                mubar = np.empty([0])
                for i in range(tottimesteps+1):

                    conc = c[i]
                    cavg = cbar[0][i]
                    mu = LiFePO4(conc, cavg, resultDir)
                    mubar = np.append(mubar, np.sum(mu)/len(mu))
                    
                ax = axes[j,k]
                ax.plot(times, mubar, label='mubar')
                ax.set_xlabel('t')
                ax.set_ylabel('mubar')
                

    return sim_output