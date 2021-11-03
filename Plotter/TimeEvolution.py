import os.path as osp
import matplotlib.pyplot as plt
import numpy as np
from numpy.testing.nose_tools.utils import jiffies
import scipy.io as sio
import os

import mpet.mod_cell as mod_cell
import mpet.geometry as geom
import mpet.props_am as props_am
import mpet.utils as utils
from muFunc import *

from mpet.config import Config, constants

def plot_c(resultDir_dic, num_pictures, volume):

    config = Config.from_dicts(resultDir_dic['sim1'])
    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    H = Npart_c * 3 
    if H > 20:
        H = 20
    L = num_pictures * 3
    if L > 20:
        L = 20

    # similarly to c_barLine, but the idea is to plot the c(x) for a selected reange of particles
    # it is already possible to select the number of frames to plot
    # the idea is to have also the possibility to plot in certain range of tim
    fig, axes = plt.subplots(Npart_c,num_pictures, sharey=True, figsize=(10, 4)) 
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        td = config["t_ref"]
        #times of the 0,0 particles are the same for all the particles
        times = sim_output['phi_applied_times'][0]*td #in mpet1.7 exist only phi_times 
        ffvec = sim_output['ffrac_c'][0]
        config = Config.from_dicts(i)

        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        
        tottimesteps = len(times)-1
        timestep_vec = np.around(np.linspace(10,tottimesteps,num=num_pictures),0)
        timestep_vec = timestep_vec.astype(np.int32)

        if Npart_c == 1:
            partStr = 'partTrodecvol'+str(volume)+'part0_c'
            c = sim_output[partStr]
            xaxis = np.arange(len(c[0]))
            # fig, axes = plt.subplots(1,num_pictures, sharey=True, figsize=(10, 4))
            j = 0
            for i in timestep_vec:

                c_of_t = c[i]
                
                ax = axes[j]
                ax.plot(xaxis, c_of_t, label='c')
                ax.set_ylabel('c(x)')
                ax.set_xlabel('x' + str(i))
                j = j + 1

        else:
            for k in range(Npart_c):
                partStr = 'partTrodecvol'+str(volume)+'part'+str(k)+'_c'
                c = sim_output[partStr]
                xaxis = np.arange(len(c[0]))
                # fig, axes = plt.subplots(1,num_pictures, sharey=True, figsize=(10, 4))
                j = 0
                for i in timestep_vec:
            
                    c_of_t = c[i]
                    
                    ax = axes[k,j]
                    ax.plot(xaxis, c_of_t, label='c')
                    ax.set_ylabel('c(x)')
                    ax.set_xlabel('x' + str(i))
                    j = j + 1

    return sim_output

def plot_Rxn(resultDir_dic, num_pictures, volume):

    config = Config.from_dicts(resultDir_dic['sim1'])
    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    H = Npart_c * 3 
    if H > 20:
        H = 20
    L = num_pictures * 3
    if L > 20:
        L = 20

    # similarly to c_barLine, but the idea is to plot the c(x) for a selected reange of particles
    # it is already possible to select the number of frames to plot
    # the idea is to have also the possibility to plot in certain range of tim
    fig, axes = plt.subplots(Npart_c,num_pictures, sharey=True, figsize=(10, 4)) 
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        td = config["t_ref"]
        
        #times of the 0,0 particles are the same for all the particles
        times = sim_output['phi_applied_times'][0]*td #in mpet1.7 exist only phi_times 
        ffvec = sim_output['ffrac_c'][0]
        config = Config.from_dicts(i)

        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        
        tottimesteps = len(times)-1

        timestep_vec = np.around(np.linspace(0,tottimesteps,num=num_pictures),0)
        timestep_vec = timestep_vec.astype(np.int32)

        if Npart_c == 1:
            partStr = 'partTrodecvol'+str(volume)+'part0_Rxn'
            Rxn = sim_output[partStr]

            xaxis = np.arange(len(Rxn[0]))
            # fig, axes = plt.subplots(1,num_pictures, sharey=True, figsize=(10, 4))
            j = 0
            for i in timestep_vec:
        
                Rxn_of_t = Rxn[i]
                
                ax = axes[j]
                ax.plot(xaxis, Rxn_of_t, label='Rxn')
                # ax.set_ylabel('Rxn(x)')
                ax.set_xlabel('x')
                j = j + 1
        else:
            for k in range(Npart_c):

                partStr = 'partTrodecvol'+str(volume)+'part'+str(k)+'_Rxn'
                Rxn = sim_output[partStr]

                # c is an array times x meters 


                xaxis = np.arange(len(Rxn[0]))
                # fig, axes = plt.subplots(1,num_pictures, sharey=True, figsize=(10, 4))
                j = 0
                for i in timestep_vec:
            
                    Rxn_of_t = Rxn[i]
                    
                    ax = axes[k,j]
                    ax.plot(xaxis, Rxn_of_t, label='Rxn')
                    # ax.set_ylabel('Rxn(x)')
                    ax.set_xlabel('x')
                    j = j + 1

    return sim_output

def plot_mu(resultDir_dic, num_pictures, volume):

    config = Config.from_dicts(resultDir_dic['sim1'])
    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    H = Npart_c * 3 
    if H > 20:
        H = 20
    L = num_pictures * 3
    if L > 20:
        L = 20

    fig, axes = plt.subplots(Npart_c,num_pictures, sharey=True, figsize=(10, 4)) 
    for i in resultDir_dic.values():
        resultDir = i
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(i)
        td = config["t_ref"]
        #times of the 0,0 particles are the same for all the particles
        times = sim_output['phi_applied_times'][0]*td #in mpet1.7 exist only phi_times 
        ffvec = sim_output['ffrac_c'][0]
        
        

        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        
        tottimesteps = len(times)-1

        timestep_vec = np.around(np.linspace(10,tottimesteps,num=num_pictures),0)
        timestep_vec = timestep_vec.astype(np.int32)

        if Npart_c ==1:
            partStr = 'partTrodecvol'+str(volume)+'part0_c'
            c = sim_output[partStr]

            partStr_bar = 'partTrodecvol'+str(volume)+'part0_cbar'
            cbar = sim_output[partStr_bar]

            xaxis = np.arange(len(c[0]))
            j = 0
            for i in timestep_vec:
        
                conc = c[i]
                cavg = cbar[0][i]
                mu = LiFePO4(conc, cavg, resultDir)

                ax = axes[j]
                ax.plot(xaxis, mu, label='mu')
                ax.set_ylabel('mu(x)')
                ax.set_xlabel('x' + str(i))

                j = j + 1
        else:
            for k in range(Npart_c):

                partStr = 'partTrodecvol'+str(volume)+'part'+str(k)+'_c'
                c = sim_output[partStr]

                partStr_bar = 'partTrodecvol'+str(volume)+'part'+str(k)+'_cbar'
                cbar = sim_output[partStr_bar]

                xaxis = np.arange(len(c[0]))
                j = 0
                for i in timestep_vec:
            
                    conc = c[i]
                    cavg = cbar[0][i]
                    mu = LiFePO4(conc, cavg, resultDir)

                    ax = axes[k,j]
                    ax.plot(xaxis, mu, label='mu')
                    ax.set_ylabel('mu(x)')
                    ax.set_xlabel('x' + str(i))                   
                    j = j + 1

    return sim_output

    