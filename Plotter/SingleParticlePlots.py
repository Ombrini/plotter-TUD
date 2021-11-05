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


    fig, ax = plt.subplots(Npart_c,Nvol_c, sharey=True, figsize=(L, H), squeeze=False)

    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(i)
        td = config["t_ref"]

        times = sim_output['phi_applied_times'][0]*td
        ffvec = sim_output['ffrac_c'][0]
        if xaxis == 'time':
            xax = times
            xlabel = 'time'
        elif xaxis == 'ffrac':
            xax = ffvec
            xlabel = 'ffrac'
        else:
            print("!!! Third parameter must be 'time' or 'ffrac', , default = 'time' ")
            xax = times
            xlabel = 'time'

        for k in range(Nvol_c):
            for j in range(Npart_c):
                partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_cbar'
                cbar = sim_output[partStr][0]

                ax[j,k].plot(xax, cbar, label='cbar')
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

def plot_dcbardt(resultDir_dic, xaxis):
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
        td = config["t_ref"]

        times = sim_output['phi_applied_times'][0]*td
        ffvec = sim_output['ffrac_c'][0]
        if xaxis == 'time':
            xax = times
            xlabel = 'time'
        elif xaxis == 'ffrac':
            xax = ffvec
            xlabel = 'ffrac'
        else:
            print("!!! Third parameter must be 'time' or 'ffrac', default = 'time' ")
            xax = times
            xlabel = 'time'


        for k in range(Nvol_c):
            for j in range(Npart_c):
                partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_dcbardt'
                dcbardt = sim_output[partStr][0]

                if Npart_c == 1:
                    ax[k].plot(xax, dcbardt, label='dcbardt')
                    ax[k].set_ylabel('dLi/dt')
                    ax[k].set_xlabel(xlabel)
                    ax[k].set_title('volume: ' +str(k+1)+' particle: '+str(j+1))
                elif Nvol_c == 1:
                    ax[j].plot(xax, dcbardt, label='dcbardt')
                    ax[j].set_ylabel('dLi/dt')
                    ax[j].set_xlabel(xlabel)
                    ax[j].set_title('volume: ' +str(k+1)+' particle: '+str(j+1))
                else:
                    ax[j,k].plot(xax, dcbardt, label='dcbardt')
                    ax[j,k].set_ylabel('dLi/dt')
                    ax[j,k].set_xlabel(xlabel)
                    ax[j,k].set_title('volume: ' +str(k+1)+' particle: '+str(j+1))

    return sim_output

def plot_mubar(resultDir_dic, xaxis):

    config = Config.from_dicts(list(resultDir_dic.values())[0])
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
        td = config["t_ref"]
        tsteps = config["tsteps"]

        #times of the 0,0 particles are the same for all the particles
        times = sim_output['phi_applied_times'][0]*td #in mpet1.7 exist only phi_times
        ffvec = sim_output['ffrac_c'][0]
        if xaxis == 'time':
            xax = times
            xlabel = 'time'
        elif xaxis == 'ffrac':
            xax = ffvec
            xlabel = 'ffrac'
        else:
            print("!!! Third parameter must be 'time' or 'ffrac', default = 'time' ")
            xax = times
            xlabel = 'time'
        config = Config.from_dicts(i)

        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']

        tottimesteps = len(times)-1

        for k in range(Nvol_c):
            for j in range(Npart_c):

                partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_c'
                c = sim_output[partStr]

                partStr_bar = 'partTrodecvol'+str(k)+'part'+str(j)+'_cbar'
                cbar = sim_output[partStr_bar][0]

                mubar = np.empty([0])
                for t in range(tottimesteps+1):

                    conc = c[t]
                    cavg = cbar[t]
                    mu = LiFePO4(conc, cavg, resultDir)
                    mubar = np.append(mubar, np.sum(mu)/len(mu))
                if Npart_c == 1:
                    ax = axes[k]
                elif Nvol_c == 1:
                    ax = axes[j]
                else:
                    ax = axes[j,k]

                ax.plot(xax, mubar, label='mubar')
                ax.set_xlabel(xlabel)
                ax.set_ylabel('mubar')


    return sim_output

def plot_mubar_vs_cbar(resultDir_dic):

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
        td = config["t_ref"]

        #times of the 0,0 particles are the same for all the particles
        times = sim_output['phi_applied_times'][0]*td #in mpet1.7 exist only phi_times
        config = Config.from_dicts(i)

        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']

        tottimesteps = len(times)-1

        for k in range(Nvol_c):
            for j in range(Npart_c):

                partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_c'
                c = sim_output[partStr]

                partStr_bar = 'partTrodecvol'+str(k)+'part'+str(j)+'_cbar'
                cbar = sim_output[partStr_bar][0]

                mubar = np.empty([0])
                for t in range(tottimesteps+1):

                    conc = c[t]
                    cavg = cbar[t]
                    mu = LiFePO4(conc, cavg, resultDir)
                    mubar = np.append(mubar, np.sum(mu)/len(mu))
                if Npart_c == 1:
                    ax = axes[k]
                elif Nvol_c == 1:
                    ax = axes[j]
                else:
                    ax = axes[j,k]

                ax.plot(cbar, mubar, label='mubar')
                ax.set_xlabel('cbar')
                ax.set_ylabel('mubar')
                ax.set_title('volume: ' +str(k+1)+' particle: '+str(j+1))


    return sim_output

def plot_mubar_singlePart(resultDir_dic, vol, part):

    config = Config.from_dicts(list(resultDir_dic.values())[0])

    fig, axes = plt.subplots(1,2, sharey=True, figsize=(10, 4))
    for i in resultDir_dic.values():
        resultDir = i
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        td = config["t_ref"]

        #times of the 0,0 particles are the same for all the particles
        times = sim_output['phi_applied_times'][0]*td #in mpet1.7 exist only phi_times
        config = Config.from_dicts(i)

        tottimesteps = len(times)-1

        partStr = 'partTrodecvol'+str(vol)+'part'+str(part)+'_c'
        c = sim_output[partStr]

        partStr_bar = 'partTrodecvol'+str(vol)+'part'+str(part)+'_cbar'
        cbar = sim_output[partStr_bar][0]

        mubar = np.empty([0])
        for t in range(tottimesteps+1):

            conc = c[t]
            cavg = cbar[t]
            mu = LiFePO4(conc, cavg, resultDir)
            mubar = np.append(mubar, np.sum(mu)/len(mu))

        ax = axes[0]
        ax.plot(cbar, mubar, label='mubar')
        ax.set_xlabel('cbar')
        ax.set_ylabel('mubar')
        ax.set_title('volume: ' +str(vol+1)+' particle: '+str(part+1))

        ax = axes[1]
        ax.plot(times, mubar, label='mubar')
        ax.set_xlabel('time')
        ax.set_ylabel('mubar')


    return sim_output