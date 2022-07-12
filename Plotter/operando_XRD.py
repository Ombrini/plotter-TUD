import os.path as osp
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as manim
import numpy as np
import scipy.io as sio
import os

import mpet.mod_cell as mod_cell
import mpet.geometry as geom
import mpet.props_am as props_am
import mpet.utils as utils
from muFunc import *

def plot_xrd(resultDir_dic, folder_string_len):
    for resultDir in resultDir_dic.values():
        labels = '$' + str(resultDir[folder_string_len:]) + '$'
        fig, ax = plt.subplots(1,1)
        matfile = osp.join(resultDir, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(resultDir)
        td = config["t_ref"]
        if config["c","type"] == "ACR2":
            stech_1 = config["c","stech_1"]
        #times of the 0,0 particles are the same for all the particles
        times = sim_output['phi_applied_times'][0]*td 
        # ffvec = sim_output['ffrac_c'][0]
        
        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        bins = np.arange(start=0,stop=1,step= 0.01)
        xrd_matrix = np.empty([np.size(times),np.size(bins)-1])
        for t in np.arange(np.size(times)):
        # for t in [int(np.size(times)/2)]:
            conc_vec_at_t = np.array([])
            for i in np.arange(Nvol_c):
                for j in np.arange(Npart_c):
                    if config["c","type"] == "ACR2":
                        partStr1 = 'partTrodecvol'+str(i)+'part'+str(j)+'_c1'
                        c1 = sim_output[partStr1][t]
                        partStr2 = 'partTrodecvol'+str(i)+'part'+str(j)+'_c2'
                        c2 = sim_output[partStr2][t]
                        c_of_t = (c1*stech_1+c2*(1-stech_1))
                    else:
                        partStr = 'partTrodecvol'+str(i)+'part'+str(j)+'_c'
                        c = sim_output[partStr][t]
                        c_of_t = c

                    conc_vec_at_t = np.append(conc_vec_at_t, c_of_t)

            conc_hist_at_t, bin_edges  = np.histogram(conc_vec_at_t, bins = bins, density=True)
            xrd_matrix[t,:] = conc_hist_at_t
        img = ax.imshow(xrd_matrix, cmap = 'jet')
        colorbar = fig.colorbar(img)
        # ax.plot(conc_hist_at_t)
        ax.set_title(labels)

def plot_xrd_smooth(resultDir_dic, folder_string_len):
    for resultDir in resultDir_dic.values():
        labels = '$' + str(resultDir[folder_string_len:]) + '$'
        fig, ax = plt.subplots(1,1)
        matfile = osp.join(resultDir, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(resultDir)
        td = config["t_ref"]
        if config["c","type"] == "ACR2":
            stech_1 = config["c","stech_1"]
        #times of the 0,0 particles are the same for all the particles
        times = sim_output['phi_applied_times'][0]*td 
        # ffvec = sim_output['ffrac_c'][0]
        
        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        bins = np.arange(start=0,stop=1,step= 0.01)
        xrd_matrix = np.empty([np.size(times),(np.size(bins)-1)*2])
        for t in np.arange(np.size(times)):
        # for t in [int(np.size(times)/2)]:
            # conc_vec_at_t = np.array([])
            tot_xrd_at_time_t = np.zeros((np.size(bins)-1)*2)
            for i in np.arange(Nvol_c):
                for j in np.arange(Npart_c):
                    if config["c","type"] == "ACR2":
                        partStr1 = 'partTrodecvol'+str(i)+'part'+str(j)+'_c1'
                        c1 = sim_output[partStr1][t]
                        partStr2 = 'partTrodecvol'+str(i)+'part'+str(j)+'_c2'
                        c2 = sim_output[partStr2][t]
                        c_of_t = (c1*stech_1+c2*(1-stech_1))
                        c_hist_part, bin_edg_part = np.histogram(c_of_t, bins=bins, density=True)
                        f_of_c = np.zeros(np.size(c_hist_part)*2)
                        # print(f_of_c)
                        x = np.linspace(-np.size(c_hist_part)/2,1.5*np.size(c_hist_part),num = np.size(c_hist_part)*2)
                        # print(x)
                        for h in np.arange(np.size(c_hist_part)):
                            f_of_h = c_hist_part[h]
                            f_of_ch = f_of_h*np.exp(-0.5*((x - h)/10)**2)
                            f_of_c = f_of_c + f_of_ch
                        tot_xrd_at_time_t = tot_xrd_at_time_t + f_of_c
                    else:
                        partStr = 'partTrodecvol'+str(i)+'part'+str(j)+'_c'
                        c = sim_output[partStr][t]
                        c_of_t = c
                        c_hist_part, bin_edg_part = np.histogram(c_of_t, bins=bins, density=True)
                        f_of_c = np.zeros(np.size(c_hist_part)*2)
                        # print(f_of_c)
                        x = np.linspace(-np.size(c_hist_part)/2,1.5*np.size(c_hist_part),num = np.size(c_hist_part)*2)
                        # print(x)
                        for h in np.arange(np.size(c_hist_part)):
                            f_of_h = c_hist_part[h]
                            f_of_ch = f_of_h*np.exp(-0.5*((x - h)/10)**2)
                            f_of_c = f_of_c + f_of_ch
                        tot_xrd_at_time_t = tot_xrd_at_time_t + f_of_c
                    # conc_vec_at_t = np.append(conc_vec_at_t, c_of_t)

            # conc_hist_at_t, bin_edges  = np.histogram(conc_vec_at_t, bins = bins, density=True)
            xrd_matrix[t,:] = tot_xrd_at_time_t
        img = ax.imshow(xrd_matrix, cmap = 'jet')
        colorbar = fig.colorbar(img)
        # ax.plot(tot_xrd_at_time_t)
        ax.set_title(labels)