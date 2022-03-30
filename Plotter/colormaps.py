from cProfile import label
import os.path as osp
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
# from mpl_toolkits import mplot3d
from scipy.interpolate import interp2d
from matplotlib.ticker import LinearLocator
from matplotlib import cm

import mpet.mod_cell as mod_cell
import mpet.geometry as geom
import mpet.props_am as props_am
import mpet.utils as utils

from mpet.config import Config, constants

# this simply plots v vs ff and v vs t in folder 
import os.path as osp
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
# import os
from matplotlib.ticker import FormatStrFormatter

from plot_utils import *

from mpet.config import Config, constants


def plot_Vmax_vs_Thick(resultDir_dic,input):

    lastV_vec = np.array([])
    thickness_vec = np.array([])
    bulk_vec = np.array([])

    for i in resultDir_dic.values():
        config = Config.from_dicts(i)
        segments = config["segments"]
        initial_Crate = segments[0][0]*config['t_ref']/60

        lastV = get_last_V(i)
        lastV_vec = np.append(lastV_vec,lastV)

        thickness = get_thickness(i,'c')
        thickness_vec = np.append(thickness_vec,thickness)

        bulk = get_bulkCon(i,'c')
        bulk_vec = np.append(bulk_vec,bulk)
        log_bulk = np.log(bulk_vec)

    if input == 'Vvst':
        plt.scatter(thickness_vec, lastV_vec, s=1000, c=np.log(bulk_vec), cmap='Greens')
        plt.colorbar(label='ln bulk cond')
        plt.xlabel('thickness')
        plt.ylabel('OCV')

    elif input == 'Vvsb':
        plt.scatter(np.log(bulk_vec), lastV_vec, s=1000, c=thickness_vec, cmap='cool')
        plt.colorbar(label='thickness')
        plt.xlabel('bulk_cond')
        plt.ylabel('OCV')
    else:
        print('options are: Vvsb or Vvst')

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.scatter3D(thickness_vec, np.log(bulk_vec), lastV_vec, s = 1000,  c=lastV_vec,cmap = 'plasma')
    # ax.set_xlabel('thickess')
    # ax.set_ylabel('bulk conductivity')
    # ax.set_zlabel('OCV')

def plot_MaxActivePart_vsCrate_and_bulk(resultDir_dic):
    thickness_vec = np.array([])
    bulk_vec = np.array([])
    crate_vec = np.array([])

    for i in resultDir_dic.values():

        crate = get_C_rate(i)
        crate = crate[0]
        crate_vec = np.append(crate_vec,crate)

        thickness = get_thickness(i,'c')
        thickness_vec = np.append(thickness_vec,thickness)

        bulk = get_bulkCon(i,'c')
        bulk_vec = np.append(bulk_vec,bulk)

    crate_sort = np.unique(np.sort(crate_vec))
    bulk_sort = np.unique(np.sort(bulk_vec))

    mat = np.zeros((np.size(bulk_sort),np.size(crate_sort)))

    for i in resultDir_dic.values():
        config = Config.from_dicts(i)
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        tot_particles = Nvol_c*Npart_c
        ffvec = sim_output['ffrac_c'][0]

        crate_i = get_C_rate(i)
        crate_i = crate_i[0]

        thickness_i = get_thickness(i,'c')

        bulk_i = get_bulkCon(i,'c')

        num_active_vec = np.array([])

        for t in range(np.size(ffvec)):
            num_active = 0
            for k in range(Nvol_c):
                for j in range(Npart_c):
                        partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_cbar'
                        cbar = sim_output[partStr][0][t]
                        if cbar > 0.15 and cbar < 0.85:
                            num_active = num_active + 1

            num_active_vec = np.append(num_active_vec,num_active)
           
        percent_active = (num_active_vec/tot_particles)*100
        percent_max = np.max(percent_active)

        index_crate = 0 
        for i in crate_sort:
            if i == crate_i:
                break
            index_crate += 1

        index_bulk = 0 
        for i in bulk_sort:
            if i == bulk_i:
                break
            index_bulk += 1
        
        mat[index_bulk,index_crate] = int(percent_max)
        # if mat[index_bulk,index_crate] == 0:
        #     if mat[index_bulk + 1,index_crate] and mat[index_bulk - 1,index_crate] != 0:
        #         mat[index_bulk,index_crate] = 0.5*(mat[index_bulk + 1,index_crate]+mat[index_bulk - 1,index_crate])
        #     else:
        #         mat[index_bulk,index_crate] = 0.5*(mat[index_bulk,index_crate - 1]+mat[index_bulk,index_crate + 1])
    print(mat)
    print('-------------------')
    y,x = np.where(mat != 0)  
    mat_inter = interp2d(x,y,mat[mat!=0],kind='linear')
    X = np.arange(len(crate_sort))
    Y = np.arange(len(bulk_sort))
    print('----------------------')
    print(mat_inter(X,Y))
    # mat = mat_inter(X,Y)
    for i in X:
        for j in Y:
            if mat[j,i] < 0:
                mat[j,i] = 0

    norm = plt.Normalize(0, 100)
    fig, ax = plt.subplots(1,1)
    img = ax.imshow(mat,extent=[0,np.size(crate_sort),0,np.size(bulk_sort)], cmap = 'jet', norm = norm)
    # img = ax.plot_surface(crate_sort,bulk_sort,mat_inter(X,Y), cmap = 'jet', linewidth = 0, antialiased = False)
    x_label_list = crate_sort
    y_label_list = bulk_sort
    xstick = np.linspace(0.5,(np.size(crate_sort)-0.5), num = (np.size(crate_sort)))
    ax.set_xticks(xstick)
    ystick = np.linspace((np.size(bulk_sort)-0.5),0.5, num = (np.size(bulk_sort)))
    ax.set_yticks(ystick)
    ax.set_xticklabels(x_label_list)
    ax.set_yticklabels(y_label_list)
    ax.set_xlabel('C rate')
    ax.set_ylabel('electrical conductivity S/m')
    colorbar = fig.colorbar(img)

    # colorbar.ax.set_ylabel('%/ active partivecles at maximum', rotation=270)

    # plt.scatter(crate_vec, bulk_vec , s=1000, c=percent_active, cmap='Green')
    # plt.colorbar(label='max active particles')
    # plt.xlabel('C rate')
    # plt.ylabel('electrical conductivity S/m')

def plot_MaxActivePart_vsCrate_1Dplot(resultDir_dic):
    thickness_vec = np.array([])
    bulk_vec = np.array([])
    crate_vec = np.array([])
    part_dim_vec = np.array([])

    for i in resultDir_dic.values():

        crate = get_C_rate(i)
        crate = crate[0]
        crate_vec = np.append(crate_vec,crate)

        thickness = get_thickness(i,'c')
        thickness_vec = np.append(thickness_vec,thickness)

        bulk = get_bulkCon(i,'c')
        bulk_vec = np.append(bulk_vec,bulk)

        part_dim = i[-27:-24]
        part_dim_vec = np.append(part_dim_vec,part_dim)

    crate_sort = np.unique(np.sort(crate_vec))
    part_dim_sort = np.unique(np.sort(part_dim_vec))
    print(part_dim_sort)
    # bulk_sort = np.unique(np.sort(bulk_vec))

    mat = np.zeros((np.size(part_dim_sort),np.size(crate_sort)))

    for i in resultDir_dic.values():
        config = Config.from_dicts(i)
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        tot_particles = Nvol_c*Npart_c
        ffvec = sim_output['ffrac_c'][0]

        crate_i = get_C_rate(i)
        crate_i = crate_i[0]

        part_dim_i = i[-27:-24]

        # thickness_i = get_thickness(i,'c')

        # bulk_i = get_bulkCon(i,'c')

        num_active_vec = np.array([])

        for t in range(np.size(ffvec)):
            num_active = 0
            for k in range(Nvol_c):
                for j in range(Npart_c):
                        partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_cbar'
                        cbar = sim_output[partStr][0][t]
                        if cbar > 0.15 and cbar < 0.85:
                            num_active = num_active + 1

            num_active_vec = np.append(num_active_vec,num_active)
           
        percent_active = (num_active_vec/tot_particles)*100
        percent_max = np.max(percent_active)

        index_crate = 0 
        for i in crate_sort:
            if i == crate_i:
                break
            index_crate += 1

        index_part = 0 
        for i in part_dim_sort:
            if i == part_dim_i:
                break
            index_part += 1

    #     index_bulk = 0 
    #     for i in bulk_sort:
    #         if i == bulk_i:
    #             break
    #         index_bulk += 1
        
        mat[index_part,index_crate] = int(percent_max)

    fig, ax = plt.subplots(1,1)
    for i in range(len(part_dim_sort)):
        ax.plot(crate_sort, mat[i,:], label = ('Avg part dim (nm) = ' + part_dim_sort[i]))
        ax.set_xlabel('C rate')
        ax.set_ylabel('Max active particles')
        ax.set_title('Electrical conductivity = 0.1 S/m')
        ax.legend()
        # if mat[index_bulk,index_crate] == 0:
    #     #     if mat[index_bulk + 1,index_crate] and mat[index_bulk - 1,index_crate] != 0:
    #     #         mat[index_bulk,index_crate] = 0.5*(mat[index_bulk + 1,index_crate]+mat[index_bulk - 1,index_crate])
    #     #     else:
    #     #         mat[index_bulk,index_crate] = 0.5*(mat[index_bulk,index_crate - 1]+mat[index_bulk,index_crate + 1])
    # print(mat)
    # print('-------------------')
    # y,x = np.where(mat != 0)  
    # mat_inter = interp2d(x,y,mat[mat!=0],kind='linear')
    # X = np.arange(len(crate_sort))
    # Y = np.arange(len(bulk_sort))
    # print('----------------------')
    # print(mat_inter(X,Y))
    # # mat = mat_inter(X,Y)
    # for i in X:
    #     for j in Y:
    #         if mat[j,i] < 0:
    #             mat[j,i] = 0

    # norm = plt.Normalize(0, 100)
    # fig, ax = plt.subplots(1,1)
    # img = ax.imshow(mat,extent=[0,np.size(crate_sort),0,np.size(bulk_sort)], cmap = 'jet', norm = norm)
    # x_label_list = crate_sort
    # y_label_list = bulk_sort
    # xstick = np.linspace(0.5,(np.size(crate_sort)-0.5), num = (np.size(crate_sort)))
    # ax.set_xticks(xstick)
    # ystick = np.linspace((np.size(bulk_sort)-0.5),0.5, num = (np.size(bulk_sort)))
    # ax.set_yticks(ystick)
    # ax.set_xticklabels(x_label_list)
    # ax.set_yticklabels(y_label_list)
    # ax.set_xlabel('C rate')
    # ax.set_ylabel('electrical conductivity S/m')
    # colorbar = fig.colorbar(img)

    # colorbar.ax.set_ylabel('%/ active partivecles at maximum', rotation=270)


    
