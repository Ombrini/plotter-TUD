# the idea is just to show the concepts of the plotting we want
# we know that the final part of data analysis will be done by ourself, but it's useful to
# have a tool that fastly permits to expand the base plots and to compare different setups


# it would olso be nice to expand the sim_output keys adding for example
# partTrode' 'vol' 'part' '_ mu (chemical potential) wich is a function of c, x (or r) and t

# other ideas:
# improving the movie in which only three colors rapresent the
# concentration in the particle towards a set of 10 or more shades from yellow to red
#
# Creating an analogous movie with the interface outside the particles that changes colors
# with the avg conc inside
#
# Also having a shade of colors in the backgroud of the particle that represent the concentration
# in the electrolyte could help visualize how the Li moves with different electrolyte

#ps: I'm sorry for the mess of files, I don't know the "best practices" of this kind of code.
# so if you just organize the code deviding the various tasks present in different functins
# we can start from there and create other plotting tools for our personal needs



import os.path as osp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as manim
from matplotlib.animation import writers
import numpy as np
import scipy.io as sio
import os
import sys

import mpet.mod_cell as mod_cell
import mpet.geometry as geom
import mpet.props_am as props_am
import mpet.utils as utils

from mpet.config import Config, constants
from elyte import *

from Voltage import *
from plot_utils import *
from SingleParticlePlots import *
from TimeEvolution import *

from OCVproject import *
# from get_values import get_C_rate, get_Nvol_c_and_Npart_c


if __name__ == '__main__':
    # os.chdir('C:\\Users\pierfrancescoo\\Documents\\Simualtions\\mpet0_1_7\\bazantgroup-mpetpath')
    try:
        simulation_folder = sys.argv[1]
    except IndexError:
        print('Provide path to simulation outputs')
        exit()
    simulations = os.listdir(simulation_folder)

    resultDir_dic = {}
    for sim in simulations:
        resultDir_dic[sim] = os.path.join(simulation_folder, sim)


    # plot_voltage(resultDir_dic)
    # plot_activeParticles(resultDir_dic)
    # plt.show()
    plot_c2D(resultDir_dic,10,0)
    # plot_c(resultDir_dic, 10 ,0)
    # plot_mu(resultDir_dic, 10,0)
    # plot_mubar_vs_cbar(resultDir_dic) #all the mubars vs the cbar of the particle
    # plot_mubar_singlePart(resultDir_dic,0,0)
    # plot_mubar_singlePart(resultDir_dic,0,0)
    # plot_mubar(resultDir_dic, 'ffrac') #all the mubars vs the cbar of the system
    # plot_cbar(resultDir_dic, 'ffrac')
    # plot_cbar(resultDir_dic, 'time')
    # plot_Crate_singleparticle(resultDir_dic, 'ffrac', 0)
    # plot_cVolume(resultDir_dic,0.49)
    # plot_Vmax_vs_Thick(resultDir_dic,'Vvsb')
    # plot_cVolume(resultDir_dic,0.5)
    # plot_cVolume(resultDir_dic,0.9)
    # elyte_c(resultDir_dic,5)
    # elyte_phi(resultDir_dic,5)

    # plt.show()
    # keys(str(resultDir_dic))

    # to test the plotting I created a bunch of differnt simulations with different Crate or
    # C segments and different Number of particles and volumes

    # of course the plottings that take plots a range of particles and volumes have to be made
    # in a folder in wich the simulations have the same number of particles and volumes

    config = Config.from_dicts(list(resultDir_dic.values())[0])
    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    fig, axes = plt.subplots(Npart_c,Nvol_c, sharey=False)
    # fig, im = plt.subplots(1,1, sharey=False)
    # fig =  plt.figure()
    norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1)
    matfile = osp.join(list(resultDir_dic.values())[0], 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    td = config["t_ref"]
    times = sim_output['phi_applied_times'][0]*td
    numtimes = np.size(times)
    ffvec = sim_output['ffrac_c'][0]
    lines = np.empty((Npart_c, Nvol_c), dtype=object)

    for k in range(Npart_c):
        for j in range(Nvol_c):
            partStr = 'partTrodecvol'+str(j)+'part'+str(k)+'_c'
            # im = axes[k,j]
            c = sim_output[partStr]
            # A = [ np.array()]
            N = np.size(c[0])
            A = np.ones((int(N/4),N))
            A[:,:] = c[0]
            line = axes[k,j].imshow(A, cmap = 'jet', norm = norm, animated = True)
            lines[k,j] = line
            # A[:,:] = c[0]
            # print(A)
            # # ax = axes[k,j]
            # ax = axes
            # im = plt.imshow(A, animated=True)

    def init():
        A[:,:] = c[50]
        lines[k,j].set_data(A)
        return [lines]
    
    def animate(t):
        for k in range(Npart_c):
            for j in range(Nvol_c):
                partStr = 'partTrodecvol'+str(j)+'part'+str(k)+'_c'
                c = sim_output[partStr]
                N = np.size(c[0])
                A = np.ones((int(N/2),N))
                toblit = []
                # a=im.get_array()
                # A = im.get_array()
                A[:,:] = c[t]
                # print(A)
                # A[:,:] = np.reshape(c[t],N)
                # a = A     
                # lines.clear()
                # im.contourf(range(N),range(N),A, cmap = 'jet', norm = norm)
                # im.imshow(A, cmap = 'jet', norm = norm)
                # print(k, ' ', j )
                lines[k,j].set_data(A)
                lines_local = lines.copy()
                toblit.append(lines_local)
                # im.set_array(A)S
        return [lines]

    anim = manim.FuncAnimation(fig, animate,
                                            frames=numtimes, interval=0.001, repeat = False, init_func= init)
                        # Writer = writers['ffmeg']
            # writer = Writer(10)
            # anim.save("mpet_anim.mp4")

    plt.show()
    exit()


