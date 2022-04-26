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



import imp
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
from cbar_c_gamma import *

from colormaps import *
from Plot_experiments import *
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

    file_name = 'ExportData_PFO_Li_LFP_170220222_1_SEE Rack 2_09_6_20220217165657.txt'
    mass = 16.4

    # plot_voltage(resultDir_dic)
    # plot_experiments(file_name, mass)
    # plot_activeParticles(resultDir_dic)
    plot_csurf_vs_cbar(resultDir_dic)
    plot_ecd_vs_cbar(resultDir_dic)
    # active_max_vs_Crate(resultDir_dic)
    # plt.show()
    # active_max_vs_bulk(resultDir_dic)
    # active_max_vs_thickness(resultDir_dic)
    # plot_MaxActivePart_vsCrate_and_bulk(resultDir_dic)
    # plot_MaxActivePart_vsCrate_1Dplot(resultDir_dic)
    # active_atSoC(resultDir_dic, 50)
    # active_atSoC(resultDir_dic, 50)
    # active_atSoC(resultDir_dic, 50)
    # plt.show()
    # plot_csld2D(resultDir_dic, save = False, directory= 0)
    # plot_c(resultDir_dic, 10 ,0)
    # plot_mu(resultDir_dic, 10,0)
    # plot_mubar_vs_cbar(resultDir_dic) #all the mubars vs the cbar of the particle
    # plot_mubar_singlePart(resultDir_dic,0,0)
    # plot_mubar_singlePart(resultDir_dic,0,0)
    # plot_mubar(resultDir_dic, 'ffrac') #all the mubars vs the cbar of the system
    plot_cbar(resultDir_dic, 'ffrac')
    # plot_cbar(resultDir_dic, 'time')
    # plot_Crate_singleparticle(resultDir_dic, 'ffrac', 0.1)
    # plot_cVolume(resultDir_dic,0.5)
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
    # in a folder in wich the simulations have the same number of particles and volume

    plt.show()
    exit()


