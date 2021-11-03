# The script is not working properly, some function are just defined and not written 
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
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import os

import mpet.mod_cell as mod_cell
import mpet.geometry as geom
import mpet.props_am as props_am
import mpet.utils as utils

from mpet.config import Config, constants

from Voltage import plot_voltage
from SingleParticlePlots import *

# os.chdir('C:\\Users\pierfrancescoo\\Documents\\Simualtions\\mpet0_1_7\\bazantgroup-mpetpath')

# here I manually create a dictionary but I need to look into the "os" package
# to see how to create it looping into a folder
resultDir_dic = {}
# resultDir_dic["sim1"] = "history\\10par_105"
# resultDir_dic["sim2"] = "history\\10par_1005"

# resultDir_dic["sim1"] = "history\\2par_505_diff"
# resultDir_dic["sim1"] = "100vol_105"
resultDir_dic["sim1"] = "Nvol5Npart3_2C"
# 
# to test the plotting I created a bunch of differnt simulations with different Crate or
# C segments and different Number of particles and volumes

num_pic = 10
volume_selected = 0
# plot_mubar(resultDir_dic)
# plot_mu(resultDir_dic,num_pic,volume_selected)
# plot_c(resultDir_dic,num_pic,volume_selected)
# plot_dcbardt(resultDir_dic)
# plot_cbar(resultDir_dic,'time')
num_pic = 10
volume_selected = 80

# plot_mubar(resultDir_dic)
# plot_mu(resultDir_dic,num_pic,volume_selected)
# plot_c(resultDir_dic,num_pic,volume_selected)

plot_cbar(resultDir_dic,'ffrac')
# plot_voltage(resultDir_dic)
plt.show()

# keys("history\\1par_1")