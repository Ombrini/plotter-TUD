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
from get_values import get_C_rate, get_Nvol_c_and_Npart_c

# os.chdir('C:\\Users\pierfrancescoo\\Documents\\Simualtions\\mpet0_1_7\\bazantgroup-mpetpath')

# here I manually create a dictionary but I need to look into the "os" package
# to see how to create it looping into a folder
resultDir_dic = {}
resultDir_dic["sim1"] = "history\\100par_105"
resultDir_dic["sim2"] = "history\\100par_1005"

# to test the plotting I created a bunch of differnt simulations with different Crate or
# C segments and different Number of particles and volumes

# of course the plottings that take plots a range of particles and volumes have to be made
# in a folder in wich the simulations have the same number of particles and volumes 

# plot_voltage(resultDir_dic)
# plt.show()

Nvol_vec, Npart_vec = get_Nvol_c_and_Npart_c("history\\100par_105")
print(Nvol_vec)

Crate = get_C_rate("history\\100par_105")
print(Crate)