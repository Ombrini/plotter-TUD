import os.path as osp
from re import L
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


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

# import mpet.mod_cell as mod_cell
# import mpet.geometry as geom
# import mpet.props_am as props_am
# import mpet.utils as utils

from mpet.config import Config, constants

# this simply plots v vs ff and v vs t in folder 
def plot_voltage(resultDir_dic, folder_string_len):
    k = constants.k                      # Boltzmann constant, J/(K Li)
    Tref = constants.T_ref               # Temp, K
    e = constants.e                      # Charge of proton, C
    # taking the mat file
    # once a dictionary is created from the folder cotaing a certain set of simulations
    # the function loop inside of it plotting all the graphs, one over the other 
    fig, ax = plt.subplots(1,2, sharey=True, figsize=(12, 6))
    max_V = np.array([])
    min_V = np.array([])
    for i in resultDir_dic.values():
        config = Config.from_dicts(i)
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        trodes = config["trodes"]
        td = config["t_ref"]
        # useful constants
        labels = '$' + str(i[folder_string_len:]) + '$'
        
        Etheta = {"a": 0.}
        for trode in trodes:
            Etheta[trode] = -(k*Tref/e) * config[trode, "phiRef"]

        data_volt = Etheta["c"] - Etheta["a"] - (k*Tref/e)*sim_output['phi_cell'][0]
        ffvec = sim_output['ffrac_c'][0]
        ffvec = sim_output['ffrac_c'][0] #+ 1 - np.max(ffvec)
        times = sim_output['phi_applied_times'][0]*td

        max_V = np.append(max_V,np.amax(data_volt))
        min_V = np.append(min_V,np.amin(data_volt))

        ax[0].plot(ffvec, data_volt, label=labels)
        ax[0].set_ylabel('Voltage (V)')
        ax[0].set_xlabel('DoD')
        ax[0].set_title('V vs DoD')

        # ax[0].yaxis.set_major_formatter(FormatStrFormatter('%g'))
        # ax[0].yaxis.set_ticks(np.arange(4.5,2.5,0.1))
        # ax[0].yaxis.set_ticks_position('left')
        # ax[0].set_ylim([2.5,4.5])
        # ax[0].set_xlim([1,0])
        ax[0].legend()

        ax[1].plot(times, data_volt, label=labels)
        ax[1].set_ylabel('Voltage (V)')
        ax[1].set_xlabel('Time (s)')
        ax[1].set_title('V vs t')
        # ax[1].yaxis.set_ticks(np.arange(4.5,2.5,0.1))
        # ax[1].yaxis.set_ticks_position('left')
        ax[1].legend()
        # ax[1].yaxis.set_major_formatter(FormatStrFormatter('%g'))
        # ax[1].yaxis.set_ticks(np.arange(3, 4, 0.05))
        # plt.savefig("employee_growth_tp.png", transparent=True)

    return sim_output





