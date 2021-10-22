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


def plot_battery(resultDir, vol, part, var, tval):
    matfile = osp.join(resultDir, 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    times = sim_output['phi_applied_times'][0]
    config = Config.from_dicts(resultDir)

    Nvol = config["Nvol"]
    trodes = config["trodes"]
    Npart = config["Npart"]
    psd_len = config["psd_len"]
    dD_e = {}
    ndD_e = {}

    # Pick out some useful calculated values
    k = constants.k                      # Boltzmann constant, J/(K Li)
    Tref = constants.T_ref               # Temp, K
    e = constants.e                      # Charge of proton, C
    F = constants.F                      # C/mol
    c_ref = constants.c_ref
    td = config["t_ref"]
    Etheta = {"a": 0.}
    for trode in trodes:
        Etheta[trode] = -(k*Tref/e) * config[trode, "phiRef"]
    Vstd = Etheta["c"] - Etheta["a"]
    
    Lfac = 1e6
    Lunit = r"$\mu$m"
    dxc = config["L"]["c"]/Nvol["c"]
    dxvec = np.array(Nvol["c"] * [dxc])
    porosvec = np.array(Nvol["c"] * [config["poros"]["c"]])
    cellsvec = dxc*np.arange(Nvol["c"]) + dxc/2.
    cellsvec_c = dxc*np.arange(Nvol["c"]) + dxc/2.
    cellsvec_c *= config["L_ref"] * Lfac
    if Nvol["s"]:
        dxs = config["L"]["s"]/Nvol["s"]
        dxvec_s = np.array(Nvol["s"] * [dxs])
        dxvec = np.hstack((dxvec_s, dxvec))
        poros_s = np.array(Nvol["s"] * [config["poros"]["s"]])
        porosvec = np.hstack((poros_s, porosvec))
        cellsvec += config["L"]["s"] / config["L"]["c"]
        cellsvec_s = dxs*np.arange(Nvol["s"]) + dxs/2.
        cellsvec = np.hstack((cellsvec_s, cellsvec))
    if "a" in trodes:
        dxa = config["L"]["a"]/Nvol["a"]
        dxvec_a = np.array(Nvol["a"] * [dxa])
        dxvec = np.hstack((dxvec_a, dxvec))
        poros_a = np.array(Nvol["a"] * [config["poros"]["a"]])
        porosvec = np.hstack((poros_a, porosvec))
        cellsvec += config["L"]["a"] / config["L"]["c"]
        cellsvec_a = dxa*np.arange(Nvol["a"]) + dxa/2.
        cellsvec = np.hstack((cellsvec_a, cellsvec))
    cellsvec_s *= config["L_ref"] * Lfac
    cellsvec *= config["L_ref"] * Lfac
    facesvec = np.insert(np.cumsum(dxvec), 0, 0.) * config["L_ref"] * Lfac
    
    # Uncomment these lines to show the variables you can plot
    # for k in sim_output.keys():
    #     print(k)
    # return

    trode = 'c'
    t = np.argmin(np.abs(times - tval))
    if tval > times[-1]:
        raise ValueError(f'Max time: {times[-1]}, got {tval}')
    
    var_names = {'c': 'concentration',
                 'phi': 'potential'}
    var_long = var_names[var]
    
    fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(8, 4))

    # elyte data elyte
    key = f'{var}_lyte_s'
    if var in 'c':
        data_elyte_s = sim_output[key]* c_ref / 1000.  # select correct elyte volume
    elif var in 'phi':
        data_elyte_s = sim_output[key]*(k*Tref/e) - Vstd
    x_elyte_s = list(range(-np.shape(data_elyte_s)[1],0))
    # elyte value at timestep chosen to plot interface region volumes
    if var in 'c':
        cval_s = data_elyte_s[t]
    elif var in 'phi':
        cval_s = data_elyte_s[t]
    data_elyte_vol_s = sim_output[key][t]  # all volumes at given t
    
    # elyte data electrode
    key = f'{var}_lyte_{trode}'
    if var in 'c':
        data_elyte_c = sim_output[key]* c_ref / 1000.  # select correct elyte volume
    elif var in 'phi':
        data_elyte_c = sim_output[key]*(k*Tref/e) - Vstd
    x_elyte_c = list(range(np.shape(data_elyte_c)[1]))
    # elyte value at timestep chosen to plot interface region volumes
    if var in 'c':
        cval_c = data_elyte_c[t]
    elif var in 'phi':
        cval_c = data_elyte_c[t]
    data_elyte_vol_c = sim_output[key][t]  # all volumes at given t
    
    # time vs variable
    ax = axes[0]
    ax.plot(times, data_elyte_c)
    ax.plot(times, data_elyte_s[:,-1], label='elyte_s last volume')

    ax.set_xlabel('time')
    ax.set_ylabel(var_long)
    ax.axvline(tval, ls='--', c='k')
#    ax.axhline(cval_c, ls='--', c='k')
    ax.legend()
    ax.set_title('Time evolution')

    #  volume vs variable
    ax = axes[1]   
#     ax.plot(x_elyte_c,cval_c, marker='o', label='elyte_c')
#     ax.plot(x_elyte_s,cval_s , label='elyte_s', c='g', marker='o')
    
#     datay = np.hstack((cval_s, cval_c))
#     ax.plot(cellsvec,datay , label='elyte_s', c='g', marker='o')
    
    ax.plot(cellsvec_s[-1]+cellsvec_c,cval_c, marker='o', label='elyte_c')
    ax.plot(cellsvec_s,cval_s , label='elyte_s', c='g', marker='o')
    
    ax.set_xlabel('Battery position [um]')
    ax.legend()
#    ax.set_ylim(0, 0.7)

    fig.suptitle(var_long)
    fig.tight_layout()
    
    return sim_output



def plot_i(resultDir, vol, part, tval):
    matfile = osp.join(resultDir, 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    times = sim_output['phi_applied_times'][0]
    config = Config.from_dicts(resultDir)

    Nvol = config["Nvol"]
    trodes = config["trodes"]
    Npart = config["Npart"]
    psd_len = config["psd_len"]
    dD_e = {}
    ndD_e = {}

    # Pick out some useful calculated values
    k = constants.k                      # Boltzmann constant, J/(K Li)
    Tref = constants.T_ref               # Temp, K
    e = constants.e                      # Charge of proton, C
    F = constants.F                      # C/mol
    c_ref = constants.c_ref
    td = config["t_ref"]
    Etheta = {"a": 0.}
    for trode in trodes:
        Etheta[trode] = -(k*Tref/e) * config[trode, "phiRef"]
    Vstd = Etheta["c"] - Etheta["a"]
    
    Lfac = 1e6
    Lunit = r"$\mu$m"
    dxc = config["L"]["c"]/Nvol["c"]
    dxvec = np.array(Nvol["c"] * [dxc])
    porosvec = np.array(Nvol["c"] * [config["poros"]["c"]])
    cellsvec = dxc*np.arange(Nvol["c"]) + dxc/2.
    cellsvec_c = dxc*np.arange(Nvol["c"]) + dxc/2.
    cellsvec_c *= config["L_ref"] * Lfac
    if Nvol["s"]:
        dxs = config["L"]["s"]/Nvol["s"]
        dxvec_s = np.array(Nvol["s"] * [dxs])
        dxvec = np.hstack((dxvec_s, dxvec))
        poros_s = np.array(Nvol["s"] * [config["poros"]["s"]])
        porosvec = np.hstack((poros_s, porosvec))
        cellsvec += config["L"]["s"] / config["L"]["c"]
        cellsvec_s = dxs*np.arange(Nvol["s"]) + dxs/2.
        cellsvec = np.hstack((cellsvec_s, cellsvec))
    if "a" in trodes:
        dxa = config["L"]["a"]/Nvol["a"]
        dxvec_a = np.array(Nvol["a"] * [dxa])
        dxvec = np.hstack((dxvec_a, dxvec))
        poros_a = np.array(Nvol["a"] * [config["poros"]["a"]])
        porosvec = np.hstack((poros_a, porosvec))
        cellsvec += config["L"]["a"] / config["L"]["c"]
        cellsvec_a = dxa*np.arange(Nvol["a"]) + dxa/2.
        cellsvec = np.hstack((cellsvec_a, cellsvec))
    cellsvec_s *= config["L_ref"] * Lfac
    cellsvec *= config["L_ref"] * Lfac
    facesvec_s = np.insert(np.cumsum(dxvec_s), 0, 0.) * config["L_ref"] * Lfac
    facesvec = np.insert(np.cumsum(dxvec), 0, 0.) * config["L_ref"] * Lfac
    
    # Uncomment these lines to show the variables you can plot
    # for k in sim_output.keys():
    #     print(k)
    # return

    trode = 'c'
    t = np.argmin(np.abs(times - tval))
    if tval > times[-1]:
        raise ValueError(f'Max time: {times[-1]}, got {tval}')
    
    datay_c = np.hstack((sim_output['c_lyte_s'], sim_output['c_lyte_c']))
    datay_p = np.hstack((sim_output['phi_lyte_s'], sim_output['phi_lyte_c']))
    
    cGP_L, pGP_L = sim_output['c_lyteGP_L'], sim_output['phi_lyteGP_L']
    cmat = np.hstack((cGP_L.reshape((-1,1)), datay_c, datay_c[:,-1].reshape((-1,1))))
    pmat = np.hstack((pGP_L.reshape((-1,1)), datay_p, datay_p[:,-1].reshape((-1,1))))
    disc = geom.get_elyte_disc(Nvol, config["L"], config["poros"], config["BruggExp"])
    numtimes = len(times)
    i_edges = np.zeros((numtimes, np.shape(cmat)[1]-1))
    
    for tInd in range(numtimes):
        i_edges[tInd, :] = mod_cell.get_lyte_internal_fluxes(cmat[tInd, :], pmat[tInd, :], disc, config)[1]
        
    ylbl = r'Current density of interface [A/m$^2$]'
    datax = range(np.shape(cmat)[1]-1)
    datay = i_edges * (F*c_ref*config["D_ref"]/config["L_ref"])
    #datay = i_edges
    
    fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(8, 4))
    
    ax = axes[0]
    ax.plot(times,datay[:,1], marker='o')
    ax.axvline(tval, ls='--', c='k')
    ax.set_xlabel('time [s]')
    ax.set_ylabel(ylbl)
    
    ax = axes[1]
    #ax.plot(datax,datay[t,:], marker='o')
    ax.plot(facesvec,datay[t,:], marker='o')
    ax.axvline(facesvec_s[-1], ls='--', c='g')
    ax.set_xlabel('Interface position [um]')
    
    #total_current_i = np.zeros(5)
    #total_current_i = datay[t,:]
    #print(total_current_i)
    
    fig.suptitle('Current density')
    fig.tight_layout()
    
    return sim_output

import matplotlib.animation
import matplotlib.pyplot as plt
import numpy as np

def animate_c(resultDir, vol, part, var, tval):
    matfile = osp.join(resultDir, 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    times = sim_output['phi_applied_times'][0]
    
    # Uncomment these lines to show the variables you can plot
    # for k in sim_output.keys():
    #     print(k)
    # return

    trode = 'c'
    t = np.argmin(np.abs(times - tval))
    if tval > times[-1]:
        raise ValueError(f'Max time: {times[-1]}, got {tval}')
    
    var_names = {'c': 'concentration',
                 'phi': 'potential'}
    var_long = var_names[var]

    # elyte data elyte
    key = f'{var}_lyte_s'
    data_elyte_s = sim_output[key]  # select correct elyte volume
    x_elyte_s = list(range(-np.shape(data_elyte_s)[1],0))
    # elyte value at timestep chosen to plot interface region volumes
    cval_s = data_elyte_s[t]
    data_elyte_vol_s = sim_output[key][t]  # all volumes at given t
    
    # elyte data electrode
    key = f'{var}_lyte_{trode}'
    data_elyte_c = sim_output[key]  # select correct elyte volume
    x_elyte_c = list(range(np.shape(data_elyte_c)[1]))
    # elyte value at timestep chosen to plot interface region volumes
    cval_c = data_elyte_c[t]
    data_elyte_vol_c = sim_output[key][t]  # all volumes at given t

    #  volume vs variable
    plt.plot(x_elyte_c,cval_c, marker='o', label='elyte_c')
    plt.plot(x_elyte_s,cval_s , label='elyte_s', c='g', marker='o')
#     plt.set_xlabel('volume')
#     plt.legend()
#     plt.set_title('Interface volumes')
#     ax.set_ylim(-2, 2)

    fig.suptitle(var_long)
    fig.tight_layout()
    
    return sim_output

plt.rcParams["animation.html"] = "jshtml"
plt.rcParams['figure.dpi'] = 150  
plt.ioff()
fig, ax = plt.subplots()

vol = 0  # volume to plot
part = 0  # particle to plot
var = "c"

def animate(tval):
    plt.cla()
    sim_output = animate_c(resultDir, vol, part, var, tval)
    plt.xlim(-10,10)

matplotlib.animation.FuncAnimation(fig, animate, frames=100)
