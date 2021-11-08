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

def interface_c(resultDir_dic, vol, part, tval):
    var = 'c'
    for i in resultDir_dic.values():
        config = Config.from_dicts(i)
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        trodes = config["trodes"]
        td = config["t_ref"]
        times = sim_output['phi_applied_times'][0]*td

        Nvol = config["Nvol"]
        trodes = config["trodes"]

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
        
        trode = 'c'
        t = np.argmin(np.abs(times - tval))
        if tval > times[-1]:
            raise ValueError(f'Max time: {times[-1]}, got {tval}')
        
        var_names = {'c': 'concentration',
                    'phi': 'potential'}
        var_long = var_names[var]
        
        fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(12, 6))
        
        key = f'interfaceTrode{trode}vol{vol}part{part}_{var}'
        nvol_i = sim_output[key].shape[1]
        
        Lfac = 1e6
        dxc = config["L"]["c"]/Nvol["c"]
        dxvec = np.array(Nvol["c"] * [dxc])
        cellsvec = dxc*np.arange(Nvol["c"]) + dxc/2.
        cellsvec_c = dxc*np.arange(Nvol["c"]) + dxc/2.
        cellsvec_c *=config["L_ref"] * Lfac
        
        dxs = config["L_i"]/nvol_i
        dxvec_i = np.array(nvol_i * [dxs])
        dxvec = np.hstack((dxvec_i, dxvec))
        cellsvec += config["L_i"] / config["L"]["c"]
        cellsvec_i = dxs*np.arange(nvol_i) + dxs/2.
        cellsvec_i *= config["L_ref"] * Lfac

        # elyte data
        key = f'{var}_lyte_{trode}'
        if var in 'c':
            data_elyte = sim_output[key][:, vol]* c_ref / 1000.  # select correct elyte volume
        elif var in 'phi':
            data_elyte = sim_output[key][:, vol]*(k*Tref/e) - Vstd  # select correct elyte volume
        # elyte value at timestep chosen to plot interface region volumes
        cval = data_elyte[t]
        
        # interface data
        key = f'interfaceTrode{trode}vol{vol}part{part}_{var}'
        # select all time steps, middle volume
        ind = -1
        if var in 'c':
            data_iface = sim_output[key][:, ind]* c_ref / 1000.  # select correct elyte volume
            data_iface_vol = sim_output[key][t]* c_ref / 1000.
        elif var in 'phi':
            data_iface = sim_output[key][:, ind]*(k*Tref/e) - Vstd  # select correct elyte volume
            data_iface_vol = sim_output[key][t]*(k*Tref/e) - Vstd  # all volumes at given t
        
        # time vs variable
        ax = axes[0]
        ax.plot(times*td, data_elyte, label='elyte')
        ax.plot(times*td, data_iface, label='iface last volume',c='m')

        ax.set_xlabel('time')
        ax.set_ylabel(var_long)
        ax.axvline(tval*td, ls='--', c='k')
        ax.axhline(cval, ls='--', c='k')
        ax.legend()
        ax.set_title('Time evolution')

        #  volume vs variable
        ax = axes[1]   
        ax.plot(cellsvec_i,data_iface_vol, marker='o', label='iface volumes',c='m')
        ax.scatter(-cellsvec_i[0], cval, label='elyte output')
        print(var_long, "elyte", cval)
        print(var_long, "iface", data_iface_vol)
        ax.set_xlabel('Interface position [um]')
        ax.legend()
        ax.set_title('Interface volumes')

        fig.suptitle(var_long)
        fig.tight_layout()
    
    return sim_output

def interface_phi(resultDir_dic, vol, part, tval):
    var = 'phi'
    for i in resultDir_dic.values():
        config = Config.from_dicts(i)
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        trodes = config["trodes"]
        td = config["t_ref"]
        times = sim_output['phi_applied_times'][0]*td

        Nvol = config["Nvol"]
        trodes = config["trodes"]

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
        
        trode = 'c'
        t = np.argmin(np.abs(times - tval))
        if tval > times[-1]:
            raise ValueError(f'Max time: {times[-1]}, got {tval}')
        
        var_names = {'c': 'concentration',
                    'phi': 'potential'}
        var_long = var_names[var]
        
        fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(12, 6))
        
        key = f'interfaceTrode{trode}vol{vol}part{part}_{var}'
        nvol_i = sim_output[key].shape[1]
        
        Lfac = 1e6
        dxc = config["L"]["c"]/Nvol["c"]
        dxvec = np.array(Nvol["c"] * [dxc])
        cellsvec = dxc*np.arange(Nvol["c"]) + dxc/2.
        cellsvec_c = dxc*np.arange(Nvol["c"]) + dxc/2.
        cellsvec_c *=config["L_ref"] * Lfac
        
        dxs = config["L_i"]/nvol_i
        dxvec_i = np.array(nvol_i * [dxs])
        dxvec = np.hstack((dxvec_i, dxvec))
        cellsvec += config["L_i"] / config["L"]["c"]
        cellsvec_i = dxs*np.arange(nvol_i) + dxs/2.
        cellsvec_i *= config["L_ref"] * Lfac

        # elyte data
        key = f'{var}_lyte_{trode}'
        if var in 'c':
            data_elyte = sim_output[key][:, vol]* c_ref / 1000.  # select correct elyte volume
        elif var in 'phi':
            data_elyte = sim_output[key][:, vol]*(k*Tref/e) - Vstd  # select correct elyte volume
        # elyte value at timestep chosen to plot interface region volumes
        cval = data_elyte[t]
        
        # interface data
        key = f'interfaceTrode{trode}vol{vol}part{part}_{var}'
        # select all time steps, middle volume
        ind = -1
        if var in 'c':
            data_iface = sim_output[key][:, ind]* c_ref / 1000.  # select correct elyte volume
            data_iface_vol = sim_output[key][t]* c_ref / 1000.
        elif var in 'phi':
            data_iface = sim_output[key][:, ind]*(k*Tref/e) - Vstd  # select correct elyte volume
            data_iface_vol = sim_output[key][t]*(k*Tref/e) - Vstd  # all volumes at given t
        
        # time vs variable
        ax = axes[0]
        ax.plot(times*td, data_elyte, label='elyte')
        ax.plot(times*td, data_iface, label='iface last volume',c='m')

        ax.set_xlabel('time')
        ax.set_ylabel(var_long)
        ax.axvline(tval*td, ls='--', c='k')
        ax.axhline(cval, ls='--', c='k')
        ax.legend()
        ax.set_title('Time evolution')

        #  volume vs variable
        ax = axes[1]   
        ax.plot(cellsvec_i,data_iface_vol, marker='o', label='iface volumes',c='m')
        ax.scatter(-cellsvec_i[0], cval, label='elyte output')
        print(var_long, "elyte", cval)
        print(var_long, "iface", data_iface_vol)
        ax.set_xlabel('Interface position [um]')
        ax.legend()
        ax.set_title('Interface volumes')

        fig.suptitle(var_long)
        fig.tight_layout()
    
    return sim_output

def interface_i(resultDir_dic, vol, part, tval):
    for i in resultDir_dic.values():
        config = Config.from_dicts(i)
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        td = config["t_ref"]
        times = sim_output['phi_applied_times'][0]*td
        numtimes = len(times)

        Nvol = config["Nvol"]
        var = 'c'
    
        # Pick out some useful calculated values
        F = constants.F                      # C/mol
        c_ref = constants.c_ref

        trode = 'c'
        t = np.argmin(np.abs(times - tval))
        if tval > times[-1]:
            raise ValueError(f'Max time: {times[-1]}, got {tval}')
        
        key = f'interfaceTrode{trode}vol{vol}part{part}_{var}'
        nvol_i = sim_output[key].shape[1]
        
        Lfac = 1e6       
        dxs = config["L_i"]/nvol_i
        dxvec_i = np.array(nvol_i * [dxs])
        facesvec_i = np.insert(np.cumsum(dxvec_i), 0, 0.) * config["L_ref"] * Lfac

        t = np.argmin(np.abs(times - tval))
        if tval > times[-1]:
            raise ValueError(f'Max time: {times[-1]}, got {tval}')

        Nvol = config["Nvol_i"]
        
        key_c = f'interfaceTrodecvol{vol}part{part}_c'
        key_phi = f'interfaceTrodecvol{vol}part{part}_phi'

        #cGP_L, pGP_L = sim_output['c_lyteGP_L'], sim_output['phi_lyteGP_L']
        c_int, phi_int = sim_output[key_c], sim_output[key_phi]
        cmat = np.hstack((c_int[:, 0].reshape((numtimes, 1)), c_int, c_int[:, -1].reshape((numtimes, 1))))
        pmat = np.hstack((phi_int[:, 0].reshape((numtimes, 1)), phi_int, phi_int[:, -1].reshape((numtimes, 1))))
        disc = geom.get_interface_disc(Nvol, config["L_i"], config["poros_i"], config["BruggExp_i"])

        i_edges = np.zeros((numtimes, np.shape(cmat)[1]-1))

        for tInd in range(numtimes):
            i_edges[tInd, :] = mod_cell.get_lyte_internal_fluxes(cmat[tInd, :], pmat[tInd, :], disc, config)[1]

        ylbl = r'Current density of interface [A/m$^2$]'
        datay = i_edges * (F*c_ref*config["D_ref"]/config["L_ref"])

        fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(12, 6))
        
        ax = axes[0]
        ax.plot(times,datay[:,1], marker='o')
        ax.axvline(tval, ls='--', c='k')
        ax.set_xlabel('time [s]')
        ax.set_ylabel(ylbl)
        
        ax = axes[1]
        ax.plot(facesvec_i,datay[t,:], marker='o')
        ax.set_xlabel('Interface position [um]')
        
        fig.suptitle('Current density')
        fig.tight_layout()
    
    return sim_output

