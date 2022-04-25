import os.path as osp
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

import mpet.mod_cell as mod_cell
import mpet.geometry as geom
import mpet.props_am as props_am
import mpet.utils as utils
from plot_utils import *

from mpet.config import Config, constants
from muFunc import *

def plot_activeParticles(resultDir_dic):
    fig, ax = plt.subplots(1,2, sharey=True, figsize=(12, 6))
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(i)
        ffvec = sim_output['ffrac_c'][0]
        td = config["t_ref"]
        times = sim_output['phi_applied_times'][0]*td
        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        tot_particles = Nvol_c*Npart_c

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
        avg_per_active = np.average(percent_active)
        avg_per_active_vec = np.empty(len(percent_active))
        avg_per_active_vec[:] = avg_per_active

        ax[0].plot(ffvec,percent_active,label = str(i[70:]))
        ax[0].set_xlabel('ffrac')
        ax[0].set_ylabel(" %/ of active particles")
        ax[0].legend()

        ax[1].plot(times,percent_active,label = str(i[70:]))
        ax[1].plot(times, avg_per_active_vec)
        ax[1].set_xlabel('t')
        ax[1].set_ylabel(" %/ of active particles")
        ax[1].legend()

def active_atSoC(resultDir_dic,SoC):
    Crate_vec = np.array([])
    percent_vec = np.array([])
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(i)
        ffvec = sim_output['ffrac_c'][0]
        td = config["t_ref"]
        times = sim_output['phi_applied_times'][0]*td
        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        tot_particles = Nvol_c*Npart_c

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
        state_of_charge = int((SoC/100)*np.size(percent_active))
        percent_at_soc = percent_active[state_of_charge]
        Crate = get_C_rate(i)[0]
        bulkcon = get_bulkCon(i, 'c')
        Crate_vec = np.append(Crate_vec,Crate)
        percent_vec = np.append(percent_vec,percent_at_soc)

    plt.plot(Crate_vec, percent_vec, 'o', color = 'green')
    plt.xlabel('Crate')
    plt.ylabel("active particles at 50%")
    plt.show()

def active_max_vs_Crate(resultDir_dic):
    Crate_vec = np.array([])
    percent_vec = np.array([])
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(i)
        ffvec = sim_output['ffrac_c'][0]
        td = config["t_ref"]
        times = sim_output['phi_applied_times'][0]*td
        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        tot_particles = Nvol_c*Npart_c

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
        Crate = get_C_rate(i)
        bulkcon = get_bulkCon(i, 'c')
        Crate_vec = np.append(Crate_vec,Crate)
        percent_vec = np.append(percent_vec,percent_max)

    plt.plot(Crate_vec, percent_vec, color = 'green')
    plt.xlabel("Crate")
    plt.ylabel("Max active particles")
    plt.show()

def active_max_vs_bulk(resultDir_dic):
    bulk_vec = np.array([])
    percent_vec = np.array([])
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(i)
        ffvec = sim_output['ffrac_c'][0]
        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        tot_particles = Nvol_c*Npart_c

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
        # Crate = get_C_rate(i)
        bulkcon = get_bulkCon(i, 'c')
        bulk_vec = np.append(bulk_vec,bulkcon)
        percent_vec = np.append(percent_vec,percent_max)

    plt.plot(bulk_vec, percent_vec, color = 'green')
    plt.xlabel("Electrical conductivity S/m")
    plt.ylabel("Max active particles")
    plt.show()

def active_max_vs_thickness(resultDir_dic):
    thick_vec = np.array([])
    percent_vec = np.array([])
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(i)
        ffvec = sim_output['ffrac_c'][0]
        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        tot_particles = Nvol_c*Npart_c
        L_c = config["L"]["c"]*config['L_ref']*1e6

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
        thick_vec = np.append(thick_vec,L_c)
        percent_vec = np.append(percent_vec,percent_max)

    plt.plot(thick_vec, percent_vec, 'o', color = 'green')
    plt.xlabel("Thickness um")
    plt.ylabel("Max active particles")
    plt.show()




def plot_cVolume(resultDir_dic, SoC):
    for i in resultDir_dic.values():
        matfile = osp.join(i, 'output_data.mat')
        sim_output = sio.loadmat(matfile)
        config = Config.from_dicts(i)
        L_c = config["L"]["c"]*config['L_ref']*1e6
        # td = config["t_ref"]
        ffvec = sim_output['ffrac_c'][0]
        Nvol_c = config["Nvol"]['c']
        Npart_c = config["Npart"]['c']
        thick_vec = np.linspace(L_c,0,num=Nvol_c)

        index = 0
        for step in ffvec:
            index = index + 1
            if step > (SoC-0.05) and step < (SoC+0.05):
                break

        c_avg_tot = np.array([])
        for k in range(Nvol_c):
            c_vol = np.array([])
            for j in range(Npart_c):
                partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_cbar'
                cbar_ffrac = sim_output[partStr][0][index]
                c_vol = np.append(c_vol,cbar_ffrac)
            c_avg_vol = np.average(c_vol)
            c_avg_tot = np.append(c_avg_tot,c_avg_vol)
        plt.plot(thick_vec,c_avg_tot,label = str(i[77:])+' - '+'SoC: '+str(SoC))
        plt.xlabel('Cathode thickness (um)')
        plt.ylabel('Normalized concentration')

        fit = np.poly1d(np.polyfit(thick_vec,c_avg_tot,deg=2),r=False)
        x_fit = np.linspace(0,L_c,num = 20)
        y_fit = np.array([])
        for t in x_fit:
            y_fit = np.append(y_fit,fit(t))
        plt.plot(x_fit,y_fit, label = str(i[50:])+' - '+'SoC: '+str(SoC))
        plt.legend()

            



def plot_cbar(resultDir_dic, xaxis):
    matfile = osp.join(list(resultDir_dic.values())[0], 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    config = Config.from_dicts(list(resultDir_dic.values())[0])

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

def plot_Crate_singleparticle(resultDir_dic, xaxis, Cthreshold):
    matfile = osp.join(list(resultDir_dic.values())[0], 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    config = Config.from_dicts(list(resultDir_dic.values())[0])

    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    print(Npart_c)
    H = Npart_c * 3
    if H > 20:
        H = 20
    L = Nvol_c * 3
    if L > 20:
        L = 20


    fig, ax = plt.subplots(Npart_c,Nvol_c, sharey=True, figsize=(L, H),squeeze=False)

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
                crate = dcbardt*3600/td

                # ax[j,k].plot(xax, crate, label=str(i[77:]) )
                ax[j,k].plot(xax, crate, label=str(i[70:]) )
                ax[j,k].set_ylabel('C-rate')
                ax[j,k].set_xlabel(xlabel)
                ax[j,k].set_title('volume: ' +str(k+1)+' particle: '+str(j+1))
                ax[j,k].legend()

                # active_Crate = np.array([])
                # index = 0
                # for c in crate:
                #     if c > Cthreshold:
                #         active_Crate = np.append(active_Crate,c)
                # avg_Crate = np.average(active_Crate)
                # for c in crate: #find position in x vec when the (de)lithiation start
                #     index = index + 1
                #     if c > Cthreshold:
                #         break
                # new_xax = np.linspace(xax[0],xax[-1],num=np.size(active_Crate))
                # avg_Crate_vec = np.array([])
                # for g in active_Crate:
                #     avg_Crate_vec = np.append(avg_Crate_vec,avg_Crate)
                # ax[j,k].plot(new_xax, avg_Crate_vec, label=str(i[77:]))
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

    fig, axes = plt.subplots(Npart_c,Nvol_c, sharey=True, figsize=(10, 4),squeeze=False)
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
                    mubar = 0.025*np.append(mubar, np.sum(mu)/len(mu))
                # if Npart_c == 1:
                #     ax = axes[k]
                # elif Nvol_c == 1:
                #     ax = axes[j]
                # else:
                ax = axes[j,k]

                ax.plot(xax, mubar, label='mubar')
                ax.set_xlabel(xlabel)
                ax.set_ylabel('mubar')


    return sim_output


def plot_csurf_vs_cbar(resultDir_dic):

    config = Config.from_dicts(list(resultDir_dic.values())[0])
    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    H = Npart_c * 3
    if H > 20:
        H = 20
    L = 1 * 3
    if L > 20:
        L = 20

    fig, axes = plt.subplots(Npart_c,Nvol_c, sharey=True, figsize=(10, 4),squeeze=False)
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

                csurf_vec = np.empty([0])
                for t in range(tottimesteps+1):

                    csurf = c[t][-1]
                    cavg = cbar[t]
                    csurf_vec = np.append(csurf_vec, (csurf-cavg))
                ax = axes[j,k]

                ax.plot(cbar, csurf_vec, label='_csurf')
                ax.set_xlabel('cbar')
                ax.set_ylabel('(c_surf - cvg)')
                ax.set_title('volume: ' +str(k+1)+' particle: '+str(j+1))


    return sim_output


def plot_ecd_vs_cbar(resultDir_dic):

    config = Config.from_dicts(list(resultDir_dic.values())[0])
    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    H = Npart_c * 3
    if H > 20:
        H = 20
    L = 1 * 3
    if L > 20:
        L = 20

    fig, axes = plt.subplots(Npart_c,Nvol_c, sharey=True, figsize=(10, 4),squeeze=False)
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
        Omega_a = config["c","Omega_a"][0][0]
        Omega_b = config["c","Omega_b"]
        Choe = config["c","B"]

        tottimesteps = len(times)

        for k in range(Nvol_c):
            for j in range(Npart_c):

                partStr = 'partTrodecvol'+str(k)+'part'+str(j)+'_c'
                c = sim_output[partStr]
      

                partStr_bar = 'partTrodecvol'+str(k)+'part'+str(j)+'_cbar'
                cbar = sim_output[partStr_bar][0]


                ecd_vec = np.empty([0])
                for t in range(tottimesteps):

                    csurf = c[t][-1]
                    cavg = cbar[t]
                    ecd1 = np.sqrt(csurf*(1-csurf))
                    ecd2 = np.exp(0.5*(Omega_a*(1-2*csurf)))
                    ecd3 = np.exp(0.5*(Omega_b*((1-2*csurf)**2-2*csurf*(1-csurf))))
                    ecd4 = np.exp(0.5*Choe*(csurf-cavg))
                    # ecd4 = 1
    
                    
                    ecd = ecd1*ecd2*ecd3*ecd4
                    ecd_vec = np.append(ecd_vec, ecd)
                ax = axes[j,k]

                ax.plot(cbar, ecd_vec, label='_csurf')
                ax.set_xlabel('cbar')
                ax.set_ylabel('ecd')
                ax.set_title('volume: ' +str(k+1)+' particle: '+str(j+1))


    return sim_output


def plot_mubar_vs_cbar(resultDir_dic):

    config = Config.from_dicts(list(resultDir_dic.values())[0])
    Nvol_c = config["Nvol"]['c']
    Npart_c = config["Npart"]['c']

    H = Npart_c * 3
    if H > 20:
        H = 20
    L = 1 * 3
    if L > 20:
        L = 20

    fig, axes = plt.subplots(Npart_c,Nvol_c, sharey=True, figsize=(10, 4),squeeze=False)
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
                    mubar = 0.025*np.append(mubar, np.sum(mu)/len(mu))
                # if Npart_c == 1:
                #     ax = axes[k]
                # elif Nvol_c == 1:
                #     ax = axes[j]
                # else:
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
        ax.plot(cbar, 0.025*mubar, label='mubar')
        ax.set_xlabel('cbar')
        ax.set_ylabel('mubar')
        ax.set_title('volume: ' +str(vol+1)+' particle: '+str(part+1))

        ax = axes[1]
        ax.plot(times, 0.025*mubar, label='mubar')
        ax.set_xlabel('time')
        ax.set_ylabel('mubar')


    return sim_output