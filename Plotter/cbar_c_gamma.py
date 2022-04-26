import os.path as osp
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import os

import matplotlib as mpl
import matplotlib.animation as manim
import matplotlib.collections as mcollect
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches

import mpet.mod_cell as mod_cell
import mpet.geometry as geom
import mpet.props_am as props_am
import mpet.utils as utils

from mpet.config import Config, constants

def cbar_c(resultDir_dic):
    pfx = 'mpet.'
    sStr = "_"
    ttl_fmt = "% = {perc:2.1f}"
    matfile = osp.join(resultDir_dic["sim1"], 'output_data.mat')
    sim_output = sio.loadmat(matfile)
    data = sim_output
    config = Config.from_dicts(resultDir_dic["sim1"])
    trodes = config["trodes"]

    Nvol = config["Nvol"]
    Npart = config["Npart"]
    td = config["t_ref"]
    psd_len = config["psd_len"]
    gamma_cont = config["gamma_contact"]

    times = sim_output['phi_applied_times'][0]*td
    numtimes = len(times)
    tmin = np.min(times)
    tmax = np.max(times)

    plot_type = "cbar_c"

    if plot_type in ["cbar_c", "cbar_a", "cbar_full"]:
        if plot_type[-4:] == "full":
            trvec = ["a", "c"]
        elif plot_type[-1] == "a":
            trvec = ["a"]
        else:
            trvec = ["c"]
        dataCbar = {}
        for trode in trodes:
            dataCbar[trode] = np.zeros((numtimes, Nvol[trode], Npart[trode]))
            for tInd in range(numtimes):
                for vInd in range(Nvol[trode]):
                    for pInd in range(Npart[trode]):
                        # dataStr = (
                        #     pfx
                        #     + "partTrode{t}vol{vInd}part{pInd}".format(
                        #         t=trode, vInd=vInd, pInd=pInd)
                        #     + sStr + "cbar")
                        # dataCbar[trode][tInd,vInd,pInd] = (
                        #     np.squeeze(utils.get_dict_key(data, dataStr))[tInd])
                        dataStr = 'partTrodecvol'+str(vInd)+'part'+str(pInd)+'_cbar'
                        dataCbar[trode][tInd,vInd,pInd] = sim_output[dataStr][0][tInd]
                        
        # Set up colors.
        # Define if you want smooth or discrete color changes
        # Option: "smooth" or "discrete"
        # color_changes = "discrete"
        color_changes = "discrete"
        # Discrete color changes:
        if color_changes == "discrete":
            # Make a discrete colormap that goes from green to yellow
            # to red instantaneously
            to_yellow = 0.3
            to_red = 0.7
            cdict = {
                "red": [(0.0, 0.0, 0.0),
                        (to_yellow, 0.0, 1.0),
                        (1.0, 1.0, 1.0)],
                "green": [(0.0, 0.502, 0.502),
                          (to_yellow, 0.502, 1.0),
                          (to_red, 1.0, 0.0),
                          (1.0, 0.0, 0.0)],
                "blue": [(0.0, 0.0, 0.0),
                         (1.0, 0.0, 0.0)]
                }
            cmap = mpl.colors.LinearSegmentedColormap(
                "discrete", cdict)
        # Smooth colormap changes:
        if color_changes == "smooth":
            # generated with colormap.org
            cmaps = np.load(r"C:\Users\tkschwietert\Documents\Phase_Field\Solid_interface_6\mpet-dev_2\mpet\plot\colormaps_custom.npz")
            cmap_data = cmaps["GnYlRd_3"]
            cmap = mpl.colors.ListedColormap(cmap_data/255.)

        plt.rc('font',size = 14)
        size_frac_min = 0.10
        figsize = (6*1.2, 2*1.2)
        fig, axs = plt.subplots(1, len(trvec), squeeze=False, figsize=figsize)
        ttlx = 0.5 if len(trvec) < 2 else 1.1
        ttl = axs[0,0].text(
            ttlx, 1.05, ttl_fmt.format(perc=0),
            transform=axs[0,0].transAxes, verticalalignment="center",
            horizontalalignment="center")
        collection = np.empty(len(trvec), dtype=object)
        collection2 = np.empty(len(trvec), dtype=object)
        for indx, trode in enumerate(trvec):
            ax = axs[0,indx]
            # Get particle sizes (and max size) (length-based)
            lens = psd_len[trode]
            len_max = np.max(lens)
            len_min = np.min(lens)
            ax.patch.set_facecolor('white')
            # Don't stretch axes to fit figure -- keep 1:1 x:y ratio.
            ax.set_aspect('equal', 'box')
            # Don't show axis ticks
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())
            ax.set_xlim(0, 1.02)
            ax.set_ylim(0, float(Npart[trode])/Nvol[trode]*1.07)
            #contact
            cont = gamma_cont[trode]
            # cont = 1-gamma_cont[trode]
            #
            # Label parts of the figure
#            ylft = ax.text(-0.07, 0.5, "Separator",
#                    transform=ax.transAxes, rotation=90,
#                    verticalalignment="center",
#                    horizontalalignment="center")
#            yrht = ax.text(1.09, 0.5, "Current Collector",
#                    transform=ax.transAxes, rotation=90,
#                    verticalalignment="center",
#                    horizontalalignment="center")
#            xbtm = ax.text(.50, -0.05, "Electrode Depth -->",
#                    transform=ax.transAxes, rotation=0,
#                    verticalalignment="center",
#                    horizontalalignment="center")
            # Geometric parameters for placing the rectangles on the axes
            spacing = 1.0 / Nvol[trode]
            size_fracs = 0.4*np.ones((Nvol[trode], Npart[trode]))
            if len_max != len_min:
                size_fracs = (lens - len_min)/(len_max - len_min)
            sizes = (size_fracs*(1-size_frac_min) + size_frac_min) / Nvol[trode]

            # Create rectangle "patches" to add to figure axes.
            rects = np.empty((Nvol[trode], Npart[trode]), dtype=object)
            rectsd = np.empty((Nvol[trode], Npart[trode]), dtype=object)
            rectsu = np.empty((Nvol[trode], Npart[trode]), dtype=object)
            rectsl = np.empty((Nvol[trode], Npart[trode]), dtype=object)
            rectsr = np.empty((Nvol[trode], Npart[trode]), dtype=object)            
            color = 'green'  # value is irrelevant -- it will be animated
            for (vInd, pInd), c in np.ndenumerate(sizes):
                size = sizes[vInd,pInd]
                center = np.array([spacing*(vInd + 0.5), spacing*(pInd + 0.5)])
                bottom_left = center - size / 4
                rects[vInd,pInd] = plt.Rectangle(
                    bottom_left, size, size, color=color,angle=0)
                lx = bottom_left[0]
                rx = bottom_left[0] + size
                ly = bottom_left[1]
                ry = bottom_left[1]+size

                # plot the the effective contact with a line around particle
                lpt = 1.5 # linewidth
                lc = 'k' #line color

                ##

                cnt = cont[vInd,pInd]
                if cnt > 0.75:
                    rectsd[vInd,pInd] = plt.plot([rx-(rx-lx)*(cnt-0.75)*4,rx], [ly,ly], color= lc ,lw = lpt)

                if 0.5 < cnt < 0.75:
                    rectsr[vInd,pInd] = plt.plot([rx,rx], [ry+(ly-ry)*(cnt-0.5)*4,ry], color= lc ,lw = lpt)
                elif cnt >= 0.75 :
                    rectsr[vInd,pInd] = plt.plot([rx,rx], [ly,ry], color= lc ,lw = lpt)

                if 0.25 < cnt < 0.5:
                    rectsu[vInd,pInd] = plt.plot([lx,lx-(lx-rx)*(cnt-0.25)*4], [ry,ry], color= lc ,lw = lpt)
                elif cnt >= 0.5:
                    rectsu[vInd,pInd] = plt.plot([lx,rx], [ry,ry], color= lc ,lw = lpt)

                if cnt < 0.25:
                    rectsl[vInd,pInd] = plt.plot([lx,lx], [ly,ly-(ly-ry)*(cnt)*4], color= lc ,lw = lpt)
                else:
                    rectsl[vInd,pInd] = plt.plot([lx,lx], [ly,ry], color= lc ,lw = lpt)

            # Create a group of rectange "patches" from the rects array
            collection[indx] = mcollect.PatchCollection(rects.reshape(-1))

            # Put them on the axes
            ax.add_collection(collection[indx])

        # Have a "background" image of rectanges representing the
        # initial state of the system.

        # def init():
        #     for indx, trode in enumerate(trvec):
        #         cbar_mat = dataCbar[trode][0,:,:]
        #         colors = cmap(cbar_mat.reshape(-1))
        #         collection[indx].set_color(colors)
        #         ttl.set_text('')
        #     out = [collection[i] for i in range(len(collection))]
        #     out.append(ttl)
        #     out = tuple(out)
        #     return out

        def animate(tind):
            for indx, trode in enumerate(trvec):
                cbar_mat = dataCbar[trode][tind,:,:]
                colors = cmap(cbar_mat.reshape(-1))
                collection[indx].set_color(colors)
            t_current = times[tind]
            tfrac = (t_current - tmin)/(tmax - tmin) * 100
            ttl.set_text(ttl_fmt.format(perc=tfrac))
            out = [collection[i] for i in range(len(collection))]
            out.append(ttl)
            out = tuple(out)
            return out

        for i in range(numtimes):
            animate(i)
            plt.pause(0.05)
            # if i == round(numtimes*0.25):
            #     plt.savefig('25soc.pdf')
            # if i == round(numtimes*0.5):
            #     plt.savefig('50soc.pdf')
            # if i == numtimes-1:
            #     plt.savefig('100soc.pdf')