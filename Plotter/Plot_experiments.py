from cProfile import label
from cv2 import polarToCart
from nbformat import read
import numpy as np
import os
import matplotlib.pyplot as plt

def plot_experiments(file_name, mass):
    current_dir = os.getcwd()
    print(current_dir)
    os.chdir(current_dir)
    os.chdir(r"C:\Users\pierfrancescoo\Documents\Experiments\Pyton_plot")
    # file_name = 'VoltageData1.txt' 
    

    file = open(file_name,"r")
    contents = file.readlines()[2:-1]


    vol = np.array([])
    cap = np.array([])
    times = np.array([])

    voltage_vec1 = np.array([])
    cap_vec1 = np.array([])
    times1 = np.array([])

    voltage_vec2 = np.array([])
    cap_vec2 = np.array([])
    times2 = np.array([])

    voltage_vec3 = np.array([])
    cap_vec3 = np.array([])
    times3 = np.array([])

    # voltage_vec4 = np.array([])
    # times4 = np.array([])

    # voltage_vec5 = np.array([])
    # times5 = np.array([])

    # voltage_vec6 = np.array([])
    # times6 = np.array([])


    # v_dis1 = np.array([])
    # timedis1 = np.array([])

    # v_dis2 = np.array([])
    # timedis2 = np.array([])

    # v_dis3 = np.array([])
    # timedis3 = np.array([])

    for line in contents:
        reading = line.split()
        cycle = float(reading[1]) 
        step = float(reading[2])
        vol_general = float(reading[8])
        time_general = float(reading[3])
        cap_general = float(reading[5])
        vol = np.append(vol,vol_general)
        times = np.append(times,time_general)
        cap = np.append(cap,cap_general)
        if reading[9] == 'D':
            if cycle == 2:
                voltage1 = float(reading[8])
                voltage_vec1 = np.append(voltage_vec1,voltage1)
                cap1 = float(reading[5])
                cap_vec1 = np.append(cap_vec1,cap1)
                time1 = float(reading[3])
                times1 = np.append(times1,time1)
            if cycle == 4:
                voltage = float(reading[8])
                voltage_vec2 = np.append(voltage_vec2,voltage)
                cap2 = float(reading[5])
                cap_vec2 = np.append(cap_vec2,cap2)
                time = float(reading[3])
                times2 = np.append(times2,time)
            if cycle == 6:
                voltage3 = float(reading[8])
                voltage_vec3 = np.append(voltage_vec3,voltage3)
                cap3 = float(reading[5])
                cap_vec3 = np.append(cap_vec3,cap3)
                time3 = float(reading[3])
                times3 = np.append(times3,time3)

    tot_cap = 2.45
    tot_cap = 150*mass*1e-3
    plt.plot(cap_vec1*1e3/tot_cap, voltage_vec1, label = '0.1C', color = 'blue', linestyle='dashed' )
    plt.plot(cap_vec2*1e3/tot_cap, voltage_vec2, label = '2C', color = 'green', linestyle='dashed')
    plt.plot(cap_vec3*1e3/tot_cap, voltage_vec3, label = '4C', color = 'red', linestyle='dashed')   
    plt.legend()

    # ax[0].plot(cap_vec1*1e3/tot_cap, voltage_vec1, label = '0.1C' )
    # ax[0].plot(cap_vec2*1e3/tot_cap, voltage_vec2, label = '2C')
    # ax[0].plot(cap_vec3*1e3/tot_cap, voltage_vec3, label = '4C')
    # # plt.xlabel('SoC')
    # plt.xlabel('Voltage/V')
    # plt.legend()
    # plt.show()
        # if cycle == 4:
        #     voltage4 = float(reading[8])
        #     voltage_vec4 = np.append(voltage_vec4,voltage4)
        #     time4 = float(reading[3])
        #     times4 = np.append(times4,time4)
        # if cycle == 5:
        #     voltage5 = float(reading[8])
        #     voltage_vec5 = np.append(voltage_vec5,voltage5)
        #     time5 = float(reading[3])
        #     times5 = np.append(times5,time5)
        # # if cycle == 6:
        #     voltage6 = float(reading[8])
        #     voltage_vec6 = np.append(voltage_vec6,voltage6)
        #     time6 = float(reading[3])
        #     times6 = np.append(times6,time6)
        # if step == 8 or step == 9:
        #     volt = float(reading[8])
        #     v_dis1 = np.append(v_dis1,volt)
        #     t = float(reading[3])
        #     timedis1 = np.append(timedis1,t)
        # if step == 18 or step == 19:
        #     volt = float(reading[8])
        #     v_dis2 = np.append(v_dis2,volt)
        #     t = float(reading[3])
        #     timedis2 = np.append(timedis2,t)
        # if step == 28 or step == 29:
        #     volt = float(reading[8])
        #     v_dis3 = np.append(v_dis3,volt)
        #     t = float(reading[3])
        #     timedis3 = np.append(timedis3,t)


    # t = np.linspace(0,len(voltage_vec),num = len(voltage_vec))
    # init1 = int(round(2828,0))
    # end1 = int(round(2890,0))
    # diff1 = end1 - init1
    # times1 = times[init1:end1]-times[init1]
    # volt1 = voltage_vec[init1:end1]


    # init2 = int(round(4695,0))
    # end2 = int(round(4775,0))
    # diff2 = end2 - init2
    # times2 = times[init2:end2]-times[init2]
    # volt2 = voltage_vec[init2:end2]


    # init3 = int(round(6677,0))
    # end3 = int(round(6729,0))
    # diff3 = end3 - init3
    # times3 = times[init3:end3]-times[init3]
    # volt3 = voltage_vec[init3:end3]

    # plt.plot(times,voltage_vec)
    # start_dis = 40000

    # # plt.plot(times,vol)
    # plt.plot(timedis1-timedis1[0]-18000,v_dis1, label = '0.1C' )
    # plt.plot(timedis2-timedis2[0]-1800,v_dis2, label = '1C')
    # plt.plot(timedis3-timedis3[0]-445,v_dis3, label = '5C')
    # plt.plot(times2-times2[0]-start_dis-18000,voltage_vec2,label = '1')
    # plt.plot(times4-times4[0]-start_dis-1400,voltage_vec4, label = '2')
    # plt.plot(times6-times6[0]-start_dis,voltage_vec6, label = '3')
    # plt.legend()
    # plt.show()
