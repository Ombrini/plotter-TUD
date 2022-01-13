import os
from subprocess import run
import subprocess
import numpy as np
from numpy.core.fromnumeric import std
from numpy.lib.function_base import sinc

def change_stddev_c(file_name,new_stddev):
    new_stddev = str(new_stddev*1e-9)
    print("changing stddev to ", new_stddev)
    file = open(file_name,"r")
    line = file.readlines()

    index = -1
    # Loop through the file line by line
    for i_line in line:  
        index = index + 1
        # checking string is present in line or not
        if '# wasting computational effort.\n'  == i_line:
            break 
    line[index+2] = 'stddev_c = '+new_stddev+'\n'
    file_out = open(file_name,"w")
    file_out.writelines(line)
    file_out.close()
    return

def change_mean_c(file_name,new_mean_c):
    new_mean_c = str(new_mean_c*1e-9)
    print("changing mean to ", new_mean_c)
    file = open(file_name,"r")
    line = file.readlines()

    index = -1
    # Loop through the file line by line
    for i_line in line:  
        index = index + 1
        # checking string is present in line or not
        if '# wasting computational effort.\n'  == i_line:
            break 
    line[index+1] = 'mean_c = '+new_mean_c+'\n'
    file_out = open(file_name,"w")
    file_out.writelines(line)
    file_out.close()
    return

def change_thick_c(file_name,new_thick_c):
    new_thick_c = str(new_thick_c*1e-6)
    print("changing thickenss to ", new_thick_c)
    file = open(file_name,"r")
    line = file.readlines()

    index = -1
    # Loop through the file line by line
    for i_line in line:  
        index = index + 1
        # checking string is present in line or not
        if '# Thicknesses, m\n'  == i_line:
            break 
    line[index+1] = 'L_c = '+new_thick_c+'\n'
    file_out = open(file_name,"w")
    file_out.writelines(line)
    file_out.close()
    return

def change_initial_part_conc(file_name,initial_part_c):
    new_initial_part_c = str(initial_part_c)
    print("changing initial part conc to ", new_initial_part_c)
    file = open(file_name,"r")
    line = file.readlines()

    index = -1
    # Loop through the file line by line
    for i_line in line:  
        index = index + 1
        # checking string is present in line or not
        if '# (for disch, anode starts full, cathode starts empty)\n'  == i_line:
            break 
    line[index+1] = 'cs0_c = '+new_initial_part_c+'\n'
    file_out = open(file_name,"w")
    file_out.writelines(line)
    file_out.close()
    return

def change_bulk_conductivity(file_name,bulk_cond):
    new_bulk_cond = str(bulk_cond)
    print("changing bulk cond to ", new_bulk_cond)
    file = open(file_name,"r")
    line = file.readlines()

    index = -1
    # Loop through the file line by line
    for i_line in line:  
        index = index + 1
        # checking string is present in line or not
        if '# Dimensional conductivity (used if simBulkCond = true), S/m\n'  == i_line:
            break 
    line[index+1] = 'sigma_s_c = ' +new_bulk_cond+'\n'
    file_out = open(file_name,"w")
    file_out.writelines(line)
    file_out.close()
    return

def change_G_mean_c(file_name,cond):
    new_cond = str(cond)
    print("changing G mean c to ", new_cond)
    file = open(file_name,"r")
    line = file.readlines()

    index = -1
    # Loop through the file line by line
    for i_line in line:  
        index = index + 1
        # checking string is present in line or not
        if '# Conductance between particles, S = 1/Ohm\n'  == i_line:
            break 
    line[index+1] = 'G_mean_c = ' +new_cond+'\n'
    line[index+2] = 'G_stddev_c = ' +str(0.05*cond)+'\n'
    file_out = open(file_name,"w")
    file_out.writelines(line)
    file_out.close()
    return


def change_CCsegments(file_name, new_Crate):
    print("changing Crate to ", new_Crate)
    file = open(file_name,"r")
    line = file.readlines()

    index = -1
    # Loop through the file line by line
    for i_line in line:  
        index = index + 1
        # checking string is present in line or not
        if '# Note: It\'s okay to leave commented lines within the segments list\n'  == i_line:
            break 
    line[index+2] = '   ('+str(new_Crate)+','+str(abs(30/new_Crate))+'),\n'
    file_out = open(file_name,"w")
    file_out.writelines(line)
    file_out.close()
    return

def change_CCsegments_complete(file_name, WritingCrate, RestingTime, ReadingCrate, SoC_Stop):
    print('Writing at ',WritingCrate,'C , Waiting for ', RestingTime, ' minutes, Reading at ', ReadingCrate,'C ', 'Stopping at ', str(SoC_Stop),'%')
    file = open(file_name,"r")
    line = file.readlines()
    timestop = 0.6*SoC_Stop
    index = -1
    # Loop through the file line by line
    for i_line in line:  
        index = index + 1
        # checking string is present in line or not
        if '# Note: It\'s okay to leave commented lines within the segments list\n'  == i_line:
            break 
    line[index+2] = '   ('+str(WritingCrate)+','+str(abs(timestop/WritingCrate))+'),\n'
    line[index+3] = '   ('+str(0)+','+str(abs(RestingTime))+'),\n'
    line[index+4] = '   ('+str(ReadingCrate)+','+str(abs(0.8*(60-timestop)/ReadingCrate))+'),\n'
    file_out = open(file_name,"w")
    file_out.writelines(line)
    file_out.close()
    return


def run_MPET(cwd):
    os.chdir(cwd)
    subprocess.call(["python", "bin\mpetrun.py", "configs\params_system.cfg"])


currentDir = os.getcwd()
# insert the path of your params_system
os.chdir(currentDir)
os.chdir(r"C:\Users\pierfrancescoo\Documents\Phase-field\mpet-dev-feature-solid\configs")
file_name = "params_system.cfg"
Sign = 1 #1 for disch, -1 for charging
if Sign == 1:
    change_initial_part_conc(file_name,0.01)
elif Sign == -1:
    change_initial_part_conc(file_name,0.97)

os.chdir(r"C:\Users\pierfrancescoo\Documents\Phase-field\mpet-dev-feature-solid")
# put the parameters or the array of parameters you want to loop
Crates_vec = Sign*np.array([1,5,10])
thickness_vec = [1]
mean_c = [140]
stddev_vec = [10]
bulk_cond_vec = [0.25]
G_mean_c_vec = [1e-5]
Soc_stop = [99] #stopping at x% SoC
ReadingCrate = Sign*1 #using 5C to read after the stop 
restingTime_vec = [0] # resting for x minutes 

for SoC in Soc_stop:
    for thick in thickness_vec:
        for stddev in stddev_vec:
            for bulk in bulk_cond_vec:
                for G_mean in G_mean_c_vec:
                    for Crate in Crates_vec:
                        for restingTime in restingTime_vec:
                            for mean in mean_c:
                                os.chdir(r"C:\Users\pierfrancescoo\Documents\Phase-field\mpet-dev-feature-solid")
                                cwd = os.getcwd()
                                os.chdir(cwd+"\configs")
                                change_thick_c(file_name,thick)
                                change_stddev_c(file_name,stddev)
                                change_bulk_conductivity(file_name, bulk)
                                change_G_mean_c(file_name,G_mean)
                                change_CCsegments_complete(file_name,Crate,restingTime,ReadingCrate,SoC)
                                change_mean_c(file_name,mean)

                                
                                run_MPET(cwd)
                                os.chdir(cwd+"\history")

                                old_name = os.listdir(cwd+"\history")[0]
                                new_name = 'Thick_'+str(thick)+'_Stddev_'+str(stddev)+'_SoCStop_'+str(SoC)+'_Crate_'+str(Crate)+'_bulkCon_'+str(bulk)+'_Gmean_'+str(G_mean)+'_MeanC_'+str(mean)
                                os.rename(old_name,new_name)
                                print('New folder: ', os.listdir(cwd+"\history")[0])


            
