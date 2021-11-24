import os
from subprocess import run
import subprocess

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

def change_Crate(file_name, new_Crate):
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
    line[index+2] = '   ('+str(new_Crate)+','+str((30/new_Crate))+'),\n'
    file_out = open(file_name,"w")
    file_out.writelines(line)
    file_out.close()
    return

def run_MPET():
    subprocess.call(["python", "bin\mpetrun.py", "configs\params_system.cfg"])



Crates_vec = [0.01,0.05,0.1,1,5]
stddev_vec = [1,10,100,200]

file_name = "params_system.cfg"

for stddev in stddev_vec:
    for Crate in Crates_vec:
        os.chdir(r"C:\Users\pierfrancescoo\Documents\Phase-field\mpet0_1_7\configs")
        change_stddev_c(file_name,stddev)
        change_Crate(file_name,Crate)
        os.chdir(r"C:\Users\pierfrancescoo\Documents\Phase-field\mpet0_1_7")
        run_MPET()
        os.chdir(r"C:\Users\pierfrancescoo\Documents\Phase-field\mpet0_1_7\history")
        old_name = os.listdir(r"C:\Users\pierfrancescoo\Documents\Phase-field\mpet0_1_7\history")[0]
        new_name = 'Stddev_'+str(stddev)+'_Crate_'+str(Crate)
        os.rename(old_name,new_name)
        print(os.listdir(r"C:\Users\pierfrancescoo\Documents\Phase-field\mpet0_1_7\history")[0])


    
