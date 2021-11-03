#@Tammo#

#preliminary script that set different parameters in loop

import subprocess
import time
import os
import shutil

# Which parameters to change

word = "poros_s" 
values = [0.7, 0.8, 0.9]

#original parameter file:
param_file = "params_system.cfg"

#make new params file and a folder to group the results in:
new_params = "tmp_params.cfg"
#new_dir = "%s\\range_%s_%s" % (os.getcwd(), word, time.strftime("%Y%m%d_%H%M%S", time.localtime()))
#os.mkdir(new_dir) 

command = "python mpetrun.py %s" % new_params

j = 0
for i in values:
    infile = open(param_file)
    outfile = open(new_params,'w')
    #find the line to change
    for line in infile:
        line_words = line.split(' ')
        if line_words[0] == word:
            #replace original value with new value
            line_words[2] = str(i)
            outfile.write(" ".join(line_words))
            outfile.write("\n")
        else: #if line does not need to change:
            outfile.write(line)
    #close the files:
    infile.close()    
    outfile.close()

    #new_dest = "%s\\%s_%s" % (new_dir,word,i) #change the name of the results 
    #print("Simulation for", word, "with value", str(i), "started.")
    subprocess.call(command)

#script finished
os.remove(new_params)