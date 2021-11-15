#!/usr/bin/env python

import os
import argparse
import ROOT
from ROOT import TFile
import math
import array
import time
from argparse import RawTextHelpFormatter

working_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='''Script that: \n\t-Executes step1_hit_perLumi_analysis.py or step3_matching.py or step4_flowerEvents_analysis.py\n\t-Submits to HTCondor\n''',
            epilog="""Typical execution\n\t python run_step.py -r 342154 --submit True """,
            formatter_class=RawTextHelpFormatter
    )
    parser.add_argument('--run','-r', type=str,help="Comma separated list of Cosmic runs to be analyzed",required=True)
    parser.add_argument('--step','-st', type=int,help="Step to execute:\n - 1: produces mini-tree for dead/hot vfat analysis\n - 3: produces mini-tree with matched muons for efficiency\n - 4: occupancy analysis",required=True, choices=[1,4])
    parser.add_argument('--num_input', '-n', type=int, help="number of input rootfiles per job. [Default: %(default)s] ", action="store", default = 1)
    parser.add_argument('--outputpath', '-o', type=str, help="Path of output files. [Default: %(default)s]",default=working_dir+"/output/")
    parser.add_argument('--execute', '-e', type=int, help="Run code on lxplus session. Set it to -1 to analyse full ntuple, otherwise set it to nevents. [Default: %(default)s]", default=0) 
    parser.add_argument('--submit', '-s', help="Submit to HTCondor. [Default: %(default)s]", default=False) 
    parser.add_argument('--reco', '-rc', help="Express, PromptReco, MinimumBias. Ntuples are stored in different folders. [Default: %(default)s]", default="Express", choices=["Express", "PromptReco", "MinimumBias"]) 
    args = parser.parse_args()

    command = ""

    ## Get list of Cosmic runs
    runList = []
    if args.run != None:
        runList = map(str, args.run.split(','))
    else:
        print "Error: provide run number"
        exit()
    print runList

    script_name = ""
    if args.step==1:
        bash_name = "run_step1_wrap.sh" 
        script_name = "step1_hit_perLumi_analysis.py" 
    elif args.step==3:
        bash_name = "run_step3_wrap.sh"
        script_name = "step3_matching.py"
    elif args.step==4:
        bash_name = "run_step4_wrap.sh"
        script_name = "step4_flowerEvents_analysis.py"
    # for HTCondor submission, prepare file
    if(args.submit):
        with open('condor.submit', 'w') as submit_file:
            submit_file.write(
            '''
universe = vanilla
executable     = {bash_exe}
getenv = True

request_memory = 4096
request_cpus   = 1
request_disk   = 16383

error   = {path}/err_{submit_time}_.$(Process)
output  = {path}/out_{submit_time}_.$(Process)
log     = {path}/log_{submit_time}__$(Process).log
+JobFlavour = "workday"
            '''.format(bash_exe = bash_name, path=args.outputpath, submit_time=str(time.time())))

    campaign = ""
    for run in runList:
        base_path_list = [] ##allowing for the files to be distributed in different folders

        if(int(run)<342090): campaign = "MWGR3"
        elif(int(run)<342200): campaign = "MWGR4"
        elif(int(run)<346100): campaign = "CRUZET"
        elif(int(run)<346246): campaign = "CRAFT"
        else: campaign = "Collisions"

        if(args.reco=="Express"): base_path_list.append("/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEMCommonNtuples/"+campaign+"/Run_"+run+"/")
        if(args.reco=="PromptReco"): base_path_list.append("/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEMCommonNtuples/"+campaign+"/"+run+"_Prompt/")
        if(args.reco=="MinimumBias"):
            for r, d, f in os.walk("/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEM/"):
                for subdir in d:
                    if 'MinimumBias' in subdir and run in subdir:
                        base_path_list.append(os.path.join("/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEM/", r, subdir))

        #generating the list of all .root files in given directory and subdirectories
        fileList = []
        for base_path in base_path_list:
            for r, d, f in os.walk(base_path): # r=root, d=directories, f = files
                for file in f:
                    if '.root' in file:
                        fileList.append(os.path.join(r, file))

        #grouping input files by n
        n_chunk = len(fileList)//args.num_input
        for base_path in base_path_list:
           print('using ntuples in '+base_path)
        print('Number of files is {0:2d}'.format(len(fileList)))
        print('Number of jobs is {0:2d}'.format(n_chunk+1))

        for file_index in range(n_chunk+1):
            sublist = '' 
            for idx, l in enumerate(fileList):
                if idx < args.num_input*(file_index+1) and idx >= args.num_input*file_index:
                    l = l.rstrip()
                    l = l+","
                    sublist = sublist+l
            if sublist=='': continue
            command = "python "+script_name+" --input_file="+sublist
            command += " --output "+args.outputpath+"/step"+str(args.step)+"_"+run+"_"+str(file_index)+".root --nevents "+str(args.execute)
            if args.step==4: command += " --VFATOFF "+args.outputpath+"/step2_hot_dead_vfat_run"+run+".txt --run_number "+run

            if(args.execute is not 0):
                print command
                os.system(command)
            elif(args.submit):
                arg_string = "--input_file="+sublist+" --output "+args.outputpath+"/step"+str(args.step)+"_"+run+"_"+str(file_index)+".root --nevents -1"
                if args.step==4: arg_string+= " --VFATOFF "+args.outputpath+"/step2_hot_dead_vfat_run"+run+".txt --run_number "+run
                with open('condor.submit', 'a') as submit_file:
                    submit_file.write(
                    '''
arguments = {arg}
queue
                    '''.format(arg = arg_string))
                
    # HTCondor submission
    if(args.submit):
        os.system("condor_submit condor.submit")
