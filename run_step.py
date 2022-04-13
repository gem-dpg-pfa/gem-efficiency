#!/usr/bin/env python
import os
import sys
import argparse
import ROOT
from ROOT import TFile
import math
import array
import time
from argparse import RawTextHelpFormatter
base_dir = os.path.expandvars('$PFA')
lib_folder = base_dir+ "/Analyzer/lib/"
sys.path.insert(1, lib_folder)
try:
    from ROOT_Utils import *
except:
   print ("ERROR:\n\tCan't find the package ROOT_Utils in ",lib_folder,"\nEXITING...\n")

working_dir = base_dir+"/VFAT_MaskMaker/"
outputpath=working_dir+"/output/"


parser = argparse.ArgumentParser(
    description='''Script that: \n\t-Submits step1_hit_perLumi_analysis.py and step2_hit_perLumi_analysis.py to HTCondor to get VFAT masks for the run\n''',
    epilog="""Typical execution\n\t python run_step.py -r 349263_Express""",
    formatter_class=RawTextHelpFormatter
    )

parser.add_argument('--run','-r', type=str,help="Space separated list of runs to be analyzed",required=True, nargs='*')
args = parser.parse_args()
if __name__ == "__main__":

    ## Get list of Cosmic runs
    runList = []
    if args.run != None:
        runList = args.run
    else:
        print "Error: provide run number"
        exit()
    print runList
    
    ## CREATING THE EXECUTABLES, BASED ON PYTHON CODES
    bash1_name = working_dir+"/CondorFiles/run_step1_wrap.sh" 
    bash2_name = working_dir+"/CondorFiles/run_step2_wrap.sh" 
    with open(bash1_name,'w') as executable_file:
        executable_file.write(
'''#!/bin/bash
python {path}/step1_hit_perLumi_analysis.py "$@" '''.format(path=working_dir))

    with open(bash2_name,'w') as executable_file:
        executable_file.write(
'''#!/bin/bash
python {path}/step2_hot_dead_chambers.py "$@" '''.format(path=working_dir))



    for run in runList:
        fileList = files_in_folder(run)
        # prepare step1 file for HTCondor submission
        condorsubmit1_file = working_dir+"/CondorFiles/condor_step1_"+run+".submit"
        with open(condorsubmit1_file, 'w') as submit1_file:
            submit1_file.write(
            '''
universe = vanilla
executable     = {bash_exe}
getenv = True

request_memory = 4096
request_cpus   = 1
request_disk   = 16383

error   = {workdir}/CondorFiles/JobFiles/err_{submit_time}_$(Process)
output  = {workdir}/CondorFiles/JobFiles/out_{submit_time}_$(Process)
log     = {workdir}/CondorFiles/JobFiles/log_{submit_time}_$(Process)
+JobFlavour = "workday"
            '''.format(workdir=working_dir,bash_exe = bash1_name, path=outputpath, submit_time=str(time.time())))

        print('Number of jobs is {0:2d}'.format(len(fileList)))


        ## ADDING THE INPUT FILES TO THE SUBMISSION FILE
        for idx,file_path in enumerate(fileList):
            arg_string = "--input_file="+file_path+" --output "+outputpath+"/step1_"+run+"_"+str(idx)+".root --nevents -1"
            with open(condorsubmit1_file, 'a') as submit1_file:
                submit1_file.write(
                '''
arguments = {arg}
queue
                '''.format(arg = arg_string))


        # prepare step2 file for HTCondor submission        
        condorsubmit2_file = working_dir+"/CondorFiles/"+"condor_step2_"+run+".submit"
        step2_argument="-r "+ run
        with open(condorsubmit2_file, "w") as submit2_file:
            submit2_file.write(
            '''
universe = vanilla
executable     = {bash_exe}
getenv = True

request_memory = 4096
request_cpus   = 1
request_disk   = 16383

error   = {workdir}/CondorFiles/JobFiles/step2_err_$(Process)
output  = {workdir}/CondorFiles/JobFiles/step2_out_$(Process)
log     = {workdir}/CondorFiles/JobFiles/step2_log_$(Process)
+JobFlavour = "workday"
arguments = {arg}
queue
                '''.format(workdir=working_dir,bash_exe = bash2_name, path=outputpath, arg=step2_argument))