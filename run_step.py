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
lib_folder = os.path.expandvars('$PFA')
lib_folder += "/lib/"
sys.path.insert(1, lib_folder)
try:
    from ROOT_Utils import *
except:
   print ("ERROR:\n\tCan't find the package ROOT_Utils in ",lib_folder,"\nEXITING...\n")

working_dir = os.path.dirname(os.path.abspath(__file__))
outputpath=working_dir+"/output/"


parser = argparse.ArgumentParser(
    description='''Script that: \n\t-Submits step1_hit_perLumi_analysis.py to HTCondor\n''',
    epilog="""Typical execution\n\t python run_step.py -r 349263_Express""",
    formatter_class=RawTextHelpFormatter
    )

parser.add_argument('--run','-r', type=str,help="Space separated list of runs to be analyzed",required=True, nargs='*')
parser.add_argument('--submit', '-s', help="Submit to HTCondor. [Default: %(default)s]", default=True) 
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

    bash_name = working_dir+"/CondorFiles/run_step1_wrap.sh" 
    with open(bash_name,'w') as job_file:
        job_file.write(
'''#!/bin/bash
python {path}/step1_hit_perLumi_analysis.py "$@" '''.format(path=working_dir))
    script_name = "step1_hit_perLumi_analysis.py" 

    for run in runList:
        fileList = files_in_folder(run)
        # for HTCondor submission, prepare file
        if(args.submit):
            condorsubmit_file = "./CondorFiles/"+"condor_"+run+".submit"
            with open(condorsubmit_file, 'w') as submit_file:
                submit_file.write(
                '''
universe = vanilla
executable     = {bash_exe}
getenv = True

request_memory = 4096
request_cpus   = 1
request_disk   = 16383

error   = ./CondorFiles/err_{submit_time}_.$(Process)
output  = ./CondorFiles/out_{submit_time}_.$(Process)
log     = ./CondorFiles/log_{submit_time}_$(Process).log
+JobFlavour = "workday"
                '''.format(bash_exe = bash_name, path=outputpath, submit_time=str(time.time())))


        n_chunk = len(fileList)
        print('Number of jobs is {0:2d}'.format(n_chunk))

        for idx,file_path in enumerate(fileList):

            arg_string = "--input_file="+file_path+" --output "+outputpath+"/step1_"+run+"_"+str(idx)+".root --nevents -1"
            with open(condorsubmit_file, 'a') as submit_file:
                submit_file.write(
                '''
arguments = {arg}
queue
                '''.format(arg = arg_string))

    os.system("condor_submit "+condorsubmit_file)
