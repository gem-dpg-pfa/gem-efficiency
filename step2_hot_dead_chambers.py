import ROOT
import csv
import os.path
from os import mkdir
import subprocess
import numpy as np
import math
import sys
import time
import argparse
import pandas as pd
from array import array
from argparse import RawTextHelpFormatter

working_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='''Scripts that: \n\t-Reads the mini-tree produced at step1\n\t-Identifies dead and hot VFATs (this will include chambers completely off\n\t-Identifies HV trips by looking for drops in the rechits vs lumi distributions (not yet implemented)\n\t-Produces list of vfats and lumis to be excluded from further analysis''',
            epilog="""Typical execution\n\t python step2_hot_dead_chambers.py -r 342154""",
            formatter_class=RawTextHelpFormatter
    )
    
    parser.add_argument('--run','-r', type=str,help="Comma separated list of Cosmic runs to be analyzed",required=True)
    parser.add_argument('--inputfile','-inF', type=str,help="Path of input file.")
    parser.add_argument('--inputpath','-inP', type=str,help="Path of directory containing input files. [Default: %(default)s] ",default=working_dir+"/output/")
    parser.add_argument('--outputpath','-out', type=str,help="Path of output files. [Default: %(default)s] ",default=working_dir+"/output/")
    parser.add_argument('--hot_rate','-hr', type=float,help="Max value of hits per lumi per VFAT. [Default: %(default)s] ",default=0.7)
    args = parser.parse_args()
    
    ROOT.gROOT.SetBatch(True)
    
    ## Get list of Cosmic runs
    runList = []
    if args.run != None:
        runList = map(str, args.run.split(','))
    else:
        print "Error: provide run number"
        exit()
    print runList
        
    base_path = args.inputpath
    for run in runList:
        #generating the list of all .root files in given directory and subdirectories
        fileList = []
        if args.inputfile != None:
            fileList.append(args.inputfile)
        else:
        #generating the list of all .root files in given directory and subdirectories
            for r, d, f in os.walk(base_path): # r=root, d=directories, f = files
                for file in f:
                    if '.root' in file and 'step1_'+run in file:
                        fileList.append(os.path.join(r, file))
    print "Input file(s):\n", fileList

    ## Check and TChain input files
    chain = ROOT.TChain("hits_perLumi")
    for f in fileList:
        if not os.path.isfile(f):
            print "File: ",f, " does not exist\n"
            continue
        else:
            chain.Add(f)

    nentries = chain.GetEntries()
    print "Analysing ",nentries," lumis"

    # prepare output file (list of hot/dead VFATs)
    fout_name = args.outputpath+"step2_hot_dead_vfat_run"+args.run.replace(",", "_")+".txt"
    fout = open(fout_name, "w")
    if os.path.isfile(fout_name):
        print "Warning: ", fout_name, " will be overwritten\n"
    else:
        print "Creating output file: ", fout_name,"\n"

    fout.write("region\tlayer\tchamber\tVFAT\treason_mask\n")
    
    endcaps = ["M", "P"]
    layers = ["1", "2"]
    for e in endcaps:
        for l in layers:
            for ch in range(36):
                for vfat in range(24):
                    chain.Draw("nhit_M1["+str(ch)+"]["+str(vfat)+"]>>htemp", "", "")
                    htemp = ROOT.gROOT.FindObject('htemp')
                    avg = htemp.GetMean()
                    if avg == 0: #dead vfat
                        fout.write(e+"\t"+l+"\t"+str(ch+1)+"\t"+str(vfat)+"\t1\n")
                    elif avg > args.hot_rate: #hot vfat
                        fout.write(e+"\t"+l+"\t"+str(ch+1)+"\t"+str(vfat)+"\t2\n")
    fout.close()
