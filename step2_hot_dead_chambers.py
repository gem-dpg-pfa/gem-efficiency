import ROOT
import os.path
import numpy as np
import argparse
import sys
from argparse import RawTextHelpFormatter
lib_folder = os.path.expandvars('$PFA')
lib_folder += "/Analyzer/lib/"
sys.path.insert(1, lib_folder)
try:
    from PFA_Analyzer_Utils import *
except:
   print ("ERROR:\n\tCan't find the package PFA_Analyzer_Utils in ",lib_folder,"\nEXITING...\n")
   sys.exit(0)


base_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__))+"/../")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='''Scripts that: \n\t-Reads the mini-tree produced at step1\n\t-Identifies dead and hot VFATs (this will include chambers completely off\n\t-Identifies HV trips by looking for drops in the rechits vs lumi distributions (not yet implemented)\n\t-Produces list of vfats and lumis to be excluded from further analysis''',
            epilog="""Typical execution\n\t python step2_hot_dead_chambers.py -r 349758_Express""",
            formatter_class=RawTextHelpFormatter
    )
    
    parser.add_argument('--run','-r', type=str,help="Single run tag previously analyzed with step1. Es. 349758_Express",required=True)
    parser.add_argument('--inputfile','-inF', type=str,help="Path of input file.")
    parser.add_argument('--inputpath','-inP', type=str,help="Path of directory containing input files. [Default: %(default)s] ",default=base_dir+"/VFAT_MaskMaker/output/")
    parser.add_argument('--outputpath','-out', type=str,help="Path of output files. [Default: %(default)s] ",default=base_dir+"/PFA_Analyzer/ExcludeMe/")
    args = parser.parse_args()
    
    ROOT.gROOT.SetBatch(True)
    
    runList = []
    if args.run != None:
        runList = [args.run]
    else:
        print "Error: provide run tag"
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
    fout_name = args.outputpath+"/ListOfDeadVFAT_run"+GetRunNumber(run)+".txt"
    fout = open(fout_name, "w")
    if os.path.isfile(fout_name):
        print "Warning: ", fout_name, " will be overwritten\n"
    else:
        print "Creating output file: ", fout_name,"\n"

    fout.write("region\tlayer\tchamber\tVFAT\treason_mask\n")
    
    endcaps = ["P","M"]
    layers = ["1", "2"]
    for e in endcaps:
        for l in layers:
            for ch in range(36):
                print "Processing GE11-"+e+'-%02d' % (ch+1) +"L"+str(l)
                avg_perVFAT = []
                peak_perVFAT = []
                for vfat in range(24):
                    chain.Draw("nhit_"+e+l+"["+str(ch)+"]["+str(vfat)+"]/nevt>>hrate1d_"+e+l+"_"+str(ch)+"_"+str(vfat), "", "")
                    htemp = ROOT.gROOT.FindObject("hrate1d_"+e+l+"_"+str(ch)+"_"+str(vfat))
                    avg = float(htemp.GetMean())
                    peak = float(htemp.GetXaxis().GetXmax())

                    avg_perVFAT.append(avg)
                    peak_perVFAT.append(peak)
                #remove zeros
                avg_perVFAT_notZero = [x for x in avg_perVFAT if x != 0.0]
                avg_perChamber = sum(avg_perVFAT_notZero)/len(avg_perVFAT_notZero) if len(avg_perVFAT_notZero)>0 else 0.0
                stdev_perChamber = np.std(avg_perVFAT_notZero)
                #print ch, avg_perChamber, stdev_perChamber
                for vfat, avg in enumerate(avg_perVFAT):
                    #print vfat, ": ",avg," peak: ",peak_perVFAT[vfat]
                    #print "GE11-",e,"-",str(ch),"L",l,"\t VFAT ",vfat,"\tAVG = ",avg,"\tchamberavg = ",avg_perChamber,"\tstd = ",stdev_perChamber
                    if (avg - avg_perChamber) < -stdev_perChamber or avg == 0:
                        fout.write(e+"\t"+l+"\t"+str(ch+1)+"\t"+str(vfat)+"\t1\n")
                    #noisy chambers: average occupancy-mean>2sigmas
                    elif (avg-avg_perChamber) > 2*stdev_perChamber :        
                        fout.write(e+"\t"+l+"\t"+str(ch+1)+"\t"+str(vfat)+"\t2\n")
                        #print "removed"
                    elif peak_perVFAT[vfat]>1.0:
                        fout.write(e+"\t"+l+"\t"+str(ch+1)+"\t"+str(vfat)+"\t2\n")
                        #print "  removed"
    fout.close()
