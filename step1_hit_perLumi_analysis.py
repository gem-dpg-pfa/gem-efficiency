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

def iEta_iPhi2VFAT(iEta,iPhi):
    try:
        etaP = int(iEta)
        phi = int(iPhi)
    except:
        print "Provided iEta and/or iPhi are not numbers.\nExiting..."
        sys.exit(0)
    
    if iEta <1 or iEta>8 or iPhi <0 or iPhi>2:
        return -1
    
    VFAT = iPhi*8 + (8-iEta)
    return VFAT

def stripHit2VFAT(firstStrip,lastStrip,etaP):
    if firstStrip <0 or firstStrip>383:
        print "Invalid Strip provided position.\nExiting..."
        sys.exit(0)

    iPhi = -1
    if (firstStrip<127 and lastStrip<127): iPhi = 0
    elif (firstStrip>255 and lastStrip>255): iPhi = 2
    elif (firstStrip>127 and lastStrip>127 and firstStrip<255 and lastStrip<255): iPhi = 1

    VFAT = iEta_iPhi2VFAT(etaP,iPhi)

    return VFAT

#### Let's add some more from different folder
#lib_folder = os.path.expandvars('$myLIB')
#sys.path.insert(1, lib_folder)
#try:
#    from ROOT_Utils import *
#    
#
#except:
#  print("ERROR:\n\tCan't find the package CMS_lumi and tdrstlye\n\tPlease verify that this file are placed in the path $myLIB/ROOT_Utils/ \n\tAdditionally keep in mind to export the environmental variable $myLIB\nEXITING...\n") 
#  sys.exit(0)
#try:
#    from PFA_Analyzer_Utils import *
#except:
#    print ("ERROR:\n\tCan't find the package PFA_Analyzer_Utils\nEXITING...\n")
#    sys.exit(0)

 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='''Scripts that: \n\t-Reads the GEMMuonNtuple\n\t-Store #events per lumi section\n\t-Store #rechits per lumi section per each endcap/layer/chamber/vfat\nProduces plain ntuples for further analysis''',
            epilog="""Typical execution\n\t python step1_hit_perLumi_analysis.py --input_file=/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEMCommonNtuples/MWGR4/Run_342154/MuDPGNtuple_2021_MWGR4_5.root, --output output/step1_342154_5.root --nevents -1 """,
            formatter_class=RawTextHelpFormatter
    )
    
    parser.add_argument('--input_file','-in', type=str,help="Comma separated list of input files to be analyzed",required=True)
    parser.add_argument('--outputname','-out', type=str,help="Path and name of output file",default="output/step1_output.root")
    parser.add_argument('--nevents','-n', type=int,help="Maximum number of events to be analysed",default=-1)
    args = parser.parse_args()
    
    ROOT.gROOT.SetBatch(True)
   
    ## Get list of input files
    files = []
    if args.input_file != None:
        files = map(str, args.input_file.split(','))
    else:
        print "Error: provide input files!"
        exit()
    
    ## Prepare output rootfile
    if os.path.isfile(args.outputname):
        print "Warning: ", args.outputname, " will be overwritten\n"
    else:
        print "Creating output file: ", args.outputname,"\n"

    outFile = ROOT.TFile(args.outputname,"RECREATE")
    outTree1 = ROOT.TTree("hits_perLumi", "hits_perLumi")
    outTree2 = ROOT.TTree("hits_perEvent", "hits_perEvent")
    
    ##branches in output file
    nevt = array('i', [0])
    run  = array('i', [0])
    lumi = array('i', [0])
    nhit_M1 = np.zeros(shape=(36,24), dtype=np.uint32)
    nhit_M2 = np.zeros(shape=(36,24), dtype=np.uint32)
    nhit_P1 = np.zeros(shape=(36,24), dtype=np.uint32)
    nhit_P2 = np.zeros(shape=(36,24), dtype=np.uint32)
    
    outTree1.Branch("nevt", nevt, "nevt/I")
    outTree1.Branch("run", run, "run/I")
    outTree1.Branch("lumi", lumi, "lumi/I")
    outTree1.Branch("nhit_M1", nhit_M1,"nhit_M1[36][24]/I")
    outTree1.Branch("nhit_M2", nhit_M2,"nhit_M2[36][24]/I")
    outTree1.Branch("nhit_P1", nhit_P1,"nhit_P1[36][24]/I")
    outTree1.Branch("nhit_P2", nhit_P2,"nhit_P2[36][24]/I")
    
    ##branches in output file
    evt_2 = array('i', [0])
    run_2  = array('i', [0])
    lumi_2 = array('i', [0])
    nhit_M1_2 = np.zeros(shape=(36), dtype=np.uint32)
    nhit_M2_2 = np.zeros(shape=(36), dtype=np.uint32)
    nhit_P1_2 = np.zeros(shape=(36), dtype=np.uint32)
    nhit_P2_2 = np.zeros(shape=(36), dtype=np.uint32)

    outTree2.Branch("evt", evt_2, "evt/I")
    outTree2.Branch("run", run_2, "run/I")
    outTree2.Branch("lumi", lumi_2, "lumi/I")
    outTree2.Branch("nhit_M1", nhit_M1_2,"nhit_M1[36]/I")
    outTree2.Branch("nhit_M2", nhit_M2_2,"nhit_M2[36]/I")
    outTree2.Branch("nhit_P1", nhit_P1_2,"nhit_P1[36]/I")
    outTree2.Branch("nhit_P2", nhit_P2_2,"nhit_P2[36]/I")

    ## Check and TChain input files
    chain = ROOT.TChain("muNtupleProducer/MuDPGTree")
    for f in files:
        if not os.path.isfile(f):
            print "File: ",f, " does not exist\n"
            continue
        else:
            chain.Add(f)
 
    nentries = chain.GetEntries()
    print "Analysing ",nentries," events"

    #disable all branches
    chain.SetBranchStatus("*",0);     
    # List of branches to be used from input ntuples
    branchList=["event_runNumber", "event_lumiBlock", "event_eventNumber", "gemRecHit_region", "gemRecHit_layer", "gemRecHit_chamber", "gemRecHit_etaPartition", "gemRecHit_firstClusterStrip", "gemRecHit_cluster_size"]
    for b in branchList:
        chain.SetBranchStatus(b,1);
    
    nevt[0]=0
    nhit_M1.fill(0)
    nhit_M2.fill(0)
    nhit_P1.fill(0)
    nhit_P2.fill(0)
    lumi_check = 0
    
    # Loop on events
    for chain_index,evt in enumerate(chain):
        #set lumi block for the first entry
        if(chain_index==0): 
            lumi_check=int(evt.event_lumiBlock)
            #print "first lumi", int(evt.event_lumiBlock)
    
        if (args.nevents>0 and chain_index>args.nevents): break
        if (chain_index%1000==0): print "Processing ",chain_index,"/",nentries," ..."
    
        #number of rechits in the event
        nrechits = len(evt.gemRecHit_chamber)
    
        if int(evt.event_lumiBlock)==lumi_check:
            nevt[0] += 1
            run[0]  = int(evt.event_runNumber)
            lumi[0] = int(evt.event_lumiBlock)
            evt_2[0] = int(evt.event_eventNumber)
            run_2[0] = int(evt.event_runNumber)
            lumi_2[0] = int(evt.event_lumiBlock)
            ## Loop on rechits
            for j in range(nrechits):
                chamber = evt.gemRecHit_chamber[j] - 1
                VFAT    = stripHit2VFAT(evt.gemRecHit_firstClusterStrip[j],evt.gemRecHit_firstClusterStrip[j]+evt.gemRecHit_cluster_size[j], evt.gemRecHit_etaPartition[j])
                if (VFAT==-1): continue
                if(evt.gemRecHit_region[j]<0 and evt.gemRecHit_layer[j]==1): nhit_M1[chamber][VFAT]+= 1
                if(evt.gemRecHit_region[j]<0 and evt.gemRecHit_layer[j]==2): nhit_M2[chamber][VFAT]+= 1
                if(evt.gemRecHit_region[j]>0 and evt.gemRecHit_layer[j]==1): nhit_P1[chamber][VFAT]+= 1
                if(evt.gemRecHit_region[j]>0 and evt.gemRecHit_layer[j]==2): nhit_P2[chamber][VFAT]+= 1

                if(evt.gemRecHit_region[j]<0 and evt.gemRecHit_layer[j]==1): nhit_M1_2[chamber]+= 1
                if(evt.gemRecHit_region[j]<0 and evt.gemRecHit_layer[j]==2): nhit_M2_2[chamber]+= 1
                if(evt.gemRecHit_region[j]>0 and evt.gemRecHit_layer[j]==1): nhit_P1_2[chamber]+= 1
                if(evt.gemRecHit_region[j]>0 and evt.gemRecHit_layer[j]==2): nhit_P2_2[chamber]+= 1
            outTree2.Fill()
            nhit_M1_2.fill(0)
            nhit_M2_2.fill(0)
            nhit_P1_2.fill(0)
            nhit_P2_2.fill(0)
        else:
            # fill branches and reset counters
            #print "filling branches with run ",run," nevt ",nevt," lumi ",lumi
            outTree1.Fill()
            nevt[0]=0
            nhit_M1.fill(0)
            nhit_M2.fill(0)
            nhit_P1.fill(0)
            nhit_P2.fill(0)
    
        # update lumi before going to next event
        lumi_check=int(evt.event_lumiBlock)
    
    outFile.cd()
    outTree1.Write("", ROOT.TObject.kOverwrite)
    outTree2.Write("", ROOT.TObject.kOverwrite)
    outFile.Close()
