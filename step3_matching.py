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
        print "Invalid iEta and/or iPhi provided position.\nExiting..."
        sys.exit(0)
    
    VFAT = iPhi*8 + (8-iEta)
    return VFAT

## Associates a propagated hit to a VFAT
## It is based on default GE11 geometry
## Might need refinements after alignment 
def propHit2VFAT(glb_r,loc_x,etaP):
    ## Magic numbers taken from strip geometry
    #  They are loc_x/glb_r (basically cos phi) for strip 127 and 255 respectively
    LeftMost_iPhi_boundary = -0.029824804
    RightMost_iPhi_boundary = 0.029362374

    prophit_cosine = loc_x/glb_r
    iPhi =  0 if  prophit_cosine<=LeftMost_iPhi_boundary else (2 if prophit_cosine>RightMost_iPhi_boundary else 1)

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
            description='''Script that: \n\t-Reads the GEMMuonNtuple\n\t-Loops on propagated muon tracks to the GE11 surface\n\t-Loops on the GEM rechHits and matches them to the propagated tracks\n\t-Reads the text files produced at step 2 to exclude "bad" VFATs\nProduces plain ntuples for further analysis''',
            epilog="""Typical execution\n\t python step3_matching.py --input_file=/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEMCommonNtuples/MWGR4/Run_342154/MuDPGNtuple_2021_MWGR4_5.root, --output output/step3_342154_5.root --nevents -1 """,
            formatter_class=RawTextHelpFormatter
    )
    
    parser.add_argument('--input_file','-in', type=str,help="Comma separated list of input files to be analyzed",required=True)
    parser.add_argument('--outputname','-out', type=str,help="Path and name of output file",default="output/step3_output.root")
    parser.add_argument('--nevents','-n', type=int,help="Maximum number of events to be analysed",default=-1)
    parser.add_argument('--isMC','-mc', type=str,help="Analysing MC sample Zmumu or Cosmics",choices=["Zmumu", "Cosmics"])
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
    outTree = ROOT.TTree("hit_matching", "hit_matching")
    
    ##branches in output file
    evtn = array('i', [0])
    run  = array('i', [0])
    lumi = array('i', [0])
    prop_region  = array('f', [0.0])
    prop_layer   = array('f', [0.0])
    prop_chamber = array('f', [0.0])
    prop_etaP    = array('f', [0.0])
    prop_pt      = array('f', [0.0])
    prop_eta     = array('f', [0.0])
    prop_errR    = array('f', [0.0])
    prop_errPhi  = array('f', [0.0])
    rec_cls      = array('f', [0.0])
    rec_bx       = array('f', [0.0])
    matched      = array('i', [0])
    double_matched = array('i', [0])
    #exclude = array('i', [0])

    outTree.Branch("evtn", evtn, "evtn/I")
    outTree.Branch("run", run, "run/I")
    outTree.Branch("lumi", lumi, "lumi/I")
    outTree.Branch("prop_region" , prop_region , "prop_region/F") 
    outTree.Branch("prop_layer"  , prop_layer  , "prop_layer/F")  
    outTree.Branch("prop_chamber", prop_chamber, "prop_chamber/F")
    outTree.Branch("prop_etaP"   , prop_etaP   , "prop_etaP/F")   
    outTree.Branch("prop_pt"     , prop_pt     , "prop_pt/F")     
    outTree.Branch("prop_eta"    , prop_eta    , "prop_eta/F")    
    outTree.Branch("prop_errR"   , prop_errR   , "prop_errR/F")   
    outTree.Branch("prop_errPhi" , prop_errPhi , "prop_errPhi/F") 
    outTree.Branch("rec_cls" , rec_cls , "rec_cls/F") 
    outTree.Branch("rec_bx"  , rec_bx  , "rec_bx/F") 
    outTree.Branch("matched"     , matched     , "matched/I")     
    outTree.Branch("double_matched", double_matched, "double_matched/I")

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
    branchList=["event_runNumber", "event_lumiBlock", "event_eventNumber",\
                "gemRecHit_region", "gemRecHit_layer", "gemRecHit_chamber", "gemRecHit_etaPartition",\
                "gemRecHit_loc_x", "gemRecHit_loc_y", "gemRecHit_loc_z", "gemRecHit_loc_r", "gemRecHit_loc_phi",\
                "gemRecHit_g_x", "gemRecHit_g_y", "gemRecHit_g_z", "gemRecHit_g_r", "gemRecHit_g_phi",\
                "gemRecHit_nRecHits", "gemRecHit_cluster_size", "gemRecHit_bx",\
                "mu_propagated_region", "mu_propagated_layer", "mu_propagated_chamber", "mu_propagated_etaP", "mu_propagated_Outermost_z",\
                "mu_propagated_pt", "mu_propagated_eta", "mu_propagated_phi", "mu_propagated_isME11",\
                "mu_propagatedGlb_r", "mu_propagatedGlb_phi", "mu_propagatedGlb_errR", "mu_propagatedGlb_errPhi",\
                "mu_propagatedGlb_x", "mu_propagatedGlb_y", "mu_propagatedGlb_z", "mu_propagatedGlb_errX", "mu_propagatedGlb_errY",\
                "mu_propagated_EtaPartition_phiMin", "mu_propagated_EtaPartition_phiMax",\
                "mu_propagated_EtaPartition_rMin", "mu_propagated_EtaPartition_rMax",\
                "mu_propagated_EtaPartition_centerX", "mu_propagated_EtaPartition_centerY",\
                "mu_pt", "mu_eta", "mu_phi", "mu_charge", "mu_isCSC", "mu_isME11",\
                ]
    for b in branchList:
        chain.SetBranchStatus(b,1);
    #enable branches for MC studies
    if(args.isMC is not None):
        branchList_MC=["genParticle_Pt", "genParticle_Eta", "genParticle_Phi",\
                       "genParticle_PdgId", "genParticle_MotherPdgId", "genParticle_GrandMotherPdgId"]
        for b in branchList_MC:
            chain.SetBranchStatus(b,1);

    # Loop on events
    for chain_index,evt in enumerate(chain):
        if (args.nevents>0 and chain_index>args.nevents): break
        if (chain_index%1000==0): print "Processing ",chain_index,"/",nentries," ..."
 
        run[0]  = int(evt.event_runNumber)
        lumi[0] = int(evt.event_lumiBlock)
        evtn[0] = int(evt.event_eventNumber)

        #number of propagated mu in the event
        nprop = len(evt.mu_propagated_pt)
            
        ## Loop on propagated mu
        for i in range(nprop):
            #bug fix propagation from opposite endcap!
            if evt.mu_propagated_Outermost_z[i]*evt.mu_propagated_region[i]<0: continue
            #requirements on fiducial region: r
            if  abs(evt.mu_propagatedGlb_r[i] - evt.mu_propagated_EtaPartition_rMin[i]) < 1: continue
            #requirements on fiducial region: phi
            phi = evt.mu_propagatedGlb_phi[i] if evt.mu_propagatedGlb_phi[i]>0 else 2*np.pi+evt.mu_propagatedGlb_phi[i]
            phiMin = evt.mu_propagated_EtaPartition_phiMin[i]
            phiMax = evt.mu_propagated_EtaPartition_phiMax[i]
            if phiMin>phiMax: phiMax = 2*np.pi + phiMax
            if abs(phi-phiMin)<0.005 or abs(phi-phiMax)<0.005: continue

            prop_region[0]  = evt.mu_propagated_region[i]
            prop_layer[0]   = evt.mu_propagated_layer[i]
            prop_chamber[0] = evt.mu_propagated_chamber[i]
            prop_etaP[0]    = evt.mu_propagated_etaP[i]
            prop_pt[0]      = evt.mu_propagated_pt[i]
            prop_eta[0]     = evt.mu_propagated_eta[i]
            prop_errR[0]    = evt.mu_propagatedGlb_errR[i]
            prop_errPhi[0]  = evt.mu_propagatedGlb_errPhi[i]
            rec_cls[0] = -99.
            rec_bx[0]  = -99.
            matched[0] = 0
            double_matched[0] = 0
            #exclude[0] = 0

            #number of propagated mu in the event
            nrechits = len(evt.gemRecHit_region)
            min_deltaphi = 999.
            min_recIndex = 0
            ## Loop on rechits
            for j in range(nrechits):
                #chamber = evt.gemRecHit_chamber[j] - 1
                #VFAT    = propHit2VFAT(evt.gemRecHit_g_r[j],evt.gemRecHit_loc_x[j], evt.gemRecHit_etaPartition[j])
                sameRegion   = bool(evt.gemRecHit_region[j] == evt.mu_propagated_region[i])
                sameLayer    = bool(evt.gemRecHit_layer[j] == evt.mu_propagated_layer[i])
                sameChamber  = bool(evt.gemRecHit_chamber[j] == evt.mu_propagated_chamber[i])
                sameEtaP     = bool(abs(evt.gemRecHit_etaPartition[j] - evt.mu_propagated_etaP[i]) < 2)
                if not sameRegion:  continue
                if not sameLayer:   continue
                if not sameChamber: continue
                if not sameEtaP:    continue

                ## Looking for min. distance
                deltaphi = abs(evt.mu_propagatedGlb_phi[i]-evt.gemRecHit_g_phi[j])
                if deltaphi<min_deltaphi:
                    min_deltaphi=deltaphi
                    min_recIndex=j
            ## Store information of the closest rechit
            rec_cls[0] = evt.gemRecHit_cluster_size[min_recIndex]
            rec_bx[0]  = evt.gemRecHit_bx[min_recIndex]

            ## Matching condition: r*dphi 1 cm
            if evt.mu_propagatedGlb_r[i]*abs(min_deltaphi)<1:
                matched[0] = 1
                ## if matched==1, look for match in other layer
                min_deltaphi_2 = 999.
                min_recIndex_2 = 0
                ## Loop on rechits
                for j in range(nrechits):
                    sameRegion   = bool(evt.gemRecHit_region[j] == evt.mu_propagated_region[i])
                    sameLayer    = bool(evt.gemRecHit_layer[j] == evt.mu_propagated_layer[i])
                    sameChamber  = bool(evt.gemRecHit_chamber[j] == evt.mu_propagated_chamber[i])
                    sameEtaP     = bool(abs(evt.gemRecHit_etaPartition[j] - evt.mu_propagated_etaP[i]) < 2)
                    if not sameRegion:  continue
                    if sameLayer:   continue #in this case we require different layer
                    if not sameChamber: continue
                    if not sameEtaP:    continue

                    ## Looking for min. distance
                    deltaphi = abs(evt.mu_propagatedGlb_phi[i]-evt.gemRecHit_g_phi[j])
                    if deltaphi<min_deltaphi_2:
                        min_deltaphi_2=deltaphi
                        min_recIndex_2=j
                ## Matching condition: r*dphi 1 cm
                if evt.mu_propagatedGlb_r[i]*abs(min_deltaphi_2)<1:
                    double_matched[0] = 1
            # For each propagated, fill branches
            outTree.Fill()
    
    outFile.cd()
    outTree.Write("", ROOT.TObject.kOverwrite)
    outFile.Close()
