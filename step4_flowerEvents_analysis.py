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
from PFA_Analyzer_Utils import *

 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='''Scripts that: \n\t-Reads the GEMMuonNtuple\n\t-Store #events per lumi section\n\t-Store #rechits per lumi section per each endcap/layer/chamber/vfat\nProduces plain ntuples for further analysis''',
            epilog="""Typical execution\n\t python step4_flowerEvent_analysis.py --input_file=/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEMCommonNtuples/MWGR4/Run_342154/MuDPGNtuple_2021_MWGR4_5.root, --output output/step4_342154_5.root --nevents -1 """,
            formatter_class=RawTextHelpFormatter
    )
    
    parser.add_argument('--input_file','-in', type=str,help="Comma separated list of input files to be analyzed",required=True)
    parser.add_argument('--outputname','-out', type=str,help="Path and name of output file",default="output/step4_output.root")
    parser.add_argument('--run_number','-run', type=str,help="Run number",default="999")
    parser.add_argument('--nevents','-n', type=int,help="Maximum number of events to be analysed",default=-1)
    parser.add_argument('--VFATOFF', type=str , help="file_path to the file containing a list of the VFAT you want to exclude from efficiency evaluation. The file must be tab separated with \tregion,layer,chamber,VFAT,reason_mask",required=False,nargs='*')
    args = parser.parse_args()
    
    ROOT.gROOT.SetBatch(True)
    
    ##Prepare list of OFF and Noisy VFATs
    VFATOFFDict = {} if args.VFATOFF is None else importOFFVFAT(args.VFATOFF, 1)   
    VFATNoisyDict = {} if args.VFATOFF is None else importOFFVFAT(args.VFATOFF, 2)   
    #print 'list of OFF vfats: '
    #print VFATOFFDict
    #print 'list of Noisy VFATs: '
    #print VFATNoisyDict
    
    ## Get list of input files
    files = []
    if args.input_file != None:
        files = map(str, args.input_file.split(','))
    else:
        print "Error: provide input files!"
        exit()
    
    ## Check and TChain input files
    chain = ROOT.TChain("muNtupleProducer/MuDPGTree")
    for f in files:
        if not os.path.isfile(f):
            print "File: ",f, " does not exist\n"
            continue
        else:
            chain.Add(f)

    ## Prepare output rootfile
    if os.path.isfile(args.outputname):
        print "Warning: ", args.outputname, " will be overwritten\n"
    else:
        print "Creating output file: ", args.outputname,"\n"
    
    outFile = ROOT.TFile(args.outputname,"RECREATE")
    outFile.cd()

    ## Book output histograms
    hrate1d_P1 = ROOT.TH1I("hrate1d_P1", "hrate1d_P1", 5000, 0, 5000)
    hrate1d_P2 = ROOT.TH1I("hrate1d_P2", "hrate1d_P2", 5000, 0, 5000)
    hrate1d_M1 = ROOT.TH1I("hrate1d_M1", "hrate1d_M1", 5000, 0, 5000)
    hrate1d_M2 = ROOT.TH1I("hrate1d_M2", "hrate1d_M2", 5000, 0, 5000)    
    hCLS_P1 = ROOT.TH1I("hCLS_P1", "hCLS_P1", 128, 0, 128)
    hCLS_P2 = ROOT.TH1I("hCLS_P2", "hCLS_P2", 128, 0, 128)
    hCLS_M1 = ROOT.TH1I("hCLS_M1", "hCLS_M1", 128, 0, 128)
    hCLS_M2 = ROOT.TH1I("hCLS_M2", "hCLS_M2", 128, 0, 128)    
  
    nentries = chain.GetEntries()
    print "Analysing ",nentries," events"
    
    #disable all branches
    chain.SetBranchStatus("*",0);     
    # List of branches to be used from input ntuples
    branchList=["event_runNumber", "event_lumiBlock", "event_eventNumber", "gemRecHit_region", "gemRecHit_layer", "gemRecHit_chamber", "gemRecHit_etaPartition", "gemRecHit_firstClusterStrip", "gemRecHit_cluster_size", "gemRecHit_g_r", "gemRecHit_loc_x"]
    for b in branchList:
        chain.SetBranchStatus(b,1);
    
    count_500hits_P1 = 0
    count_500hits_P2 = 0
    count_500hits_M1 = 0
    count_500hits_M2 = 0
    count_500hits_AND = 0
    count_1000hits_P1 = 0
    count_1000hits_P2 = 0
    count_1000hits_M1 = 0
    count_1000hits_M2 = 0
    count_1000hits_AND = 0

    hits_vector_P1 = []
    hits_vector_P2 = []
    hits_vector_M1 = []
    hits_vector_M2 = []
    avgCls_vector_P1 = []
    avgCls_vector_P2 = []
    avgCls_vector_M1 = []
    avgCls_vector_M2 = []
    evt_vector  = []
    lumi_vector = []

    # Loop on events
    for chain_index,evt in enumerate(chain):
        if (args.nevents>0 and chain_index>args.nevents): break
        if (chain_index%1000==0): print "Processing ",chain_index,"/",nentries," ..."

        evtn = int(evt.event_eventNumber)
        run  = int(evt.event_runNumber)
        lumi = int(evt.event_lumiBlock)

        hits_P1 = 0
        hits_P2 = 0
        hits_M1 = 0
        hits_M2 = 0
        cls_vector_P1 = []
        cls_vector_P2 = []
        cls_vector_M1 = []
        cls_vector_M2 = []
    
        #number of rechits in the event
        nrechits = len(evt.gemRecHit_chamber)
    
        ## Loop on rechits
        for j in range(nrechits):
            region = evt.gemRecHit_region[j]
            chamber = evt.gemRecHit_chamber[j]
            layer = evt.gemRecHit_layer[j]
            etaP = evt.gemRecHit_etaPartition[j]
            rec_glb_r = evt.gemRecHit_g_r[j]
            rec_loc_x = evt.gemRecHit_loc_x[j]
            cls = evt.gemRecHit_cluster_size[j]
            RecHitEtaPartitionID =  region*(100*chamber+10*layer+etaP)
            VFAT = propHit2VFAT(rec_glb_r, rec_loc_x, etaP)
            if RecHitEtaPartitionID in VFATOFFDict:
                #Skip OFF VFATs
                if VFAT in VFATOFFDict[RecHitEtaPartitionID]:
                   continue
                #Skip Noisy VFATs
                if VFAT in VFATNoisyDict[RecHitEtaPartitionID]:
                   continue
            if region>0 and layer==1: 
                hits_P1 += 1
                hCLS_P1.Fill(cls)
                cls_vector_P1.append(cls)
            if region>0 and layer==2: 
                hits_P2 += 1 
                hCLS_P2.Fill(cls)
                cls_vector_P2.append(cls)
            if region<0 and layer==1:
                hits_M1 += 1
                hCLS_M1.Fill(cls) 
                cls_vector_M1.append(cls)
            if region<0 and layer==2: 
                hits_M2 += 1
                hCLS_M2.Fill(cls) 
                cls_vector_M2.append(cls)

        #avg cls
        avgCls_vector_P1.append(sum(cls_vector_P1)/float(len(cls_vector_P1))) if hits_P1>0 else avgCls_vector_P1.append(0.0)
        avgCls_vector_P2.append(sum(cls_vector_P2)/float(len(cls_vector_P2))) if hits_P2>0 else avgCls_vector_P2.append(0.0)
        avgCls_vector_M1.append(sum(cls_vector_M1)/float(len(cls_vector_M1))) if hits_M1>0 else avgCls_vector_M1.append(0.0)
        avgCls_vector_M2.append(sum(cls_vector_M2)/float(len(cls_vector_M2))) if hits_M2>0 else avgCls_vector_M2.append(0.0)
        hrate1d_P1.Fill(hits_P1)
        hrate1d_P2.Fill(hits_P2)
        hrate1d_M1.Fill(hits_M1)
        hrate1d_M2.Fill(hits_M2)   

        hits_vector_P1.append(hits_P1) 
        hits_vector_P2.append(hits_P2) 
        hits_vector_M1.append(hits_M1) 
        hits_vector_M2.append(hits_M2) 
        evt_vector.append(evtn)
        if(hits_P1>500): count_500hits_P1 += 1
        if(hits_P2>500): count_500hits_P2 += 1
        if(hits_M1>500): count_500hits_M1 += 1
        if(hits_M2>500): count_500hits_M2 += 1
        if(hits_P1>500 and hits_P2>500 and hits_M1>500 and hits_M2>500): count_500hits_AND += 1
        if(hits_P1>1000): count_1000hits_P1 += 1
        if(hits_P2>1000): count_1000hits_P2 += 1
        if(hits_M1>1000): count_1000hits_M1 += 1
        if(hits_M2>1000): count_1000hits_M2 += 1
        if(hits_P1>1000 and hits_P2>1000 and hits_M1>1000 and hits_M2>1000): count_1000hits_AND += 1

    dim=len(evt_vector)
    evt_vector_P1=evt_vector
    evt_vector_P2=evt_vector
    evt_vector_M1=evt_vector
    evt_vector_M2=evt_vector
    evt_vector_P1, hits_vector_P1 = (list(t) for t in zip(*sorted(zip(evt_vector_P1, hits_vector_P1))))
    evt_vector_P2, hits_vector_P2 = (list(t) for t in zip(*sorted(zip(evt_vector_P2, hits_vector_P2))))
    evt_vector_M1, hits_vector_M1 = (list(t) for t in zip(*sorted(zip(evt_vector_M1, hits_vector_M1))))
    evt_vector_M2, hits_vector_M2 = (list(t) for t in zip(*sorted(zip(evt_vector_M2, hits_vector_M2))))

    evt_vector_P1_2=evt_vector
    evt_vector_P2_2=evt_vector
    evt_vector_M1_2=evt_vector
    evt_vector_M2_2=evt_vector
    evt_vector_P1_2, avgCls_vector_P1 = (list(t) for t in zip(*sorted(zip(evt_vector_P1_2, avgCls_vector_P1))))
    evt_vector_P2_2, avgCls_vector_P2 = (list(t) for t in zip(*sorted(zip(evt_vector_P2_2, avgCls_vector_P2))))
    evt_vector_M1_2, avgCls_vector_M1 = (list(t) for t in zip(*sorted(zip(evt_vector_M1_2, avgCls_vector_M1))))
    evt_vector_M2_2, avgCls_vector_M2 = (list(t) for t in zip(*sorted(zip(evt_vector_M2_2, avgCls_vector_M2))))

    g_P1 = ROOT.TGraph(dim, np.array(evt_vector_P1, 'd'), np.array(hits_vector_P1, 'd'))
    g_P2 = ROOT.TGraph(dim, np.array(evt_vector_P2, 'd'), np.array(hits_vector_P2, 'd'))
    g_M1 = ROOT.TGraph(dim, np.array(evt_vector_M1, 'd'), np.array(hits_vector_M1, 'd'))
    g_M2 = ROOT.TGraph(dim, np.array(evt_vector_M2, 'd'), np.array(hits_vector_M2, 'd'))
 
    g_avgCls_P1 = ROOT.TGraph(dim, np.array(evt_vector_P1, 'd'), np.array(avgCls_vector_P1, 'd'))
    g_avgCls_P2 = ROOT.TGraph(dim, np.array(evt_vector_P2, 'd'), np.array(avgCls_vector_P2, 'd'))
    g_avgCls_M1 = ROOT.TGraph(dim, np.array(evt_vector_M1, 'd'), np.array(avgCls_vector_M1, 'd'))
    g_avgCls_M2 = ROOT.TGraph(dim, np.array(evt_vector_M2, 'd'), np.array(avgCls_vector_M2, 'd'))
 
    outFile.cd()
    g_P1.Write()
    g_P2.Write()
    g_M1.Write()
    g_M2.Write()
    g_avgCls_P1.Write()
    g_avgCls_P2.Write()
    g_avgCls_M1.Write()
    g_avgCls_M2.Write()
    hrate1d_P1.Write()
    hrate1d_P2.Write()
    hrate1d_M1.Write()
    hrate1d_M2.Write()  

    #prepare canvases
    c2 = ROOT.TCanvas("c1","c1",150,10,800,800)
    ROOT.gStyle.SetOptStat(0)
    c2.SetLogy();
    hrate1d_P1.SetLineColor(ROOT.kBlack);
    hrate1d_P2.SetLineColor(ROOT.kRed);
    hrate1d_M1.SetLineColor(ROOT.kBlue);
    hrate1d_M2.SetLineColor(ROOT.kMagenta);

    hrate1d_P1.SetTitle("run "+args.run_number);
    hrate1d_P1.GetYaxis().SetRangeUser(1, 130000);
    hrate1d_P1.GetXaxis().SetTitle("n. recHits per event");

    hrate1d_P1.Draw("hist");
    hrate1d_P2.Draw("hist same");
    hrate1d_M1.Draw("hist same");
    hrate1d_M2.Draw("hist same");
    c2.Update();

    leg = ROOT.TLegend(0.6,0.7,0.9,0.9);
    leg.AddEntry(hrate1d_P1,"GE+1/1 - Layer 1","f");
    leg.AddEntry(hrate1d_P2,"GE+1/1 - Layer 2","f");
    leg.AddEntry(hrate1d_M1,"GE-1/1 - Layer 1","f");
    leg.AddEntry(hrate1d_M2,"GE-1/1 - Layer 2","f");
    leg.Draw();
    #c2.SaveAs("run"+args.run_number+"_hit_event_v2.png");
    c2.Write()

    c3 = ROOT.TCanvas("c3","c3",150,10,800,800);
    ROOT.gStyle.SetOptStat(0);
    #c3->SetLogy();
    g_P1.SetLineColor(ROOT.kBlack);
    g_P2.SetLineColor(ROOT.kRed);
    g_M1.SetLineColor(ROOT.kBlue);
    g_M2.SetLineColor(ROOT.kMagenta);

    g_P1.SetMarkerColor(ROOT.kBlack);
    g_P2.SetMarkerColor(ROOT.kRed);
    g_M1.SetMarkerColor(ROOT.kBlue);
    g_M2.SetMarkerColor(ROOT.kMagenta);

    g_P1.SetMarkerStyle(20)
    g_P2.SetMarkerStyle(21)
    g_M1.SetMarkerStyle(22)
    g_M2.SetMarkerStyle(23)

    g_P1.SetTitle("run "+args.run_number);
    g_P1.GetYaxis().SetRangeUser(1, 5000);
    g_P1.GetXaxis().SetTitle("event number");
    g_P1.GetYaxis().SetTitle("n. recHits");

    g_P1.Draw("AP");
    g_P2.Draw("P same");
    g_M1.Draw("P same");
    g_M2.Draw("P same");
    c3.Update();

    leg3 = ROOT.TLegend(0.6,0.7,0.9,0.9);
    leg3.AddEntry(g_P1,"GE+1/1 - Layer 1","f");
    leg3.AddEntry(g_P2,"GE+1/1 - Layer 2","f");
    leg3.AddEntry(g_M1,"GE-1/1 - Layer 1","f");
    leg3.AddEntry(g_M2,"GE-1/1 - Layer 2","f");
    leg3.Draw();
    #c3.SaveAs("run"+args.run_number+"_trend_event.png");
    c3.Write()

    with open('output/step4_summary_run'+args.run_number+'.txt', 'w') as f:
        f.write('count_500hits_P1\t'+str(count_500hits_P1)+'\n')
        f.write('count_500hits_P2\t'+str(count_500hits_P2)+'\n')
        f.write('count_500hits_M1\t'+str(count_500hits_M1)+'\n')
        f.write('count_500hits_M2\t'+str(count_500hits_M2)+'\n')
        f.write('count_500hits_AND\t'+str(count_500hits_AND)+'\n')
        f.write('count_1000hits_P1\t'+str(count_1000hits_P1)+'\n')
        f.write('count_1000hits_P2\t'+str(count_1000hits_P2)+'\n')
        f.write('count_1000hits_M1\t'+str(count_1000hits_M1)+'\n')
        f.write('count_1000hits_M2\t'+str(count_1000hits_M2)+'\n')
        f.write('count_1000hits_AND\t'+str(count_1000hits_AND)+'\n')

    outFile.Close()
