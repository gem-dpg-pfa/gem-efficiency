# gem-efficiency
## Install the code
```sh
ssh -XtY my.name@lxplus.cern.ch
cd some/working/dir
git clone git@github.com:gem-dpg-pfa/gem-efficiency.git
cd gem_efficiency
source setup.sh
```

## How does the code work
1.  Script `step1_hit_perLumi_analysis.py`: 
 - Reads the GEMMuonNtuple
 - Stores #events per lumi section
 - Stores #rechits per lumi section per each endcap/layer/chamber/vfat
 - Produces plain ntuples for further analysis

  Typical execution:
```sh
python step1_hit_perLumi_analysis.py --input_file=/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEMCommonNtuples/MWGR4/Run_342154/MuDPGNtuple_2021_MWGR4_5.root, --output output/step1_342154_5.root --nevents -1
```

2.  Script `step2_hot_dead_chambers.py`: 
 - Reads the mini-tree produced at step1 (will chain all files for a given run (multiple runs also supported)
 - Identifies dead and hot VFATs (this will include chambers completely off
 - Produces list of vfats and lumis to be excluded from further analysis

  Typical execution:
```sh
python step2_hot_dead_chambers.py -r 342154
```

3.  Script `step3_matching.py`: 
 - Reads the GEMMuonNtuple
 - Loops on propagated muon tracks to the GE11 surface
 - Loops on the GEM rechHits and matches them to the propagated tracks
 - Reads the text files produced at step 2 to exclude "bad" VFATs (not yet implemented)
 - Produces plain ntuples for further analysis and plotting

  Typical execution:
```sh
python step3_matching.py --input_file=/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEMCommonNtuples/MWGR4/Run_342154/MuDPGNtuple_2021_MWGR4_5.root, --output output/step3_342154_5.root --nevents -1
```

## How to run on multiple files (and multiple data taking runs)
The script `run_step.py` is a wrap script which executes step1 and step3 or submit jobs to HTCondor, while step2 is intented to be run on terminal.
### Examples
 - `python run_step.py -r 342154 --step 1 --execute 10000 --o output/`
   - Will run `step1_hit_perLumi_analysis.py `
   - Will analyse the first 10000 events from ntuples of run 342154 running on terminal. Input files will be analysed one by one
   - Plain ntuples for noise monitoring will be produced in the output path given, with name `step1_342154_0.root
   - If number of events exceeds content of first ntuple file found on `\eos\`, multiple output rootfiles will be produced

 - `python run_step.py -r 342154 --step 1 --execute -1 -n 2`
   - Will run `step1_hit_perLumi_analysis.py `
   - Will analyse all events from ntuples of run 342154 running on terminal
   - Input files will be analysed in groups of 2
   - Plain ntuples for noise monitoring will be produced in the default output path, with names `step1_342154_*.root

 - `python run_step.py -r 342154 --step 1 --submit True`
   - Will produce a `condor.submit` file, containing commands for executing `step1_hit_perLumi_analysis.py `
   - Each command will be intended to analyse one input file from run 342154 ntuples
   - Plain ntuples for noise monitoring will be produced in the default output path, with names `step1_342154_*.root
   - err, out, log HTCondor files will be produced in the default output path

 - `python run_step.py -r 342810,342966,343034,343082 --step 1 --submit True`
   - Will produce a `condor.submit` file, containing commands for executing `step1_hit_perLumi_analysis.py `
   - Each command will be intended to analyse one input file
   - All files for the runs 342810, 342966, 343034 and 343082 will be analysed
   - Plain ntuples for noise monitoring will be produced in the default output path, with names `step1_$run_*.root
   - err, out, log HTCondor files will be produced in the default output path
