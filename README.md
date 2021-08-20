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
1.  step1_hit_perLumi_analysis.py
Script that: 
 - Reads the GEMMuonNtuple
 - Stores #events per lumi section
 - Stores #rechits per lumi section per each endcap/layer/chamber/vfat
 - Produces plain ntuples for further analysis
Typical execution:
```sh
python step1_hit_perLumi_analysis.py --input_file=/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEMCommonNtuples/MWGR4/Run_342154/MuDPGNtuple_2021_MWGR4_5.root, --output output/step1_342154_5.root --nevents -1
```
