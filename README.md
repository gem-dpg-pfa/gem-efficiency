# VFAT_MaskMaker
## Purpose
Analyse the GEM Common muon ntuples to get list of VFATs OFF/Noisy

## Details 
Needs the `PFA_Analyzer_Utils.py` library 

Step1 -- takes as input the path of GEM Common muon ntuples file and produces a mini tree containing useful info

Step2 -- takes as input the output of step1 and produces a list of bad/noisy VFAT as txt file directly into `$DOC2_PFA/Analyzer/ExcludeMe`

## Condor submit
Condor can ease your life, parallelizing the jobs. 

By using `run_step.py` condor submit jobs are crated in the folder CondorFiles, ready to be submitted to the cluster.

Once you know which tag do your GEM Common Ntuples have, let's say 349834_Express, you just need to run
```
python run_step.py -r 349263_Express
cd CondorFiles/
condor_submit condor_step1_349263_Express.submit
Then, when the previous ends, execute:
condor_submit condor_step2_349263_Express.submit
```
