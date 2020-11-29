# Ponia-OniaPhoton

This package is mean to be run using AOD for HI analysis (it may work on MINIAOD as well)

* Setup: (it has being tested on 8_0_x should run in any of the recent cmssw releases)

```
export SCRAM_ARCH=slc7_amd64_gcc530
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_8_0_35
cd CMSSW_8_0_35/src/
cmsenv
git clone https://github.com/alberto-sanchez/chi-analysis-miniaod.git -b hi_analysis Ponia/OniaPhoton
scram b

```

* Run: (use your favorite input sample)

```
cmsRun Ponia/OniaPhoton/test/run_chic_hi.py (for chic reconstruction)
```

In test directory you may find other examples to run over data and mc, as well as some crab cfgs
