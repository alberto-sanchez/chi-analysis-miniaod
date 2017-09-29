# Ponia-OniaPhoton

This package is mean to be run using MINIAOD (2017 version and up, only)

* Setup: (it has being tested on 9_2_x should run in any of the recent cmssw releases)

```
export SCRAM_ARCH=slc6_amd64_gcc530
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_9_2_10
cd CMSSW_9_2_10/src/
cmsenv
git clone https://github.com/alberto-sanchez/chi-analysis-miniaod.git Ponia/OniaPhoton
scram b

```

* Run: (use your favorite input sample)

```
cmsRun Ponia/OniaPhoton/test/run-chib-miniaod.py (for chib reconstruction)
```

In test directory you can find other examples to run over data and mc.

Note: Be aware that data/mc previous to 2017, can be converted to MINIAOD
following instructions in cmssw pr #19304 (https://github.com/cms-sw/cmssw/pull/19304)

