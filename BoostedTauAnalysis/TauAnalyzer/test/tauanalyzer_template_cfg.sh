#!/bin/bash

jobDir="/afs/cern.ch/user/y/yohay/CMSSW_5_2_4_patch3_20Sep12/src/BoostedTauAnalysis/TauAnalyzer/test/DIR"
fileNamePrefix="tauanalyzer_template_cfg_JOB"
outputFilePrefix="DIR_tau_analysis_JOB"
outputFile="${outputFilePrefix}.root"

cd $jobDir
eval `scramv1 runtime -sh`
cd -
cp $jobDir/$fileNamePrefix.py .
cmsRun $fileNamePrefix.py
cmsStage -f $outputFile /store/user/yohay/
rm $outputFile $outputFilePrefix.txt

exit 0
