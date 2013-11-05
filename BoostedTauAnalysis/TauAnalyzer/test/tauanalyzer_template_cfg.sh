#!/bin/bash

jobDir="/afs/cern.ch/user/f/friccita/myOtherOther533Area/CMSSW_5_3_3/src/BoostedTauAnalysis/TauAnalyzer/test/DIR"
fileNamePrefix="tauanalyzer_template_cfg_JOB"
outputFilePrefix="DIR_tau_analysis_JOB"
outputFile="${outputFilePrefix}.root"

cd $jobDir
eval `scramv1 runtime -sh`
cd -
cp $jobDir/$fileNamePrefix.py .
cmsRun $fileNamePrefix.py
cmsStage -f $outputFile /store/user/friccita/TauAna
rm $outputFile $outputFilePrefix.txt

exit 0
