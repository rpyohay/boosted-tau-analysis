#!/bin/bash

jobDir="${CMSSW_BASE}/src/BoostedTauAnalysis/TauAnalyzer/test/DIR"
fileNamePrefix="tauanalyzer_template_cfg_JOB"
outputFilePrefix="DIR_tau_analysis_JOB"
outputFile="${outputFilePrefix}.root"

cd $jobDir
eval `scramv1 runtime -sh`
cd -
cp $jobDir/$fileNamePrefix.py .
cmsRun $fileNamePrefix.py
cmsStage -f $outputFile /store/user/`whoami`
rm $outputFile $outputFilePrefix.txt

exit 0
