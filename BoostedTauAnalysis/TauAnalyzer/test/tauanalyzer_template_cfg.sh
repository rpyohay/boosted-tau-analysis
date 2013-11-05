#!/bin/bash

<<<<<<< HEAD
jobDir="/afs/cern.ch/user/f/friccita/myOtherOther533Area/CMSSW_5_3_3/src/BoostedTauAnalysis/TauAnalyzer/test/DIR"
=======
jobDir="/afs/cern.ch/user/y/yohay/CMSSW_5_2_4_patch3_20Sep12/src/BoostedTauAnalysis/TauAnalyzer/test/DIR"
>>>>>>> 46247d24eb20eea62a68284a461b2dc6bfa58c65
fileNamePrefix="tauanalyzer_template_cfg_JOB"
outputFilePrefix="DIR_tau_analysis_JOB"
outputFile="${outputFilePrefix}.root"

cd $jobDir
eval `scramv1 runtime -sh`
cd -
cp $jobDir/$fileNamePrefix.py .
cmsRun $fileNamePrefix.py
<<<<<<< HEAD
cmsStage -f $outputFile /store/user/friccita/TauAna
=======
cmsStage -f $outputFile /store/user/yohay/
>>>>>>> 46247d24eb20eea62a68284a461b2dc6bfa58c65
rm $outputFile $outputFilePrefix.txt

exit 0
