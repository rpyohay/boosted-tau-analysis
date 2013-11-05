#!/bin/bash

jobDir="/afs/cern.ch/user/y/yohay/CMSSW_5_3_2_patch4/src/BoostedTauAnalysis/TauSkimmer/test/DIR"
cfgFileName="CFG_FILE_NAME"
outputFileName="OUTPUT_FILE_NAME"

cd $jobDir
eval `scramv1 runtime -sh`
cd -
cp $jobDir/$cfgFileName .
#rm $jobDir/$cfgFileName*
cmsRun $cfgFileName
cmsStage -f $outputFileName /store/user/yohay/
rm $cfgFileName* $outputFileName

exit 0
