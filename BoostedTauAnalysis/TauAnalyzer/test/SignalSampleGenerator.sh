#!/bin/tcsh

setenv SCRAM_ARCH slc5_amd64_gcc462
cd /afs/cern.ch/user/y/yohay/CMSSW_5_2_4_patch3/src
eval `scramv1 runtime -csh`
cd -
source /afs/cern.ch/cms/caf/setup.csh
mv /afs/cern.ch/user/y/yohay/CMSSW_5_2_4_patch3/src/BoostedTauAnalysis/TauAnalyzer/test/DIRNAME/SCRIPTNAME_JOBNUM.* .
cmsRun SCRIPTNAME_JOBNUM.py
cmsStage -f Summer12_NMSSMHiggs_JOBNUM.root /store/user/yohay/DIRNAME/
cmsStage -f RandomEngineState_JOBNUM.log /store/user/yohay/DIRNAME/
rm SCRIPTNAME_JOBNUM.* Summer12_NMSSMHiggs_JOBNUM.root* RandomEngineState_JOBNUM.log* histProbFunction_JOBNUM.root*

exit 0
