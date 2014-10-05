#!/bin/bash

export SCRAM_ARCH=slc5_amd64_gcc462
cd /afs/cern.ch/user/y/yohay/CMSSW_5_3_3_Git/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
mv /afs/cern.ch/user/y/yohay/CMSSW_5_3_3_Git/src/BoostedTauAnalysis/TauAnalyzer/test/DIRNAME/SCRIPTNAME_JOBNUM.* .
cmsRun SCRIPTNAME_JOBNUM.py
cmsStage -f Summer12_DR53X_NMSSMHiggsDATATIER_JOBNUM.root /store/user/yohay/DIRNAME/
rm SCRIPTNAME_JOBNUM.* Summer12_DR53X_NMSSMHiggsDATATIER_JOBNUM.root*
if [ -e RandomEngineState_JOBNUM.log ]
    then
    rm RandomEngineState_JOBNUM.log*
fi
if [ -e histProbFunction_JOBNUM.root ]
    then
    rm histProbFunction_JOBNUM.root*
fi

exit 0
