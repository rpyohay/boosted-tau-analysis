#!/bin/bash

#usage
if [ $# -ne 2 ]
    then
    echo "Usage: ./runMassPoints.sh <version_narrow> <version_normal> "
    exit 0
fi

#version
versionNarrow=$1
versionNormal=$2

#setup
eval `scramv1 runtime -sh`

#MT bin
iBeg=0
MTBin=( "_highMT" "_lowMT" )
aMass=( "_a5" "_a7" "_a9" "_a11" "_a13" "_a15" )
nMTBins=${#MTBin[@]}
iEndMTBin=`expr $nMTBins - 1`
nMass=${#aMass[@]}
iEndMass=`expr $nMass - 1`

#run--order matters here!
for iMTBin in `seq $iBeg $iEndMTBin`
  do
  for iMass in `seq $iBeg $iEndMass`
    do
    root -l -b -q 'formatSigPlots.C("'${versionNormal}'", "'${versionNormal}'", "'${versionNarrow}'", "", "'${aMass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'
    done
done

exit 0