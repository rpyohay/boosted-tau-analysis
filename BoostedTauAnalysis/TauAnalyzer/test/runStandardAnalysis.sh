#!/bin/bash

#usage
if [ $# -ne 2 ]
    then
    echo "Usage: ./runStandardAnalysis.sh <version_normal> <HLT_path>"
    exit 0
fi

#version
versionNormal=$1
HLTPath=$2

#setup
eval `scramv1 runtime -sh`

#MT bin
iBeg=0
MTBin=( "_highMT" "_lowMT" )
#MTBin=( "_highMT" )
#MTBin=( "_lowMT" )
nMTBins=${#MTBin[@]}
iEndMTBin=`expr $nMTBins - 1`

#run--order matters here!
for iMTBin in `seq $iBeg $iEndMTBin`
  do
  root -l -b -q 'formatDataBkgPlots.C("'${versionNormal}'", "'${versionNormal}'", "'${MTBin[${iMTBin}]}'", 5)'
  root -l -b -q 'formatSigPlots.C("'${versionNormal}'", "'${versionNormal}'", "", "_a9", "'${MTBin[${iMTBin}]}'", 5)'
#  root -l -b -q 'savePlots.C("'${versionNormal}'", "_a9", "'${MTBin[${iMTBin}]}'", false, "'${HLTPath}'")'
done

exit 0