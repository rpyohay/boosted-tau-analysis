#!/bin/bash

#usage
if [ $# -ne 1 ]
    then
    echo "Usage: ./runStandardAnalysis.sh <version>"
    exit 0
fi

#version
version=$1

#setup
eval `scramv1 runtime -sh`

#MT bin
iBeg=0
MTBin=( "_highMT" "_lowMT" )
nMTBins=${#MTBin[@]}
iEndMTBin=`expr $nMTBins - 1`

#run
for iMTBin in `seq $iBeg $iEndMTBin`
  do
  root -l -b -q 'formatDataBkgPlots.C("'${version}'", "'${version}'", "'${MTBin[${iMTBin}]}'")'
  root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "", "_a9", "'${MTBin[${iMTBin}]}'")'
done

exit 0