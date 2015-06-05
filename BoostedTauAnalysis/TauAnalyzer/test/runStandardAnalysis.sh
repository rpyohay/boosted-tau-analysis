#!/bin/bash

#usage
if [ $# -ne 3 ]
    then
    echo "Usage: ./runStandardAnalysis.sh <version_narrow> <version_normal> <HLT_path>"
    exit 0
fi

#version
versionNarrow=$1
versionNormal=$2
HLTPath=$3

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
  root -l -b -q 'formatBoostedResonanceDataBkgPlots.C("'${versionNarrow}'", "'${versionNarrow}'", "'${MTBin[${iMTBin}]}'", 9, "'${HLTPath}'")'
  root -l -b -q 'savePlots.C("'${versionNarrow}'", "_a9", "'${MTBin[${iMTBin}]}'", true, "'${HLTPath}'")'
  root -l -b -q 'formatDataBkgPlots.C("'${versionNormal}'", "'${versionNormal}'", "'${versionNarrow}'", "'${MTBin[${iMTBin}]}'", 3)'
  root -l -b -q 'formatSigPlots.C("'${versionNormal}'", "'${versionNormal}'", "'${versionNarrow}'", "", "_a9", "'${MTBin[${iMTBin}]}'", 3)'
  root -l -b -q 'savePlots.C("'${versionNormal}'", "_a9", "'${MTBin[${iMTBin}]}'", false, "'${HLTPath}'")'
done

exit 0