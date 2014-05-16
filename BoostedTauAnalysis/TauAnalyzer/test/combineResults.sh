#!/bin/bash

#usage
if [ $# -ne 4 ]
    then
    echo "Usage: ./combineResults.sh minJob maxJob exclString dataDir"
    exit 0
fi

#arguments
minJob=$1
maxJob=$2
exclString=$3
dataDir=$4

#macros
filePrefix="${dataDir}/res/CMSSW_"
suffix=".stdout"
isPhotonDataset=`echo $dataDir | grep Photon`

#totals
nProcessedTot=0
nWMuNuTot=0
nHLTTot=0
nWMuonPTTot=0
nWMuonOrPhotonIsoTot=0
nJetTot=0
nTauMuonPTTot=0
nTauMuonSoftTot=0
nMuHadIsoTot=0
nMuHadTot=0

#jobs to exclude
excludeList=`echo $exclString | sed -e "s%--exclude=\(.*\)%\1%"`
iExclude=1
go=1
while [ $go -eq 1 ]
  do
  jobToExclude=`echo $excludeList | sed -e "s%\([0-9]*\).*%\1%"`
  excludeArray[$iExclude]=$jobToExclude
  preList="$excludeList"
  excludeList=`echo $excludeList | sed -e "s%[0-9]*,\(.*\)%\1%"`
  postList="$excludeList"
  iExclude=`expr $iExclude + 1`
  if [ "$preList" = "$postList" ]
      then
      go=0
  fi
done

#loop over job outputs
for iJob in `seq $minJob $maxJob`
  do

  #determine if this job should be skipped in the sum
  skip=0
  iSkip=1
  while [ $skip -eq 0 ] && [ $iSkip -le ${#excludeArray[@]} ]
    do
    if [ $iJob -eq ${excludeArray[$iSkip]} ]
	then
	skip=1
    fi
    iSkip=`expr $iSkip + 1`
  done

  #get the numbers of events passing each cut
  if [ $skip -eq 0 ]
      then
      if [ "$isPhotonDataset" = "" ]
	  then
	  nProcessed=`grep IsoMu24eta2p1Selector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*\([0-9]*\)[ ]*[0-9]*.*%\1%"`
	  nHLT=`grep IsoMu24eta2p1Selector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
	  nWMuonPT=`grep WMuonPTSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
	  nWMuonOrPhotonIso=`grep WIsoMuonSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
      else
	  nProcessed=`grep singlePhotonHLTSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*\([0-9]*\)[ ]*[0-9]*.*%\1%"`
	  nHLT=`grep singlePhotonHLTSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
	  nWMuonOrPhotonIso=`grep photonSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
      fi
      #nJet=`grep jetSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
      nTauMuonPT=`grep tauMuonPTSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
      nTauMuonSoft=`grep tauMuonSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
      #nMuHadIso=`grep muHadIsoTauSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
      nMuHad=`grep muHadTauSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`

      #increment totals
      nProcessedTot=`expr $nProcessedTot + $nProcessed`
      nHLTTot=`expr $nHLTTot + $nHLT`
      if [ "$isPhotonDataset" = "" ]
	  then
	  nWMuonPTTot=`expr $nWMuonPTTot + $nWMuonPT`
      fi
      nWMuonOrPhotonIsoTot=`expr $nWMuonOrPhotonIsoTot + $nWMuonOrPhotonIso`
      #nJetTot=`expr $nJetTot + $nJet`
      nTauMuonPTTot=`expr $nTauMuonPTTot + $nTauMuonPT`
      nTauMuonSoftTot=`expr $nTauMuonSoftTot + $nTauMuonSoft`
      #nMuHadIsoTot=`expr $nMuHadIsoTot + $nMuHadIso`
      nMuHadTot=`expr $nMuHadTot + $nMuHad`
  fi
done

#print totals
echo "nProcessedTot = $nProcessedTot"
echo "nHLTTot = $nHLTTot"
if [ "$isPhotonDataset" = "" ]
    then
    echo "nWMuonPTTot = $nWMuonPTTot"
fi
echo "nWMuonOrPhotonIsoTot = $nWMuonOrPhotonIsoTot"
#echo "nJetTot = $nJetTot"
echo "nTauMuonPTTot = $nTauMuonPTTot"
echo "nTauMuonSoftTot = $nTauMuonSoftTot"
#echo "nMuHadIsoTot = $nMuHadIsoTot"
echo "nMuHadTot = $nMuHadTot"

exit 0