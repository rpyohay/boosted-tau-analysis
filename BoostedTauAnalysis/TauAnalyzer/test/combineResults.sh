#!/bin/bash

#arguments
minJob=$1
maxJob=$2
exclString=$3

#macros
filePrefix="../../TauSkimmer/test/W4JetsToLNu_TuneZ2Star_8TeV-madgraph-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v2/res/CMSSW_"
suffix=".stdout"

#totals
nProcessedTot=0
nWMuNuTot=0
nHLTTot=0
nWMuonPTTot=0
nWMuonIsoTot=0
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

  #get the numbers of events passing each cut
#  nProcessed=`grep genWMuNuSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*\([0-9]*\)[ ]*[0-9]*.*%\1%"`
  nProcessed=`grep IsoMu24eta2p1Selector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*\([0-9]*\)[ ]*[0-9]*.*%\1%"`
#  nWMuNu=`grep genWMuNuSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nHLT=`grep IsoMu24eta2p1Selector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nWMuonPT=`grep WMuonPTSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nWMuonIso=`grep WIsoMuonSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
#  nJet=`grep jetSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nTauMuonPT=`grep tauMuonPTSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nTauMuonSoft=`grep tauMuonSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
#  nMuHadIso=`grep muHadIsoTauSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nMuHad=`grep muHadTauSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`

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

  #increment totals
  if [ $skip -eq 0 ]
      then
      nProcessedTot=`expr $nProcessedTot + $nProcessed`
#      nWMuNuTot=`expr $nWMuNuTot + $nWMuNu`
      nHLTTot=`expr $nHLTTot + $nHLT`
      nWMuonPTTot=`expr $nWMuonPTTot + $nWMuonPT`
      nWMuonIsoTot=`expr $nWMuonIsoTot + $nWMuonIso`
#  nJetTot=`expr $nJetTot + $nJet`
      nTauMuonPTTot=`expr $nTauMuonPTTot + $nTauMuonPT`
      nTauMuonSoftTot=`expr $nTauMuonSoftTot + $nTauMuonSoft`
#  nMuHadIsoTot=`expr $nMuHadIsoTot + $nMuHadIso`
      nMuHadTot=`expr $nMuHadTot + $nMuHad`
  fi
done

#print totals
echo "nProcessedTot = $nProcessedTot"
#echo "nWMuNuTot = $nWMuNuTot"
echo "nHLTTot = $nHLTTot"
echo "nWMuonPTTot = $nWMuonPTTot"
echo "nWMuonIsoTot = $nWMuonIsoTot"
#echo "nJetTot = $nJetTot"
echo "nTauMuonPTTot = $nTauMuonPTTot"
echo "nTauMuonSoftTot = $nTauMuonSoftTot"
#echo "nMuHadIsoTot = $nMuHadIsoTot"
echo "nMuHadTot = $nMuHadTot"

exit 0