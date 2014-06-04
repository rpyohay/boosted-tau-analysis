#!/bin/bash

#arguments
discriminator=$1
sample=$2
minJob=1
maxJob=202
#filePrefix=$5
#suffix=$6

#macros
#filePrefix="${sample}_${discriminator}/analyzeSelectionTemplate_${discriminator}_"
filePrefix="ZZ_TuneZ2star_8TeV_pythia6_tauola-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v3/res/CMSSW_"
suffix=".stdout"
#suffix=".txt"

#totals
nProcessedTot=0
nWMuNuTot=0
nHLTTot=0
nWMuonPTTot=0
nWMuonIsoTot=0
nJetTot=0
#nNsubjTot=0
nTauMuonPTTot=0
nTauMuonSoftTot=0
nMuHadTot=0
#nSecondJetTot=0

#loop over job outputs
for iJob in `seq $minJob $maxJob`
  do
  echo "Job #$iJob"

  #get the numbers of events passing each cut
  nProcessed=`grep IsoMu24eta2p1Selector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*\([0-9]*\)[ ]*[0-9]*.*%\1%"`
#  nWMuNu=`grep genWMuNuSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nHLT=`grep IsoMu24eta2p1Selector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nWMuonPT=`grep WMuonPTSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nWMuonIso=`grep WIsoMuonSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
#  nJet=`grep jetSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  #nNsubj=`grep NsubjFilter ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nTauMuonPT=`grep tauMuonPTSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nTauMuonSoft=`grep tauMuonSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nMuHad=`grep muHadTauSelector ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
#  nSecondJet=`grep SecondJetCut ${filePrefix}${iJob}${suffix} | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`

  #if [ $iJob == 54 ]
  #    then
  #    nProcessed=58117
  #    nWMuNu=19275
  #fi
  #if [ $iJob == 71 ]
  #    then
  #    nProcessed=58117
  #    nWMuNu=19413
  #fi
  #if [ $iJob == 81 ]
  #    then
  #    nProcessed=58117
  #    nWMuNu=19239
  #fi



  #increment totals
  nProcessedTot=`expr $nProcessedTot + $nProcessed`
#  nWMuNuTot=`expr $nWMuNuTot + $nWMuNu`
  nHLTTot=`expr $nHLTTot + $nHLT`
  nWMuonPTTot=`expr $nWMuonPTTot + $nWMuonPT`
  nWMuonIsoTot=`expr $nWMuonIsoTot + $nWMuonIso`
#  nJetTot=`expr $nJetTot + $nJet`
  #nNsubjTot=`expr $nNsubjTot + $nNsubj`
  nTauMuonPTTot=`expr $nTauMuonPTTot + $nTauMuonPT`
  nTauMuonSoftTot=`expr $nTauMuonSoftTot + $nTauMuonSoft`
  nMuHadTot=`expr $nMuHadTot + $nMuHad`
#  nSecondJetTot=`expr $nSecondJetTot + $nSecondJet`
done

#print totals
echo "nProcessedTot = $nProcessedTot"
#echo "nWMuNuTot = $nWMuNuTot"
echo "nHLTTot = $nHLTTot"
echo "nWMuonPTTot = $nWMuonPTTot"
echo "nWMuonIsoTot = $nWMuonIsoTot"
#echo "nJetTot = $nJetTot"
#echo "nNsubjTot = $nNsubjTot"
echo "nTauMuonPTTot = $nTauMuonPTTot"
echo "nTauMuonSoftTot = $nTauMuonSoftTot"
echo "nMuHadTot = $nMuHadTot"
#echo "nSecondJetTot = $nSecondJetTot"

exit 0
