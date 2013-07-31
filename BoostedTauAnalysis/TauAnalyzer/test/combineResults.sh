#!/bin/bash

#arguments
discriminator=$1
sample=$2
minJob=$3
maxJob=$4
#filePrefix=$5
#suffix=$6

#macros
#filePrefix="${sample}_${discriminator}/analyzeSelectionTemplate_${discriminator}_"
filePrefix="../../TauSkimmer/test/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v2/res/CMSSW_"
suffix=".stdout"
#suffix=".txt"

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

#loop over job outputs
for iJob in `seq $minJob $maxJob`
  do
  echo "Job #$iJob"

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

  #increment totals
  #jobs 710, 843, and 850 are status 0 in CRAB but no file can be found for them
  #job 173 is done but wasn't included in run over files in my EOS area
#  if [ "$nWMuNu" != "" ] && [ "$nHLT" != "" ] && [ "$nWMuonPT" != "" ] && [ "$nWMuonIso" != "" ] && [ $iJob -ne 70 ] && [ $iJob -ne 135 ] && [ $iJob -ne 173 ] && [ $iJob -ne 251 ] && [ $iJob -ne 326 ] && [ $iJob -ne 390 ] && [ $iJob -ne 501 ] && [ $iJob -ne 546 ] && [ $iJob -ne 568 ] && [ $iJob -ne 710 ] && [ $iJob -ne 769 ] && [ $iJob -ne 835 ] && [ $iJob -ne 836 ] && [ $iJob -ne 843 ] && [ $iJob -ne 846 ] && [ $iJob -ne 850 ] && [ $iJob -ne 892 ] && [ $iJob -ne 1131 ]
#      then
      nProcessedTot=`expr $nProcessedTot + $nProcessed`
#      nWMuNuTot=`expr $nWMuNuTot + $nWMuNu`
      nHLTTot=`expr $nHLTTot + $nHLT`
      nWMuonPTTot=`expr $nWMuonPTTot + $nWMuonPT`
      nWMuonIsoTot=`expr $nWMuonIsoTot + $nWMuonIso`
#  fi
#  nJetTot=`expr $nJetTot + $nJet`
  nTauMuonPTTot=`expr $nTauMuonPTTot + $nTauMuonPT`
  nTauMuonSoftTot=`expr $nTauMuonSoftTot + $nTauMuonSoft`
#  nMuHadIsoTot=`expr $nMuHadIsoTot + $nMuHadIso`
  nMuHadTot=`expr $nMuHadTot + $nMuHad`
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