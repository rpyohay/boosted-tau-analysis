#!/bin/bash

#command line arguments and macros
if [ $# -ne 1 ]
    then
    echo "Usage: ./getRawCutFlow.sh <version>"
    exit 0
fi
version=$1

#vector of samples
samples=( "DYJetsToLL_M-10To50" "DYJetsToLL_M-50" "data" "TTJets" "W1JetsToLNu" "W2JetsToLNu" "W3JetsToLNu" "W4JetsToLNu" "Tbar_s-channel" "Tbar_t-channel" "T_s-channel" "T_t-channel" "WZ" "WW" "ZZ" "nonIsoWData" )
nSamples=${#samples[@]}
iEndSample=`expr $nSamples - 1`

#loop over regions
for iRegion in "Iso" "NonIso" ""
  do

  #prepare name of b veto filter module for sample with no isolation cut
  bVetoFilterModuleName="${iRegion}BVetoFilter"
  if [ -z $iRegion ]
      then
      bVetoFilterModuleName="AllBVetoFilter"
  fi

  #choose the correct instance of TriggerObjectFilter and MTFilter and the right region
  TriggerObjectFilterLine=1
  MTFilterLine=1
  signFilterLine=1
  trigTauFilterLine=1
  region=A
  if [ -z $iRegion ]
      then
      TriggerObjectFilterLine=5
      MTFilterLine=3
      signFilterLine=5
      trigTauFilterLine=5
      region="A + B (e.g. no isolation cut)"
  elif [ $iRegion = "NonIso" ]
      then
      TriggerObjectFilterLine=3
      MTFilterLine=2
      trigTauFilterLine=3
      region=B
  fi

  #running totals
  nTauTot=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
  nBVetoTauTot=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
  nTrgMatchedTot=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
  nSameSgnMuTot=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
  nOppSgnTauTot=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
  nTrigLepFilterTot=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
  nHighMTTot=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
  nLowMTTot=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )

  #loop over all job files in a version directory
  for iJob in `ls -alh ${version} | grep LSFJOB | awk '{ print $9 }'`
    do
    file=${version}/${iJob}/STDOUT
    #get the name of the sample
    sample=`grep Initiating $file | sed -e "s%.*root.*user/[a-z]*/\(.*\)/data_no_selection.*%\1%" | head -n 1`

    #get the cut flow
    nTau=`grep "muHad${iRegion}TauSelector" $file | head -n 1 | awk '{ print $5 }'`
    nBVetoTau=`grep $bVetoFilterModuleName $file | head -n 1 | awk '{ print $5 }'`
    nTrgMatched=`grep TriggerObjectFilter $file | head -n $TriggerObjectFilterLine | tail -n 1 | awk '{ print $5 }'`
    nSameSgnMu=`grep "OSSFFilter${iRegion}" $file | head -n $signFilterLine | tail -n 1 | awk '{ print $5 }'`
    nOppSgnTau=`grep "SSSFFilter${iRegion}" $file | head -n $signFilterLine | tail -n 1 | awk '{ print $5 }'`
    nTrigTauFilter=`grep "trigMuonTauFilter" $file | head -n $trigTauFilterLine | tail -n 1 | awk '{ print $5 }'`
    nHighMT=`grep "highMTFilter" $file | head -n $MTFilterLine | tail -n 1 | awk '{ print $5 }'`
    nLowMT=`grep "lowMTFilter" $file | head -n $MTFilterLine | tail -n 1 | awk '{ print $5 }'`

    #update running totals
    data=`echo $sample | grep 22Jan2013-v1_AOD_skim_v3`
    nonIsoWData=`echo $sample | grep 22Jan2013-v1_AOD_skim_v4`
    if [ ${#data} -gt 0 ]
	then
	iSample=2
    elif [ ${#nonIsoWData} -gt 0 ]
	then
	iSample=$iEndSample
    else
	foundMCSample=`echo $sample | grep ${samples[0]}`
	iSample=0
	while [ $iSample -lt $nSamples ] && [ -z $foundMCSample ]
	  do
	  iSample=`expr $iSample + 1`
	  foundMCSample=`echo $sample | grep ${samples[${iSample}]}`
	done
    fi
    nTauTot[${iSample}]=`expr ${nTauTot[${iSample}]} + $nTau`
    nBVetoTauTot[${iSample}]=`expr ${nBVetoTauTot[${iSample}]} + $nBVetoTau`
    nTrgMatchedTot[${iSample}]=`expr ${nTrgMatchedTot[${iSample}]} + $nTrgMatched`
    nSameSgnMuTot[${iSample}]=`expr ${nSameSgnMuTot[${iSample}]} + $nSameSgnMu`
    nOppSgnTauTot[${iSample}]=`expr ${nOppSgnTauTot[${iSample}]} + $nOppSgnTau`
    nTrigLepFilterTot[${iSample}]=`expr ${nTrigLepFilterTot[${iSample}]} + $nTrigTauFilter`
    nHighMTTot[${iSample}]=`expr ${nHighMTTot[${iSample}]} + $nHighMT`
    nLowMTTot[${iSample}]=`expr ${nLowMTTot[${iSample}]} + $nLowMT`
  done

  #print the cut flow
  echo "Region: $region"
  echo ""
  for iSample in `seq 0 $iEndSample`
    do
    if [ ${samples[${iSample}]} != "data" ] #blinded
	then
	echo "Sample: ${samples[${iSample}]}"
	echo "HPS isolation: ${nTauTot[${iSample}]}"
	echo "b veto: ${nBVetoTauTot[${iSample}]}"
	echo "HLT matching: ${nTrgMatchedTot[${iSample}]}"
	echo "W muon, soft muon same sign: ${nSameSgnMuTot[${iSample}]}"
	echo "Soft muon, tau opposite sign: ${nOppSgnTauTot[${iSample}]}"
	echo "Lepton filter: ${nTrigLepFilterTot[${iSample}]}"
	echo "High MT: ${nHighMTTot[${iSample}]}"
	echo "Low MT: ${nLowMTTot[${iSample}]}"
	echo ""
    fi
  done
done

exit 0
