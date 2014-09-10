#!/bin/bash

#usage
if [ $# -ne 1 ]
    then
    echo "Usage: ./calculateBVetoUncertainties.sh <version>"
    exit 0
fi

#version
version=$1

#a1 mass
a1Mass=( "_a5" "_a7" "_a9" "_a11" "_a13" "_a15" )
nMasses=${#a1Mass[@]}
iBeg=0
iEndMass=`expr $nMasses - 1`

#MT bin
MTBin=( "_highMT" "_lowMT" )
nMTBins=${#MTBin[@]}
iEndMTBin=`expr $nMTBins - 1`

#loop over a1 mass
for iMass in `seq $iBeg $iEndMass`
do
  echo "a1 mass: ${a1Mass[${iMass}]}"

  #loop over MT bin
  for iMTBin in `seq $iBeg $iEndMTBin`
  do
    if [ "${MTBin[${iMTBin}]}" = "_lowMT" ]
	then
	echo "MT <= 50 GeV"
    elif [ "${MTBin[${iMTBin}]}" = "_highMT" ]
	then
	echo "MT > 50 GeV"
    fi

    #get the nominal expected signal
    outputNom=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    nWh1Nom=`echo $outputNom | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    echo "No. Wh1: $nWh1Nom"

    #get the +/-1sigma b veto data/MC scale expected signal
    outputMinus1SigmaBVetoScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaBVetoScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    outputPlus1SigmaBVetoScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaBVetoScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

    #calculate the uncertainty due to b veto data/MC scale
    nWh1Minus1SigmaBVetoScale=`echo $outputMinus1SigmaBVetoScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    nWh1Plus1SigmaBVetoScale=`echo $outputPlus1SigmaBVetoScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    Wh1ErrLowBVetoScale=$(echo "scale=3; $nWh1Minus1SigmaBVetoScale - $nWh1Nom" | bc -q 2>/dev/null)
    Wh1ErrHighBVetoScale=$(echo "scale=3; $nWh1Nom - $nWh1Plus1SigmaBVetoScale" | bc -q 2>/dev/null)
    Wh1FracErrLowBVetoScale=$(echo "scale=3; $Wh1ErrLowBVetoScale / $nWh1Nom" | bc -q 2>/dev/null)
    Wh1FracErrHighBVetoScale=$(echo "scale=3; $Wh1ErrHighBVetoScale / $nWh1Nom" | bc -q 2>/dev/null)
    echo "- $Wh1ErrLowBVetoScale($Wh1FracErrLowBVetoScale%) + $Wh1ErrHighBVetoScale($Wh1FracErrHighBVetoScale%) (b veto data/MC scale)"

    #get the +/-1sigma e/gamma data/MC scale expected signal
    outputMinus1SigmaEGScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaEGScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    outputPlus1SigmaEGScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaEGScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

    #calculate the uncertainty due to e/g data/MC scale
    nWh1Minus1SigmaEGScale=`echo $outputMinus1SigmaEGScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    nWh1Plus1SigmaEGScale=`echo $outputPlus1SigmaEGScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    Wh1ErrLowEGScale=$(echo "scale=3; $nWh1Minus1SigmaEGScale - $nWh1Nom" | bc -q 2>/dev/null)
    Wh1ErrHighEGScale=$(echo "scale=3; $nWh1Nom - $nWh1Plus1SigmaEGScale" | bc -q 2>/dev/null)
    Wh1FracErrLowEGScale=$(echo "scale=3; $Wh1ErrLowEGScale / $nWh1Nom" | bc -q 2>/dev/null)
    Wh1FracErrHighEGScale=$(echo "scale=3; $Wh1ErrHighEGScale / $nWh1Nom" | bc -q 2>/dev/null)
    echo "- $Wh1ErrLowEGScale($Wh1FracErrLowEGScale%) + $Wh1ErrHighEGScale($Wh1FracErrHighEGScale%) (e/gamma data/MC scale)"
    echo "-1sigma: $nWh1Minus1SigmaEGScale, +1sigma: $nWh1Plus1SigmaEGScale"

    #get the +/-1sigma jet energy data/MC scale expected signal
    outputMinus1SigmaJES=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaJES", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    outputPlus1SigmaJES=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaJES", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

    #calculate the uncertainty due to jet energy data/MC scale
    nWh1Minus1SigmaJES=`echo $outputMinus1SigmaJES | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    nWh1Plus1SigmaJES=`echo $outputPlus1SigmaJES | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    Wh1ErrLowJES=$(echo "scale=3; $nWh1Minus1SigmaJES - $nWh1Nom" | bc -q 2>/dev/null)
    Wh1ErrHighJES=$(echo "scale=3; $nWh1Nom - $nWh1Plus1SigmaJES" | bc -q 2>/dev/null)
    Wh1FracErrLowJES=$(echo "scale=3; $Wh1ErrLowJES / $nWh1Nom" | bc -q 2>/dev/null)
    Wh1FracErrHighJES=$(echo "scale=3; $Wh1ErrHighJES / $nWh1Nom" | bc -q 2>/dev/null)
    echo "- $Wh1ErrLowJES($Wh1FracErrLowJES%) + $Wh1ErrHighJES($Wh1FracErrHighJES%) (jet energy data/MC scale)"
    echo "-1sigma: $nWh1Minus1SigmaJES, +1sigma: $nWh1Plus1SigmaJES"

    #get the +/-1sigma MC jet energy resolution smeared expected signal
    outputMinus1SigmaJER=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaJER", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    outputPlus1SigmaJER=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaJER", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

    #calculate the uncertainty due to MC jet energy resolution smearing
    nWh1Minus1SigmaJER=`echo $outputMinus1SigmaJER | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    nWh1Plus1SigmaJER=`echo $outputPlus1SigmaJER | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    Wh1ErrLowJER=$(echo "scale=3; $nWh1Minus1SigmaJER - $nWh1Nom" | bc -q 2>/dev/null)
    Wh1ErrHighJER=$(echo "scale=3; $nWh1Nom - $nWh1Plus1SigmaJER" | bc -q 2>/dev/null)
    Wh1FracErrLowJER=$(echo "scale=3; $Wh1ErrLowJER / $nWh1Nom" | bc -q 2>/dev/null)
    Wh1FracErrHighJER=$(echo "scale=3; $Wh1ErrHighJER / $nWh1Nom" | bc -q 2>/dev/null)
    echo "- $Wh1ErrLowJER($Wh1FracErrLowJER%) + $Wh1ErrHighJER($Wh1FracErrHighJER%) (jet energy MC resolution smeared)"
    echo "-1sigma: $nWh1Minus1SigmaJER, +1sigma: $nWh1Plus1SigmaJER"

    #get the +/-1sigma muon energy data/MC scale expected signal
    outputMinus1SigmaMuEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaMuEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    outputPlus1SigmaMuEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaMuEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

    #calculate the uncertainty due to muon energy data/MC scale
    nWh1Minus1SigmaMuEnergyScale=`echo $outputMinus1SigmaMuEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    nWh1Plus1SigmaMuEnergyScale=`echo $outputPlus1SigmaMuEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    Wh1ErrLowMuEnergyScale=$(echo "scale=3; $nWh1Minus1SigmaMuEnergyScale - $nWh1Nom" | bc -q 2>/dev/null)
    Wh1ErrHighMuEnergyScale=$(echo "scale=3; $nWh1Nom - $nWh1Plus1SigmaMuEnergyScale" | bc -q 2>/dev/null)
    Wh1FracErrLowMuEnergyScale=$(echo "scale=3; $Wh1ErrLowMuEnergyScale / $nWh1Nom" | bc -q 2>/dev/null)
    Wh1FracErrHighMuEnergyScale=$(echo "scale=3; $Wh1ErrHighMuEnergyScale / $nWh1Nom" | bc -q 2>/dev/null)
    echo "- $Wh1ErrLowMuEnergyScale($Wh1FracErrLowMuEnergyScale%) + $Wh1ErrHighMuEnergyScale($Wh1FracErrHighMuEnergyScale%) (muon energy data/MC scale)"
    echo "-1sigma: $nWh1Minus1SigmaMuEnergyScale, +1sigma: $nWh1Plus1SigmaMuEnergyScale"

    #get the +/-1sigma tau energy data/MC scale expected signal
    outputMinus1SigmaTauEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaTauEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    outputPlus1SigmaTauEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaTauEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

    #calculate the uncertainty due to tau energy data/MC scale
    nWh1Minus1SigmaTauEnergyScale=`echo $outputMinus1SigmaTauEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    nWh1Plus1SigmaTauEnergyScale=`echo $outputPlus1SigmaTauEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    Wh1ErrLowTauEnergyScale=$(echo "scale=3; $nWh1Minus1SigmaTauEnergyScale - $nWh1Nom" | bc -q 2>/dev/null)
    Wh1ErrHighTauEnergyScale=$(echo "scale=3; $nWh1Nom - $nWh1Plus1SigmaTauEnergyScale" | bc -q 2>/dev/null)
    Wh1FracErrLowTauEnergyScale=$(echo "scale=3; $Wh1ErrLowTauEnergyScale / $nWh1Nom" | bc -q 2>/dev/null)
    Wh1FracErrHighTauEnergyScale=$(echo "scale=3; $Wh1ErrHighTauEnergyScale / $nWh1Nom" | bc -q 2>/dev/null)
    echo "- $Wh1ErrLowTauEnergyScale($Wh1FracErrLowTauEnergyScale%) + $Wh1ErrHighTauEnergyScale($Wh1FracErrHighTauEnergyScale%) (tau energy data/MC scale)"
    echo "-1sigma: $nWh1Minus1SigmaTauEnergyScale, +1sigma: $nWh1Plus1SigmaTauEnergyScale"

    #get the +/-1sigma unclustered energy data/MC scale expected signal
    outputMinus1SigmaUnclusteredEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaUnclusteredEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    outputPlus1SigmaUnclusteredEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaUnclusteredEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

    #calculate the uncertainty due to unclustered energy data/MC scale
    nWh1Minus1SigmaUnclusteredEnergyScale=`echo $outputMinus1SigmaUnclusteredEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    nWh1Plus1SigmaUnclusteredEnergyScale=`echo $outputPlus1SigmaUnclusteredEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
    Wh1ErrLowUnclusteredEnergyScale=$(echo "scale=3; $nWh1Minus1SigmaUnclusteredEnergyScale - $nWh1Nom" | bc -q 2>/dev/null)
    Wh1ErrHighUnclusteredEnergyScale=$(echo "scale=3; $nWh1Nom - $nWh1Plus1SigmaUnclusteredEnergyScale" | bc -q 2>/dev/null)
    Wh1FracErrLowUnclusteredEnergyScale=$(echo "scale=3; $Wh1ErrLowUnclusteredEnergyScale / $nWh1Nom" | bc -q 2>/dev/null)
    Wh1FracErrHighUnclusteredEnergyScale=$(echo "scale=3; $Wh1ErrHighUnclusteredEnergyScale / $nWh1Nom" | bc -q 2>/dev/null)
    echo "- $Wh1ErrLowUnclusteredEnergyScale($Wh1FracErrLowUnclusteredEnergyScale%) + $Wh1ErrHighUnclusteredEnergyScale($Wh1FracErrHighUnclusteredEnergyScale%) (unclustered energy data/MC scale)"
    echo "-1sigma: $nWh1Minus1SigmaUnclusteredEnergyScale, +1sigma: $nWh1Plus1SigmaUnclusteredEnergyScale"

    #add energy scale and resolution uncertainties in quadrature
    #(wait on this until the individual uncertainties are understood)

    #temporary until gg fusion skims are done
    if [ "${a1Mass[${iMass}]}" = "_a9" ]
	then
	nggNom=`echo $outputNom | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
	echo "No. gg: $nggNom"

        #calculate the uncertainty due to b veto data/MC scale
	nggMinus1SigmaBVetoScale=`echo $outputMinus1SigmaBVetoScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
	nggPlus1SigmaBVetoScale=`echo $outputPlus1SigmaBVetoScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
	ggErrLowBVetoScale=$(echo "scale=3; $nggMinus1SigmaBVetoScale - $nggNom" | bc -q 2>/dev/null)
	ggErrHighBVetoScale=$(echo "scale=3; $nggNom - $nggPlus1SigmaBVetoScale" | bc -q 2>/dev/null)
	ggFracErrLowBVetoScale=$(echo "scale=3; $ggErrLowBVetoScale / $nggNom" | bc -q 2>/dev/null)
	ggFracErrHighBVetoScale=$(echo "scale=3; $ggErrHighBVetoScale / $nggNom" | bc -q 2>/dev/null)
	echo "- $ggErrLowBVetoScale($ggFracErrLowBVetoScale%) + $ggErrHighBVetoScale($ggFracErrHighBVetoScale%) (b veto data/MC scale)"

        #calculate the uncertainty due to e/g data/MC scale
	nggMinus1SigmaEGScale=`echo $outputMinus1SigmaEGScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
	nggPlus1SigmaEGScale=`echo $outputPlus1SigmaEGScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
	ggErrLowEGScale=$(echo "scale=3; $nggMinus1SigmaEGScale - $nggNom" | bc -q 2>/dev/null)
	ggErrHighEGScale=$(echo "scale=3; $nggNom - $nggPlus1SigmaEGScale" | bc -q 2>/dev/null)
	ggFracErrLowEGScale=$(echo "scale=3; $ggErrLowEGScale / $nggNom" | bc -q 2>/dev/null)
	ggFracErrHighEGScale=$(echo "scale=3; $ggErrHighEGScale / $nggNom" | bc -q 2>/dev/null)
	echo "- $ggErrLowEGScale($ggFracErrLowEGScale%) + $ggErrHighEGScale($ggFracErrHighEGScale%) (e/gamma data/MC scale)"
	echo "-1sigma: $nggMinus1SigmaEGScale, +1sigma: $nggPlus1SigmaEGScale"


        #get the +/-1sigma jet energy data/MC scale expected signal
        outputMinus1SigmaJES=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaJES", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
        outputPlus1SigmaJES=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaJES", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

        #calculate the uncertainty due to jet energy data/MC scale
        nggMinus1SigmaJES=`echo $outputMinus1SigmaJES | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
        nggPlus1SigmaJES=`echo $outputPlus1SigmaJES | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
        ggErrLowJES=$(echo "scale=3; $nggMinus1SigmaJES - $nggNom" | bc -q 2>/dev/null)
        ggErrHighJES=$(echo "scale=3; $nggNom - $nggPlus1SigmaJES" | bc -q 2>/dev/null)
        ggFracErrLowJES=$(echo "scale=3; $ggErrLowJES / $nggNom" | bc -q 2>/dev/null)
        ggFracErrHighJES=$(echo "scale=3; $ggErrHighJES / $nggNom" | bc -q 2>/dev/null)
        echo "- $ggErrLowJES($ggFracErrLowJES%) + $ggErrHighJES($ggFracErrHighJES%) (jet energy data/MC scale)"
	echo "-1sigma: $nggMinus1SigmaJES, +1sigma: $nggPlus1SigmaJES"

        #get the +/-1sigma MC jet energy resolution smeared expected signal
        outputMinus1SigmaJER=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaJER", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
        outputPlus1SigmaJER=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaJER", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

        #calculate the uncertainty due to MC jet energy resolution smearing
        nggMinus1SigmaJER=`echo $outputMinus1SigmaJER | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
        nggPlus1SigmaJER=`echo $outputPlus1SigmaJER | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
        ggErrLowJER=$(echo "scale=3; $nggMinus1SigmaJER - $nggNom" | bc -q 2>/dev/null)
        ggErrHighJER=$(echo "scale=3; $nggNom - $nggPlus1SigmaJER" | bc -q 2>/dev/null)
        ggFracErrLowJER=$(echo "scale=3; $ggErrLowJER / $nggNom" | bc -q 2>/dev/null)
        ggFracErrHighJER=$(echo "scale=3; $ggErrHighJER / $nggNom" | bc -q 2>/dev/null)
        echo "- $ggErrLowJER($ggFracErrLowJER%) + $ggErrHighJER($ggFracErrHighJER%) (jet energy MC resolution smeared)"
	echo "-1sigma: $nggMinus1SigmaJER, +1sigma: $nggPlus1SigmaJER"

        #get the +/-1sigma muon energy data/MC scale expected signal
        outputMinus1SigmaMuEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaMuEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
        outputPlus1SigmaMuEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaMuEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

        #calculate the uncertainty due to muon energy data/MC scale
        nggMinus1SigmaMuEnergyScale=`echo $outputMinus1SigmaMuEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
        nggPlus1SigmaMuEnergyScale=`echo $outputPlus1SigmaMuEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
        ggErrLowMuEnergyScale=$(echo "scale=3; $nggMinus1SigmaMuEnergyScale - $nggNom" | bc -q 2>/dev/null)
        ggErrHighMuEnergyScale=$(echo "scale=3; $nggNom - $nggPlus1SigmaMuEnergyScale" | bc -q 2>/dev/null)
        ggFracErrLowMuEnergyScale=$(echo "scale=3; $ggErrLowMuEnergyScale / $nggNom" | bc -q 2>/dev/null)
        ggFracErrHighMuEnergyScale=$(echo "scale=3; $ggErrHighMuEnergyScale / $nggNom" | bc -q 2>/dev/null)
        echo "- $ggErrLowMuEnergyScale($ggFracErrLowMuEnergyScale%) + $ggErrHighMuEnergyScale($ggFracErrHighMuEnergyScale%) (muon energy data/MC scale)"
	echo "-1sigma: $nggMinus1SigmaMuEnergyScale, +1sigma: $nggPlus1SigmaMuEnergyScale"

        #get the +/-1sigma tau energy data/MC scale expected signal
        outputMinus1SigmaTauEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaTauEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
        outputPlus1SigmaTauEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaTauEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

        #calculate the uncertainty due to tau energy data/MC scale
        nggMinus1SigmaTauEnergyScale=`echo $outputMinus1SigmaTauEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
        nggPlus1SigmaTauEnergyScale=`echo $outputPlus1SigmaTauEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
        ggErrLowTauEnergyScale=$(echo "scale=3; $nggMinus1SigmaTauEnergyScale - $nggNom" | bc -q 2>/dev/null)
        ggErrHighTauEnergyScale=$(echo "scale=3; $nggNom - $nggPlus1SigmaTauEnergyScale" | bc -q 2>/dev/null)
        ggFracErrLowTauEnergyScale=$(echo "scale=3; $ggErrLowTauEnergyScale / $nggNom" | bc -q 2>/dev/null)
        ggFracErrHighTauEnergyScale=$(echo "scale=3; $ggErrHighTauEnergyScale / $nggNom" | bc -q 2>/dev/null)
        echo "- $ggErrLowTauEnergyScale($ggFracErrLowTauEnergyScale%) + $ggErrHighTauEnergyScale($ggFracErrHighTauEnergyScale%) (tau energy data/MC scale)"
	echo "-1sigma: $nggMinus1SigmaTauEnergyScale, +1sigma: $nggPlus1SigmaTauEnergyScale"

        #get the +/-1sigma unclustered energy data/MC scale expected signal
        outputMinus1SigmaUnclusteredEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaUnclusteredEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
        outputPlus1SigmaUnclusteredEnergyScale=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaUnclusteredEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`

        #calculate the uncertainty due to unclustered energy data/MC scale
        nggMinus1SigmaUnclusteredEnergyScale=`echo $outputMinus1SigmaUnclusteredEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
        nggPlus1SigmaUnclusteredEnergyScale=`echo $outputPlus1SigmaUnclusteredEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
        ggErrLowUnclusteredEnergyScale=$(echo "scale=3; $nggMinus1SigmaUnclusteredEnergyScale - $nggNom" | bc -q 2>/dev/null)
        ggErrHighUnclusteredEnergyScale=$(echo "scale=3; $nggNom - $nggPlus1SigmaUnclusteredEnergyScale" | bc -q 2>/dev/null)
        ggFracErrLowUnclusteredEnergyScale=$(echo "scale=3; $ggErrLowUnclusteredEnergyScale / $nggNom" | bc -q 2>/dev/null)
        ggFracErrHighUnclusteredEnergyScale=$(echo "scale=3; $ggErrHighUnclusteredEnergyScale / $nggNom" | bc -q 2>/dev/null)
        echo "- $ggErrLowUnclusteredEnergyScale($ggFracErrLowUnclusteredEnergyScale%) + $ggErrHighUnclusteredEnergyScale($ggFracErrHighUnclusteredEnergyScale%) (unclustered energy data/MC scale)"
	echo "-1sigma: $nggMinus1SigmaUnclusteredEnergyScale, +1sigma: $nggPlus1SigmaUnclusteredEnergyScale"
    fi
  done
  echo ""
done

exit 0