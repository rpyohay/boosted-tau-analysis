#!/bin/bash

#usage
if [ $# -ne 2 ]
    then
    echo "Usage: ./calculateBVetoUncertainties.sh <version> <narrowVersion>"
    exit 0
fi

#version
version=$1

#narrowbins version
narrowVersion=$2

#a1 mass
a1Mass=( "_a5" "_a7" "_a9" "_a11" "_a13" "_a15" ) #ggH
#a1Mass=( "_a5" "_a7" "_a9" "_a11" "_a13" ) #ggH
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

    ### WH ###

    #get the nominal expected signal
    outputNom=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    nWh1Nom=`echo $outputNom | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    echo "No. Wh1 (nominal): $nWh1Nom"

    #get the +/-1sigma b veto data/MC scale expected signal
    outputMinus1SigmaBVetoScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaBVetoScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    outputPlus1SigmaBVetoScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaBVetoScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

    #initialize fractional errors to 0--only set nonzero if nWh1Nom != 0
    Wh1FracErrLowBVetoScale=0
    Wh1FracErrHighBVetoScale=0

    #calculate the uncertainty due to b veto data/MC scale
    nWh1Minus1SigmaBVetoScale=`echo $outputMinus1SigmaBVetoScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    nWh1Plus1SigmaBVetoScale=`echo $outputPlus1SigmaBVetoScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    Wh1ErrLowBVetoScale=$(echo "scale=4; $nWh1Nom - $nWh1Minus1SigmaBVetoScale" | bc -q 2>/dev/null)
    Wh1ErrHighBVetoScale=$(echo "scale=4; $nWh1Nom - $nWh1Plus1SigmaBVetoScale" | bc -q 2>/dev/null)
    if [ "$nWh1Nom" != "0" ]
	then
	Wh1FracErrLowBVetoScale=$(echo "scale=4; ($Wh1ErrLowBVetoScale / $nWh1Nom) * 100" | bc -q 2>/dev/null)
	Wh1FracErrHighBVetoScale=$(echo "scale=4; ($Wh1ErrHighBVetoScale / $nWh1Nom) * 100" | bc -q 2>/dev/null)
    fi
    echo "$Wh1ErrLowBVetoScale($Wh1FracErrLowBVetoScale%)/$Wh1ErrHighBVetoScale($Wh1FracErrHighBVetoScale%) (b veto data/MC scale)"

    #get the central expected signal for MET uncertainties
    outputCen=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_central", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    nWh1Cen=`echo $outputCen | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    echo "No. Wh1 (central value for MET uncertainties): $nWh1Cen"

    #get the +/-1sigma e/gamma data/MC scale expected signal
    outputMinus1SigmaEGScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaEGScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    outputPlus1SigmaEGScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaEGScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

    #initialize fractional errors to 0--only set nonzero if nWh1Cen != 0
    Wh1FracErrLowEGScale=0
    Wh1FracErrHighEGScale=0
    Wh1FracErrLowJES=0
    Wh1FracErrHighJES=0
    Wh1FracErrLowJER=0
    Wh1FracErrHighJER=0
    Wh1FracErrLowMuEnergyScale=0
    Wh1FracErrHighMuEnergyScale=0
    Wh1FracErrLowTauEnergyScale=0
    Wh1FracErrHighTauEnergyScale=0
    Wh1FracErrLowUnclusteredEnergyScale=0
    Wh1FracErrHighUnclusteredEnergyScale=0

    #calculate the uncertainty due to e/g data/MC scale
    nWh1Minus1SigmaEGScale=`echo $outputMinus1SigmaEGScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    nWh1Plus1SigmaEGScale=`echo $outputPlus1SigmaEGScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    Wh1ErrLowEGScale=$(echo "scale=4; $nWh1Cen - $nWh1Minus1SigmaEGScale" | bc -q 2>/dev/null)
    Wh1ErrHighEGScale=$(echo "scale=4; $nWh1Cen - $nWh1Plus1SigmaEGScale" | bc -q 2>/dev/null)
    if [ "$nWh1Cen" != "0" ]
	then
	Wh1FracErrLowEGScale=$(echo "scale=4; ($Wh1ErrLowEGScale / $nWh1Cen) * 100" | bc -q 2>/dev/null)
	Wh1FracErrHighEGScale=$(echo "scale=4; ($Wh1ErrHighEGScale / $nWh1Cen) * 100" | bc -q 2>/dev/null)
    fi
    echo "$Wh1ErrLowEGScale($Wh1FracErrLowEGScale%)/$Wh1ErrHighEGScale($Wh1FracErrHighEGScale%) (e/gamma data/MC scale)"
    echo "-1sigma: $nWh1Minus1SigmaEGScale, +1sigma: $nWh1Plus1SigmaEGScale"

    #get the +/-1sigma jet energy data/MC scale expected signal
    outputMinus1SigmaJES=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaJES", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    outputPlus1SigmaJES=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaJES", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

    #calculate the uncertainty due to jet energy data/MC scale
    nWh1Minus1SigmaJES=`echo $outputMinus1SigmaJES | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    nWh1Plus1SigmaJES=`echo $outputPlus1SigmaJES | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    Wh1ErrLowJES=$(echo "scale=4; $nWh1Cen - $nWh1Minus1SigmaJES" | bc -q 2>/dev/null)
    Wh1ErrHighJES=$(echo "scale=4; $nWh1Cen - $nWh1Plus1SigmaJES" | bc -q 2>/dev/null)
    if [ "$nWh1Cen" != "0" ]
	then
	Wh1FracErrLowJES=$(echo "scale=4; ($Wh1ErrLowJES / $nWh1Cen) * 100" | bc -q 2>/dev/null)
	Wh1FracErrHighJES=$(echo "scale=4; ($Wh1ErrHighJES / $nWh1Cen) * 100" | bc -q 2>/dev/null)
    fi
    echo "$Wh1ErrLowJES($Wh1FracErrLowJES%)/$Wh1ErrHighJES($Wh1FracErrHighJES%) (jet energy data/MC scale)"
    echo "-1sigma: $nWh1Minus1SigmaJES, +1sigma: $nWh1Plus1SigmaJES"

#v6 WH and v7 ggH skims don't apply JER smearing
#    #get the +/-1sigma MC jet energy resolution smeared expected signal
#    outputMinus1SigmaJER=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaJER", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
#    outputPlus1SigmaJER=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaJER", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

#    #calculate the uncertainty due to MC jet energy resolution smearing
#    nWh1Minus1SigmaJER=`echo $outputMinus1SigmaJER | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
#    nWh1Plus1SigmaJER=`echo $outputPlus1SigmaJER | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
#    Wh1ErrLowJER=$(echo "scale=4; $nWh1Cen - $nWh1Minus1SigmaJER" | bc -q 2>/dev/null)
#    Wh1ErrHighJER=$(echo "scale=4; $nWh1Cen - $nWh1Plus1SigmaJER" | bc -q 2>/dev/null)
#    if [ "$nWh1Cen" != "0" ]
#	then
#	Wh1FracErrLowJER=$(echo "scale=4; ($Wh1ErrLowJER / $nWh1Cen) * 100" | bc -q 2>/dev/null)
#	Wh1FracErrHighJER=$(echo "scale=4; ($Wh1ErrHighJER / $nWh1Cen) * 100" | bc -q 2>/dev/null)
#    fi
#    echo "$Wh1ErrLowJER($Wh1FracErrLowJER%)/$Wh1ErrHighJER($Wh1FracErrHighJER%) (jet energy MC resolution smeared)"
#    echo "-1sigma: $nWh1Minus1SigmaJER, +1sigma: $nWh1Plus1SigmaJER"

    #get the +/-1sigma muon energy data/MC scale expected signal
    outputMinus1SigmaMuEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaMuEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    outputPlus1SigmaMuEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaMuEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

    #calculate the uncertainty due to muon energy data/MC scale
    nWh1Minus1SigmaMuEnergyScale=`echo $outputMinus1SigmaMuEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    nWh1Plus1SigmaMuEnergyScale=`echo $outputPlus1SigmaMuEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    Wh1ErrLowMuEnergyScale=$(echo "scale=4; $nWh1Cen - $nWh1Minus1SigmaMuEnergyScale" | bc -q 2>/dev/null)
    Wh1ErrHighMuEnergyScale=$(echo "scale=4; $nWh1Cen - $nWh1Plus1SigmaMuEnergyScale" | bc -q 2>/dev/null)
    if [ "$nWh1Cen" != "0" ]
	then
	Wh1FracErrLowMuEnergyScale=$(echo "scale=4; ($Wh1ErrLowMuEnergyScale / $nWh1Cen) * 100" | bc -q 2>/dev/null)
	Wh1FracErrHighMuEnergyScale=$(echo "scale=4; ($Wh1ErrHighMuEnergyScale / $nWh1Cen) * 100" | bc -q 2>/dev/null)
    fi
    echo "$Wh1ErrLowMuEnergyScale($Wh1FracErrLowMuEnergyScale%)/$Wh1ErrHighMuEnergyScale($Wh1FracErrHighMuEnergyScale%) (muon energy data/MC scale)"
    echo "-1sigma: $nWh1Minus1SigmaMuEnergyScale, +1sigma: $nWh1Plus1SigmaMuEnergyScale"

    #get the +/-1sigma tau energy data/MC scale expected signal
    outputMinus1SigmaTauEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaTauEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    outputPlus1SigmaTauEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaTauEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

    #calculate the uncertainty due to tau energy data/MC scale
    nWh1Minus1SigmaTauEnergyScale=`echo $outputMinus1SigmaTauEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    nWh1Plus1SigmaTauEnergyScale=`echo $outputPlus1SigmaTauEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    Wh1ErrLowTauEnergyScale=$(echo "scale=4; $nWh1Cen - $nWh1Minus1SigmaTauEnergyScale" | bc -q 2>/dev/null)
    Wh1ErrHighTauEnergyScale=$(echo "scale=4; $nWh1Cen - $nWh1Plus1SigmaTauEnergyScale" | bc -q 2>/dev/null)
    if [ "$nWh1Cen" != "0" ]
	then
	Wh1FracErrLowTauEnergyScale=$(echo "scale=4; ($Wh1ErrLowTauEnergyScale / $nWh1Cen) * 100" | bc -q 2>/dev/null)
	Wh1FracErrHighTauEnergyScale=$(echo "scale=4; ($Wh1ErrHighTauEnergyScale / $nWh1Cen) * 100" | bc -q 2>/dev/null)
    fi
    echo "$Wh1ErrLowTauEnergyScale($Wh1FracErrLowTauEnergyScale%)/$Wh1ErrHighTauEnergyScale($Wh1FracErrHighTauEnergyScale%) (tau energy data/MC scale)"
    echo "-1sigma: $nWh1Minus1SigmaTauEnergyScale, +1sigma: $nWh1Plus1SigmaTauEnergyScale"

    #get the +/-1sigma unclustered energy data/MC scale expected signal
    outputMinus1SigmaUnclusteredEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaUnclusteredEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    outputPlus1SigmaUnclusteredEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaUnclusteredEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

    #calculate the uncertainty due to unclustered energy data/MC scale
    nWh1Minus1SigmaUnclusteredEnergyScale=`echo $outputMinus1SigmaUnclusteredEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    nWh1Plus1SigmaUnclusteredEnergyScale=`echo $outputPlus1SigmaUnclusteredEnergyScale | sed -e "s%.*Region A signal 0, m > 4: \([0-9\.]*\).*%\1%"`
    Wh1ErrLowUnclusteredEnergyScale=$(echo "scale=4; $nWh1Cen - $nWh1Minus1SigmaUnclusteredEnergyScale" | bc -q 2>/dev/null)
    Wh1ErrHighUnclusteredEnergyScale=$(echo "scale=4; $nWh1Cen - $nWh1Plus1SigmaUnclusteredEnergyScale" | bc -q 2>/dev/null)
    if [ "$nWh1Cen" != "0" ]
	then
	Wh1FracErrLowUnclusteredEnergyScale=$(echo "scale=4; ($Wh1ErrLowUnclusteredEnergyScale / $nWh1Cen) * 100" | bc -q 2>/dev/null)
	Wh1FracErrHighUnclusteredEnergyScale=$(echo "scale=4; ($Wh1ErrHighUnclusteredEnergyScale / $nWh1Cen) * 100" | bc -q 2>/dev/null)
    fi
    echo "$Wh1ErrLowUnclusteredEnergyScale($Wh1FracErrLowUnclusteredEnergyScale%)/$Wh1ErrHighUnclusteredEnergyScale($Wh1FracErrHighUnclusteredEnergyScale%) (unclustered energy data/MC scale)"
    echo "-1sigma: $nWh1Minus1SigmaUnclusteredEnergyScale, +1sigma: $nWh1Plus1SigmaUnclusteredEnergyScale"

    #get the larger of the +/-1sigma e/g scale shifts
    Wh1FracErrEGScale=$Wh1FracErrLowEGScale
    isLargerWh1EGScale=`echo ${Wh1FracErrHighEGScale#-}'>'${Wh1FracErrLowEGScale#-} | bc -l`
    echo $Wh1FracErrEGScale
    echo $isLargerWh1EGScale
    if [ "$isLargerWh1EGScale" = "1" ]
    then
	Wh1FracErrEGScale=$Wh1FracErrHighEGScale
    fi

    #get the larger of the +/-1sigma JES shifts
    Wh1FracErrJES=$Wh1FracErrLowJES
    isLargerWh1JES=`echo ${Wh1FracErrHighJES#-}'>'${Wh1FracErrLowJES#-} | bc -l`
    if [ "$isLargerWh1JES" = "1" ]
    then
	Wh1FracErrJES=$Wh1FracErrHighJES
    fi

    #get the larger of the +/-1sigma JER shifts
    Wh1FracErrJER=$Wh1FracErrLowJER
    isLargerWh1JER=`echo ${Wh1FracErrHighJER#-}'>'${Wh1FracErrLowJER#-} | bc -l`
    if [ "$isLargerWh1JER" = "1" ]
    then
	Wh1FracErrJER=$Wh1FracErrHighJER
    fi

    #get the larger of the +/-1sigma muon energy scale shifts
    Wh1FracErrMuEnergyScale=$Wh1FracErrLowMuEnergyScale
    isLargerWh1MuEnergyScale=`echo ${Wh1FracErrHighMuEnergyScale#-}'>'${Wh1FracErrLowMuEnergyScale#-} | bc -l`
    if [ "$isLargerWh1MuEnergyScale" = "1" ]
    then
	Wh1FracErrMuEnergyScale=$Wh1FracErrHighMuEnergyScale
    fi

    #get the larger of the +/-1sigma tau energy scale shifts
    Wh1FracErrTauEnergyScale=$Wh1FracErrLowTauEnergyScale
    isLargerWh1TauEnergyScale=`echo ${Wh1FracErrHighTauEnergyScale#-}'>'${Wh1FracErrLowTauEnergyScale#-} | bc -l`
    if [ "$isLargerWh1TauEnergyScale" = "1" ]
    then
	Wh1FracErrTauEnergyScale=$Wh1FracErrHighTauEnergyScale
    fi

    #get the larger of the +/-1sigma unclustered energy scale shifts
    Wh1FracErrUnclusteredEnergyScale=$Wh1FracErrLowUnclusteredEnergyScale
    isLargerWh1UnclusteredEnergyScale=`echo ${Wh1FracErrHighUnclusteredEnergyScale#-}'>'${Wh1FracErrLowUnclusteredEnergyScale#-} | bc -l`
    if [ "$isLargerWh1UnclusteredEnergyScale" = "1" ]
    then
	Wh1FracErrUnclusteredEnergyScale=$Wh1FracErrHighUnclusteredEnergyScale
    fi

    #add energy scale and resolution uncertainties in quadrature
    #use the larger of the +/-1sigma shifts as the symmetric error
    Wh1METFracErr=$(echo "scale=4; sqrt(($Wh1FracErrEGScale * $Wh1FracErrEGScale) + ($Wh1FracErrJES * $Wh1FracErrJES) + ($Wh1FracErrJER * $Wh1FracErrJER) + ($Wh1FracErrMuEnergyScale * $Wh1FracErrMuEnergyScale) + ($Wh1FracErrTauEnergyScale * $Wh1FracErrTauEnergyScale) + ($Wh1FracErrUnclusteredEnergyScale * $Wh1FracErrUnclusteredEnergyScale))" | bc -q 2>/dev/null)
    echo "+/-$Wh1METFracErr% (total MET scale)"

    ### ggH ###

    #get the nominal expected signal
    nggNom=`echo $outputNom | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    echo "No. gg (nominal): $nggNom"

    #initialize fractional errors to 0--only set nonzero if nggNom != 0
    ggFracErrLowBVetoScale=0
    ggFracErrHighBVetoScale=0

    #calculate the uncertainty due to b veto data/MC scale
    nggMinus1SigmaBVetoScale=`echo $outputMinus1SigmaBVetoScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    nggPlus1SigmaBVetoScale=`echo $outputPlus1SigmaBVetoScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    ggErrLowBVetoScale=$(echo "scale=4; $nggNom - $nggMinus1SigmaBVetoScale" | bc -q 2>/dev/null)
    ggErrHighBVetoScale=$(echo "scale=4; $nggNom - $nggPlus1SigmaBVetoScale" | bc -q 2>/dev/null)
    if [ "$nggNom" != "0" ]
	then
	ggFracErrLowBVetoScale=$(echo "scale=4; ($ggErrLowBVetoScale / $nggNom) * 100" | bc -q 2>/dev/null)
	ggFracErrHighBVetoScale=$(echo "scale=4; ($ggErrHighBVetoScale / $nggNom) * 100" | bc -q 2>/dev/null)
    fi
    echo "$ggErrLowBVetoScale($ggFracErrLowBVetoScale%)/$ggErrHighBVetoScale($ggFracErrHighBVetoScale%) (b veto data/MC scale)"

    #get the central expected signal for MET uncertainties
    nggCen=`echo $outputCen | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    echo "No. gg (central value for MET uncertainties): $nggCen"

    #initialize fractional errors to 0--only set nonzero if nggCen != 0
    ggFracErrLowEGScale=0
    ggFracErrHighEGScale=0
    ggFracErrLowJES=0
    ggFracErrHighJES=0
    ggFracErrLowJER=0
    ggFracErrHighJER=0
    ggFracErrLowMuEnergyScale=0
    ggFracErrHighMuEnergyScale=0
    ggFracErrLowTauEnergyScale=0
    ggFracErrHighTauEnergyScale=0
    ggFracErrLowUnclusteredEnergyScale=0
    ggFracErrHighUnclusteredEnergyScale=0

    #calculate the uncertainty due to e/g data/MC scale
    nggMinus1SigmaEGScale=`echo $outputMinus1SigmaEGScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    nggPlus1SigmaEGScale=`echo $outputPlus1SigmaEGScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    ggErrLowEGScale=$(echo "scale=4; $nggCen - $nggMinus1SigmaEGScale" | bc -q 2>/dev/null)
    ggErrHighEGScale=$(echo "scale=4; $nggCen - $nggPlus1SigmaEGScale" | bc -q 2>/dev/null)
    if [ "$nggCen" != "0" ]
	then
	ggFracErrLowEGScale=$(echo "scale=4; ($ggErrLowEGScale / $nggCen) * 100" | bc -q 2>/dev/null)
	ggFracErrHighEGScale=$(echo "scale=4; ($ggErrHighEGScale / $nggCen) * 100" | bc -q 2>/dev/null)
    fi
    echo "$ggErrLowEGScale($ggFracErrLowEGScale%)/$ggErrHighEGScale($ggFracErrHighEGScale%) (e/gamma data/MC scale)"
    echo "-1sigma: $nggMinus1SigmaEGScale, +1sigma: $nggPlus1SigmaEGScale"
	
    #get the +/-1sigma jet energy data/MC scale expected signal
    outputMinus1SigmaJES=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaJES", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    outputPlus1SigmaJES=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaJES", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

    #calculate the uncertainty due to jet energy data/MC scale
    nggMinus1SigmaJES=`echo $outputMinus1SigmaJES | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    nggPlus1SigmaJES=`echo $outputPlus1SigmaJES | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    ggErrLowJES=$(echo "scale=4; $nggCen - $nggMinus1SigmaJES" | bc -q 2>/dev/null)
    ggErrHighJES=$(echo "scale=4; $nggCen - $nggPlus1SigmaJES" | bc -q 2>/dev/null)
    if [ "$nggCen" != "0" ]
	then
	ggFracErrLowJES=$(echo "scale=4; ($ggErrLowJES / $nggCen) * 100" | bc -q 2>/dev/null)
	ggFracErrHighJES=$(echo "scale=4; ($ggErrHighJES / $nggCen) * 100" | bc -q 2>/dev/null)
    fi
    echo "$ggErrLowJES($ggFracErrLowJES%)/$ggErrHighJES($ggFracErrHighJES%) (jet energy data/MC scale)"
    echo "-1sigma: $nggMinus1SigmaJES, +1sigma: $nggPlus1SigmaJES"

#v6 WH and v7 ggH skims don't apply JER smearing
#    #get the +/-1sigma MC jet energy resolution smeared expected signal
#    outputMinus1SigmaJER=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaJER", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
#    outputPlus1SigmaJER=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaJER", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

#    #calculate the uncertainty due to MC jet energy resolution smearing
#    nggMinus1SigmaJER=`echo $outputMinus1SigmaJER | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
#    nggPlus1SigmaJER=`echo $outputPlus1SigmaJER | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
#    ggErrLowJER=$(echo "scale=4; $nggCen - $nggMinus1SigmaJER" | bc -q 2>/dev/null)
#    ggErrHighJER=$(echo "scale=4; $nggCen - $nggPlus1SigmaJER" | bc -q 2>/dev/null)
#    if [ "$nggCen" != "0" ]
#	then
#	ggFracErrLowJER=$(echo "scale=4; ($ggErrLowJER / $nggCen) * 100" | bc -q 2>/dev/null)
#	ggFracErrHighJER=$(echo "scale=4; ($ggErrHighJER / $nggCen) * 100" | bc -q 2>/dev/null)
#    fi
#    echo "$ggErrLowJER($ggFracErrLowJER%)/$ggErrHighJER($ggFracErrHighJER%) (jet energy MC resolution smeared)"
#    echo "-1sigma: $nggMinus1SigmaJER, +1sigma: $nggPlus1SigmaJER"

    #get the +/-1sigma muon energy data/MC scale expected signal
    outputMinus1SigmaMuEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaMuEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    outputPlus1SigmaMuEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaMuEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

    #calculate the uncertainty due to muon energy data/MC scale
    nggMinus1SigmaMuEnergyScale=`echo $outputMinus1SigmaMuEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    nggPlus1SigmaMuEnergyScale=`echo $outputPlus1SigmaMuEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    ggErrLowMuEnergyScale=$(echo "scale=4; $nggCen - $nggMinus1SigmaMuEnergyScale" | bc -q 2>/dev/null)
    ggErrHighMuEnergyScale=$(echo "scale=4; $nggCen - $nggPlus1SigmaMuEnergyScale" | bc -q 2>/dev/null)
    if [ "$nggCen" != "0" ]
	then
	ggFracErrLowMuEnergyScale=$(echo "scale=4; ($ggErrLowMuEnergyScale / $nggCen) * 100" | bc -q 2>/dev/null)
	ggFracErrHighMuEnergyScale=$(echo "scale=4; ($ggErrHighMuEnergyScale / $nggCen) * 100" | bc -q 2>/dev/null)
    fi
    echo "$ggErrLowMuEnergyScale($ggFracErrLowMuEnergyScale%)/$ggErrHighMuEnergyScale($ggFracErrHighMuEnergyScale%) (muon energy data/MC scale)"
    echo "-1sigma: $nggMinus1SigmaMuEnergyScale, +1sigma: $nggPlus1SigmaMuEnergyScale"

    #get the +/-1sigma tau energy data/MC scale expected signal
    outputMinus1SigmaTauEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaTauEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    outputPlus1SigmaTauEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaTauEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

    #calculate the uncertainty due to tau energy data/MC scale
    nggMinus1SigmaTauEnergyScale=`echo $outputMinus1SigmaTauEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    nggPlus1SigmaTauEnergyScale=`echo $outputPlus1SigmaTauEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    ggErrLowTauEnergyScale=$(echo "scale=4; $nggCen - $nggMinus1SigmaTauEnergyScale" | bc -q 2>/dev/null)
    ggErrHighTauEnergyScale=$(echo "scale=4; $nggCen - $nggPlus1SigmaTauEnergyScale" | bc -q 2>/dev/null)
    if [ "$nggCen" != "0" ]
	then
	ggFracErrLowTauEnergyScale=$(echo "scale=4; ($ggErrLowTauEnergyScale / $nggCen) * 100" | bc -q 2>/dev/null)
	ggFracErrHighTauEnergyScale=$(echo "scale=4; ($ggErrHighTauEnergyScale / $nggCen) * 100" | bc -q 2>/dev/null)
    fi
    echo "$ggErrLowTauEnergyScale($ggFracErrLowTauEnergyScale%)/$ggErrHighTauEnergyScale($ggFracErrHighTauEnergyScale%) (tau energy data/MC scale)"
    echo "-1sigma: $nggMinus1SigmaTauEnergyScale, +1sigma: $nggPlus1SigmaTauEnergyScale"

    #get the +/-1sigma unclustered energy data/MC scale expected signal
    outputMinus1SigmaUnclusteredEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_minus1SigmaUnclusteredEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`
    outputPlus1SigmaUnclusteredEnergyScale=`root -l -b -q 'formatSigPlots.C("'${version}'", "'${version}'", "'${narrowVersion}'", "_plus1SigmaUnclusteredEnergyScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'", 3)'`

    #calculate the uncertainty due to unclustered energy data/MC scale
    nggMinus1SigmaUnclusteredEnergyScale=`echo $outputMinus1SigmaUnclusteredEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    nggPlus1SigmaUnclusteredEnergyScale=`echo $outputPlus1SigmaUnclusteredEnergyScale | sed -e "s%.*Region A signal 1, m > 4: \([0-9\.]*\).*%\1%"`
    ggErrLowUnclusteredEnergyScale=$(echo "scale=4; $nggCen - $nggMinus1SigmaUnclusteredEnergyScale" | bc -q 2>/dev/null)
    ggErrHighUnclusteredEnergyScale=$(echo "scale=4; $nggCen - $nggPlus1SigmaUnclusteredEnergyScale" | bc -q 2>/dev/null)
    if [ "$nggCen" != "0" ]
	then
	ggFracErrLowUnclusteredEnergyScale=$(echo "scale=4; ($ggErrLowUnclusteredEnergyScale / $nggCen) * 100" | bc -q 2>/dev/null)
	ggFracErrHighUnclusteredEnergyScale=$(echo "scale=4; ($ggErrHighUnclusteredEnergyScale / $nggCen) * 100" | bc -q 2>/dev/null)
    fi
    echo "$ggErrLowUnclusteredEnergyScale($ggFracErrLowUnclusteredEnergyScale%)/$ggErrHighUnclusteredEnergyScale($ggFracErrHighUnclusteredEnergyScale%) (unclustered energy data/MC scale)"
    echo "-1sigma: $nggMinus1SigmaUnclusteredEnergyScale, +1sigma: $nggPlus1SigmaUnclusteredEnergyScale"

    #get the larger of the +/-1sigma e/g scale shifts
    ggFracErrEGScale=$ggFracErrLowEGScale
    isLargerggEGScale=`echo ${ggFracErrHighEGScale#-}'>'${ggFracErrLowEGScale#-} | bc -l`
    if [ "$isLargerggEGScale" = "1" ]
    then
	ggFracErrEGScale=$ggFracErrHighEGScale
    fi

    #get the larger of the +/-1sigma JES shifts
    ggFracErrJES=$ggFracErrLowJES
    isLargerggJES=`echo ${ggFracErrHighJES#-}'>'${ggFracErrLowJES#-} | bc -l`
    if [ "$isLargerggJES" = "1" ]
    then
	ggFracErrJES=$ggFracErrHighJES
    fi

    #get the larger of the +/-1sigma JER shifts
    ggFracErrJER=$ggFracErrLowJER
    isLargerggJER=`echo ${ggFracErrHighJER#-}'>'${ggFracErrLowJER#-} | bc -l`
    if [ "$isLargerggJER" = "1" ]
    then
	ggFracErrJER=$ggFracErrHighJER
    fi

    #get the larger of the +/-1sigma muon energy scale shifts
    ggFracErrMuEnergyScale=$ggFracErrLowMuEnergyScale
    isLargerggMuEnergyScale=`echo ${ggFracErrHighMuEnergyScale#-}'>'${ggFracErrLowMuEnergyScale#-} | bc -l`
    if [ "$isLargerggMuEnergyScale" = "1" ]
    then
	ggFracErrMuEnergyScale=$ggFracErrHighMuEnergyScale
    fi

    #get the larger of the +/-1sigma tau energy scale shifts
    ggFracErrTauEnergyScale=$ggFracErrLowTauEnergyScale
    isLargerggTauEnergyScale=`echo ${ggFracErrHighTauEnergyScale#-}'>'${ggFracErrLowTauEnergyScale#-} | bc -l`
    if [ "$isLargerggTauEnergyScale" = "1" ]
    then
	ggFracErrTauEnergyScale=$ggFracErrHighTauEnergyScale
    fi

    #get the larger of the +/-1sigma unclustered energy scale shifts
    ggFracErrUnclusteredEnergyScale=$ggFracErrLowUnclusteredEnergyScale
    isLargerggUnclusteredEnergyScale=`echo ${ggFracErrHighUnclusteredEnergyScale#-}'>'${ggFracErrLowUnclusteredEnergyScale#-} | bc -l`
    if [ "$isLargerggUnclusteredEnergyScale" = "1" ]
    then
	ggFracErrUnclusteredEnergyScale=$ggFracErrHighUnclusteredEnergyScale
    fi

    #add energy scale and resolution uncertainties in quadrature
    #use the larger of the +/-1sigma shifts as the symmetric error
    ggMETFracErr=$(echo "scale=4; sqrt(($ggFracErrEGScale * $ggFracErrEGScale) + ($ggFracErrJES * $ggFracErrJES) + ($ggFracErrJER * $ggFracErrJER) + ($ggFracErrMuEnergyScale * $ggFracErrMuEnergyScale) + ($ggFracErrTauEnergyScale * $ggFracErrTauEnergyScale) + ($ggFracErrUnclusteredEnergyScale * $ggFracErrUnclusteredEnergyScale))" | bc -q 2>/dev/null)
    echo "+/-$ggMETFracErr% (total MET scale)"
  done
  echo ""
done

exit 0