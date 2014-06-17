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

  #loop over MT bin
  for iMTBin in `seq $iBeg $iEndMTBin`
  do

    #get the nominal expected signal
    outputNom=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    nWh1Nom=`echo $outputNom | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`

    #get the expected signal - 1sigma b veto scale
    outputMinus1Sigma=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_minus1SigmaBVetoScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    nWh1Minus1Sigma=`echo $outputMinus1Sigma | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`

    #get the expected signal + 1sigma b veto scale
    outputPlus1Sigma=`root -l -b -q 'formatPlots.C("'${version}'", "'${version}'", false, "_plus1SigmaBVetoScale", "'${a1Mass[${iMass}]}'", "'${MTBin[${iMTBin}]}'")'`
    nWh1Plus1Sigma=`echo $outputPlus1Sigma | sed -e "s%.*Region A signal 0, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`

    #calculate the absolute difference between the nominal and - 1sigma b veto scale
    Wh1ErrLow=$(echo "scale=3; $nWh1Minus1Sigma - $nWh1Nom" | bc -q 2>/dev/null)

    #calculate the absolute difference between the nominal and + 1sigma b veto scale
    Wh1ErrHigh=$(echo "scale=3; $nWh1Nom - $nWh1Plus1Sigma" | bc -q 2>/dev/null)

    #calculate the fractional difference between the nominal and - 1sigma b veto scale
    Wh1FracErrLow=$(echo "scale=3; $Wh1ErrLow / $nWh1Nom" | bc -q 2>/dev/null)

    #calculate the fractional difference between the nominal and + 1sigma b veto scale
    Wh1FracErrHigh=$(echo "scale=3; $Wh1ErrHigh / $nWh1Nom" | bc -q 2>/dev/null)

    #print the absolute and fractional errors
    echo "ma1 = ${a1Mass[${iMass}]}, ${MTBin[${iMTBin}]}"
    echo "No. Wh1: $nWh1Nom - $Wh1ErrLow($Wh1FracErrLow) + $Wh1ErrHigh($Wh1FracErrHigh)"

    #temporary until gg fusion skims are done
    if [ "${a1Mass[${iMass}]}" = "_a9" ]
	then
	nggNom=`echo $outputNom | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
	nggMinus1Sigma=`echo $outputMinus1Sigma | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
	nggPlus1Sigma=`echo $outputPlus1Sigma | sed -e "s%.*Region A signal 1, m > 4: \([0-9]*\.[0-9]*\).*%\1%"`
	ggErrLow=$(echo "scale=3; $nggMinus1Sigma - $nggNom" | bc -q 2>/dev/null)
	ggErrHigh=$(echo "scale=3; $nggNom - $nggPlus1Sigma" | bc -q 2>/dev/null)
	ggFracErrLow=$(echo "scale=3; $ggErrLow / $nggNom" | bc -q 2>/dev/null)
	ggFracErrHigh=$(echo "scale=3; $ggErrHigh / $nggNom" | bc -q 2>/dev/null)
	echo "ma1 = ${a1Mass[${iMass}]}, ${MTBin[${iMTBin}]}"
	echo "No. gg: $nggNom - $ggErrLow($ggFracErrLow) + $ggErrHigh($ggFracErrHigh)"
    fi
  done
done

exit 0