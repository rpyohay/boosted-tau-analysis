#!/bin/bash

if [ $# -ne 2 ]
    then
    echo "Usage: ./generateWh1TauAnalyzerCfg.sh <version> <template cfg>"
    exit 0
fi

####STUFF TO CONFIGURE####

#version
version=$1
templateCfg=$2
infoTag=""
dir=$version

#number of samples
nSamples=1
iBeg=0
iEnd=`expr $nSamples - 1`

#input file prefix
inputFilePrefix="file:/data1/yohay/Wh1_Medium/EDM_files/"

#CleanJets output file prefix
cleanJetsOutputFilePrefix="`pwd`/${dir}/"
#cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
tauAnalyzerOutputFilePrefix="/data1/`whoami`/Wh1_Medium/"
#tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/Wh1_Medium/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_new.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Wh1.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_Wh1_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_Wh1_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_Wh1_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}Wh1${infoTag}_${version}.root" )

#samples
samples=( "Wh1" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_cfg.py

  #generate cfg file for the sample with no isolation cut
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.tauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_all_cfg.py
done

exit 0
