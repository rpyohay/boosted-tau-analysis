#!/bin/bash

#input files
inputFileHeader="file:/data1/yohay/"
inputFileTrailer=".root"
inputFiles=( "'${inputFileHeader}NMSSMHiggs_gg_files1-250_24Sep12${inputFileTrailer}',\n    '${inputFileHeader}NMSSMHiggs_gg_files251-500_24Sep12${inputFileTrailer}',\n    '${inputFileHeader}NMSSMHiggs_gg_files501-750_24Sep12${inputFileTrailer}',\n    '${inputFileHeader}NMSSMHiggs_gg_files751-1000_24Sep12${inputFileTrailer}'" "'${inputFileHeader}Summer12DrellYan_PFTauReReco${inputFileTrailer}'" )

#sequences
sequences=( signalSequence ZTauTauSequence )
medIsoSequences=( medIsoSignalSequence medIsoZTauTauSequence )

#generate a cfg file for each input file
iLastFile=`expr ${#inputFiles[@]} - 1`
for iFile in `seq 0 $iLastFile`
  do
  #sed -e "s%INPUT_FILES%${inputFiles[${iFile}]}%" -e "s%SEQUENCE%${sequences[${iFile}]}%" tauIsoEff.py > tauIsoEff_${sequences[${iFile}]}.py
  sed -e "s%INPUT_FILES%${inputFiles[${iFile}]}%" -e "s%SEQUENCE%${medIsoSequences[${iFile}]}%" tauIsoEff.py > tauIsoEff_${medIsoSequences[${iFile}]}.py
done

exit 0

