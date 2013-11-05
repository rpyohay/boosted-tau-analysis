#!/bin/bash

#input files
#inputFiles="NMSSMHiggs_gg_skim_1000Files Summer12DrellYan"
inputFiles="Summer12DrellYan NMSSMHiggs_WH_skim_1000Files"

#sequences
tightNonIsoSequences=( tightNonIsoSignalSequence tightNonIsoZMuMuSequence )
tightDetectorIsoSequences=( tightDetectorIsoSignalSequence tightDetectorIsoZMuMuSequence )
tightPFIsoNoPUSubtractionSequences=( tightPFIsoNoPUSubtractionSignalSequence tightPFIsoNoPUSubtractionZMuMuSequence )
tightMuonIsolationSequences=( tightMuonIsolationSignalSequence tightMuonIsolationZMuMuSequence )

#generate a cfg file for each input file
i=0
for iFile in $inputFiles
  do
  #sed -e "s%FILE_NAME%$iFile%" -e "s%SEQUENCE%${tightNonIsoSequences[${i}]}%" muIsoEff.py > muIsoEff_${tightNonIsoSequences[${i}]}.py
  #sed -e "s%FILE_NAME%$iFile%" -e "s%SEQUENCE%${tightDetectorIsoSequences[${i}]}%" muIsoEff.py > muIsoEff_${tightDetectorIsoSequences[${i}]}.py
  #sed -e "s%FILE_NAME%$iFile%" -e "s%SEQUENCE%${tightPFIsoNoPUSubtractionSequences[${i}]}%" muIsoEff.py > muIsoEff_${tightPFIsoNoPUSubtractionSequences[${i}]}.py
  sed -e "s%FILE_NAME%$iFile%" -e "s%SEQUENCE%${tightMuonIsolationSequences[${i}]}%" muIsoEff.py > muIsoEff_${tightMuonIsolationSequences[${i}]}.py
  i=`expr $i + 1`
done

exit 0

