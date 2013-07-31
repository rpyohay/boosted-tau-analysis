#!/bin/bash

#input files
inputFileHeader="file:/data1/yohay/"
inputFileTrailer="_24Sep12.root"
inputFiles=( "'${inputFileHeader}NMSSMHiggs_WH_files1-250${inputFileTrailer}',\n    '${inputFileHeader}NMSSMHiggs_WH_files251-500${inputFileTrailer}',\n    '${inputFileHeader}NMSSMHiggs_WH_files501-750${inputFileTrailer}',\n    '${inputFileHeader}NMSSMHiggs_WH_files751-1000${inputFileTrailer}'" "'${inputFileHeader}Summer12DrellYan_PFTauReReco${inputFileTrailer}'" )

#mom PDG IDs (MOM_PDGID)
momPDGIDs=( A_PDGID Z_PDGID )

#sister tau decay types (SISTER_TAU_DECAY_TYPE)
sisterTauDecayTypes=( TAU_HAD TAU_MU )

#count sister (COUNT_SISTER)
countSister=( False True )

#minimum number of gen objects required to proceed (MIN_NUM_GEN_OBJECTS_TO_PASS_FILTER)
minNumGenObjectsToPassFilter=( 1 2 )

#analyzer (ANALYZER)
analyzer=( TauEfficiencyAnalyzer TauTagAndProbeEfficiencyAnalyzer )

#sample (SAMPLE)
sample=( NMSSMHiggs_WH_muHad ZTauTau_hadHad )

#tag1 (TAG1)
tag1=( denominatorTag tagFailTag )

#tag2 (TAG2)
tag2=( numeratorTag tagTagTag )

#HLT selection (HLT_SELECTION)
#HLTSelection=( "process.genWMuNuSelector*process.IsoMu24eta2p1Selector*" "" )
HLTSelection=( "process.genWMuNuSelector*" "" )

#generate a cfg file for each input file
for i in `seq 0 1`
  do
  sed -e "s%INPUT_FILES%${inputFiles[${i}]}%" -e "s%MOM_PDGID%${momPDGIDs[${i}]}%" -e "s%SISTER_TAU_DECAY_TYPE%${sisterTauDecayTypes[${i}]}%" -e "s%COUNT_SISTER%${countSister[${i}]}%" -e "s%MIN_NUM_GEN_OBJECTS_TO_PASS_FILTER%${minNumGenObjectsToPassFilter[${i}]}%" -e "s%ANALYZER%${analyzer[${i}]}%" -e "s%SAMPLE%${sample[${i}]}%" -e "s%TAG1%${tag1[${i}]}%" -e "s%TAG2%${tag2[${i}]}%" -e "s%HLT_SELECTION%${HLTSelection[${i}]}%" calcHPSTauIsoEff.py > calcHPSTauIsoEff_${sample[${i}]}.py
done

exit 0

