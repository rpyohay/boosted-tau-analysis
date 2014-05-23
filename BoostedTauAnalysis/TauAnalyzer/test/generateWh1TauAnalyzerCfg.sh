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
nSamples=6
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
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a5.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a7.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_new.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a11.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a13.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a15.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Wh1_a5.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Wh1_a7.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Wh1_a9.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Wh1_a11.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Wh1_a13.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Wh1_a15.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_Wh1_a5_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_Wh1_a7_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_Wh1_a9_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_Wh1_a11_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_Wh1_a13_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_Wh1_a15_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_Wh1_a5_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_Wh1_a7_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_Wh1_a9_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_Wh1_a11_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_Wh1_a13_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_Wh1_a15_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_Wh1_a5_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_Wh1_a7_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_Wh1_a9_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_Wh1_a11_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_Wh1_a13_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_Wh1_a15_${version}.root" )
isoTauBVetoScaleMinus1SigmaAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_minus1SigmaBVetoScale_Wh1_a5_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_minus1SigmaBVetoScale_Wh1_a7_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_minus1SigmaBVetoScale_Wh1_a9_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_minus1SigmaBVetoScale_Wh1_a11_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_minus1SigmaBVetoScale_Wh1_a13_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_minus1SigmaBVetoScale_Wh1_a15_${version}.root" )
isoTauBVetoScalePlus1SigmaAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_plus1SigmaBVetoScale_Wh1_a5_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_plus1SigmaBVetoScale_Wh1_a7_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_plus1SigmaBVetoScale_Wh1_a9_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_plus1SigmaBVetoScale_Wh1_a11_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_plus1SigmaBVetoScale_Wh1_a13_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_plus1SigmaBVetoScale_Wh1_a15_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}Wh1_a5${infoTag}_${version}.root" "${EDMOutputFilePrefix}Wh1_a7${infoTag}_${version}.root" "${EDMOutputFilePrefix}Wh1_a9${infoTag}_${version}.root" "${EDMOutputFilePrefix}Wh1_a11${infoTag}_${version}.root" "${EDMOutputFilePrefix}Wh1_a13${infoTag}_${version}.root" "${EDMOutputFilePrefix}Wh1_a15${infoTag}_${version}.root" )

#samples
samples=( "Wh1_a5" "Wh1_a7" "Wh1_a9" "Wh1_a11" "Wh1_a13" "Wh1_a15" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample, nominal b veto scale
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%True%" -e "s%BTAGSCALESHIFT%mean%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_nominal_b_veto_scale_cfg.py

  #generate cfg file for the isolated sample, nominal b veto scale - 1sigma
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauBVetoScaleMinus1SigmaAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%True%" -e "s%BTAGSCALESHIFT%min%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_minus1sigma_b_veto_scale_cfg.py

  #generate cfg file for the isolated sample, nominal b veto scale + 1sigma
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauBVetoScalePlus1SigmaAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%True%" -e "s%BTAGSCALESHIFT%max%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_plus1sigma_b_veto_scale_cfg.py

  #generate cfg file for the sample with no isolation cut, nominal b veto scale
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.tauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%True%" -e "s%BTAGSCALESHIFT%max%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_all_cfg.py
done

#generate run cfg that runs all nominal b veto scale files in the directory
cat <<EOF > runWh1IsoNominalBVetoScaleTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *Wh1_a*_iso_nominal*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file >& \$outFile &
done

exit 0
EOF
chmod a+x runWh1IsoNominalBVetoScaleTauAnalyzerCfgs.sh

#generate run cfg that runs all - 1sigma b veto scale files in the directory
cat <<EOF > runWh1IsoMinus1SigmaBVetoScaleTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *Wh1_a*_iso_minus1sigma*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file >& \$outFile &
done

exit 0
EOF
chmod a+x runWh1IsoMinus1SigmaBVetoScaleTauAnalyzerCfgs.sh

#generate run cfg that runs all + 1sigma b veto scale files in the directory
cat <<EOF > runWh1IsoPlus1SigmaBVetoScaleTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *Wh1_a*_iso_plus1sigma*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file >& \$outFile &
done

exit 0
EOF
chmod a+x runWh1IsoPlus1SigmaBVetoScaleTauAnalyzerCfgs.sh

exit 0
