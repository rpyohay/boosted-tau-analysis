#!/bin/bash

if [ $# -ne 2 ]
    then
    echo "Usage: ./generateggMETUncertaintyTauAnalyzerCfg.sh <version> <template cfg>"
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
inputFilePrefix="file:/data1/yohay/gg/EDM_files/"

#CleanJets output file prefix
cleanJetsOutputFilePrefix="`pwd`/${dir}/"
#cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
tauAnalyzerOutputFilePrefix="/data1/`whoami`/gg/"
#tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/gg/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a5_v7.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a7_v7.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a9_v7.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a11_v7.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a13_v7.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a15_v7.root'\n    ])" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}gg_a5${infoTag}_${version}.root" "${EDMOutputFilePrefix}gg_a7${infoTag}_${version}.root" "${EDMOutputFilePrefix}gg_a9${infoTag}_${version}.root" "${EDMOutputFilePrefix}gg_a11${infoTag}_${version}.root" "${EDMOutputFilePrefix}gg_a13${infoTag}_${version}.root" "${EDMOutputFilePrefix}gg_a15${infoTag}_${version}.root" )

#samples
samples=( "gg_a5" "gg_a7" "gg_a9" "gg_a11" "gg_a13" "gg_a15" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample, energy scale uncertainties shifted up(down) by +(-)1 
  #sigma
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%PREFIX%${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis%g" -e "s%VERSION%${version}%g" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%HIGGSREW%True%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%${samples[${i}]}%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_signal_MET_uncertainties_cfg.py

done

#generate run cfg that runs all signal MET uncertainty files in the directory
cat <<EOF > runggIsoSignalMETUncertaintyTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *gg_a*_iso_signal_MET_uncertainties*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file >& \$outFile
done

exit 0
EOF
chmod a+x runggIsoSignalMETUncertaintyTauAnalyzerCfgs.sh

exit 0
