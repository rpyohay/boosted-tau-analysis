#!/bin/bash

if [ $# -ne 2 ]
    then
    echo "Usage: ./generateWh1METUncertaintyTauAnalyzerCfg.sh <version> <template cfg>"
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
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a5_v6.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a7_v6.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a9_v6.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a11_v6.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a13_v6.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_a15_v6.root'\n    ])" )

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

  #generate cfg file for the isolated sample, energy and b veto efficiency data/MC scale 
  #uncertainties shifted up(down) by +(-)1 sigma
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%PREFIX%${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis%g" -e "s%VERSION%${version}%g" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%HIGGSREW%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%${samples[${i}]}%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_signal_MET_uncertainties_cfg.py

done

#generate run cfg that runs all signal MET uncertainty files in the directory
cat <<EOF > runWh1IsoSignalMETUncertaintyTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *Wh1_a*_iso_signal_MET_uncertainties*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file >& \$outFile
done

exit 0
EOF
chmod a+x runWh1IsoSignalMETUncertaintyTauAnalyzerCfgs.sh

exit 0
