#!/bin/bash

if [ $# -ne 2 ]
    then
    echo "Usage: ./generateTTJetsDrellYanAnalyzerCfgs.sh <version> <template cfg>"
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
inputFilePrefix="root://eoscms//eos/cms/store/user/yohay/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-Summer12_DR53X-PU_S10_START53_V7A-v2-AODSIM_skim_DYEnriched/"

#DrellYanAnalyzer output file prefix
#DrellYanAnalyzerOutputFilePrefix="/data1/yohay/TTJets/analysis/"
DrellYanAnalyzerOutputFilePrefix=""

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_10_1_1i6.root',\n    '${inputFilePrefix}data_no_selection_11_1_q81.root',\n    '${inputFilePrefix}data_no_selection_12_1_aSy.root',\n    '${inputFilePrefix}data_no_selection_13_1_XxE.root',\n    '${inputFilePrefix}data_no_selection_14_1_ekn.root',\n    '${inputFilePrefix}data_no_selection_1_1_4dd.root',\n    '${inputFilePrefix}data_no_selection_2_1_iwm.root',\n    '${inputFilePrefix}data_no_selection_3_1_4Ho.root',\n    '${inputFilePrefix}data_no_selection_4_1_UvP.root',\n    '${inputFilePrefix}data_no_selection_5_1_n5p.root',\n    '${inputFilePrefix}data_no_selection_6_1_6Qj.root',\n    '${inputFilePrefix}data_no_selection_7_1_Aa5.root',\n    '${inputFilePrefix}data_no_selection_8_1_tBu.root',\n    '${inputFilePrefix}data_no_selection_9_1_fyZ.root'\n    ])" )

#DrellYanAnalyzer output files
DrellYanAnalyzerOutputFiles=( "${DrellYanAnalyzerOutputFilePrefix}DrellYanAnalysis_TTJets_${version}.root" )

#samples
samples=( "TTJets" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%OUTFILE%${DrellYanAnalyzerOutputFiles[${i}]}%" -e "s%PUSCENARIO%S10%" -e "s%MCFLAG%True%" ../${templateCfg} > drellyananalyzer_${samples[${i}]}_cfg.py

  #generate job submission script for LSF
  cat <<EOF > drellyananalyzer_${samples[${i}]}_cfg.sh
#!/bin/bash

jobDir="`pwd`"
fileNamePrefix="drellyananalyzer_${samples[${i}]}"

cd \$jobDir
eval \`scramv1 runtime -sh\`
cd -
cp \$jobDir/\${fileNamePrefix}_cfg.py .
cmsRun \${fileNamePrefix}_cfg.py
cmsStage -f ${DrellYanAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
rm ${DrellYanAnalyzerOutputFiles[${i}]}

exit 0
EOF
  chmod a+x drellyananalyzer_${samples[${i}]}_cfg.sh
done

#generate run cfg that runs all files in the directory
cat <<EOF > runTTJetsDrellYanAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *TTJets*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runTTJetsDrellYanAnalyzerCfgs.sh

#generate script that submits all jobs to LSF
cat <<EOF > submitTTJetsDrellYanAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh drellyananalyzer*TTJets*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitTTJetsDrellYanAnalyzerJobs.sh

#generate script that copies all files locally from EOS
cat <<EOF > copyDrellYanEnrichedTTJetsFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  cmsStage -f /store/user/`whoami`/DrellYanAnalysis_TTJets_${version}.root /data1/`whoami`/TTJets/analysis/
  cmsRm /store/user/`whoami`/DrellYanAnalysis_TTJets_${version}.root
done

exit 0
EOF
chmod a+x copyDrellYanEnrichedTTJetsFromEOS.sh

exit 0
