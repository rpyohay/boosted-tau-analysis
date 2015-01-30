#!/bin/bash

if [ $# -gt 3 ]
    then
    echo "Usage: ./generateQCDBTauAnalyzerCfgs.sh <version> <template cfg> [reweightOnly]"
    exit 0
fi

####STUFF TO CONFIGURE####

#version
version=$1
templateCfg=$2
infoTag=""
reweightOnly=0
if [ "$3" == "reweightOnly" ]
    then
    reweightOnly=1
fi
dir=$version

#number of samples
nSamples=4
iBeg=0
iEnd=`expr $nSamples - 1`

#input file prefix
inputFilePrefix="root://eoscms//eos/cms/store/user/friccita/"

#input file suffix
inputFileSuffix="_bEnriched_TuneZ2star_8Tev-pythia6-evtgen/Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v1/"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/QCD/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/QCDB/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_10_1_sQy.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_11_1_Ir8.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_12_1_0ep.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_13_1_E4l.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_14_1_cPa.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_15_1_w76.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_16_1_8mD.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_17_1_EGj.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_18_1_knb.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_19_1_ijn.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_1_1_4Q6.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_20_1_sH0.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_2_1_s9k.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_3_1_pwt.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_4_1_J8V.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_5_1_DNx.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_6_1_tAm.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_7_1_IWV.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_8_1_UOL.root',\n    '${inputFilePrefix}QCD_Pt-15To30${inputFileSuffix}data_no_selection_9_1_Rsi.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_10_1_4Ts.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_11_1_ZKz.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_12_1_vCC.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_13_1_hu9.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_14_1_JXs.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_15_1_uB7.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_16_1_Xwe.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_17_1_Mj7.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_18_1_9fV.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_19_1_S2Q.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_1_1_3e0.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_20_1_qwb.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_21_1_4St.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_22_1_b1H.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_2_1_UDV.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_3_1_okS.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_4_1_VEV.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_5_1_ySF.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_6_1_p9L.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_7_1_rr8.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_8_1_XQ4.root',\n    '${inputFilePrefix}QCD_Pt-30To50${inputFileSuffix}data_no_selection_9_1_32Z.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_10_1_rRJ.root',\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_11_1_nhw.root',\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_1_1_3BQ.root',\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_2_1_JJs.root',\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_3_1_2qi.root',\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_4_1_Cpk.root',\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_5_1_00M.root',\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_6_1_4bE.root',\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_7_1_kZK.root',\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_8_1_Nuu.root',\n    '${inputFilePrefix}QCD_Pt-50To150${inputFileSuffix}data_no_selection_9_1_hrj.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefix}QCD_Pt-150${inputFileSuffix}data_no_selection_1_1_XN5.root',\n    '${inputFilePrefix}QCD_Pt-150${inputFileSuffix}data_no_selection_2_1_TPS.root',\n    '${inputFilePrefix}QCD_Pt-150${inputFileSuffix}data_no_selection_3_1_a5t.root',\n    ])" )
#"readFiles.extend([\n    ])"
#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_QCDB_Pt-15To30.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_QCDB_Pt-30To50.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_QCDB_Pt-50To150.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_QCDB_Pt-150.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_QCDB_Pt-15To30_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_QCDB_Pt-30To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_QCDB_Pt-50To150_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_QCDB_Pt-150_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_QCDB_Pt-15To30_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_QCDB_Pt-30To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_QCDB_Pt-50To150_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_QCDB_Pt-150_${version}.root" )
nonIsoReweightTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_QCDB_Pt-15To30_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_QCDB_Pt-30To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_QCDB_Pt-50To150_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_QCDB_Pt-150_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_QCDB_Pt-15To30_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_QCDB_Pt-30To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_QCDB_Pt-50To150_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_QCDB_Pt-150_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}QCDB_Pt-15To30${infoTag}_${version}.root" "${EDMOutputFilePrefix}QCDB_Pt-30To50${infoTag}_${version}.root" "${EDMOutputFilePrefix}QCDB_Pt-50To150${infoTag}_${version}.root" "${EDMOutputFilePrefix}QCDB_Pt-150${infoTag}_${version}.root" )

#samples
samples=( "QCDB_Pt-15To30" "QCDB_Pt-30To50" "QCDB_Pt-50To150" "QCDB_Pt-150" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%HIGGSREW%False" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_cfg.py

  #generate cfg file for the non-isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%HIGGSREW%False" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_nonIso_cfg.py

  #generate cfg file for the non-isolated, reweighted sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoReweightTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%HIGGSREW%False" -e "s%REWEIGHT%False" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_nonIsoReweight_cfg.py

  #generate cfg file for the sample with no isolation cut
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.tauAnalysisSequence%" -e "s%HIGGSREW%False" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_all_cfg.py

  #generate iso+nonIso+reweight job submission script for LSF
  cat <<EOF > tauanalyzer_${samples[${i}]}_cfg.sh
#!/bin/bash

jobDir="`pwd`"
fileNamePrefix="tauanalyzer_${samples[${i}]}"

cd \$jobDir
eval \`scramv1 runtime -sh\`
cd -
cp \$jobDir/\${fileNamePrefix}_iso_cfg.py \$jobDir/\${fileNamePrefix}_nonIso_cfg.py \$jobDir/\${fileNamePrefix}_nonIsoReweight_cfg.py .
cmsRun \${fileNamePrefix}_iso_cfg.py
if [ $reweightOnly -eq 0 ]
    then
    cmsRun \${fileNamePrefix}_nonIso_cfg.py
    cmsStage -f ${nonIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
    rm ${nonIsoTauAnalyzerOutputFiles[${i}]}
fi
cmsRun \${fileNamePrefix}_nonIsoReweight_cfg.py
cmsStage -f ${isoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${nonIsoReweightTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
rm ${isoTauAnalyzerOutputFiles[${i}]} ${nonIsoReweightTauAnalyzerOutputFiles[${i}]}

exit 0
EOF
  chmod a+x tauanalyzer_${samples[${i}]}_cfg.sh

  #generate noIsoCut job submission script for LSF
  cat <<EOF > tauanalyzer_${samples[${i}]}_all_cfg.sh
#!/bin/bash

jobDir="`pwd`"
fileNamePrefix="tauanalyzer_${samples[${i}]}"

cd \$jobDir
eval \`scramv1 runtime -sh\`
cd -
cp \$jobDir/\${fileNamePrefix}_all_cfg.py .
cmsRun \${fileNamePrefix}_all_cfg.py
cmsStage -f ${allTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
rm ${allTauAnalyzerOutputFiles[${i}]}

exit 0
EOF
  chmod a+x tauanalyzer_${samples[${i}]}_all_cfg.sh
done

#generate run cfg that runs all files in the directory
cat <<EOF > runQCDBTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *QCDB_*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runQCDBTauAnalyzerCfgs.sh

#generate script that submits all iso+nonIso+reweight jobs to LSF
cat <<EOF > submitQCDBTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*QCDB_*.sh | grep -v all | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitQCDBTauAnalyzerJobs.sh

#generate script that submits all noIsoCut jobs to LSF
cat <<EOF > submitQCDBAllTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*QCDB_*all*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitQCDBAllTauAnalyzerJobs.sh

#generate script that copies all iso+nonIso+reweight files locally from EOS
cat <<EOF > copyQCDBFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in "15To30" "30To50" "50To150" "150"
  do
  for cut in Iso NonIso NonIsoReweight
    do
    if [ "\$cut" != "NonIso" ] || [ $reweightOnly -eq 0 ]
        then
        cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_QCDB_Pt-\${sample}_${version}.root /data1/`whoami`/QCDB/analysis/
        cmsRm /store/user/`whoami`/muHad\${cut}Analysis_QCDB_Pt-\${sample}_${version}.root
    fi
  done
done

exit 0
EOF
chmod a+x copyQCDBFromEOS.sh

#generate script that copies all noIsoCut files locally from EOS
cat <<EOF > copyAllQCDBFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in "15To30" "30To50" "50To150" "150"
  do
  cmsStage -f /store/user/`whoami`/muHadAnalysis_QCDB_Pt-\${sample}_${version}.root /data1/`whoami`/QCDB/analysis/
  cmsRm /store/user/`whoami`/muHadAnalysis_QCDB_Pt-\${sample}_${version}.root
done

exit 0
EOF
chmod a+x copyAllQCDBFromEOS.sh

exit 0
