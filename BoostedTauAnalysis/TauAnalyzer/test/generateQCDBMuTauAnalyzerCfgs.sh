#!/bin/bash

if [ $# -gt 3 ]
    then
    echo "Usage: ./generateQCDBMuTauAnalyzerCfgs.sh <version> <template cfg> [reweightOnly]"
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
inputFileSuffix="_bEnriched_MuEnrichedPt14_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v1/"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/QCD/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/QCDBMu/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_10_1_77n.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_11_1_Vgr.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_12_1_6Cr.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_13_1_CVV.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_14_1_a2Y.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_15_1_NLR.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_16_1_0Hc.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_17_1_tvf.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_18_1_nt5.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_19_1_1Jl.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_1_1_VgX.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_20_1_Uh5.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_21_1_5KU.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_22_1_2Cy.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_23_1_5IS.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_24_1_syn.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_25_1_Z1W.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_2_1_Ju1.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_3_1_KUV.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_4_1_GyN.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_5_1_FPB.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_6_1_HmL.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_7_1_c6W.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_8_1_aG7.root',\n    '${inputFilePrefix}QCD_pt15to30${inputFileSuffix}data_no_selection_9_1_nbd.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_10_1_0PW.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_11_1_gw1.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_12_1_uby.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_13_1_UDP.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_14_1_hM0.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_15_1_T3L.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_16_1_jxY.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_17_1_ibR.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_18_1_832.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_19_1_EIR.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_1_1_yOz.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_20_1_FV7.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_21_1_PFJ.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_22_1_p7K.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_23_1_rag.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_2_1_POu.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_3_1_c3u.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_4_1_ws9.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_5_1_ErM.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_6_1_AnO.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_7_1_Hb0.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_8_1_gY5.root',\n    '${inputFilePrefix}QCD_pt30to50${inputFileSuffix}data_no_selection_9_1_hmz.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_10_1_yis.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_11_1_X7k.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_12_1_x8D.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_13_1_iC7.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_14_1_l9L.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_15_1_Moo.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_16_1_aCF.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_17_1_WdI.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_18_1_vwr.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_19_1_DIn.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_1_1_vX7.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_20_1_0re.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_2_1_2y9.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_3_1_hsi.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_4_1_j0l.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_5_1_ueb.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_6_1_czh.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_7_1_8wN.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_8_1_3EC.root',\n    '${inputFilePrefix}QCD_pt50to150${inputFileSuffix}data_no_selection_9_1_mQ1.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_10_1_ENt.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_11_1_X8c.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_12_1_MO2.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_13_1_LJf.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_14_1_mMa.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_15_1_Phi.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_16_1_y1L.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_17_1_l2x.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_18_1_7Lh.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_19_1_CMf.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_1_1_iCX.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_20_1_kDM.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_21_1_Mxa.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_2_1_7q8.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_3_1_NVL.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_4_1_uSU.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_5_1_qzJ.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_6_1_XQb.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_7_1_hYB.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_8_1_n3X.root',\n    '${inputFilePrefix}QCD_pt150${inputFileSuffix}data_no_selection_9_1_qak.root',\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_QCDBMu_pt15to30.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_QCDBMu_pt30to50.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_QCDBMu_pt50to150.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_QCDBMu_pt150.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_QCDBMu_pt15to30_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_QCDBMu_pt30to50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_QCDBMu_pt50to150_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_QCDBMu_pt150_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_QCDBMu_pt15to30_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_QCDBMu_pt30to50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_QCDBMu_pt50to150_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_QCDBMu_pt150_${version}.root" )
nonIsoReweightTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_QCDBMu_pt15to30_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_QCDBMu_pt30to50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_QCDBMu_pt50to150_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_QCDBMu_pt150_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_QCDBMu_pt15to30_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_QCDBMu_pt30to50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_QCDBMu_pt50to150_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_QCDBMu_pt150_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}QCDBMu_pt15to30${infoTag}_${version}.root" "${EDMOutputFilePrefix}QCDBMu_pt30to50${infoTag}_${version}.root" "${EDMOutputFilePrefix}QCDBMu_pt50to150${infoTag}_${version}.root" "${EDMOutputFilePrefix}QCDBMu_pt150${infoTag}_${version}.root" )

#samples
samples=( "QCDBMu_pt15to30" "QCDBMu_pt30to50" "QCDBMu_pt50to150" "QCDBMu_pt150" )

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
cat <<EOF > runQCDBMuTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *QCDBMu_*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runQCDBMuTauAnalyzerCfgs.sh

#generate script that submits all iso+nonIso+reweight jobs to LSF
cat <<EOF > submitQCDBMuTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*QCDBMu_*.sh | grep -v all | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitQCDBMuTauAnalyzerJobs.sh

#generate script that submits all noIsoCut jobs to LSF
cat <<EOF > submitQCDBMuAllTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*QCDBMu_*all*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitQCDBMuAllTauAnalyzerJobs.sh

#generate script that copies all iso+nonIso+reweight files locally from EOS
cat <<EOF > copyQCDBMuFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in "15to30" "30to50" "50to150" "150"
  do
  for cut in Iso NonIso NonIsoReweight
    do
    if [ "\$cut" != "NonIso" ] || [ $reweightOnly -eq 0 ]
        then
        cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_QCDBMu_pt\${sample}_${version}.root /data1/`whoami`/QCDBMu/analysis/
        cmsRm /store/user/`whoami`/muHad\${cut}Analysis_QCDBMu_pt\${sample}_${version}.root
    fi
  done
done

exit 0
EOF
chmod a+x copyQCDBMuFromEOS.sh

#generate script that copies all noIsoCut files locally from EOS
cat <<EOF > copyAllQCDBMuFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in "15to30" "30to50" "50to150" "150"
  do
  cmsStage -f /store/user/`whoami`/muHadAnalysis_QCDBMu_pt\${sample}_${version}.root /data1/`whoami`/QCDBMu/analysis/
  cmsRm /store/user/`whoami`/muHadAnalysis_QCDBMu_pt\${sample}_${version}.root
done

exit 0
EOF
chmod a+x copyAllQCDBMuFromEOS.sh

exit 0
