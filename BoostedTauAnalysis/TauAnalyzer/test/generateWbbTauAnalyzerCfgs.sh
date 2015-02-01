#!/bin/bash

if [ $# -gt 3 ]
    then
    echo "Usage: ./generateWbbTauAnalyzerCfgs.sh <version> <template cfg> [reweightOnly]"
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
nSamples=1
iBeg=0
iEnd=`expr $nSamples - 1`

#input file prefix
inputFilePrefix="root://eoscms//eos/cms/store/user/yohay/WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v2/"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/Wbb/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/Wbb/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_10_1_gbq.root',\n    '${inputFilePrefix}data_no_selection_11_1_Jo0.root',\n    '${inputFilePrefix}data_no_selection_12_1_T0m.root',\n    '${inputFilePrefix}data_no_selection_13_1_WxJ.root',\n    '${inputFilePrefix}data_no_selection_14_1_qVS.root',\n    '${inputFilePrefix}data_no_selection_15_1_lfD.root',\n    '${inputFilePrefix}data_no_selection_16_2_pWY.root',\n    '${inputFilePrefix}data_no_selection_17_1_pRf.root',\n    '${inputFilePrefix}data_no_selection_18_1_RdB.root',\n    '${inputFilePrefix}data_no_selection_19_1_dU6.root',\n    '${inputFilePrefix}data_no_selection_1_1_SHo.root',\n    '${inputFilePrefix}data_no_selection_20_1_0nl.root',\n    '${inputFilePrefix}data_no_selection_21_1_QBM.root',\n    '${inputFilePrefix}data_no_selection_22_1_gqs.root',\n    '${inputFilePrefix}data_no_selection_23_1_WqF.root',\n    '${inputFilePrefix}data_no_selection_24_1_zd0.root',\n    '${inputFilePrefix}data_no_selection_25_1_jMJ.root',\n    '${inputFilePrefix}data_no_selection_26_1_XwN.root',\n    '${inputFilePrefix}data_no_selection_27_1_VDJ.root',\n    '${inputFilePrefix}data_no_selection_28_1_h0b.root',\n    '${inputFilePrefix}data_no_selection_29_1_O0T.root',\n    '${inputFilePrefix}data_no_selection_2_1_lIb.root',\n    '${inputFilePrefix}data_no_selection_30_1_jjt.root',\n    '${inputFilePrefix}data_no_selection_31_1_DqD.root',\n    '${inputFilePrefix}data_no_selection_32_1_IF0.root',\n    '${inputFilePrefix}data_no_selection_33_1_TTg.root',\n    '${inputFilePrefix}data_no_selection_34_1_AVS.root',\n    '${inputFilePrefix}data_no_selection_35_1_Ezp.root',\n    '${inputFilePrefix}data_no_selection_36_1_6Qm.root',\n    '${inputFilePrefix}data_no_selection_37_1_JBL.root',\n    '${inputFilePrefix}data_no_selection_38_1_7Mf.root',\n    '${inputFilePrefix}data_no_selection_39_1_Eip.root',\n    '${inputFilePrefix}data_no_selection_3_1_maD.root',\n    '${inputFilePrefix}data_no_selection_40_1_afn.root',\n    '${inputFilePrefix}data_no_selection_41_1_xG3.root',\n    '${inputFilePrefix}data_no_selection_42_1_667.root',\n    '${inputFilePrefix}data_no_selection_43_1_ZHa.root',\n    '${inputFilePrefix}data_no_selection_44_1_NsM.root',\n    '${inputFilePrefix}data_no_selection_45_1_XXm.root',\n    '${inputFilePrefix}data_no_selection_46_1_ezl.root',\n    '${inputFilePrefix}data_no_selection_47_1_3pT.root',\n    '${inputFilePrefix}data_no_selection_48_1_Uqu.root',\n    '${inputFilePrefix}data_no_selection_49_1_oQp.root',\n    '${inputFilePrefix}data_no_selection_4_1_s2m.root',\n    '${inputFilePrefix}data_no_selection_50_1_PCY.root',\n    '${inputFilePrefix}data_no_selection_51_1_27Y.root',\n    '${inputFilePrefix}data_no_selection_52_1_p7F.root',\n    '${inputFilePrefix}data_no_selection_53_1_Ue6.root',\n    '${inputFilePrefix}data_no_selection_54_1_Mh5.root',\n    '${inputFilePrefix}data_no_selection_55_1_Q3C.root',\n    '${inputFilePrefix}data_no_selection_56_1_pP6.root',\n    '${inputFilePrefix}data_no_selection_57_1_Cji.root',\n    '${inputFilePrefix}data_no_selection_58_1_GyW.root',\n    '${inputFilePrefix}data_no_selection_59_1_TRC.root',\n    '${inputFilePrefix}data_no_selection_5_1_DLn.root',\n    '${inputFilePrefix}data_no_selection_60_1_yoc.root',\n    '${inputFilePrefix}data_no_selection_61_1_G6r.root',\n    '${inputFilePrefix}data_no_selection_62_1_1kK.root',\n    '${inputFilePrefix}data_no_selection_63_1_RKa.root',\n    '${inputFilePrefix}data_no_selection_64_1_gTe.root',\n    '${inputFilePrefix}data_no_selection_65_1_4VJ.root',\n    '${inputFilePrefix}data_no_selection_66_1_PY3.root',\n    '${inputFilePrefix}data_no_selection_67_1_TEG.root',\n    '${inputFilePrefix}data_no_selection_68_1_QCV.root',\n    '${inputFilePrefix}data_no_selection_69_1_y7r.root',\n    '${inputFilePrefix}data_no_selection_6_1_XU3.root',\n    '${inputFilePrefix}data_no_selection_70_1_m5I.root',\n    '${inputFilePrefix}data_no_selection_71_2_NDf.root',\n    '${inputFilePrefix}data_no_selection_72_1_vgd.root',\n    '${inputFilePrefix}data_no_selection_73_1_6Ig.root',\n    '${inputFilePrefix}data_no_selection_74_1_dTW.root',\n    '${inputFilePrefix}data_no_selection_75_1_VUX.root',\n    '${inputFilePrefix}data_no_selection_76_1_oEp.root',\n    '${inputFilePrefix}data_no_selection_77_1_Mur.root',\n    '${inputFilePrefix}data_no_selection_78_1_mQZ.root',\n    '${inputFilePrefix}data_no_selection_79_1_b4I.root',\n    '${inputFilePrefix}data_no_selection_7_1_zZX.root',\n    '${inputFilePrefix}data_no_selection_80_1_Izd.root',\n    '${inputFilePrefix}data_no_selection_81_1_chB.root',\n    '${inputFilePrefix}data_no_selection_82_1_S3O.root',\n    '${inputFilePrefix}data_no_selection_83_2_Ojy.root',\n    '${inputFilePrefix}data_no_selection_84_2_KQc.root',\n    '${inputFilePrefix}data_no_selection_8_1_rw2.root',\n    '${inputFilePrefix}data_no_selection_9_1_kgK.root',\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Wbb.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_Wbb_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_Wbb_${version}.root" )
nonIsoReweightTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_Wbb_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_Wbb_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}Wbb${infoTag}_${version}.root" )

#samples
samples=( "Wbb" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%HIGGSREW%False%" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_cfg.py

  #generate cfg file for the non-isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%HIGGSREW%False%" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_nonIso_cfg.py

  #generate cfg file for the non-isolated, reweighted sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoReweightTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%HIGGSREW%False%" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_nonIsoReweight_cfg.py

  #generate cfg file for the sample with no isolation cut
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.tauAnalysisSequence%"  -e "s%HIGGSREW%False%" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_all_cfg.py

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
cat <<EOF > runWbbTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *Wbb*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runWbbTauAnalyzerCfgs.sh

#generate script that submits all iso+nonIso+reweight jobs to LSF
cat <<EOF > submitWbbTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*Wbb*.sh | grep -v all | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 1nd -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitWbbTauAnalyzerJobs.sh

#generate script that submits all noIsoCut jobs to LSF
cat <<EOF > submitWbbAllTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*Wbb*all*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 1nd -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitWbbAllTauAnalyzerJobs.sh

#generate script that copies all iso+nonIso+reweight files locally from EOS
cat <<EOF > copyWbbFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  for cut in Iso NonIso NonIsoReweight
    do
    if [ "\$cut" != "NonIso" ] || [ $reweightOnly -eq 0 ]
        then
        cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_Wbb_${version}.root /data1/`whoami`/Wbb/analysis/
        cmsRm /store/user/`whoami`/muHad\${cut}Analysis_Wbb_${version}.root
    fi
  done
done

exit 0
EOF
chmod a+x copyWbbFromEOS.sh

#generate script that copies all noIsoCut files locally from EOS
cat <<EOF > copyAllWbbFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  cmsStage -f /store/user/`whoami`/muHadAnalysis_Wbb_${version}.root /data1/`whoami`/Wbb/analysis/
  cmsRm /store/user/`whoami`/muHadAnalysis_Wbb_${version}.root
done

exit 0
EOF
chmod a+x copyAllWbbFromEOS.sh

exit 0
