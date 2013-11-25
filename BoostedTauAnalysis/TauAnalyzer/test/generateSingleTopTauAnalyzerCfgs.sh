#!/bin/bash

if [ $# -gt 2 ]
    then
    echo "Usage: ./generateSingleTopTauAnalyzerCfgs.sh <version> [reweightOnly|]"
    exit 0
fi

####STUFF TO CONFIGURE####

#version
version=$1
infoTag=""
reweightOnly=0
if [ "$2" == "reweightOnly" ]
    then
    reweightOnly=1
fi
dir=$version

#number of samples
nSamples=4
iBeg=0
iEnd=`expr $nSamples - 1`

#input file prefix and suffix
inputFilePrefix="root://eoscms//eos/cms/store/user/yohay/T"
inputFileSuffix="-channel_TuneZ2star_8TeV-powheg-tauola-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v2/"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/SingleTop/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/SingleTop/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}_s${inputFileSuffix}data_no_selection_1_1_jWr.root',\n    '${inputFilePrefix}_s${inputFileSuffix}data_no_selection_2_1_gDH.root',\n    '${inputFilePrefix}_s${inputFileSuffix}data_no_selection_3_1_ADm.root',\n    '${inputFilePrefix}_s${inputFileSuffix}data_no_selection_4_1_8SP.root',\n    '${inputFilePrefix}_s${inputFileSuffix}data_no_selection_5_1_jrE.root',\n    '${inputFilePrefix}_s${inputFileSuffix}data_no_selection_6_1_mfi.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}bar_s${inputFileSuffix}data_no_selection_1_1_lcL.root',\n    '${inputFilePrefix}bar_s${inputFileSuffix}data_no_selection_2_1_kQc.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_11_1_qpn.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_13_1_dgn.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_14_1_rIC.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_15_1_n4E.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_16_2_iQp.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_17_1_fW4.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_18_1_CKD.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_19_1_n9P.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_1_1_nCB.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_20_1_AwY.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_21_1_7Rf.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_22_1_p9P.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_23_1_fTS.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_24_1_dEA.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_25_1_AyZ.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_26_1_1W4.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_27_1_1EA.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_28_1_KOU.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_29_1_rGL.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_2_1_UW9.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_30_1_HfT.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_31_1_F4C.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_32_1_eFW.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_34_2_CqN.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_35_1_b5D.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_36_1_vMy.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_37_1_8TV.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_39_1_Igs.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_3_1_cfn.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_40_1_sQw.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_41_1_yB7.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_42_1_fx7.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_43_1_baR.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_44_1_z87.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_45_1_SiT.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_46_1_IVa.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_48_1_us0.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_49_1_WO2.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_4_1_VUo.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_50_1_RP3.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_51_1_UEN.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_53_1_b3G.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_54_1_0yk.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_55_2_QQu.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_56_1_ZM2.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_57_1_EQw.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_58_1_Ikx.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_59_2_4Xp.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_5_2_Sd7.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_60_1_Rxt.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_61_1_qkG.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_62_1_isY.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_63_1_LIa.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_64_1_pSj.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_65_1_8qR.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_66_1_mR0.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_67_1_gbf.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_68_2_OsF.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_69_1_ttn.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_6_1_CmV.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_70_1_T7f.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_71_1_8Xh.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_72_1_Fi6.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_73_1_eIN.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_74_1_rV4.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_75_1_lCi.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_76_1_Q70.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_77_1_3jW.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_7_1_KZb.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_8_1_4nn.root',\n    '${inputFilePrefix}_t${inputFileSuffix}data_no_selection_9_1_suA.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_10_2_Yx8.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_11_2_HZb.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_12_2_Iyq.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_13_2_a1l.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_14_2_hdc.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_15_2_F9i.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_16_2_N88.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_17_2_LYN.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_18_2_ohC.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_19_2_M7f.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_1_2_4yc.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_20_3_Y6b.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_21_2_s5M.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_22_2_Chp.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_23_2_Pvn.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_24_2_uRR.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_25_3_v6g.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_26_3_gAj.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_27_2_qtp.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_28_2_zHw.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_29_2_tbD.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_2_2_myu.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_30_2_u8x.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_31_2_Ffv.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_32_2_QEE.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_33_2_7wZ.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_34_2_wDy.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_35_2_yHq.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_36_2_Dth.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_37_2_lvp.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_38_2_Jmz.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_39_2_r82.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_3_2_jfk.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_4_2_jht.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_5_2_2XM.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_6_2_s56.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_7_2_vxX.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_8_2_07G.root',\n    '${inputFilePrefix}bar_t${inputFileSuffix}data_no_selection_9_2_Zoe.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_T_s-channel.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Tbar_s-channel.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_T_t-channel.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Tbar_t-channel.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_T_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_Tbar_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_T_t-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_Tbar_t-channel_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_T_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_Tbar_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_T_t-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_Tbar_t-channel_${version}.root" )
nonIsoReweightTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_T_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_Tbar_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_T_t-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_Tbar_t-channel_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_T_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_Tbar_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_T_t-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_Tbar_t-channel_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}T_s-channel${infoTag}_${version}.root" "${EDMOutputFilePrefix}Tbar_s-channel${infoTag}_${version}.root" "${EDMOutputFilePrefix}T_t-channel${infoTag}_${version}.root" "${EDMOutputFilePrefix}Tbar_t-channel${infoTag}_${version}.root" )

#samples
samples=( "T_s-channel" "Tbar_s-channel" "T_t-channel" "Tbar_t-channel" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" ../tauanalyzer_WNJetsToLNu_Wh1_template_cfg.py > tauanalyzer_${samples[${i}]}_iso_cfg.py

  #generate cfg file for the non-isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" ../tauanalyzer_WNJetsToLNu_Wh1_template_cfg.py > tauanalyzer_${samples[${i}]}_nonIso_cfg.py

  #generate cfg file for the non-isolated, reweighted sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoReweightTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%True%" ../tauanalyzer_WNJetsToLNu_Wh1_template_cfg.py > tauanalyzer_${samples[${i}]}_nonIsoReweight_cfg.py

  #generate cfg file for the sample with no isolation cut
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.tauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" ../tauanalyzer_WNJetsToLNu_Wh1_template_cfg.py > tauanalyzer_${samples[${i}]}_all_cfg.py

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
cat <<EOF > runSingleTopTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *T*channel*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runSingleTopTauAnalyzerCfgs.sh

#generate script that submits all iso+nonIso+reweight jobs to LSF
cat <<EOF > submitSingleTopTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*T*channel*.sh | grep -v all | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitSingleTopTauAnalyzerJobs.sh

#generate script that submits all noIsoCut jobs to LSF
cat <<EOF > submitSingleTopAllTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*T*channel*all*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitSingleTopAllTauAnalyzerJobs.sh

#generate script that copies all iso+nonIso+reweight files locally from EOS
cat <<EOF > copySingleTopFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in "_s" "bar_s" "_t" "bar_t"
  do
  for cut in Iso NonIso NonIsoReweight
    do
    if [ "\$cut" != "NonIso" ] || [ $reweightOnly -eq 0 ]
        then
        cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_T\${sample}-channel_${version}.root /data1/`whoami`/SingleTop/analysis/
        cmsRm /store/user/`whoami`/muHad\${cut}Analysis_T\${sample}-channel_${version}.root
    fi
  done
done

exit 0
EOF
chmod a+x copySingleTopFromEOS.sh

#generate script that copies all noIsoCut files locally from EOS
cat <<EOF > copyAllSingleTopFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in "_s" "bar_s" "_t" "bar_t"
  do
  cmsStage -f /store/user/`whoami`/muHadAnalysis_T\${sample}-channel_${version}.root /data1/`whoami`/SingleTop/analysis/
  cmsRm /store/user/`whoami`/muHadAnalysis_T\${sample}-channel_${version}.root
done

exit 0
EOF
chmod a+x copyAllSingleTopFromEOS.sh

exit 0
