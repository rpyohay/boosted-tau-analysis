#!/bin/bash

if [ $# -gt 3 ]
    then
    echo "Usage: ./generateTTJetsTauAnalyzerCfgs.sh <version> <template cfg> [reweightOnly]"
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
inputFilePrefix="root://eoscms//eos/cms/store/user/yohay/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-Summer12_DR53X-PU_S10_START53_V7A-v2-AODSIM_skim_v2/"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/TTJets/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/TTJets/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_100_1_6HW.root',\n    '${inputFilePrefix}data_no_selection_101_1_tJp.root',\n    '${inputFilePrefix}data_no_selection_102_1_FTA.root',\n    '${inputFilePrefix}data_no_selection_103_1_ool.root',\n    '${inputFilePrefix}data_no_selection_104_1_ahq.root',\n    '${inputFilePrefix}data_no_selection_105_1_G2C.root',\n    '${inputFilePrefix}data_no_selection_106_1_Xr8.root',\n    '${inputFilePrefix}data_no_selection_107_1_9oR.root',\n    '${inputFilePrefix}data_no_selection_108_1_6Ff.root',\n    '${inputFilePrefix}data_no_selection_109_1_u4q.root',\n    '${inputFilePrefix}data_no_selection_10_1_CVv.root',\n    '${inputFilePrefix}data_no_selection_110_1_UN8.root',\n    '${inputFilePrefix}data_no_selection_111_1_Xgo.root',\n    '${inputFilePrefix}data_no_selection_112_1_WE4.root',\n    '${inputFilePrefix}data_no_selection_113_1_PlU.root',\n    '${inputFilePrefix}data_no_selection_114_1_V21.root',\n    '${inputFilePrefix}data_no_selection_115_1_zeL.root',\n    '${inputFilePrefix}data_no_selection_116_1_CwE.root',\n    '${inputFilePrefix}data_no_selection_117_1_NZd.root',\n    '${inputFilePrefix}data_no_selection_118_1_5KA.root',\n    '${inputFilePrefix}data_no_selection_119_1_URP.root',\n    '${inputFilePrefix}data_no_selection_11_1_fri.root',\n    '${inputFilePrefix}data_no_selection_120_1_ctA.root',\n    '${inputFilePrefix}data_no_selection_121_1_bUM.root',\n    '${inputFilePrefix}data_no_selection_122_1_OnI.root',\n    '${inputFilePrefix}data_no_selection_123_1_tDs.root',\n    '${inputFilePrefix}data_no_selection_124_1_uoi.root',\n    '${inputFilePrefix}data_no_selection_125_1_xeV.root',\n    '${inputFilePrefix}data_no_selection_126_1_ioa.root',\n    '${inputFilePrefix}data_no_selection_127_1_U4r.root',\n    '${inputFilePrefix}data_no_selection_128_1_EfX.root',\n    '${inputFilePrefix}data_no_selection_129_1_mce.root',\n    '${inputFilePrefix}data_no_selection_12_1_tXm.root',\n    '${inputFilePrefix}data_no_selection_130_1_qrx.root',\n    '${inputFilePrefix}data_no_selection_131_1_vrW.root',\n    '${inputFilePrefix}data_no_selection_132_1_clw.root',\n    '${inputFilePrefix}data_no_selection_133_1_rz0.root',\n    '${inputFilePrefix}data_no_selection_134_1_RmT.root',\n    '${inputFilePrefix}data_no_selection_135_1_G9F.root',\n    '${inputFilePrefix}data_no_selection_136_1_tmF.root',\n    '${inputFilePrefix}data_no_selection_137_1_hy9.root',\n    '${inputFilePrefix}data_no_selection_13_1_c54.root',\n    '${inputFilePrefix}data_no_selection_14_1_uvU.root',\n    '${inputFilePrefix}data_no_selection_15_1_1e3.root',\n    '${inputFilePrefix}data_no_selection_16_1_FWz.root',\n    '${inputFilePrefix}data_no_selection_17_1_QWQ.root',\n    '${inputFilePrefix}data_no_selection_18_1_JvZ.root',\n    '${inputFilePrefix}data_no_selection_19_1_mpq.root',\n    '${inputFilePrefix}data_no_selection_1_1_E7J.root',\n    '${inputFilePrefix}data_no_selection_20_1_o3B.root',\n    '${inputFilePrefix}data_no_selection_21_1_cEg.root',\n    '${inputFilePrefix}data_no_selection_22_1_Qkq.root',\n    '${inputFilePrefix}data_no_selection_23_1_pWl.root',\n    '${inputFilePrefix}data_no_selection_24_1_cpn.root',\n    '${inputFilePrefix}data_no_selection_25_1_stP.root',\n    '${inputFilePrefix}data_no_selection_26_1_wbR.root',\n    '${inputFilePrefix}data_no_selection_27_1_vfL.root',\n    '${inputFilePrefix}data_no_selection_28_1_Sg1.root',\n    '${inputFilePrefix}data_no_selection_29_1_ACC.root',\n    '${inputFilePrefix}data_no_selection_2_1_Hpw.root',\n    '${inputFilePrefix}data_no_selection_30_1_CHP.root',\n    '${inputFilePrefix}data_no_selection_31_1_sg6.root',\n    '${inputFilePrefix}data_no_selection_32_1_YpO.root',\n    '${inputFilePrefix}data_no_selection_33_1_Nva.root',\n    '${inputFilePrefix}data_no_selection_34_1_RYu.root',\n    '${inputFilePrefix}data_no_selection_35_1_g0B.root',\n    '${inputFilePrefix}data_no_selection_36_1_kLL.root',\n    '${inputFilePrefix}data_no_selection_37_1_WUn.root',\n    '${inputFilePrefix}data_no_selection_38_1_gTc.root',\n    '${inputFilePrefix}data_no_selection_39_1_9iG.root',\n    '${inputFilePrefix}data_no_selection_3_1_Mqo.root',\n    '${inputFilePrefix}data_no_selection_40_1_EUp.root',\n    '${inputFilePrefix}data_no_selection_41_1_zY9.root',\n    '${inputFilePrefix}data_no_selection_42_1_QRs.root',\n    '${inputFilePrefix}data_no_selection_43_1_82E.root',\n    '${inputFilePrefix}data_no_selection_44_1_xIU.root',\n    '${inputFilePrefix}data_no_selection_45_1_vmz.root',\n    '${inputFilePrefix}data_no_selection_46_1_ohz.root',\n    '${inputFilePrefix}data_no_selection_47_1_e24.root',\n    '${inputFilePrefix}data_no_selection_48_1_6IR.root',\n    '${inputFilePrefix}data_no_selection_49_1_wyu.root',\n    '${inputFilePrefix}data_no_selection_4_1_exB.root',\n    '${inputFilePrefix}data_no_selection_50_1_zcA.root',\n    '${inputFilePrefix}data_no_selection_51_1_e5w.root',\n    '${inputFilePrefix}data_no_selection_52_1_Kfv.root',\n    '${inputFilePrefix}data_no_selection_53_1_6vQ.root',\n    '${inputFilePrefix}data_no_selection_54_1_lBl.root',\n    '${inputFilePrefix}data_no_selection_55_1_x9L.root',\n    '${inputFilePrefix}data_no_selection_56_1_2tp.root',\n    '${inputFilePrefix}data_no_selection_57_1_fkh.root',\n    '${inputFilePrefix}data_no_selection_58_1_ZyX.root',\n    '${inputFilePrefix}data_no_selection_59_1_X2c.root',\n    '${inputFilePrefix}data_no_selection_5_1_6Tn.root',\n    '${inputFilePrefix}data_no_selection_60_1_WCg.root',\n    '${inputFilePrefix}data_no_selection_61_1_2dn.root',\n    '${inputFilePrefix}data_no_selection_62_1_Zmq.root',\n    '${inputFilePrefix}data_no_selection_63_1_Ij3.root',\n    '${inputFilePrefix}data_no_selection_64_1_g2l.root',\n    '${inputFilePrefix}data_no_selection_65_1_8CN.root',\n    '${inputFilePrefix}data_no_selection_66_1_9R9.root',\n    '${inputFilePrefix}data_no_selection_67_1_8Rj.root',\n    '${inputFilePrefix}data_no_selection_68_1_X1F.root',\n    '${inputFilePrefix}data_no_selection_69_1_Epu.root',\n    '${inputFilePrefix}data_no_selection_6_1_595.root',\n    '${inputFilePrefix}data_no_selection_70_1_3VC.root',\n    '${inputFilePrefix}data_no_selection_71_1_YO8.root',\n    '${inputFilePrefix}data_no_selection_72_1_y3F.root',\n    '${inputFilePrefix}data_no_selection_73_1_Owv.root',\n    '${inputFilePrefix}data_no_selection_74_1_mc4.root',\n    '${inputFilePrefix}data_no_selection_75_1_GgW.root',\n    '${inputFilePrefix}data_no_selection_76_1_K2C.root',\n    '${inputFilePrefix}data_no_selection_77_1_kEN.root',\n    '${inputFilePrefix}data_no_selection_78_1_xax.root',\n    '${inputFilePrefix}data_no_selection_79_1_RfC.root',\n    '${inputFilePrefix}data_no_selection_7_1_r6v.root',\n    '${inputFilePrefix}data_no_selection_80_1_Fdi.root',\n    '${inputFilePrefix}data_no_selection_81_1_J76.root',\n    '${inputFilePrefix}data_no_selection_82_1_qck.root',\n    '${inputFilePrefix}data_no_selection_83_1_T5H.root',\n    '${inputFilePrefix}data_no_selection_84_1_1HF.root',\n    '${inputFilePrefix}data_no_selection_85_1_wqo.root',\n    '${inputFilePrefix}data_no_selection_86_1_L2k.root',\n    '${inputFilePrefix}data_no_selection_87_1_C57.root',\n    '${inputFilePrefix}data_no_selection_88_1_j6H.root',\n    '${inputFilePrefix}data_no_selection_89_1_cDd.root',\n    '${inputFilePrefix}data_no_selection_8_1_LId.root',\n    '${inputFilePrefix}data_no_selection_90_1_PMj.root',\n    '${inputFilePrefix}data_no_selection_91_1_nsP.root',\n    '${inputFilePrefix}data_no_selection_92_1_AF8.root',\n    '${inputFilePrefix}data_no_selection_93_1_Wvp.root',\n    '${inputFilePrefix}data_no_selection_94_1_NFU.root',\n    '${inputFilePrefix}data_no_selection_95_1_gE3.root',\n    '${inputFilePrefix}data_no_selection_96_1_S8x.root',\n    '${inputFilePrefix}data_no_selection_97_1_VQP.root',\n    '${inputFilePrefix}data_no_selection_98_1_1y7.root',\n    '${inputFilePrefix}data_no_selection_99_1_o5p.root',\n    '${inputFilePrefix}data_no_selection_9_1_iiR.root',\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_TTJets.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_TTJets_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_TTJets_${version}.root" )
nonIsoReweightTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_TTJets_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_TTJets_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}TTJets${infoTag}_${version}.root" )

#samples
samples=( "TTJets" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" -e "s%BTAGSCALESHIFT%mean%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_cfg.py

  #generate cfg file for the non-isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" -e "s%BTAGSCALESHIFT%mean%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_nonIso_cfg.py

  #generate cfg file for the non-isolated, reweighted sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoReweightTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%True%" -e "s%BTAGSCALESHIFT%mean%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_nonIsoReweight_cfg.py

  #generate cfg file for the sample with no isolation cut
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.tauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" -e "s%BTAGSCALESHIFT%mean%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_all_cfg.py

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
#cmsRun \${fileNamePrefix}_nonIsoReweight_cfg.py
cmsStage -f ${isoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
#cmsStage -f ${nonIsoReweightTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
rm ${isoTauAnalyzerOutputFiles[${i}]}
#rm ${nonIsoReweightTauAnalyzerOutputFiles[${i}]} 

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
cat <<EOF > runTTJetsTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *TTJets*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runTTJetsTauAnalyzerCfgs.sh

#generate script that submits all iso+nonIso+reweight jobs to LSF
cat <<EOF > submitTTJetsTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*TTJets*.sh | grep -v all | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitTTJetsTauAnalyzerJobs.sh

#generate script that submits all noIsoCut jobs to LSF
cat <<EOF > submitTTJetsAllTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*TTJets*all*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitTTJetsAllTauAnalyzerJobs.sh

#generate script that copies all iso+nonIso+reweight files locally from EOS
cat <<EOF > copyTTJetsFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  #for cut in Iso NonIso NonIsoReweight
  for cut in Iso NonIso
    do
    if [ "\$cut" != "NonIso" ] || [ $reweightOnly -eq 0 ]
        then
        cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_TTJets_${version}.root /data1/`whoami`/TTJets/analysis/
        cmsRm /store/user/`whoami`/muHad\${cut}Analysis_TTJets_${version}.root
    fi
  done
done

exit 0
EOF
chmod a+x copyTTJetsFromEOS.sh

#generate script that copies all noIsoCut files locally from EOS
cat <<EOF > copyAllTTJetsFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  cmsStage -f /store/user/`whoami`/muHadAnalysis_TTJets_${version}.root /data1/`whoami`/TTJets/analysis/
  cmsRm /store/user/`whoami`/muHadAnalysis_TTJets_${version}.root
done

exit 0
EOF
chmod a+x copyAllTTJetsFromEOS.sh

exit 0
