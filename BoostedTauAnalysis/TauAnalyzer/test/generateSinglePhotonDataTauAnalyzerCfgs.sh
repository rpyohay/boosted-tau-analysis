#!/bin/bash

if [ $# -gt 3 ]
    then
    echo "Usage: ./generateSinglePhotonDataTauAnalyzerCfgs.sh <version> <template cfg> [reweightOnly]"
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

#number of jobs
nJobs=12
iBeg=0
iEndJob=`expr $nJobs - 1`

#samples
samples=( "SinglePhotonParked_Run2012D" )
nSamples=${#samples[@]}
iEndSample=`expr $nSamples - 1`

#input file prefix
inputFilePrefix="root://eoscms//eos/cms/store/user/yohay/SinglePhotonParked_Run2012D-22Jan2013-v1_AOD_skim_v2/"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/data/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/data/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_100_1_DrU.root',\n    '${inputFilePrefix}data_no_selection_101_1_ZUk.root',\n    '${inputFilePrefix}data_no_selection_102_1_VNc.root',\n    '${inputFilePrefix}data_no_selection_103_1_gBZ.root',\n    '${inputFilePrefix}data_no_selection_104_2_Gv7.root',\n    '${inputFilePrefix}data_no_selection_105_2_hPP.root',\n    '${inputFilePrefix}data_no_selection_106_2_HuQ.root',\n    '${inputFilePrefix}data_no_selection_107_2_OPM.root',\n    '${inputFilePrefix}data_no_selection_108_1_c8r.root',\n    '${inputFilePrefix}data_no_selection_109_1_gN7.root',\n    '${inputFilePrefix}data_no_selection_110_1_i8j.root',\n    '${inputFilePrefix}data_no_selection_111_2_OTr.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_112_1_uiM.root',\n    '${inputFilePrefix}data_no_selection_113_1_BE2.root',\n    '${inputFilePrefix}data_no_selection_114_1_ysW.root',\n    '${inputFilePrefix}data_no_selection_115_3_y5B.root',\n    '${inputFilePrefix}data_no_selection_116_1_b1m.root',\n    '${inputFilePrefix}data_no_selection_117_2_WTS.root',\n    '${inputFilePrefix}data_no_selection_118_1_aSV.root',\n    '${inputFilePrefix}data_no_selection_119_1_Bsq.root',\n    '${inputFilePrefix}data_no_selection_11_2_1tR.root',\n    '${inputFilePrefix}data_no_selection_120_1_dXW.root',\n    '${inputFilePrefix}data_no_selection_121_2_gQv.root',\n    '${inputFilePrefix}data_no_selection_122_2_zwP.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_123_1_DPH.root',\n    '${inputFilePrefix}data_no_selection_124_1_otq.root',\n    '${inputFilePrefix}data_no_selection_125_1_imY.root',\n    '${inputFilePrefix}data_no_selection_126_1_Abu.root',\n    '${inputFilePrefix}data_no_selection_127_1_grj.root',\n    '${inputFilePrefix}data_no_selection_128_2_6i5.root',\n    '${inputFilePrefix}data_no_selection_129_2_g1Q.root',\n    '${inputFilePrefix}data_no_selection_12_1_p9m.root',\n    '${inputFilePrefix}data_no_selection_130_1_tcf.root',\n    '${inputFilePrefix}data_no_selection_131_1_gXy.root',\n    '${inputFilePrefix}data_no_selection_132_2_krp.root',\n    '${inputFilePrefix}data_no_selection_133_2_uSl.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_134_2_EjD.root',\n    '${inputFilePrefix}data_no_selection_135_2_oIS.root',\n    '${inputFilePrefix}data_no_selection_136_1_0oB.root',\n    '${inputFilePrefix}data_no_selection_137_1_VeC.root',\n    '${inputFilePrefix}data_no_selection_138_2_A74.root',\n    '${inputFilePrefix}data_no_selection_139_1_rGN.root',\n    '${inputFilePrefix}data_no_selection_13_1_SHe.root',\n    '${inputFilePrefix}data_no_selection_140_1_C6W.root',\n    '${inputFilePrefix}data_no_selection_141_1_qtK.root',\n    '${inputFilePrefix}data_no_selection_142_1_y21.root',\n    '${inputFilePrefix}data_no_selection_143_1_joJ.root',\n    '${inputFilePrefix}data_no_selection_144_1_xQi.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_145_1_P85.root',\n    '${inputFilePrefix}data_no_selection_146_1_1bh.root',\n    '${inputFilePrefix}data_no_selection_147_1_B5H.root',\n    '${inputFilePrefix}data_no_selection_148_1_1CW.root',\n    '${inputFilePrefix}data_no_selection_149_1_Ubi.root',\n    '${inputFilePrefix}data_no_selection_14_1_G6H.root',\n    '${inputFilePrefix}data_no_selection_15_2_94w.root',\n    '${inputFilePrefix}data_no_selection_16_2_wI2.root',\n    '${inputFilePrefix}data_no_selection_17_1_G2Z.root',\n    '${inputFilePrefix}data_no_selection_18_2_seS.root',\n    '${inputFilePrefix}data_no_selection_19_1_z2F.root',\n    '${inputFilePrefix}data_no_selection_1_1_H8x.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_20_1_eVr.root',\n    '${inputFilePrefix}data_no_selection_21_2_UIq.root',\n    '${inputFilePrefix}data_no_selection_22_1_c2R.root',\n    '${inputFilePrefix}data_no_selection_23_2_wQy.root',\n    '${inputFilePrefix}data_no_selection_24_1_Sq2.root',\n    '${inputFilePrefix}data_no_selection_25_1_IIW.root',\n    '${inputFilePrefix}data_no_selection_26_1_Ilt.root',\n    '${inputFilePrefix}data_no_selection_27_1_dAE.root',\n    '${inputFilePrefix}data_no_selection_28_1_YUY.root',\n    '${inputFilePrefix}data_no_selection_29_1_HZH.root',\n    '${inputFilePrefix}data_no_selection_2_1_IOH.root',\n    '${inputFilePrefix}data_no_selection_31_2_gcs.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_32_1_cmL.root',\n    '${inputFilePrefix}data_no_selection_34_2_1wG.root',\n    '${inputFilePrefix}data_no_selection_35_1_6zL.root',\n    '${inputFilePrefix}data_no_selection_36_1_Zht.root',\n    '${inputFilePrefix}data_no_selection_37_1_mOh.root',\n    '${inputFilePrefix}data_no_selection_38_1_zoj.root',\n    '${inputFilePrefix}data_no_selection_39_1_6cO.root',\n    '${inputFilePrefix}data_no_selection_3_1_j2Z.root',\n    '${inputFilePrefix}data_no_selection_40_2_JHr.root',\n    '${inputFilePrefix}data_no_selection_41_2_OnS.root',\n    '${inputFilePrefix}data_no_selection_42_1_lox.root',\n    '${inputFilePrefix}data_no_selection_44_1_Qtu.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_45_1_X6C.root',\n    '${inputFilePrefix}data_no_selection_46_1_fvu.root',\n    '${inputFilePrefix}data_no_selection_47_2_aD0.root',\n    '${inputFilePrefix}data_no_selection_48_2_r2T.root',\n    '${inputFilePrefix}data_no_selection_49_1_yG1.root',\n    '${inputFilePrefix}data_no_selection_4_1_ley.root',\n    '${inputFilePrefix}data_no_selection_50_1_Wzj.root',\n    '${inputFilePrefix}data_no_selection_51_1_pY2.root',\n    '${inputFilePrefix}data_no_selection_52_3_z27.root',\n    '${inputFilePrefix}data_no_selection_53_2_N4j.root',\n    '${inputFilePrefix}data_no_selection_54_2_1NV.root',\n    '${inputFilePrefix}data_no_selection_55_2_5ES.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_56_2_WGo.root',\n    '${inputFilePrefix}data_no_selection_58_2_V3d.root',\n    '${inputFilePrefix}data_no_selection_59_2_ych.root',\n    '${inputFilePrefix}data_no_selection_5_1_bwB.root',\n    '${inputFilePrefix}data_no_selection_60_2_XSx.root',\n    '${inputFilePrefix}data_no_selection_61_2_12I.root',\n    '${inputFilePrefix}data_no_selection_62_2_PWr.root',\n    '${inputFilePrefix}data_no_selection_63_1_li7.root',\n    '${inputFilePrefix}data_no_selection_64_2_76a.root',\n    '${inputFilePrefix}data_no_selection_65_2_pQV.root',\n    '${inputFilePrefix}data_no_selection_66_1_wqx.root',\n    '${inputFilePrefix}data_no_selection_67_1_N8F.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_68_2_il6.root',\n    '${inputFilePrefix}data_no_selection_69_1_wjg.root',\n    '${inputFilePrefix}data_no_selection_6_2_v5J.root',\n    '${inputFilePrefix}data_no_selection_70_2_7Ev.root',\n    '${inputFilePrefix}data_no_selection_71_2_rUI.root',\n    '${inputFilePrefix}data_no_selection_72_3_mpd.root',\n    '${inputFilePrefix}data_no_selection_73_1_4jt.root',\n    '${inputFilePrefix}data_no_selection_74_2_Tp6.root',\n    '${inputFilePrefix}data_no_selection_75_1_3AU.root',\n    '${inputFilePrefix}data_no_selection_76_1_ISy.root',\n    '${inputFilePrefix}data_no_selection_77_1_Fd7.root',\n    '${inputFilePrefix}data_no_selection_78_1_BLV.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_79_1_1zJ.root',\n    '${inputFilePrefix}data_no_selection_7_2_TfR.root',\n    '${inputFilePrefix}data_no_selection_80_1_Oh7.root',\n    '${inputFilePrefix}data_no_selection_81_1_hnD.root',\n    '${inputFilePrefix}data_no_selection_82_1_D0L.root',\n    '${inputFilePrefix}data_no_selection_83_1_g1A.root',\n    '${inputFilePrefix}data_no_selection_84_1_ti4.root',\n    '${inputFilePrefix}data_no_selection_85_2_Qdy.root',\n    '${inputFilePrefix}data_no_selection_86_2_DOu.root',\n    '${inputFilePrefix}data_no_selection_87_1_De3.root',\n    '${inputFilePrefix}data_no_selection_88_1_iBP.root',\n    '${inputFilePrefix}data_no_selection_89_1_WIM.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_8_1_pSf.root',\n    '${inputFilePrefix}data_no_selection_90_1_vGG.root',\n    '${inputFilePrefix}data_no_selection_91_1_6YT.root',\n    '${inputFilePrefix}data_no_selection_92_1_y9S.root',\n    '${inputFilePrefix}data_no_selection_93_1_aKc.root',\n    '${inputFilePrefix}data_no_selection_94_1_Ktf.root',\n    '${inputFilePrefix}data_no_selection_95_1_2Uu.root',\n    '${inputFilePrefix}data_no_selection_96_1_bSs.root',\n    '${inputFilePrefix}data_no_selection_97_1_tmr.root',\n    '${inputFilePrefix}data_no_selection_98_1_mA2.root',\n    '${inputFilePrefix}data_no_selection_99_1_ctV.root',\n    '${inputFilePrefix}data_no_selection_9_1_Y9h.root'\n    ])" )

#arrays of output files
cleanJetsOutFiles=( )
isoTauAnalyzerOutputFiles=( )
nonIsoTauAnalyzerOutputFiles=( )
nonIsoReweightTauAnalyzerOutputFiles=( )
allTauAnalyzerOutputFiles=( )
EDMOutputFiles=( )

#loop over number of samples
for iSample in `seq $iBeg $iEndSample`
  do

  #loop over number of jobs
  for iJob in `seq $iBeg $iEndJob`
    do
    hashedIndex=`expr ${iSample} \* ${nJobs} + ${iJob}`

    #CleanJets output file
    cleanJetsOutFiles[$hashedIndex]="${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_${samples[${iSample}]}_${iJob}.root"

    #TauAnalyzer output files
    isoTauAnalyzerOutputFiles[$hashedIndex]="${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_${samples[${iSample}]}_${iJob}_${version}.root"
    nonIsoTauAnalyzerOutputFiles[$hashedIndex]="${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_${samples[${iSample}]}_${iJob}_${version}.root"
    nonIsoReweightTauAnalyzerOutputFiles[$hashedIndex]="${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_${samples[${iSample}]}_${iJob}_${version}.root"
    allTauAnalyzerOutputFiles[$hashedIndex]="${tauAnalyzerOutputFilePrefix}muHadAnalysis_${samples[${iSample}]}_${iJob}_${version}.root"

    #EDM output files
    EDMOutputFiles[$hashedIndex]="${EDMOutputFilePrefix}${samples[${iSample}]}${infoTag}_${iJob}_${version}.root"
  done
done

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for iSample in `seq $iBeg $iEndSample`
  do

  #loop over number of jobs
  for iJob in `seq $iBeg $iEndJob`
    do
    hashedIndex=`expr ${iSample} \* ${nJobs} + ${iJob}`

    #generate cfg file for the isolated sample
    sed -e "s%FILES%${inputFileBlocks[${hashedIndex}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${hashedIndex}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${hashedIndex}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%CUSTOMTAUSELECTOR%CustomTauSepFromPhotonSelector%" -e "s%OVERLAPCANDTAG%photonSelector%g" -e "s%MUONORPHOTONTAG%photonTag = cms.InputTag('photonSelector')%" -e "s%TRIGGEROBJECTFILTER%photonTriggerObjectFilter%g" -e "s%OSSFFILTERISO%%g" -e "s%OSSFFILTERNONISO%%g" -e "s%OSSFFILTER%%g" ../${templateCfg} > tauanalyzer_${samples[${iSample}]}_${iJob}_iso_cfg.py

    #generate cfg file for the non-isolated sample
    sed -e "s%FILES%${inputFileBlocks[${hashedIndex}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${hashedIndex}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${hashedIndex}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%CUSTOMTAUSELECTOR%CustomTauSepFromPhotonSelector%" -e "s%OVERLAPCANDTAG%photonSelector%g" -e "s%MUONORPHOTONTAG%photonTag = cms.InputTag('photonSelector')%" -e "s%TRIGGEROBJECTFILTER%photonTriggerObjectFilter%g" -e "s%OSSFFILTERISO%%g" -e "s%OSSFFILTERNONISO%%g" -e "s%OSSFFILTER%%g" ../${templateCfg} > tauanalyzer_${samples[${iSample}]}_${iJob}_nonIso_cfg.py

    #generate cfg file for the non-isolated, reweighted sample
    sed -e "s%FILES%${inputFileBlocks[${hashedIndex}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${hashedIndex}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoReweightTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${hashedIndex}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%CUSTOMTAUSELECTOR%CustomTauSepFromPhotonSelector%" -e "s%OVERLAPCANDTAG%photonSelector%g" -e "s%MUONORPHOTONTAG%photonTag = cms.InputTag('photonSelector')%" -e "s%TRIGGEROBJECTFILTER%photonTriggerObjectFilter%g" -e "s%OSSFFILTERISO%%g" -e "s%OSSFFILTERNONISO%%g" -e "s%OSSFFILTER%%g" ../${templateCfg} > tauanalyzer_${samples[${iSample}]}_${iJob}_nonIsoReweight_cfg.py

    #generate cfg file for the sample with no isolation cut
    sed -e "s%FILES%${inputFileBlocks[${hashedIndex}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${hashedIndex}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${hashedIndex}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${hashedIndex}]}%" -e "s%SEQUENCE%process.tauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%CUSTOMTAUSELECTOR%CustomTauSepFromPhotonSelector%" -e "s%OVERLAPCANDTAG%photonSelector%g" -e "s%MUONORPHOTONTAG%photonTag = cms.InputTag('photonSelector')%" -e "s%TRIGGEROBJECTFILTER%photonTriggerObjectFilter%g" -e "s%OSSFFILTERISO%%g" -e "s%OSSFFILTERNONISO%%g" -e "s%OSSFFILTER%%g" ../${templateCfg} > tauanalyzer_${samples[${iSample}]}_${iJob}_all_cfg.py

    #generate iso+nonIso+reweight job submission script for LSF
    cat <<EOF > tauanalyzer_${samples[${iSample}]}_${iJob}_cfg.sh
#!/bin/bash

jobDir="`pwd`"
fileNamePrefix="tauanalyzer_${samples[${iSample}]}_${iJob}"
cd \$jobDir
eval \`scramv1 runtime -sh\`
cd -
cp \$jobDir/\${fileNamePrefix}_iso_cfg.py \$jobDir/\${fileNamePrefix}_nonIso_cfg.py \$jobDir/\${fileNamePrefix}_nonIsoReweight_cfg.py .
cmsRun \${fileNamePrefix}_iso_cfg.py
if [ $reweightOnly -eq 0 ]
    then
    cmsRun \${fileNamePrefix}_nonIso_cfg.py
    cmsStage -f ${nonIsoTauAnalyzerOutputFiles[${hashedIndex}]} /store/user/`whoami`/
    rm ${nonIsoTauAnalyzerOutputFiles[${hashedIndex}]}
fi
cmsRun \${fileNamePrefix}_nonIsoReweight_cfg.py
cmsStage -f ${isoTauAnalyzerOutputFiles[${hashedIndex}]} /store/user/`whoami`/
cmsStage -f ${nonIsoReweightTauAnalyzerOutputFiles[${hashedIndex}]} /store/user/`whoami`/
rm ${isoTauAnalyzerOutputFiles[${hashedIndex}]}
rm ${nonIsoReweightTauAnalyzerOutputFiles[${hashedIndex}]}

exit 0
EOF
    chmod a+x tauanalyzer_${samples[${iSample}]}_${iJob}_cfg.sh

    #generate noIsoCut job submission script for LSF
    cat <<EOF > tauanalyzer_${samples[${iSample}]}_${iJob}_all_cfg.sh
#!/bin/bash

jobDir="`pwd`"
fileNamePrefix="tauanalyzer_${samples[${iSample}]}_${iJob}"
cd \$jobDir
eval \`scramv1 runtime -sh\`
cd -
cp \$jobDir/\${fileNamePrefix}_all_cfg.py .
cmsRun \${fileNamePrefix}_all_cfg.py
cmsStage -f ${allTauAnalyzerOutputFiles[${hashedIndex}]} /store/user/`whoami`/
rm ${allTauAnalyzerOutputFiles[${hashedIndex}]}

exit 0
EOF
    chmod a+x tauanalyzer_${samples[${iSample}]}_${iJob}_all_cfg.sh
  done
done

#generate run cfg that runs all files in the directory
cat <<EOF > runSinglePhotonDataTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *SinglePhotonParked*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  isIso=\`echo \$file | sed -e "s%.*\(iso\).*%\1%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runSinglePhotonDataTauAnalyzerCfgs.sh

#generate script that submits all iso+nonIso+reweight jobs to LSF
cat <<EOF > submitSinglePhotonDataTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*SinglePhotonParked*.sh | grep -v all | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitSinglePhotonDataTauAnalyzerJobs.sh

#generate script that submits all noIsoCut jobs to LSF
cat <<EOF > submitSinglePhotonDataAllTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*SinglePhotonParked*all*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitSinglePhotonDataAllTauAnalyzerJobs.sh

#generate script that copies all iso+nonIso+reweight files locally from EOS
cat <<EOF > copySinglePhotonDataFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in ${samples}
  do
  for iJob in \`seq ${iBeg} ${iEndJob}\`
    do
    for cut in Iso NonIso NonIsoReweight
      do
      if [ "\$cut" != "NonIso" ] || [ $reweightOnly -eq 0 ]
          then
          cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_\${sample}_\${iJob}_${version}.root /data1/`whoami`/SinglePhotonParkedData/analysis/
          cmsRm /store/user/`whoami`/muHad\${cut}Analysis_\${sample}_\${iJob}_${version}.root
      fi
    done
  done
done
exit 0
EOF
chmod a+x copySinglePhotonDataFromEOS.sh

#generate script that copies all noIsoCut files locally from EOS
cat <<EOF > copyAllSinglePhotonDataFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in ${samples}
  do
  for iJob in \`seq ${iBeg} ${iEndJob}\`
    do
    cmsStage -f /store/user/`whoami`/muHadAnalysis_\${sample}_\${iJob}_${version}.root /data1/`whoami`/SinglePhotonParkedData/analysis/
    cmsRm /store/user/`whoami`/muHadAnalysis_\${sample}_\${iJob}_${version}.root
  done
done
exit 0
EOF
chmod a+x copyAllSinglePhotonDataFromEOS.sh

exit 0
