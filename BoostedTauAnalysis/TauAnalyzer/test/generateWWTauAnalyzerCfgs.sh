#!/bin/bash

if [ $# -gt 3 ]
    then
    echo "Usage: ./generateWWTauAnalyzerCfgs.sh <version> <template cfg> [reweightOnly]"
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
inputFilePrefix="root://eoscms//eos/cms/store/user/yohay/WW_TuneZ2star_8TeV_pythia6_tauola-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v2/"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/WW/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/WW/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_100_1_DRZ.root',\n    '${inputFilePrefix}data_no_selection_101_1_ctI.root',\n    '${inputFilePrefix}data_no_selection_102_1_rUm.root',\n    '${inputFilePrefix}data_no_selection_103_1_AEb.root',\n    '${inputFilePrefix}data_no_selection_104_1_BAr.root',\n    '${inputFilePrefix}data_no_selection_105_1_tJc.root',\n    '${inputFilePrefix}data_no_selection_106_1_NSA.root',\n    '${inputFilePrefix}data_no_selection_107_1_noF.root',\n    '${inputFilePrefix}data_no_selection_108_1_7kf.root',\n    '${inputFilePrefix}data_no_selection_10_1_4RA.root',\n    '${inputFilePrefix}data_no_selection_111_1_QCU.root',\n    '${inputFilePrefix}data_no_selection_112_1_ow5.root',\n    '${inputFilePrefix}data_no_selection_113_1_tYY.root',\n    '${inputFilePrefix}data_no_selection_114_1_ADb.root',\n    '${inputFilePrefix}data_no_selection_115_1_ie1.root',\n    '${inputFilePrefix}data_no_selection_116_1_bxF.root',\n    '${inputFilePrefix}data_no_selection_117_1_ee3.root',\n    '${inputFilePrefix}data_no_selection_118_1_8Fd.root',\n    '${inputFilePrefix}data_no_selection_119_1_pwr.root',\n    '${inputFilePrefix}data_no_selection_11_1_vwY.root',\n    '${inputFilePrefix}data_no_selection_120_1_eud.root',\n    '${inputFilePrefix}data_no_selection_121_1_ksy.root',\n    '${inputFilePrefix}data_no_selection_122_1_gjQ.root',\n    '${inputFilePrefix}data_no_selection_123_1_3Nu.root',\n    '${inputFilePrefix}data_no_selection_124_1_Hl7.root',\n    '${inputFilePrefix}data_no_selection_125_1_H6e.root',\n    '${inputFilePrefix}data_no_selection_128_1_0qW.root',\n    '${inputFilePrefix}data_no_selection_12_1_ZFL.root',\n    '${inputFilePrefix}data_no_selection_130_1_b43.root',\n    '${inputFilePrefix}data_no_selection_131_1_bxv.root',\n    '${inputFilePrefix}data_no_selection_133_1_ZE5.root',\n    '${inputFilePrefix}data_no_selection_134_1_Liz.root',\n    '${inputFilePrefix}data_no_selection_135_1_4Nr.root',\n    '${inputFilePrefix}data_no_selection_136_1_kLf.root',\n    '${inputFilePrefix}data_no_selection_137_1_wXE.root',\n    '${inputFilePrefix}data_no_selection_138_1_BJ8.root',\n    '${inputFilePrefix}data_no_selection_139_1_UHJ.root',\n    '${inputFilePrefix}data_no_selection_13_1_tWP.root',\n    '${inputFilePrefix}data_no_selection_141_1_AXt.root',\n    '${inputFilePrefix}data_no_selection_142_1_aFe.root',\n    '${inputFilePrefix}data_no_selection_145_1_gUH.root',\n    '${inputFilePrefix}data_no_selection_146_1_rK8.root',\n    '${inputFilePrefix}data_no_selection_147_1_Zrf.root',\n    '${inputFilePrefix}data_no_selection_148_1_TXe.root',\n    '${inputFilePrefix}data_no_selection_149_1_4Ig.root',\n    '${inputFilePrefix}data_no_selection_151_1_cA7.root',\n    '${inputFilePrefix}data_no_selection_152_1_yVD.root',\n    '${inputFilePrefix}data_no_selection_154_1_csM.root',\n    '${inputFilePrefix}data_no_selection_155_1_AG0.root',\n    '${inputFilePrefix}data_no_selection_156_1_Kkp.root',\n    '${inputFilePrefix}data_no_selection_157_1_Gl1.root',\n    '${inputFilePrefix}data_no_selection_158_1_iid.root',\n    '${inputFilePrefix}data_no_selection_159_1_w8z.root',\n    '${inputFilePrefix}data_no_selection_160_1_a7l.root',\n    '${inputFilePrefix}data_no_selection_161_1_fdI.root',\n    '${inputFilePrefix}data_no_selection_162_1_O2J.root',\n    '${inputFilePrefix}data_no_selection_163_1_AJo.root',\n    '${inputFilePrefix}data_no_selection_165_1_dcM.root',\n    '${inputFilePrefix}data_no_selection_166_1_ep2.root',\n    '${inputFilePrefix}data_no_selection_167_1_GBI.root',\n    '${inputFilePrefix}data_no_selection_16_1_YQj.root',\n    '${inputFilePrefix}data_no_selection_171_1_PCl.root',\n    '${inputFilePrefix}data_no_selection_172_1_uiE.root',\n    '${inputFilePrefix}data_no_selection_173_1_P6a.root',\n    '${inputFilePrefix}data_no_selection_174_1_egK.root',\n    '${inputFilePrefix}data_no_selection_175_1_1xZ.root',\n    '${inputFilePrefix}data_no_selection_176_1_yn3.root',\n    '${inputFilePrefix}data_no_selection_177_1_RqC.root',\n    '${inputFilePrefix}data_no_selection_178_1_MXd.root',\n    '${inputFilePrefix}data_no_selection_179_1_pK4.root',\n    '${inputFilePrefix}data_no_selection_17_1_FUR.root',\n    '${inputFilePrefix}data_no_selection_180_1_4mD.root',\n    '${inputFilePrefix}data_no_selection_181_1_4fb.root',\n    '${inputFilePrefix}data_no_selection_182_1_Yr0.root',\n    '${inputFilePrefix}data_no_selection_183_1_rUn.root',\n    '${inputFilePrefix}data_no_selection_184_1_LDw.root',\n    '${inputFilePrefix}data_no_selection_185_1_HSu.root',\n    '${inputFilePrefix}data_no_selection_186_1_4Q2.root',\n    '${inputFilePrefix}data_no_selection_187_1_tYs.root',\n    '${inputFilePrefix}data_no_selection_188_1_0X4.root',\n    '${inputFilePrefix}data_no_selection_189_1_Eby.root',\n    '${inputFilePrefix}data_no_selection_190_1_eW2.root',\n    '${inputFilePrefix}data_no_selection_191_1_4kY.root',\n    '${inputFilePrefix}data_no_selection_192_1_YXn.root',\n    '${inputFilePrefix}data_no_selection_193_1_o1F.root',\n    '${inputFilePrefix}data_no_selection_194_1_YmJ.root',\n    '${inputFilePrefix}data_no_selection_195_1_Ety.root',\n    '${inputFilePrefix}data_no_selection_196_1_4AU.root',\n    '${inputFilePrefix}data_no_selection_197_1_Fyj.root',\n    '${inputFilePrefix}data_no_selection_198_1_zp0.root',\n    '${inputFilePrefix}data_no_selection_199_1_ila.root',\n    '${inputFilePrefix}data_no_selection_19_1_Vgf.root',\n    '${inputFilePrefix}data_no_selection_1_1_K6P.root',\n    '${inputFilePrefix}data_no_selection_200_1_p65.root',\n    '${inputFilePrefix}data_no_selection_201_1_MKl.root',\n    '${inputFilePrefix}data_no_selection_202_1_HeF.root',\n    '${inputFilePrefix}data_no_selection_21_1_Jci.root',\n    '${inputFilePrefix}data_no_selection_22_1_vye.root',\n    '${inputFilePrefix}data_no_selection_23_1_ydq.root',\n    '${inputFilePrefix}data_no_selection_24_1_BkI.root',\n    '${inputFilePrefix}data_no_selection_25_1_0Pg.root',\n    '${inputFilePrefix}data_no_selection_26_1_nmb.root',\n    '${inputFilePrefix}data_no_selection_27_1_IV6.root',\n    '${inputFilePrefix}data_no_selection_28_1_Gmk.root',\n    '${inputFilePrefix}data_no_selection_29_1_FHL.root',\n    '${inputFilePrefix}data_no_selection_2_1_zeH.root',\n    '${inputFilePrefix}data_no_selection_30_1_ZD4.root',\n    '${inputFilePrefix}data_no_selection_31_1_KoR.root',\n    '${inputFilePrefix}data_no_selection_33_1_JhR.root',\n    '${inputFilePrefix}data_no_selection_34_1_yYc.root',\n    '${inputFilePrefix}data_no_selection_36_1_MCu.root',\n    '${inputFilePrefix}data_no_selection_38_1_5r1.root',\n    '${inputFilePrefix}data_no_selection_39_1_Hek.root',\n    '${inputFilePrefix}data_no_selection_3_1_q0w.root',\n    '${inputFilePrefix}data_no_selection_40_1_6vC.root',\n    '${inputFilePrefix}data_no_selection_41_1_TSi.root',\n    '${inputFilePrefix}data_no_selection_42_1_Vcp.root',\n    '${inputFilePrefix}data_no_selection_43_1_SpG.root',\n    '${inputFilePrefix}data_no_selection_44_1_ZsZ.root',\n    '${inputFilePrefix}data_no_selection_45_1_Wa5.root',\n    '${inputFilePrefix}data_no_selection_46_1_WW1.root',\n    '${inputFilePrefix}data_no_selection_47_1_ntP.root',\n    '${inputFilePrefix}data_no_selection_48_1_kNi.root',\n    '${inputFilePrefix}data_no_selection_49_1_pBF.root',\n    '${inputFilePrefix}data_no_selection_4_1_nYN.root',\n    '${inputFilePrefix}data_no_selection_50_1_h9K.root',\n    '${inputFilePrefix}data_no_selection_51_1_dET.root',\n    '${inputFilePrefix}data_no_selection_53_1_xKm.root',\n    '${inputFilePrefix}data_no_selection_54_1_3y7.root',\n    '${inputFilePrefix}data_no_selection_55_1_3BE.root',\n    '${inputFilePrefix}data_no_selection_56_1_yQ1.root',\n    '${inputFilePrefix}data_no_selection_57_1_h94.root',\n    '${inputFilePrefix}data_no_selection_58_1_Bx6.root',\n    '${inputFilePrefix}data_no_selection_59_1_35P.root',\n    '${inputFilePrefix}data_no_selection_5_1_99u.root',\n    '${inputFilePrefix}data_no_selection_60_1_XZg.root',\n    '${inputFilePrefix}data_no_selection_61_1_i9c.root',\n    '${inputFilePrefix}data_no_selection_62_1_CEo.root',\n    '${inputFilePrefix}data_no_selection_63_1_fLk.root',\n    '${inputFilePrefix}data_no_selection_65_1_5pe.root',\n    '${inputFilePrefix}data_no_selection_66_1_Wqf.root',\n    '${inputFilePrefix}data_no_selection_67_1_smk.root',\n    '${inputFilePrefix}data_no_selection_68_1_m1G.root',\n    '${inputFilePrefix}data_no_selection_69_1_gMK.root',\n    '${inputFilePrefix}data_no_selection_6_1_sZ4.root',\n    '${inputFilePrefix}data_no_selection_70_1_yNy.root',\n    '${inputFilePrefix}data_no_selection_71_1_JrI.root',\n    '${inputFilePrefix}data_no_selection_72_1_1en.root',\n    '${inputFilePrefix}data_no_selection_75_1_eGZ.root',\n    '${inputFilePrefix}data_no_selection_76_1_Dlj.root',\n    '${inputFilePrefix}data_no_selection_77_1_T1B.root',\n    '${inputFilePrefix}data_no_selection_78_1_RnB.root',\n    '${inputFilePrefix}data_no_selection_79_1_2B6.root',\n    '${inputFilePrefix}data_no_selection_7_1_7vC.root',\n    '${inputFilePrefix}data_no_selection_80_1_8On.root',\n    '${inputFilePrefix}data_no_selection_81_1_wtB.root',\n    '${inputFilePrefix}data_no_selection_82_1_dgj.root',\n    '${inputFilePrefix}data_no_selection_83_1_Ugc.root',\n    '${inputFilePrefix}data_no_selection_84_1_1Li.root',\n    '${inputFilePrefix}data_no_selection_85_1_VIV.root',\n    '${inputFilePrefix}data_no_selection_86_1_eqj.root',\n    '${inputFilePrefix}data_no_selection_87_1_1WY.root',\n    '${inputFilePrefix}data_no_selection_88_1_QTD.root',\n    '${inputFilePrefix}data_no_selection_89_1_35i.root',\n    '${inputFilePrefix}data_no_selection_8_1_FMK.root',\n    '${inputFilePrefix}data_no_selection_90_1_1eY.root',\n    '${inputFilePrefix}data_no_selection_91_1_EqF.root',\n    '${inputFilePrefix}data_no_selection_92_1_GdU.root',\n    '${inputFilePrefix}data_no_selection_93_1_GYP.root',\n    '${inputFilePrefix}data_no_selection_94_1_mbf.root',\n    '${inputFilePrefix}data_no_selection_95_1_Ljq.root',\n    '${inputFilePrefix}data_no_selection_96_1_dHp.root',\n    '${inputFilePrefix}data_no_selection_97_1_ukS.root',\n    '${inputFilePrefix}data_no_selection_98_1_nNc.root',\n    '${inputFilePrefix}data_no_selection_99_1_MMP.root',\n    '${inputFilePrefix}data_no_selection_9_1_NjG.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_WW.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_WW_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_WW_${version}.root" )
nonIsoReweightTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_WW_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_WW_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}WW${infoTag}_${version}.root" )

#samples
samples=( "WW" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_cfg.py

  #generate cfg file for the non-isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_nonIso_cfg.py

  #generate cfg file for the non-isolated, reweighted sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoReweightTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%True%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_nonIsoReweight_cfg.py

  #generate cfg file for the sample with no isolation cut
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.tauAnalysisSequence%" -e "s%PUSCENARIO%S10%" -e "s%REWEIGHT%False%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_all_cfg.py

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
cat <<EOF > runWWTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *WW*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runWWTauAnalyzerCfgs.sh

#generate script that submits all iso+nonIso+reweight jobs to LSF
cat <<EOF > submitWWTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*WW*.sh | grep -v all | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitWWTauAnalyzerJobs.sh

#generate script that submits all noIsoCut jobs to LSF
cat <<EOF > submitWWAllTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*WW*all*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitWWAllTauAnalyzerJobs.sh

#generate script that copies all iso+nonIso+reweight files locally from EOS
cat <<EOF > copyWWFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  for cut in Iso NonIso NonIsoReweight
    do
    if [ "\$cut" != "NonIso" ] || [ $reweightOnly -eq 0 ]
        then
        cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_WW_${version}.root /data1/`whoami`/WW/analysis/
        cmsRm /store/user/`whoami`/muHad\${cut}Analysis_WW_${version}.root
    fi
  done
done

exit 0
EOF
chmod a+x copyWWFromEOS.sh

#generate script that copies all noIsoCut files locally from EOS
cat <<EOF > copyAllWWFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  cmsStage -f /store/user/`whoami`/muHadAnalysis_WW_${version}.root /data1/`whoami`/WW/analysis/
  cmsRm /store/user/`whoami`/muHadAnalysis_WW_${version}.root
done

exit 0
EOF
chmod a+x copyAllWWFromEOS.sh

exit 0
