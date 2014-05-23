#!/bin/bash

if [ $# -gt 3 ]
    then
    echo "Usage: ./generateZZTauAnalyzerCfgs.sh <version> <template cfg> [reweightOnly]"
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
inputFilePrefix="root://eoscms//eos/cms/store/user/yohay/ZZ_TuneZ2star_8TeV_pythia6_tauola-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v2/"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/ZZ/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/ZZ/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_100_3_rXa.root',\n    '${inputFilePrefix}data_no_selection_101_3_yOF.root',\n    '${inputFilePrefix}data_no_selection_102_3_2nB.root',\n    '${inputFilePrefix}data_no_selection_103_3_xAY.root',\n    '${inputFilePrefix}data_no_selection_104_3_3Ek.root',\n    '${inputFilePrefix}data_no_selection_105_3_wIE.root',\n    '${inputFilePrefix}data_no_selection_106_3_qb0.root',\n    '${inputFilePrefix}data_no_selection_107_3_Tvf.root',\n    '${inputFilePrefix}data_no_selection_108_3_tCS.root',\n    '${inputFilePrefix}data_no_selection_109_3_jf1.root',\n    '${inputFilePrefix}data_no_selection_10_1_fX9.root',\n    '${inputFilePrefix}data_no_selection_110_3_v4P.root',\n    '${inputFilePrefix}data_no_selection_111_3_t19.root',\n    '${inputFilePrefix}data_no_selection_112_3_YTA.root',\n    '${inputFilePrefix}data_no_selection_113_3_K7q.root',\n    '${inputFilePrefix}data_no_selection_114_3_KdL.root',\n    '${inputFilePrefix}data_no_selection_115_3_26H.root',\n    '${inputFilePrefix}data_no_selection_116_3_DDJ.root',\n    '${inputFilePrefix}data_no_selection_117_3_2FE.root',\n    '${inputFilePrefix}data_no_selection_118_3_kxZ.root',\n    '${inputFilePrefix}data_no_selection_119_3_0ld.root',\n    '${inputFilePrefix}data_no_selection_11_1_ivB.root',\n    '${inputFilePrefix}data_no_selection_120_4_8NK.root',\n    '${inputFilePrefix}data_no_selection_122_3_AfI.root',\n    '${inputFilePrefix}data_no_selection_123_3_JVw.root',\n    '${inputFilePrefix}data_no_selection_124_3_Hin.root',\n    '${inputFilePrefix}data_no_selection_125_3_Bh2.root',\n    '${inputFilePrefix}data_no_selection_126_3_eiv.root',\n    '${inputFilePrefix}data_no_selection_127_3_C2T.root',\n    '${inputFilePrefix}data_no_selection_128_3_br0.root',\n    '${inputFilePrefix}data_no_selection_129_3_98u.root',\n    '${inputFilePrefix}data_no_selection_12_1_QoV.root',\n    '${inputFilePrefix}data_no_selection_130_3_Ns8.root',\n    '${inputFilePrefix}data_no_selection_131_3_oRw.root',\n    '${inputFilePrefix}data_no_selection_132_4_uY4.root',\n    '${inputFilePrefix}data_no_selection_133_3_4zV.root',\n    '${inputFilePrefix}data_no_selection_134_3_3U0.root',\n    '${inputFilePrefix}data_no_selection_135_3_C2H.root',\n    '${inputFilePrefix}data_no_selection_136_3_n8g.root',\n    '${inputFilePrefix}data_no_selection_137_3_yLd.root',\n    '${inputFilePrefix}data_no_selection_138_3_pkG.root',\n    '${inputFilePrefix}data_no_selection_139_3_loR.root',\n    '${inputFilePrefix}data_no_selection_13_1_Z9p.root',\n    '${inputFilePrefix}data_no_selection_140_3_PvS.root',\n    '${inputFilePrefix}data_no_selection_141_3_xiq.root',\n    '${inputFilePrefix}data_no_selection_142_3_Une.root',\n    '${inputFilePrefix}data_no_selection_143_3_PfL.root',\n    '${inputFilePrefix}data_no_selection_144_3_Nxs.root',\n    '${inputFilePrefix}data_no_selection_145_3_CNu.root',\n    '${inputFilePrefix}data_no_selection_146_1_fbY.root',\n    '${inputFilePrefix}data_no_selection_147_1_sw2.root',\n    '${inputFilePrefix}data_no_selection_148_1_ljI.root',\n    '${inputFilePrefix}data_no_selection_149_2_tyG.root',\n    '${inputFilePrefix}data_no_selection_150_1_llY.root',\n    '${inputFilePrefix}data_no_selection_151_1_H04.root',\n    '${inputFilePrefix}data_no_selection_152_1_ZcO.root',\n    '${inputFilePrefix}data_no_selection_153_1_8i0.root',\n    '${inputFilePrefix}data_no_selection_154_1_0DT.root',\n    '${inputFilePrefix}data_no_selection_155_1_Ifz.root',\n    '${inputFilePrefix}data_no_selection_157_1_nNy.root',\n    '${inputFilePrefix}data_no_selection_158_1_3mL.root',\n    '${inputFilePrefix}data_no_selection_159_1_8PR.root',\n    '${inputFilePrefix}data_no_selection_15_1_bE4.root',\n    '${inputFilePrefix}data_no_selection_160_1_DWV.root',\n    '${inputFilePrefix}data_no_selection_161_1_wWw.root',\n    '${inputFilePrefix}data_no_selection_162_1_Ybc.root',\n    '${inputFilePrefix}data_no_selection_163_2_22Z.root',\n    '${inputFilePrefix}data_no_selection_164_1_L91.root',\n    '${inputFilePrefix}data_no_selection_165_2_QYm.root',\n    '${inputFilePrefix}data_no_selection_166_1_MgC.root',\n    '${inputFilePrefix}data_no_selection_167_2_aUe.root',\n    '${inputFilePrefix}data_no_selection_168_1_Bfi.root',\n    '${inputFilePrefix}data_no_selection_169_1_pXn.root',\n    '${inputFilePrefix}data_no_selection_16_1_3M5.root',\n    '${inputFilePrefix}data_no_selection_170_1_9Xl.root',\n    '${inputFilePrefix}data_no_selection_171_1_KXv.root',\n    '${inputFilePrefix}data_no_selection_172_1_kIM.root',\n    '${inputFilePrefix}data_no_selection_173_2_HFl.root',\n    '${inputFilePrefix}data_no_selection_174_1_NSg.root',\n    '${inputFilePrefix}data_no_selection_175_1_Qc1.root',\n    '${inputFilePrefix}data_no_selection_176_1_34n.root',\n    '${inputFilePrefix}data_no_selection_177_1_1cB.root',\n    '${inputFilePrefix}data_no_selection_178_2_SEj.root',\n    '${inputFilePrefix}data_no_selection_179_1_jIi.root',\n    '${inputFilePrefix}data_no_selection_17_2_FrN.root',\n    '${inputFilePrefix}data_no_selection_180_1_gtW.root',\n    '${inputFilePrefix}data_no_selection_181_1_jo1.root',\n    '${inputFilePrefix}data_no_selection_182_1_WEA.root',\n    '${inputFilePrefix}data_no_selection_183_1_Ajk.root',\n    '${inputFilePrefix}data_no_selection_184_1_LKy.root',\n    '${inputFilePrefix}data_no_selection_185_1_bCI.root',\n    '${inputFilePrefix}data_no_selection_186_1_sSX.root',\n    '${inputFilePrefix}data_no_selection_187_1_nam.root',\n    '${inputFilePrefix}data_no_selection_188_1_D0x.root',\n    '${inputFilePrefix}data_no_selection_189_2_9vj.root',\n    '${inputFilePrefix}data_no_selection_18_1_K6a.root',\n    '${inputFilePrefix}data_no_selection_190_1_pJz.root',\n    '${inputFilePrefix}data_no_selection_191_1_a0d.root',\n    '${inputFilePrefix}data_no_selection_192_1_8ob.root',\n    '${inputFilePrefix}data_no_selection_193_1_RGQ.root',\n    '${inputFilePrefix}data_no_selection_194_1_0SH.root',\n    '${inputFilePrefix}data_no_selection_195_1_UhG.root',\n    '${inputFilePrefix}data_no_selection_196_1_NGV.root',\n    '${inputFilePrefix}data_no_selection_197_1_qLA.root',\n    '${inputFilePrefix}data_no_selection_198_1_jHy.root',\n    '${inputFilePrefix}data_no_selection_199_1_vDG.root',\n    '${inputFilePrefix}data_no_selection_19_1_xuY.root',\n    '${inputFilePrefix}data_no_selection_1_1_ubM.root',\n    '${inputFilePrefix}data_no_selection_200_1_Yxe.root',\n    '${inputFilePrefix}data_no_selection_201_1_TH8.root',\n    '${inputFilePrefix}data_no_selection_202_1_i4B.root',\n    '${inputFilePrefix}data_no_selection_20_1_Dm0.root',\n    '${inputFilePrefix}data_no_selection_21_1_1uS.root',\n    '${inputFilePrefix}data_no_selection_22_1_7Jz.root',\n    '${inputFilePrefix}data_no_selection_23_1_Nct.root',\n    '${inputFilePrefix}data_no_selection_24_1_fwG.root',\n    '${inputFilePrefix}data_no_selection_25_1_xDK.root',\n    '${inputFilePrefix}data_no_selection_26_1_sK0.root',\n    '${inputFilePrefix}data_no_selection_27_1_Qrq.root',\n    '${inputFilePrefix}data_no_selection_28_1_Ojn.root',\n    '${inputFilePrefix}data_no_selection_29_1_RLB.root',\n    '${inputFilePrefix}data_no_selection_2_1_WWT.root',\n    '${inputFilePrefix}data_no_selection_30_1_Yib.root',\n    '${inputFilePrefix}data_no_selection_32_1_j5t.root',\n    '${inputFilePrefix}data_no_selection_33_2_iHh.root',\n    '${inputFilePrefix}data_no_selection_34_1_Z2x.root',\n    '${inputFilePrefix}data_no_selection_35_1_tOr.root',\n    '${inputFilePrefix}data_no_selection_36_1_bWt.root',\n    '${inputFilePrefix}data_no_selection_37_1_UWh.root',\n    '${inputFilePrefix}data_no_selection_38_1_7ck.root',\n    '${inputFilePrefix}data_no_selection_39_1_nlK.root',\n    '${inputFilePrefix}data_no_selection_3_1_5wE.root',\n    '${inputFilePrefix}data_no_selection_40_1_tyV.root',\n    '${inputFilePrefix}data_no_selection_42_1_30L.root',\n    '${inputFilePrefix}data_no_selection_43_1_NJQ.root',\n    '${inputFilePrefix}data_no_selection_44_1_2oa.root',\n    '${inputFilePrefix}data_no_selection_45_1_DIp.root',\n    '${inputFilePrefix}data_no_selection_46_1_eae.root',\n    '${inputFilePrefix}data_no_selection_47_1_Dfg.root',\n    '${inputFilePrefix}data_no_selection_48_1_idF.root',\n    '${inputFilePrefix}data_no_selection_49_1_nCP.root',\n    '${inputFilePrefix}data_no_selection_4_1_HVR.root',\n    '${inputFilePrefix}data_no_selection_50_1_BOQ.root',\n    '${inputFilePrefix}data_no_selection_51_3_ci6.root',\n    '${inputFilePrefix}data_no_selection_53_3_Cw0.root',\n    '${inputFilePrefix}data_no_selection_54_3_0sJ.root',\n    '${inputFilePrefix}data_no_selection_55_3_tb8.root',\n    '${inputFilePrefix}data_no_selection_56_3_tLb.root',\n    '${inputFilePrefix}data_no_selection_57_3_cjv.root',\n    '${inputFilePrefix}data_no_selection_59_3_woO.root',\n    '${inputFilePrefix}data_no_selection_5_1_OPd.root',\n    '${inputFilePrefix}data_no_selection_60_3_s9b.root',\n    '${inputFilePrefix}data_no_selection_61_3_lqs.root',\n    '${inputFilePrefix}data_no_selection_62_3_T1m.root',\n    '${inputFilePrefix}data_no_selection_63_3_h6i.root',\n    '${inputFilePrefix}data_no_selection_64_3_XAM.root',\n    '${inputFilePrefix}data_no_selection_65_4_N4L.root',\n    '${inputFilePrefix}data_no_selection_66_3_orq.root',\n    '${inputFilePrefix}data_no_selection_67_4_fzm.root',\n    '${inputFilePrefix}data_no_selection_68_4_Aro.root',\n    '${inputFilePrefix}data_no_selection_6_1_VM8.root',\n    '${inputFilePrefix}data_no_selection_70_3_lGL.root',\n    '${inputFilePrefix}data_no_selection_73_3_30M.root',\n    '${inputFilePrefix}data_no_selection_74_3_68N.root',\n    '${inputFilePrefix}data_no_selection_75_3_Snw.root',\n    '${inputFilePrefix}data_no_selection_76_3_7vp.root',\n    '${inputFilePrefix}data_no_selection_77_3_Yop.root',\n    '${inputFilePrefix}data_no_selection_78_3_K3m.root',\n    '${inputFilePrefix}data_no_selection_79_3_uu3.root',\n    '${inputFilePrefix}data_no_selection_7_1_CsI.root',\n    '${inputFilePrefix}data_no_selection_80_3_Dac.root',\n    '${inputFilePrefix}data_no_selection_82_3_RxY.root',\n    '${inputFilePrefix}data_no_selection_83_3_0yu.root',\n    '${inputFilePrefix}data_no_selection_84_3_lMX.root',\n    '${inputFilePrefix}data_no_selection_85_3_Hqy.root',\n    '${inputFilePrefix}data_no_selection_86_3_TvI.root',\n    '${inputFilePrefix}data_no_selection_87_3_ynH.root',\n    '${inputFilePrefix}data_no_selection_88_3_lCl.root',\n    '${inputFilePrefix}data_no_selection_89_3_OtI.root',\n    '${inputFilePrefix}data_no_selection_8_1_2Sb.root',\n    '${inputFilePrefix}data_no_selection_90_3_nKo.root',\n    '${inputFilePrefix}data_no_selection_91_3_vVw.root',\n    '${inputFilePrefix}data_no_selection_92_3_Fcv.root',\n    '${inputFilePrefix}data_no_selection_93_3_HiE.root',\n    '${inputFilePrefix}data_no_selection_95_3_Lu8.root',\n    '${inputFilePrefix}data_no_selection_96_3_ntj.root',\n    '${inputFilePrefix}data_no_selection_97_3_Kim.root',\n    '${inputFilePrefix}data_no_selection_98_3_nPT.root',\n    '${inputFilePrefix}data_no_selection_99_3_6vg.root',\n    '${inputFilePrefix}data_no_selection_9_1_E1i.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_ZZ.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_ZZ_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_ZZ_${version}.root" )
nonIsoReweightTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_ZZ_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_ZZ_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}ZZ${infoTag}_${version}.root" )

#samples
samples=( "ZZ" )

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
cat <<EOF > runZZTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *ZZ*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runZZTauAnalyzerCfgs.sh

#generate script that submits all iso+nonIso+reweight jobs to LSF
cat <<EOF > submitZZTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*ZZ*.sh | grep -v all | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitZZTauAnalyzerJobs.sh

#generate script that submits all noIsoCut jobs to LSF
cat <<EOF > submitZZAllTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*ZZ*all*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitZZAllTauAnalyzerJobs.sh

#generate script that copies all iso+nonIso+reweight files locally from EOS
cat <<EOF > copyZZFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  #for cut in Iso NonIso NonIsoReweight
  for cut in Iso NonIso
    do
    if [ "\$cut" != "NonIso" ] || [ $reweightOnly -eq 0 ]
        then
        cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_ZZ_${version}.root /data1/`whoami`/ZZ/analysis/
        cmsRm /store/user/`whoami`/muHad\${cut}Analysis_ZZ_${version}.root
    fi
  done
done

exit 0
EOF
chmod a+x copyZZFromEOS.sh

#generate script that copies all noIsoCut files locally from EOS
cat <<EOF > copyAllZZFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  cmsStage -f /store/user/`whoami`/muHadAnalysis_ZZ_${version}.root /data1/`whoami`/ZZ/analysis/
  cmsRm /store/user/`whoami`/muHadAnalysis_ZZ_${version}.root
done

exit 0
EOF
chmod a+x copyAllZZFromEOS.sh

exit 0
