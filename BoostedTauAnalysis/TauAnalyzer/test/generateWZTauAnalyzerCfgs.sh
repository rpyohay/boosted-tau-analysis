#!/bin/bash

####STUFF TO CONFIGURE####

#version
<<<<<<< HEAD
version="v9"
=======
version="v40"
>>>>>>> 46247d24eb20eea62a68284a461b2dc6bfa58c65
infoTag=""
dir=$version

#number of samples
nSamples=1
iBeg=0
iEnd=`expr $nSamples - 1`

#input file prefix
inputFilePrefix="root://eoscms//eos/cms/store/user/yohay/WZ_TuneZ2star_8TeV_pythia6_tauola-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v2/"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/WZ/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/WZ/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_100_1_fTn.root',\n    '${inputFilePrefix}data_no_selection_101_1_DTw.root',\n    '${inputFilePrefix}data_no_selection_102_1_f6J.root',\n    '${inputFilePrefix}data_no_selection_103_1_e36.root',\n    '${inputFilePrefix}data_no_selection_104_1_CTU.root',\n    '${inputFilePrefix}data_no_selection_105_1_rBR.root',\n    '${inputFilePrefix}data_no_selection_106_1_K0h.root',\n    '${inputFilePrefix}data_no_selection_107_1_Wme.root',\n    '${inputFilePrefix}data_no_selection_108_1_gFh.root',\n    '${inputFilePrefix}data_no_selection_109_1_m4Y.root',\n    '${inputFilePrefix}data_no_selection_10_2_Isy.root',\n    '${inputFilePrefix}data_no_selection_110_1_W9K.root',\n    '${inputFilePrefix}data_no_selection_111_1_8WO.root',\n    '${inputFilePrefix}data_no_selection_112_1_wrg.root',\n    '${inputFilePrefix}data_no_selection_113_1_yRR.root',\n    '${inputFilePrefix}data_no_selection_114_1_xt5.root',\n    '${inputFilePrefix}data_no_selection_115_1_Sc4.root',\n    '${inputFilePrefix}data_no_selection_116_1_ss9.root',\n    '${inputFilePrefix}data_no_selection_117_1_HBr.root',\n    '${inputFilePrefix}data_no_selection_118_1_wZv.root',\n    '${inputFilePrefix}data_no_selection_119_1_mVR.root',\n    '${inputFilePrefix}data_no_selection_11_2_Y4c.root',\n    '${inputFilePrefix}data_no_selection_120_1_LbQ.root',\n    '${inputFilePrefix}data_no_selection_121_1_Gba.root',\n    '${inputFilePrefix}data_no_selection_122_1_1Mx.root',\n    '${inputFilePrefix}data_no_selection_123_1_3PZ.root',\n    '${inputFilePrefix}data_no_selection_124_1_TD3.root',\n    '${inputFilePrefix}data_no_selection_125_1_Uac.root',\n    '${inputFilePrefix}data_no_selection_126_1_Ffk.root',\n    '${inputFilePrefix}data_no_selection_127_1_Hdh.root',\n    '${inputFilePrefix}data_no_selection_128_1_6RP.root',\n    '${inputFilePrefix}data_no_selection_129_1_grd.root',\n    '${inputFilePrefix}data_no_selection_12_2_cF0.root',\n    '${inputFilePrefix}data_no_selection_130_1_VeN.root',\n    '${inputFilePrefix}data_no_selection_131_1_QZe.root',\n    '${inputFilePrefix}data_no_selection_132_1_RLs.root',\n    '${inputFilePrefix}data_no_selection_133_1_FJP.root',\n    '${inputFilePrefix}data_no_selection_134_1_k3i.root',\n    '${inputFilePrefix}data_no_selection_135_2_HBI.root',\n    '${inputFilePrefix}data_no_selection_136_1_Gd2.root',\n    '${inputFilePrefix}data_no_selection_137_1_tx2.root',\n    '${inputFilePrefix}data_no_selection_138_1_cq6.root',\n    '${inputFilePrefix}data_no_selection_139_1_Js8.root',\n    '${inputFilePrefix}data_no_selection_13_2_qtG.root',\n    '${inputFilePrefix}data_no_selection_140_1_8H2.root',\n    '${inputFilePrefix}data_no_selection_141_1_aJt.root',\n    '${inputFilePrefix}data_no_selection_142_1_NO4.root',\n    '${inputFilePrefix}data_no_selection_143_1_uEn.root',\n    '${inputFilePrefix}data_no_selection_144_1_g6n.root',\n    '${inputFilePrefix}data_no_selection_145_1_5rX.root',\n    '${inputFilePrefix}data_no_selection_146_1_yRQ.root',\n    '${inputFilePrefix}data_no_selection_147_1_SKV.root',\n    '${inputFilePrefix}data_no_selection_148_1_Xgx.root',\n    '${inputFilePrefix}data_no_selection_149_1_Oyy.root',\n    '${inputFilePrefix}data_no_selection_14_2_a5C.root',\n    '${inputFilePrefix}data_no_selection_150_1_cRa.root',\n    '${inputFilePrefix}data_no_selection_151_1_4i9.root',\n    '${inputFilePrefix}data_no_selection_152_1_fDH.root',\n    '${inputFilePrefix}data_no_selection_153_1_dPI.root',\n    '${inputFilePrefix}data_no_selection_154_1_3lr.root',\n    '${inputFilePrefix}data_no_selection_155_1_ASk.root',\n    '${inputFilePrefix}data_no_selection_156_1_5nM.root',\n    '${inputFilePrefix}data_no_selection_157_1_mIu.root',\n    '${inputFilePrefix}data_no_selection_159_1_jxT.root',\n    '${inputFilePrefix}data_no_selection_15_2_bxN.root',\n    '${inputFilePrefix}data_no_selection_160_1_Xio.root',\n    '${inputFilePrefix}data_no_selection_161_1_lBj.root',\n    '${inputFilePrefix}data_no_selection_162_1_l6A.root',\n    '${inputFilePrefix}data_no_selection_163_1_H2a.root',\n    '${inputFilePrefix}data_no_selection_164_1_HGL.root',\n    '${inputFilePrefix}data_no_selection_165_1_JCX.root',\n    '${inputFilePrefix}data_no_selection_166_1_12t.root',\n    '${inputFilePrefix}data_no_selection_167_1_wSB.root',\n    '${inputFilePrefix}data_no_selection_168_1_nV4.root',\n    '${inputFilePrefix}data_no_selection_169_1_ISp.root',\n    '${inputFilePrefix}data_no_selection_16_2_fhC.root',\n    '${inputFilePrefix}data_no_selection_170_2_5OC.root',\n    '${inputFilePrefix}data_no_selection_171_1_hls.root',\n    '${inputFilePrefix}data_no_selection_172_1_9my.root',\n    '${inputFilePrefix}data_no_selection_173_1_OHs.root',\n    '${inputFilePrefix}data_no_selection_174_1_kFI.root',\n    '${inputFilePrefix}data_no_selection_175_1_sOM.root',\n    '${inputFilePrefix}data_no_selection_176_1_kJ8.root',\n    '${inputFilePrefix}data_no_selection_177_1_o5c.root',\n    '${inputFilePrefix}data_no_selection_178_1_XQR.root',\n    '${inputFilePrefix}data_no_selection_179_2_by4.root',\n    '${inputFilePrefix}data_no_selection_17_2_dug.root',\n    '${inputFilePrefix}data_no_selection_180_1_Oow.root',\n    '${inputFilePrefix}data_no_selection_181_1_6UA.root',\n    '${inputFilePrefix}data_no_selection_182_1_3iZ.root',\n    '${inputFilePrefix}data_no_selection_183_1_8Ja.root',\n    '${inputFilePrefix}data_no_selection_185_1_QiQ.root',\n    '${inputFilePrefix}data_no_selection_186_1_RiW.root',\n    '${inputFilePrefix}data_no_selection_187_1_YDe.root',\n    '${inputFilePrefix}data_no_selection_188_1_7sW.root',\n    '${inputFilePrefix}data_no_selection_189_1_gVc.root',\n    '${inputFilePrefix}data_no_selection_18_2_RJw.root',\n    '${inputFilePrefix}data_no_selection_190_1_OMt.root',\n    '${inputFilePrefix}data_no_selection_192_1_ukq.root',\n    '${inputFilePrefix}data_no_selection_193_1_Jkk.root',\n    '${inputFilePrefix}data_no_selection_194_1_CQ7.root',\n    '${inputFilePrefix}data_no_selection_195_1_Ztq.root',\n    '${inputFilePrefix}data_no_selection_196_1_Zm8.root',\n    '${inputFilePrefix}data_no_selection_197_1_DNF.root',\n    '${inputFilePrefix}data_no_selection_198_1_e7j.root',\n    '${inputFilePrefix}data_no_selection_199_1_Knu.root',\n    '${inputFilePrefix}data_no_selection_19_2_Dna.root',\n    '${inputFilePrefix}data_no_selection_1_2_Vjl.root',\n    '${inputFilePrefix}data_no_selection_200_1_ArY.root',\n    '${inputFilePrefix}data_no_selection_201_1_Oul.root',\n    '${inputFilePrefix}data_no_selection_202_1_74n.root',\n    '${inputFilePrefix}data_no_selection_203_1_Cg5.root',\n    '${inputFilePrefix}data_no_selection_20_2_FNQ.root',\n    '${inputFilePrefix}data_no_selection_21_2_BYs.root',\n    '${inputFilePrefix}data_no_selection_22_2_p3E.root',\n    '${inputFilePrefix}data_no_selection_23_2_B5v.root',\n    '${inputFilePrefix}data_no_selection_24_2_wV3.root',\n    '${inputFilePrefix}data_no_selection_25_3_mqo.root',\n    '${inputFilePrefix}data_no_selection_26_3_xzw.root',\n    '${inputFilePrefix}data_no_selection_27_2_SH4.root',\n    '${inputFilePrefix}data_no_selection_28_2_V3A.root',\n    '${inputFilePrefix}data_no_selection_29_2_52C.root',\n    '${inputFilePrefix}data_no_selection_2_3_tHZ.root',\n    '${inputFilePrefix}data_no_selection_30_2_Rzn.root',\n    '${inputFilePrefix}data_no_selection_31_2_N1C.root',\n    '${inputFilePrefix}data_no_selection_32_2_2vL.root',\n    '${inputFilePrefix}data_no_selection_33_3_rsA.root',\n    '${inputFilePrefix}data_no_selection_34_2_nTv.root',\n    '${inputFilePrefix}data_no_selection_35_2_JAQ.root',\n    '${inputFilePrefix}data_no_selection_36_2_Zc0.root',\n    '${inputFilePrefix}data_no_selection_37_2_5h0.root',\n    '${inputFilePrefix}data_no_selection_38_2_Vqr.root',\n    '${inputFilePrefix}data_no_selection_39_2_OeD.root',\n    '${inputFilePrefix}data_no_selection_3_3_xZC.root',\n    '${inputFilePrefix}data_no_selection_40_2_zpQ.root',\n    '${inputFilePrefix}data_no_selection_41_2_wH7.root',\n    '${inputFilePrefix}data_no_selection_42_2_K8h.root',\n    '${inputFilePrefix}data_no_selection_43_2_Uye.root',\n    '${inputFilePrefix}data_no_selection_44_2_TTg.root',\n    '${inputFilePrefix}data_no_selection_45_2_7Rj.root',\n    '${inputFilePrefix}data_no_selection_46_2_Pti.root',\n    '${inputFilePrefix}data_no_selection_47_2_FhV.root',\n    '${inputFilePrefix}data_no_selection_48_2_2MJ.root',\n    '${inputFilePrefix}data_no_selection_49_2_BC2.root',\n    '${inputFilePrefix}data_no_selection_4_2_0AQ.root',\n    '${inputFilePrefix}data_no_selection_50_2_M36.root',\n    '${inputFilePrefix}data_no_selection_51_2_9tk.root',\n    '${inputFilePrefix}data_no_selection_52_2_6m8.root',\n    '${inputFilePrefix}data_no_selection_53_2_ZFZ.root',\n    '${inputFilePrefix}data_no_selection_54_2_EHi.root',\n    '${inputFilePrefix}data_no_selection_55_2_rtZ.root',\n    '${inputFilePrefix}data_no_selection_57_2_Poh.root',\n    '${inputFilePrefix}data_no_selection_58_2_h3v.root',\n    '${inputFilePrefix}data_no_selection_59_3_4u0.root',\n    '${inputFilePrefix}data_no_selection_5_3_cx6.root',\n    '${inputFilePrefix}data_no_selection_60_2_m9U.root',\n    '${inputFilePrefix}data_no_selection_61_2_ZUu.root',\n    '${inputFilePrefix}data_no_selection_62_2_zHI.root',\n    '${inputFilePrefix}data_no_selection_63_2_EnF.root',\n    '${inputFilePrefix}data_no_selection_64_2_TAS.root',\n    '${inputFilePrefix}data_no_selection_65_2_BfX.root',\n    '${inputFilePrefix}data_no_selection_66_2_tnf.root',\n    '${inputFilePrefix}data_no_selection_67_2_hL3.root',\n    '${inputFilePrefix}data_no_selection_68_2_kO3.root',\n    '${inputFilePrefix}data_no_selection_69_2_b3I.root',\n    '${inputFilePrefix}data_no_selection_6_2_qPA.root',\n    '${inputFilePrefix}data_no_selection_70_2_shf.root',\n    '${inputFilePrefix}data_no_selection_71_2_1zB.root',\n    '${inputFilePrefix}data_no_selection_72_2_Taq.root',\n    '${inputFilePrefix}data_no_selection_73_2_zKx.root',\n    '${inputFilePrefix}data_no_selection_74_2_0wy.root',\n    '${inputFilePrefix}data_no_selection_75_2_Xkj.root',\n    '${inputFilePrefix}data_no_selection_76_2_wAW.root',\n    '${inputFilePrefix}data_no_selection_77_2_Laq.root',\n    '${inputFilePrefix}data_no_selection_78_2_U3F.root',\n    '${inputFilePrefix}data_no_selection_79_2_SDq.root',\n    '${inputFilePrefix}data_no_selection_7_2_ZT4.root',\n    '${inputFilePrefix}data_no_selection_80_2_yKT.root',\n    '${inputFilePrefix}data_no_selection_81_2_mNv.root',\n    '${inputFilePrefix}data_no_selection_82_2_2cr.root',\n    '${inputFilePrefix}data_no_selection_83_2_m1Q.root',\n    '${inputFilePrefix}data_no_selection_84_2_OIt.root',\n    '${inputFilePrefix}data_no_selection_85_1_tRf.root',\n    '${inputFilePrefix}data_no_selection_86_1_P0I.root',\n    '${inputFilePrefix}data_no_selection_87_1_bKr.root',\n    '${inputFilePrefix}data_no_selection_88_1_6iF.root',\n    '${inputFilePrefix}data_no_selection_89_1_bzH.root',\n    '${inputFilePrefix}data_no_selection_8_2_u35.root',\n    '${inputFilePrefix}data_no_selection_90_1_bcc.root',\n    '${inputFilePrefix}data_no_selection_91_1_Cza.root',\n    '${inputFilePrefix}data_no_selection_92_1_NOu.root',\n    '${inputFilePrefix}data_no_selection_93_1_8Wt.root',\n    '${inputFilePrefix}data_no_selection_94_1_cuq.root',\n    '${inputFilePrefix}data_no_selection_95_1_R3D.root',\n    '${inputFilePrefix}data_no_selection_96_1_Bdi.root',\n    '${inputFilePrefix}data_no_selection_97_1_KWp.root',\n    '${inputFilePrefix}data_no_selection_98_1_rkB.root',\n    '${inputFilePrefix}data_no_selection_99_1_p4P.root',\n    '${inputFilePrefix}data_no_selection_9_2_LtI.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_WZ.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_WZ_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_WZ_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_WZ_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}WZ${infoTag}_${version}.root" )

#samples
samples=( "WZ" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" ../tauanalyzer_WNJetsToLNu_Wh1_template_cfg.py > tauanalyzer_${samples[${i}]}_iso_cfg.py

  #generate cfg file for the non-isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" ../tauanalyzer_WNJetsToLNu_Wh1_template_cfg.py > tauanalyzer_${samples[${i}]}_nonIso_cfg.py

  #generate job submission script for LSF
  cat <<EOF > tauanalyzer_${samples[${i}]}_cfg.sh
#!/bin/bash

jobDir="`pwd`"
fileNamePrefix="tauanalyzer_${samples[${i}]}"

cd \$jobDir
eval \`scramv1 runtime -sh\`
cd -
cp \$jobDir/\${fileNamePrefix}_iso_cfg.py \$jobDir/\${fileNamePrefix}_nonIso_cfg.py .
cmsRun \${fileNamePrefix}_iso_cfg.py
cmsRun \${fileNamePrefix}_nonIso_cfg.py
cmsStage -f ${isoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${nonIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
rm ${isoTauAnalyzerOutputFiles[${i}]} ${nonIsoTauAnalyzerOutputFiles[${i}]} 

exit 0
EOF
  chmod a+x tauanalyzer_${samples[${i}]}_cfg.sh
done

#generate run cfg that runs all files in the directory
cat <<EOF > runWZTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *WZ*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runWZTauAnalyzerCfgs.sh

#generate script that submits all jobs to LSF
cat <<EOF > submitWZTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*WZ*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitWZTauAnalyzerJobs.sh

#generate script that copies all files locally from EOS
cat <<EOF > copyWZFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  for cut in Iso NonIso
    do
    cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_WZ_${version}.root /data1/`whoami`/WZ/analysis/
  done
done

exit 0
EOF
chmod a+x copyWZFromEOS.sh

exit 0
