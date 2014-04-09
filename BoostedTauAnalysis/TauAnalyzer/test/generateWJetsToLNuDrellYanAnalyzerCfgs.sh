#!/bin/bash

if [ $# -ne 2 ]
    then
    echo "Usage: ./generateWJetsToLNuDrellYanAnalyzerCfgs.sh <version> <template cfg>"
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
inputFilePrefix="root://eoscms//eos/cms/store/user/yohay/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_DYEnriched/"

#DrellYanAnalyzer output file prefix
#DrellYanAnalyzerOutputFilePrefix="/data1/yohay/WJetsToLNu/analysis/"
DrellYanAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/WJetsToLNu/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}data_no_selection_100_1_BQc.root',\n    '${inputFilePrefix}data_no_selection_101_1_bjW.root',\n    '${inputFilePrefix}data_no_selection_102_1_7tO.root',\n    '${inputFilePrefix}data_no_selection_103_1_u9K.root',\n    '${inputFilePrefix}data_no_selection_104_1_aWz.root',\n    '${inputFilePrefix}data_no_selection_105_1_xHV.root',\n    '${inputFilePrefix}data_no_selection_106_1_xlf.root',\n    '${inputFilePrefix}data_no_selection_107_1_wZn.root',\n    '${inputFilePrefix}data_no_selection_108_1_6Vy.root',\n    '${inputFilePrefix}data_no_selection_109_1_h35.root',\n    '${inputFilePrefix}data_no_selection_10_1_wCA.root',\n    '${inputFilePrefix}data_no_selection_110_1_jyn.root',\n    '${inputFilePrefix}data_no_selection_111_1_56i.root',\n    '${inputFilePrefix}data_no_selection_112_1_iOi.root',\n    '${inputFilePrefix}data_no_selection_113_1_PRJ.root',\n    '${inputFilePrefix}data_no_selection_114_1_rEX.root',\n    '${inputFilePrefix}data_no_selection_115_1_s7B.root',\n    '${inputFilePrefix}data_no_selection_116_1_TWd.root',\n    '${inputFilePrefix}data_no_selection_117_1_S4c.root',\n    '${inputFilePrefix}data_no_selection_118_1_ZKe.root',\n    '${inputFilePrefix}data_no_selection_119_1_LyO.root',\n    '${inputFilePrefix}data_no_selection_11_1_SzQ.root',\n    '${inputFilePrefix}data_no_selection_120_1_Wcv.root',\n    '${inputFilePrefix}data_no_selection_121_1_6kT.root',\n    '${inputFilePrefix}data_no_selection_122_1_Tz0.root',\n    '${inputFilePrefix}data_no_selection_123_1_BjO.root',\n    '${inputFilePrefix}data_no_selection_124_1_QSE.root',\n    '${inputFilePrefix}data_no_selection_125_1_p4p.root',\n    '${inputFilePrefix}data_no_selection_126_1_DW8.root',\n    '${inputFilePrefix}data_no_selection_127_1_y7X.root',\n    '${inputFilePrefix}data_no_selection_128_1_QYd.root',\n    '${inputFilePrefix}data_no_selection_129_1_Trn.root',\n    '${inputFilePrefix}data_no_selection_12_1_EPK.root',\n    '${inputFilePrefix}data_no_selection_130_1_rL6.root',\n    '${inputFilePrefix}data_no_selection_131_1_Mws.root',\n    '${inputFilePrefix}data_no_selection_132_1_H8w.root',\n    '${inputFilePrefix}data_no_selection_133_1_whb.root',\n    '${inputFilePrefix}data_no_selection_134_1_maJ.root',\n    '${inputFilePrefix}data_no_selection_135_1_xyR.root',\n    '${inputFilePrefix}data_no_selection_136_1_BPM.root',\n    '${inputFilePrefix}data_no_selection_137_1_pst.root',\n    '${inputFilePrefix}data_no_selection_138_1_tkj.root',\n    '${inputFilePrefix}data_no_selection_139_1_rvI.root',\n    '${inputFilePrefix}data_no_selection_13_1_Y5n.root',\n    '${inputFilePrefix}data_no_selection_140_1_jBc.root',\n    '${inputFilePrefix}data_no_selection_141_1_wRw.root',\n    '${inputFilePrefix}data_no_selection_142_1_OmC.root',\n    '${inputFilePrefix}data_no_selection_143_1_0BU.root',\n    '${inputFilePrefix}data_no_selection_144_1_IAq.root',\n    '${inputFilePrefix}data_no_selection_145_1_JJy.root',\n    '${inputFilePrefix}data_no_selection_146_1_shX.root',\n    '${inputFilePrefix}data_no_selection_147_1_Jkz.root',\n    '${inputFilePrefix}data_no_selection_148_1_Am8.root',\n    '${inputFilePrefix}data_no_selection_149_1_j4O.root',\n    '${inputFilePrefix}data_no_selection_14_1_azh.root',\n    '${inputFilePrefix}data_no_selection_150_1_Wtf.root',\n    '${inputFilePrefix}data_no_selection_151_1_0nK.root',\n    '${inputFilePrefix}data_no_selection_152_1_Izo.root',\n    '${inputFilePrefix}data_no_selection_153_1_9wA.root',\n    '${inputFilePrefix}data_no_selection_154_1_kYW.root',\n    '${inputFilePrefix}data_no_selection_155_1_nVD.root',\n    '${inputFilePrefix}data_no_selection_156_1_uSh.root',\n    '${inputFilePrefix}data_no_selection_157_1_vUl.root',\n    '${inputFilePrefix}data_no_selection_158_1_dVg.root',\n    '${inputFilePrefix}data_no_selection_159_1_NdY.root',\n    '${inputFilePrefix}data_no_selection_15_1_AOM.root',\n    '${inputFilePrefix}data_no_selection_160_1_TPJ.root',\n    '${inputFilePrefix}data_no_selection_161_1_76W.root',\n    '${inputFilePrefix}data_no_selection_162_1_Gxu.root',\n    '${inputFilePrefix}data_no_selection_163_1_WmR.root',\n    '${inputFilePrefix}data_no_selection_164_1_6iJ.root',\n    '${inputFilePrefix}data_no_selection_165_1_ARz.root',\n    '${inputFilePrefix}data_no_selection_166_1_vEX.root',\n    '${inputFilePrefix}data_no_selection_167_1_dL0.root',\n    '${inputFilePrefix}data_no_selection_168_1_gDR.root',\n    '${inputFilePrefix}data_no_selection_169_1_u0s.root',\n    '${inputFilePrefix}data_no_selection_16_1_wfl.root',\n    '${inputFilePrefix}data_no_selection_170_1_0qj.root',\n    '${inputFilePrefix}data_no_selection_171_1_FcK.root',\n    '${inputFilePrefix}data_no_selection_172_1_0vv.root',\n    '${inputFilePrefix}data_no_selection_173_1_fZL.root',\n    '${inputFilePrefix}data_no_selection_174_1_aBK.root',\n    '${inputFilePrefix}data_no_selection_175_1_Utc.root',\n    '${inputFilePrefix}data_no_selection_176_1_7cz.root',\n    '${inputFilePrefix}data_no_selection_177_1_XXr.root',\n    '${inputFilePrefix}data_no_selection_178_1_aHd.root',\n    '${inputFilePrefix}data_no_selection_179_1_tf7.root',\n    '${inputFilePrefix}data_no_selection_17_1_Jfe.root',\n    '${inputFilePrefix}data_no_selection_180_1_Tl8.root',\n    '${inputFilePrefix}data_no_selection_181_1_1w2.root',\n    '${inputFilePrefix}data_no_selection_182_1_PnM.root',\n    '${inputFilePrefix}data_no_selection_183_1_yL5.root',\n    '${inputFilePrefix}data_no_selection_184_1_Trz.root',\n    '${inputFilePrefix}data_no_selection_185_1_uDT.root',\n    '${inputFilePrefix}data_no_selection_186_1_xFQ.root',\n    '${inputFilePrefix}data_no_selection_187_1_p77.root',\n    '${inputFilePrefix}data_no_selection_188_1_PqU.root',\n    '${inputFilePrefix}data_no_selection_18_1_74M.root',\n    '${inputFilePrefix}data_no_selection_19_1_I79.root',\n    '${inputFilePrefix}data_no_selection_1_1_Qnx.root',\n    '${inputFilePrefix}data_no_selection_20_1_jx8.root',\n    '${inputFilePrefix}data_no_selection_21_1_KyV.root',\n    '${inputFilePrefix}data_no_selection_22_1_Knm.root',\n    '${inputFilePrefix}data_no_selection_23_1_oPf.root',\n    '${inputFilePrefix}data_no_selection_24_1_atk.root',\n    '${inputFilePrefix}data_no_selection_25_1_1vP.root',\n    '${inputFilePrefix}data_no_selection_26_1_NbJ.root',\n    '${inputFilePrefix}data_no_selection_27_1_gze.root',\n    '${inputFilePrefix}data_no_selection_28_1_NXU.root',\n    '${inputFilePrefix}data_no_selection_29_1_NLc.root',\n    '${inputFilePrefix}data_no_selection_2_1_rvY.root',\n    '${inputFilePrefix}data_no_selection_30_1_c55.root',\n    '${inputFilePrefix}data_no_selection_31_1_cJk.root',\n    '${inputFilePrefix}data_no_selection_32_1_jB4.root',\n    '${inputFilePrefix}data_no_selection_33_1_byh.root',\n    '${inputFilePrefix}data_no_selection_34_1_x7D.root',\n    '${inputFilePrefix}data_no_selection_35_1_xcI.root',\n    '${inputFilePrefix}data_no_selection_36_1_JVW.root',\n    '${inputFilePrefix}data_no_selection_37_1_Wqj.root',\n    '${inputFilePrefix}data_no_selection_38_1_nPM.root',\n    '${inputFilePrefix}data_no_selection_39_1_ZcM.root',\n    '${inputFilePrefix}data_no_selection_3_1_FSz.root',\n    '${inputFilePrefix}data_no_selection_40_1_cyO.root',\n    '${inputFilePrefix}data_no_selection_41_1_ou0.root',\n    '${inputFilePrefix}data_no_selection_42_1_bVf.root',\n    '${inputFilePrefix}data_no_selection_43_1_1Uq.root',\n    '${inputFilePrefix}data_no_selection_44_1_5sV.root',\n    '${inputFilePrefix}data_no_selection_45_1_QFC.root',\n    '${inputFilePrefix}data_no_selection_46_1_XA8.root',\n    '${inputFilePrefix}data_no_selection_47_1_aVw.root',\n    '${inputFilePrefix}data_no_selection_48_1_wgk.root',\n    '${inputFilePrefix}data_no_selection_49_1_E8G.root',\n    '${inputFilePrefix}data_no_selection_4_1_dPb.root',\n    '${inputFilePrefix}data_no_selection_50_1_soH.root',\n    '${inputFilePrefix}data_no_selection_51_1_53E.root',\n    '${inputFilePrefix}data_no_selection_52_1_JTP.root',\n    '${inputFilePrefix}data_no_selection_53_1_EiC.root',\n    '${inputFilePrefix}data_no_selection_54_1_l9p.root',\n    '${inputFilePrefix}data_no_selection_55_1_CyK.root',\n    '${inputFilePrefix}data_no_selection_56_1_vLr.root',\n    '${inputFilePrefix}data_no_selection_57_1_ZLz.root',\n    '${inputFilePrefix}data_no_selection_58_1_o7R.root',\n    '${inputFilePrefix}data_no_selection_59_1_wAD.root',\n    '${inputFilePrefix}data_no_selection_5_1_Znv.root',\n    '${inputFilePrefix}data_no_selection_60_1_ECS.root',\n    '${inputFilePrefix}data_no_selection_61_1_RfL.root',\n    '${inputFilePrefix}data_no_selection_62_1_ijy.root',\n    '${inputFilePrefix}data_no_selection_63_1_okh.root',\n    '${inputFilePrefix}data_no_selection_64_1_rw3.root',\n    '${inputFilePrefix}data_no_selection_65_1_Ltj.root',\n    '${inputFilePrefix}data_no_selection_66_1_QNq.root',\n    '${inputFilePrefix}data_no_selection_67_1_sjB.root',\n    '${inputFilePrefix}data_no_selection_68_1_9Pp.root',\n    '${inputFilePrefix}data_no_selection_69_1_7zQ.root',\n    '${inputFilePrefix}data_no_selection_6_1_dB3.root',\n    '${inputFilePrefix}data_no_selection_70_1_Ng6.root',\n    '${inputFilePrefix}data_no_selection_71_1_ISQ.root',\n    '${inputFilePrefix}data_no_selection_72_1_dIl.root',\n    '${inputFilePrefix}data_no_selection_73_1_TYR.root',\n    '${inputFilePrefix}data_no_selection_74_1_02r.root',\n    '${inputFilePrefix}data_no_selection_75_1_OqH.root',\n    '${inputFilePrefix}data_no_selection_76_1_BF9.root',\n    '${inputFilePrefix}data_no_selection_77_1_JzM.root',\n    '${inputFilePrefix}data_no_selection_78_1_W1t.root',\n    '${inputFilePrefix}data_no_selection_79_1_s3B.root',\n    '${inputFilePrefix}data_no_selection_7_1_Ju6.root',\n    '${inputFilePrefix}data_no_selection_80_1_Rqr.root',\n    '${inputFilePrefix}data_no_selection_81_1_fph.root',\n    '${inputFilePrefix}data_no_selection_82_1_UNJ.root',\n    '${inputFilePrefix}data_no_selection_83_1_m5A.root',\n    '${inputFilePrefix}data_no_selection_84_1_EV2.root',\n    '${inputFilePrefix}data_no_selection_85_1_wQt.root',\n    '${inputFilePrefix}data_no_selection_86_1_y0C.root',\n    '${inputFilePrefix}data_no_selection_87_1_IcP.root',\n    '${inputFilePrefix}data_no_selection_88_1_AlL.root',\n    '${inputFilePrefix}data_no_selection_89_1_vJH.root',\n    '${inputFilePrefix}data_no_selection_8_1_hmS.root',\n    '${inputFilePrefix}data_no_selection_90_1_Mjz.root',\n    '${inputFilePrefix}data_no_selection_91_1_cjg.root',\n    '${inputFilePrefix}data_no_selection_92_1_8t8.root',\n    '${inputFilePrefix}data_no_selection_93_1_XHH.root',\n    '${inputFilePrefix}data_no_selection_94_1_VoZ.root',\n    '${inputFilePrefix}data_no_selection_95_1_trF.root',\n    '${inputFilePrefix}data_no_selection_96_1_bON.root',\n    '${inputFilePrefix}data_no_selection_97_1_PUE.root',\n    '${inputFilePrefix}data_no_selection_98_1_d5r.root',\n    '${inputFilePrefix}data_no_selection_99_1_agm.root',\n    '${inputFilePrefix}data_no_selection_9_1_RO2.root'\n    ])" )

#DrellYanAnalyzer output files
DrellYanAnalyzerOutputFiles=( "${DrellYanAnalyzerOutputFilePrefix}DrellYanAnalysis_WJetsToLNu_${version}.root" )

#samples
samples=( "WJetsToLNu" )

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
cat <<EOF > runWJetsToLNuDrellYanAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *WJetsToLNu*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runWJetsToLNuDrellYanAnalyzerCfgs.sh

#generate script that submits all jobs to LSF
cat <<EOF > submitWJetsToLNuDrellYanAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh drellyananalyzer*WJetsToLNu*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitWJetsToLNuDrellYanAnalyzerJobs.sh

#generate script that copies all files locally from EOS
cat <<EOF > copyDrellYanEnrichedWJetsToLNuFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  cmsStage -f /store/user/`whoami`/DrellYanAnalysis_WJetsToLNu_${version}.root /data1/`whoami`/WJetsToLNu/analysis/
  cmsRm /store/user/`whoami`/DrellYanAnalysis_WJetsToLNu_${version}.root
done

exit 0
EOF
chmod a+x copyDrellYanEnrichedWJetsToLNuFromEOS.sh

exit 0
