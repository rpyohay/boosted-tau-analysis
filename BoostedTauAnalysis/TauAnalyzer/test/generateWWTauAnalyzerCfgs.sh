#!/bin/bash

if [ $# -gt 2 ]
    then
    echo "Usage: ./generateWWTauAnalyzerCfgs.sh <version> <template cfg>"
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
inputFilePrefix="root://eoscms//eos/cms/store/user/friccita/WW_TuneZ2star_8TeV_pythia6_tauola-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v3/data_no_selection_"

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
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}100_1_VQn.root',\n    '${inputFilePrefix}101_1_ta4.root',\n    '${inputFilePrefix}102_1_jMb.root',\n    '${inputFilePrefix}103_1_8uj.root',\n    '${inputFilePrefix}104_1_uCO.root',\n    '${inputFilePrefix}105_1_AXv.root',\n    '${inputFilePrefix}106_1_o1Q.root',\n    '${inputFilePrefix}107_1_2zX.root',\n    '${inputFilePrefix}108_1_sn5.root',\n    '${inputFilePrefix}109_1_Rl5.root',\n    '${inputFilePrefix}10_1_gMG.root',\n    '${inputFilePrefix}110_1_LYe.root',\n    '${inputFilePrefix}111_1_Ln2.root',\n    '${inputFilePrefix}112_1_ayy.root',\n    '${inputFilePrefix}113_1_tkZ.root',\n    '${inputFilePrefix}114_1_UvK.root',\n    '${inputFilePrefix}115_1_sRn.root',\n    '${inputFilePrefix}116_1_GbQ.root',\n    '${inputFilePrefix}117_1_7fr.root',\n    '${inputFilePrefix}118_1_7Dn.root',\n    '${inputFilePrefix}119_1_THe.root',\n    '${inputFilePrefix}11_1_BCN.root',\n    '${inputFilePrefix}120_1_8Xy.root',\n    '${inputFilePrefix}121_1_AQI.root',\n    '${inputFilePrefix}122_1_DGG.root',\n    '${inputFilePrefix}123_1_248.root',\n    '${inputFilePrefix}124_1_ZyA.root',\n    '${inputFilePrefix}125_1_OtU.root',\n    '${inputFilePrefix}126_1_Udk.root',\n    '${inputFilePrefix}127_1_haD.root',\n    '${inputFilePrefix}128_1_4Cb.root',\n    '${inputFilePrefix}129_1_eW7.root',\n    '${inputFilePrefix}12_1_u8U.root',\n    '${inputFilePrefix}130_1_4ia.root',\n    '${inputFilePrefix}131_1_FLO.root',\n    '${inputFilePrefix}132_1_We7.root',\n    '${inputFilePrefix}133_1_nWe.root',\n    '${inputFilePrefix}134_1_tHt.root',\n    '${inputFilePrefix}135_1_y8v.root',\n    '${inputFilePrefix}136_1_30m.root',\n    '${inputFilePrefix}137_1_C0w.root',\n    '${inputFilePrefix}138_1_Fth.root',\n    '${inputFilePrefix}139_1_tPO.root',\n    '${inputFilePrefix}13_1_bOY.root',\n    '${inputFilePrefix}140_1_Txj.root',\n    '${inputFilePrefix}141_1_uOh.root',\n    '${inputFilePrefix}142_1_H65.root',\n    '${inputFilePrefix}143_1_X0l.root',\n    '${inputFilePrefix}144_1_iRO.root',\n    '${inputFilePrefix}145_1_y7d.root',\n    '${inputFilePrefix}146_2_imL.root',\n    '${inputFilePrefix}147_1_35U.root',\n    '${inputFilePrefix}148_1_VgY.root',\n    '${inputFilePrefix}149_1_s4r.root',\n    '${inputFilePrefix}14_1_LAh.root',\n    '${inputFilePrefix}150_1_gsZ.root',\n    '${inputFilePrefix}151_1_DLV.root',\n    '${inputFilePrefix}152_1_sGR.root',\n    '${inputFilePrefix}153_1_L31.root',\n    '${inputFilePrefix}154_1_GrW.root',\n    '${inputFilePrefix}155_1_Qdm.root',\n    '${inputFilePrefix}156_1_lwQ.root',\n    '${inputFilePrefix}157_1_Qtv.root',\n    '${inputFilePrefix}158_1_o2K.root',\n    '${inputFilePrefix}159_1_KVK.root',\n    '${inputFilePrefix}15_1_ZtZ.root',\n    '${inputFilePrefix}160_1_R0R.root',\n    '${inputFilePrefix}161_1_Do2.root',\n    '${inputFilePrefix}162_1_KcF.root',\n    '${inputFilePrefix}163_1_EDh.root',\n    '${inputFilePrefix}164_1_ufA.root',\n    '${inputFilePrefix}165_1_23Z.root',\n    '${inputFilePrefix}166_1_VMI.root',\n    '${inputFilePrefix}167_1_05y.root',\n    '${inputFilePrefix}168_1_GVV.root',\n    '${inputFilePrefix}169_1_tDy.root',\n    '${inputFilePrefix}16_1_IxP.root',\n    '${inputFilePrefix}170_1_rX9.root',\n    '${inputFilePrefix}171_1_MUs.root',\n    '${inputFilePrefix}172_1_ABj.root',\n    '${inputFilePrefix}173_1_yzz.root',\n    '${inputFilePrefix}174_1_5pp.root',\n    '${inputFilePrefix}175_1_Vhw.root',\n    '${inputFilePrefix}176_1_ETE.root',\n    '${inputFilePrefix}177_1_7Rz.root',\n    '${inputFilePrefix}178_1_aA4.root',\n    '${inputFilePrefix}179_1_VXO.root',\n    '${inputFilePrefix}17_1_o8a.root',\n    '${inputFilePrefix}180_1_AXe.root',\n    '${inputFilePrefix}181_1_3jX.root',\n    '${inputFilePrefix}182_1_gzz.root',\n    '${inputFilePrefix}183_1_S8j.root',\n    '${inputFilePrefix}184_2_M2P.root',\n    '${inputFilePrefix}185_1_LKd.root',\n    '${inputFilePrefix}186_1_WaX.root',\n    '${inputFilePrefix}187_1_0i1.root',\n    '${inputFilePrefix}188_1_KOC.root',\n    '${inputFilePrefix}189_1_wJi.root',\n    '${inputFilePrefix}18_1_R2f.root',\n    '${inputFilePrefix}190_1_CJE.root',\n    '${inputFilePrefix}191_1_cTd.root',\n    '${inputFilePrefix}192_1_tD2.root',\n    '${inputFilePrefix}193_1_qJG.root',\n    '${inputFilePrefix}194_1_tUW.root',\n    '${inputFilePrefix}195_1_xn0.root',\n    '${inputFilePrefix}196_1_c34.root',\n    '${inputFilePrefix}197_1_k8q.root',\n    '${inputFilePrefix}198_1_xyc.root',\n    '${inputFilePrefix}199_1_lK2.root',\n    '${inputFilePrefix}19_1_Gjj.root',\n    '${inputFilePrefix}1_1_O7A.root',\n    '${inputFilePrefix}200_1_5aF.root',\n    '${inputFilePrefix}201_1_qBp.root',\n    '${inputFilePrefix}202_1_wIX.root',\n    '${inputFilePrefix}20_1_oW9.root',\n    '${inputFilePrefix}21_1_5EK.root',\n    '${inputFilePrefix}22_1_Z6w.root',\n    '${inputFilePrefix}23_1_jDm.root',\n    '${inputFilePrefix}24_1_UB8.root',\n    '${inputFilePrefix}25_1_x7Q.root',\n    '${inputFilePrefix}26_1_7Oi.root',\n    '${inputFilePrefix}27_1_TC6.root',\n    '${inputFilePrefix}28_1_Xva.root',\n    '${inputFilePrefix}29_1_xSO.root',\n    '${inputFilePrefix}2_1_Ocu.root',\n    '${inputFilePrefix}30_1_B1X.root',\n    '${inputFilePrefix}31_1_HZi.root',\n    '${inputFilePrefix}32_1_6rd.root',\n    '${inputFilePrefix}33_1_Wcb.root',\n    '${inputFilePrefix}34_1_tOV.root',\n    '${inputFilePrefix}35_1_Oka.root',\n    '${inputFilePrefix}36_1_B0R.root',\n    '${inputFilePrefix}37_1_Rb8.root',\n    '${inputFilePrefix}38_1_kr7.root',\n    '${inputFilePrefix}39_1_Jc5.root',\n    '${inputFilePrefix}3_1_ehl.root',\n    '${inputFilePrefix}40_1_CoZ.root',\n    '${inputFilePrefix}41_1_k9m.root',\n    '${inputFilePrefix}42_1_3Mb.root',\n    '${inputFilePrefix}43_1_5BV.root',\n    '${inputFilePrefix}44_1_0u9.root',\n    '${inputFilePrefix}45_1_QZN.root',\n    '${inputFilePrefix}46_1_sBN.root',\n    '${inputFilePrefix}47_1_OOy.root',\n    '${inputFilePrefix}48_1_R7s.root',\n    '${inputFilePrefix}49_1_COu.root',\n    '${inputFilePrefix}4_1_feo.root',\n    '${inputFilePrefix}50_1_JqR.root',\n    '${inputFilePrefix}51_1_OBe.root',\n    '${inputFilePrefix}52_1_xtz.root',\n    '${inputFilePrefix}53_1_0Vi.root',\n    '${inputFilePrefix}54_1_nSu.root',\n    '${inputFilePrefix}55_1_jhP.root',\n    '${inputFilePrefix}56_1_qNZ.root',\n    '${inputFilePrefix}57_1_gNR.root',\n    '${inputFilePrefix}58_1_1BU.root',\n    '${inputFilePrefix}59_1_Fh5.root',\n    '${inputFilePrefix}5_1_OhB.root',\n    '${inputFilePrefix}60_1_oFB.root',\n    '${inputFilePrefix}61_1_rDQ.root',\n    '${inputFilePrefix}62_1_B1q.root',\n    '${inputFilePrefix}63_1_0Q3.root',\n    '${inputFilePrefix}64_1_65e.root',\n    '${inputFilePrefix}65_1_rzo.root',\n    '${inputFilePrefix}66_1_RvX.root',\n    '${inputFilePrefix}67_1_lkQ.root',\n    '${inputFilePrefix}68_1_Ewb.root',\n    '${inputFilePrefix}69_1_yzZ.root',\n    '${inputFilePrefix}6_1_PEY.root',\n    '${inputFilePrefix}70_1_Flr.root',\n    '${inputFilePrefix}71_1_ra4.root',\n    '${inputFilePrefix}72_1_9kP.root',\n    '${inputFilePrefix}73_1_gwv.root',\n    '${inputFilePrefix}74_1_HDz.root',\n    '${inputFilePrefix}75_1_EIR.root',\n    '${inputFilePrefix}76_1_In1.root',\n    '${inputFilePrefix}77_1_Wyp.root',\n    '${inputFilePrefix}78_1_3sc.root',\n    '${inputFilePrefix}79_1_ByV.root',\n    '${inputFilePrefix}7_1_IyP.root',\n    '${inputFilePrefix}80_1_kBu.root',\n    '${inputFilePrefix}81_1_Hbf.root',\n    '${inputFilePrefix}82_1_bos.root',\n    '${inputFilePrefix}83_1_7Xc.root',\n    '${inputFilePrefix}84_1_pf9.root',\n    '${inputFilePrefix}85_1_Oc5.root',\n    '${inputFilePrefix}86_1_7SU.root',\n    '${inputFilePrefix}87_1_2eM.root',\n    '${inputFilePrefix}88_1_v1r.root',\n    '${inputFilePrefix}89_1_Tea.root',\n    '${inputFilePrefix}8_1_X0z.root',\n    '${inputFilePrefix}90_1_ZFn.root',\n    '${inputFilePrefix}91_1_0qe.root',\n    '${inputFilePrefix}92_1_iMV.root',\n    '${inputFilePrefix}93_1_YKF.root',\n    '${inputFilePrefix}94_1_ok1.root',\n    '${inputFilePrefix}95_1_gHd.root',\n    '${inputFilePrefix}96_1_u8s.root',\n    '${inputFilePrefix}97_1_wuW.root',\n    '${inputFilePrefix}98_1_Z0K.root',\n    '${inputFilePrefix}99_1_GPu.root',\n    '${inputFilePrefix}9_1_Pfo.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_WW.root" )

#TauAnalyzer output files
highMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_WW_${version}.root" )
lowMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_WW_${version}.root" )
highMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_WW_${version}.root" )
lowMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_WW_${version}.root" )
highMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_WW_${version}.root" )
lowMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_WW_${version}.root" )

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

  #generate cfg file
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%HIGHMTNONISOTAUANALYZEROUTFILE%${highMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTALLTAUANALYZEROUTFILE%${highMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTISOTAUANALYZEROUTFILE%${highMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTNONISOTAUANALYZEROUTFILE%${lowMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTALLTAUANALYZEROUTFILE%${lowMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTISOTAUANALYZEROUTFILE%${lowMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%HIGGSREW%False%" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_cfg.py

  #generate job submission script for LSF
  cat <<EOF > tauanalyzer_${samples[${i}]}_cfg.sh
#!/bin/bash

jobDir="`pwd`"
fileNamePrefix="tauanalyzer_${samples[${i}]}"

cd \$jobDir
eval \`scramv1 runtime -sh\`
cd -
cp \$jobDir/\${fileNamePrefix}_cfg.py .
cmsRun \${fileNamePrefix}_cfg.py
cmsStage -f ${highMTNonIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${lowMTNonIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${highMTIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${lowMTIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${highMTAllTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${lowMTAllTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
rm ${highMTNonIsoTauAnalyzerOutputFiles[${i}]} ${lowMTNonIsoTauAnalyzerOutputFiles[${i}]} ${highMTIsoTauAnalyzerOutputFiles[${i}]} ${lowMTIsoTauAnalyzerOutputFiles[${i}]} ${highMTAllTauAnalyzerOutputFiles[${i}]} ${lowMTAllTauAnalyzerOutputFiles[${i}]}

exit 0
EOF
  chmod a+x tauanalyzer_${samples[${i}]}_cfg.sh

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

#generate script that submits all jobs to LSF
cat <<EOF > submitWWTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*WW*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitWWTauAnalyzerJobs.sh

#generate script that copies all files locally from EOS
cat <<EOF > copyWWFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  for cut in "Iso" "NonIso" ""
    do
    for MTBin in high low
      do
      cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_WW_${version}.root /data1/`whoami`/WW/analysis/
      cmsRm /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_WW_${version}.root
    done
  done
done

exit 0
EOF
chmod a+x copyWWFromEOS.sh

exit 0
