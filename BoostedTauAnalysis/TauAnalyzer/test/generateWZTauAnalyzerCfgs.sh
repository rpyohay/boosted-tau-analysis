#!/bin/bash

if [ $# -gt 2 ]
    then
    echo "Usage: ./generateWZTauAnalyzerCfgs.sh <version> <template cfg>"
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
inputFilePrefix="root://eoscms//eos/cms/store/user/friccita/WZ_TuneZ2star_8TeV_pythia6_tauola-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v3/data_no_selection_"

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
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}100_2_IPb.root',\n    '${inputFilePrefix}101_2_iYs.root',\n    '${inputFilePrefix}102_2_jnD.root',\n    '${inputFilePrefix}103_2_p0A.root',\n    '${inputFilePrefix}104_2_rHK.root',\n    '${inputFilePrefix}105_2_dxl.root',\n    '${inputFilePrefix}106_2_bqA.root',\n    '${inputFilePrefix}107_1_Voj.root',\n    '${inputFilePrefix}108_1_gPw.root',\n    '${inputFilePrefix}109_1_1P7.root',\n    '${inputFilePrefix}10_1_0cq.root',\n    '${inputFilePrefix}110_1_uX8.root',\n    '${inputFilePrefix}111_1_W5p.root',\n    '${inputFilePrefix}112_2_Feg.root',\n    '${inputFilePrefix}113_1_vbf.root',\n    '${inputFilePrefix}114_2_kdE.root',\n    '${inputFilePrefix}115_1_Xa2.root',\n    '${inputFilePrefix}116_2_zqe.root',\n    '${inputFilePrefix}117_2_63q.root',\n    '${inputFilePrefix}118_2_pv3.root',\n    '${inputFilePrefix}119_2_mk3.root',\n    '${inputFilePrefix}11_1_7hY.root',\n    '${inputFilePrefix}120_1_xoM.root',\n    '${inputFilePrefix}121_1_y2N.root',\n    '${inputFilePrefix}122_1_haA.root',\n    '${inputFilePrefix}123_1_vqy.root',\n    '${inputFilePrefix}124_2_QtV.root',\n    '${inputFilePrefix}125_2_Xix.root',\n    '${inputFilePrefix}126_2_wRK.root',\n    '${inputFilePrefix}127_2_KUY.root',\n    '${inputFilePrefix}128_2_85Q.root',\n    '${inputFilePrefix}129_2_CT5.root',\n    '${inputFilePrefix}12_1_rNS.root',\n    '${inputFilePrefix}130_1_ehp.root',\n    '${inputFilePrefix}131_1_raU.root',\n    '${inputFilePrefix}132_2_j0X.root',\n    '${inputFilePrefix}133_2_Wxz.root',\n    '${inputFilePrefix}134_1_m8t.root',\n    '${inputFilePrefix}135_2_hXc.root',\n    '${inputFilePrefix}136_1_kiM.root',\n    '${inputFilePrefix}137_2_wwQ.root',\n    '${inputFilePrefix}138_2_m88.root',\n    '${inputFilePrefix}139_2_H4W.root',\n    '${inputFilePrefix}13_1_OGe.root',\n    '${inputFilePrefix}140_1_Pfe.root',\n    '${inputFilePrefix}141_1_MMb.root',\n    '${inputFilePrefix}142_1_EMn.root',\n    '${inputFilePrefix}143_1_YIX.root',\n    '${inputFilePrefix}144_2_euM.root',\n    '${inputFilePrefix}145_1_XPq.root',\n    '${inputFilePrefix}146_1_1bg.root',\n    '${inputFilePrefix}147_1_ePa.root',\n    '${inputFilePrefix}148_2_BtO.root',\n    '${inputFilePrefix}149_2_ihu.root',\n    '${inputFilePrefix}14_1_WIE.root',\n    '${inputFilePrefix}150_2_ZdG.root',\n    '${inputFilePrefix}151_1_wVI.root',\n    '${inputFilePrefix}152_1_Rkm.root',\n    '${inputFilePrefix}153_1_pRG.root',\n    '${inputFilePrefix}154_1_7CE.root',\n    '${inputFilePrefix}155_1_huy.root',\n    '${inputFilePrefix}156_1_eAc.root',\n    '${inputFilePrefix}157_1_GJS.root',\n    '${inputFilePrefix}158_1_j3v.root',\n    '${inputFilePrefix}159_1_xzX.root',\n    '${inputFilePrefix}15_1_IbV.root',\n    '${inputFilePrefix}160_1_ol3.root',\n    '${inputFilePrefix}161_1_Y7K.root',\n    '${inputFilePrefix}162_1_mdp.root',\n    '${inputFilePrefix}163_1_UMv.root',\n    '${inputFilePrefix}164_1_FoI.root',\n    '${inputFilePrefix}165_1_TNL.root',\n    '${inputFilePrefix}166_1_hMv.root',\n    '${inputFilePrefix}167_1_VGO.root',\n    '${inputFilePrefix}168_1_Wyn.root',\n    '${inputFilePrefix}169_1_HfW.root',\n    '${inputFilePrefix}16_1_6TR.root',\n    '${inputFilePrefix}170_1_UI1.root',\n    '${inputFilePrefix}171_1_hPG.root',\n    '${inputFilePrefix}172_1_tOD.root',\n    '${inputFilePrefix}173_1_4Rn.root',\n    '${inputFilePrefix}174_2_Zb3.root',\n    '${inputFilePrefix}175_1_tYR.root',\n    '${inputFilePrefix}176_1_EUJ.root',\n    '${inputFilePrefix}177_1_95W.root',\n    '${inputFilePrefix}178_2_75j.root',\n    '${inputFilePrefix}179_2_R6s.root',\n    '${inputFilePrefix}17_1_cIk.root',\n    '${inputFilePrefix}180_1_wdL.root',\n    '${inputFilePrefix}181_2_o1p.root',\n    '${inputFilePrefix}182_2_Fyb.root',\n    '${inputFilePrefix}183_1_dP7.root',\n    '${inputFilePrefix}184_1_UQX.root',\n    '${inputFilePrefix}185_1_eap.root',\n    '${inputFilePrefix}186_1_w8I.root',\n    '${inputFilePrefix}187_2_CBx.root',\n    '${inputFilePrefix}188_1_KFm.root',\n    '${inputFilePrefix}189_1_XC5.root',\n    '${inputFilePrefix}18_1_Ojx.root',\n    '${inputFilePrefix}190_1_HOv.root',\n    '${inputFilePrefix}191_1_ewT.root',\n    '${inputFilePrefix}192_2_Tch.root',\n    '${inputFilePrefix}193_1_7rK.root',\n    '${inputFilePrefix}194_1_c7j.root',\n    '${inputFilePrefix}195_1_uYS.root',\n    '${inputFilePrefix}196_2_HLd.root',\n    '${inputFilePrefix}197_1_XXP.root',\n    '${inputFilePrefix}198_2_1tQ.root',\n    '${inputFilePrefix}199_1_h8Q.root',\n    '${inputFilePrefix}19_1_CNK.root',\n    '${inputFilePrefix}1_1_bRw.root',\n    '${inputFilePrefix}200_1_dPz.root',\n    '${inputFilePrefix}201_1_Z2b.root',\n    '${inputFilePrefix}202_1_Xy1.root',\n    '${inputFilePrefix}203_2_AsV.root',\n    '${inputFilePrefix}20_1_Xqe.root',\n    '${inputFilePrefix}21_1_gmS.root',\n    '${inputFilePrefix}22_1_I7i.root',\n    '${inputFilePrefix}23_1_HU0.root',\n    '${inputFilePrefix}24_1_baO.root',\n    '${inputFilePrefix}25_1_B0s.root',\n    '${inputFilePrefix}26_1_4t9.root',\n    '${inputFilePrefix}27_1_fou.root',\n    '${inputFilePrefix}28_1_jUV.root',\n    '${inputFilePrefix}29_1_xLn.root',\n    '${inputFilePrefix}2_1_TZv.root',\n    '${inputFilePrefix}30_1_0re.root',\n    '${inputFilePrefix}31_1_n8K.root',\n    '${inputFilePrefix}32_1_AWc.root',\n    '${inputFilePrefix}33_1_CMX.root',\n    '${inputFilePrefix}34_1_wUP.root',\n    '${inputFilePrefix}35_1_gRY.root',\n    '${inputFilePrefix}36_1_Bnk.root',\n    '${inputFilePrefix}37_1_Bcq.root',\n    '${inputFilePrefix}38_1_6hO.root',\n    '${inputFilePrefix}39_1_T6a.root',\n    '${inputFilePrefix}3_1_IiP.root',\n    '${inputFilePrefix}40_1_zMS.root',\n    '${inputFilePrefix}41_1_0oQ.root',\n    '${inputFilePrefix}42_1_Lz8.root',\n    '${inputFilePrefix}43_1_keZ.root',\n    '${inputFilePrefix}44_1_5r4.root',\n    '${inputFilePrefix}45_1_mxE.root',\n    '${inputFilePrefix}46_1_vFC.root',\n    '${inputFilePrefix}47_1_TCJ.root',\n    '${inputFilePrefix}48_1_nmY.root',\n    '${inputFilePrefix}49_1_mlH.root',\n    '${inputFilePrefix}4_1_rir.root',\n    '${inputFilePrefix}50_1_LgE.root',\n    '${inputFilePrefix}51_1_CZw.root',\n    '${inputFilePrefix}52_1_dOv.root',\n    '${inputFilePrefix}53_1_SVi.root',\n    '${inputFilePrefix}54_1_pme.root',\n    '${inputFilePrefix}55_1_eXu.root',\n    '${inputFilePrefix}56_1_Q6W.root',\n    '${inputFilePrefix}57_2_mlw.root',\n    '${inputFilePrefix}58_1_nvz.root',\n    '${inputFilePrefix}59_1_a92.root',\n    '${inputFilePrefix}5_1_5ID.root',\n    '${inputFilePrefix}60_1_HnB.root',\n    '${inputFilePrefix}61_1_8kR.root',\n    '${inputFilePrefix}62_1_tns.root',\n    '${inputFilePrefix}63_1_wPa.root',\n    '${inputFilePrefix}64_1_bxO.root',\n    '${inputFilePrefix}65_1_Jqv.root',\n    '${inputFilePrefix}66_1_wOM.root',\n    '${inputFilePrefix}67_1_lwp.root',\n    '${inputFilePrefix}68_1_riD.root',\n    '${inputFilePrefix}69_1_Hjd.root',\n    '${inputFilePrefix}6_1_abD.root',\n    '${inputFilePrefix}70_1_HUb.root',\n    '${inputFilePrefix}71_1_oSa.root',\n    '${inputFilePrefix}72_1_Geo.root',\n    '${inputFilePrefix}73_1_6Wd.root',\n    '${inputFilePrefix}74_1_YYi.root',\n    '${inputFilePrefix}75_1_Aiz.root',\n    '${inputFilePrefix}76_1_rco.root',\n    '${inputFilePrefix}77_1_rmR.root',\n    '${inputFilePrefix}78_1_ZdU.root',\n    '${inputFilePrefix}79_1_voH.root',\n    '${inputFilePrefix}7_1_uOY.root',\n    '${inputFilePrefix}80_1_6aa.root',\n    '${inputFilePrefix}81_1_dHG.root',\n    '${inputFilePrefix}82_1_CAs.root',\n    '${inputFilePrefix}83_2_0yK.root',\n    '${inputFilePrefix}84_1_r3F.root',\n    '${inputFilePrefix}85_1_vQc.root',\n    '${inputFilePrefix}86_1_CXA.root',\n    '${inputFilePrefix}87_1_6Za.root',\n    '${inputFilePrefix}88_1_Q4a.root',\n    '${inputFilePrefix}89_1_UYQ.root',\n    '${inputFilePrefix}8_1_lRe.root',\n    '${inputFilePrefix}90_1_j97.root',\n    '${inputFilePrefix}91_1_uYP.root',\n    '${inputFilePrefix}92_1_Szz.root',\n    '${inputFilePrefix}93_1_4w7.root',\n    '${inputFilePrefix}94_1_4KV.root',\n    '${inputFilePrefix}95_2_fpn.root',\n    '${inputFilePrefix}96_2_Oiy.root',\n    '${inputFilePrefix}97_2_6PR.root',\n    '${inputFilePrefix}98_2_f2I.root',\n    '${inputFilePrefix}99_2_hy6.root',\n    '${inputFilePrefix}9_1_lmk.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_WZ.root" )

#TauAnalyzer output files
highMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_WZ_${version}.root" )
lowMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_WZ_${version}.root" )
highMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_WZ_${version}.root" )
lowMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_WZ_${version}.root" )
highMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_WZ_${version}.root" )
lowMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_WZ_${version}.root" )

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
  for cut in "Iso" "NonIso" ""
    do
    for MTBin in high low
      do
      cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_WZ_${version}.root /data1/`whoami`/WZ/analysis/
      cmsRm /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_WZ_${version}.root
    done
  done
done

exit 0
EOF
chmod a+x copyWZFromEOS.sh

exit 0
