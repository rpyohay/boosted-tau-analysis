#!/bin/bash

if [ $# -gt 2 ]
    then
    echo "Usage: ./generateSingleTopTauAnalyzerCfgs.sh <version> <template cfg>"
    exit 0
fi

####STUFF TO CONFIGURE####

#version
version=$1
templateCfg=$2
infoTag=""
dir=$version

#number of samples
nSamples=4
iBeg=0
iEnd=`expr $nSamples - 1`

#input file prefix and suffix
begInputFilePrefix="root://eoscms//eos/cms/store/user/friccita/T"
endInputFilePrefix="-channel_TuneZ2star_8TeV-powheg-tauola-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v3/data_no_selection_"
inputFilePrefixTSChannel="${begInputFilePrefix}_s${endInputFilePrefix}"
inputFilePrefixTbarSChannel="${begInputFilePrefix}bar_s${endInputFilePrefix}"
inputFilePrefixTTChannel="${begInputFilePrefix}_t${endInputFilePrefix}"
inputFilePrefixTbarTChannel="${begInputFilePrefix}bar_t${endInputFilePrefix}"

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
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefixTSChannel}1_1_osz.root',\n    '${inputFilePrefixTSChannel}2_1_iCr.root',\n    '${inputFilePrefixTSChannel}3_1_zsv.root',\n    '${inputFilePrefixTSChannel}4_1_o4R.root',\n    '${inputFilePrefixTSChannel}5_1_xfa.root',\n    '${inputFilePrefixTSChannel}6_1_xqp.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefixTbarSChannel}1_1_iAN.root',\n    '${inputFilePrefixTbarSChannel}2_1_QD4.root',\n    '${inputFilePrefixTbarSChannel}3_1_Egf.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefixTTChannel}10_1_wO8.root',\n    '${inputFilePrefixTTChannel}11_1_ysJ.root',\n    '${inputFilePrefixTTChannel}12_1_AIb.root',\n    '${inputFilePrefixTTChannel}13_1_ptm.root',\n    '${inputFilePrefixTTChannel}14_1_UZk.root',\n    '${inputFilePrefixTTChannel}15_1_S19.root',\n    '${inputFilePrefixTTChannel}16_1_0kc.root',\n    '${inputFilePrefixTTChannel}17_1_61V.root',\n    '${inputFilePrefixTTChannel}18_1_qx4.root',\n    '${inputFilePrefixTTChannel}19_1_Xsk.root',\n    '${inputFilePrefixTTChannel}1_1_yZW.root',\n    '${inputFilePrefixTTChannel}20_1_Vme.root',\n    '${inputFilePrefixTTChannel}21_1_9S9.root',\n    '${inputFilePrefixTTChannel}22_1_oVU.root',\n    '${inputFilePrefixTTChannel}23_1_utZ.root',\n    '${inputFilePrefixTTChannel}24_1_fP3.root',\n    '${inputFilePrefixTTChannel}25_1_wPi.root',\n    '${inputFilePrefixTTChannel}26_1_yui.root',\n    '${inputFilePrefixTTChannel}27_1_ZH9.root',\n    '${inputFilePrefixTTChannel}28_1_GSa.root',\n    '${inputFilePrefixTTChannel}29_1_PwD.root',\n    '${inputFilePrefixTTChannel}2_1_MkE.root',\n    '${inputFilePrefixTTChannel}30_1_G81.root',\n    '${inputFilePrefixTTChannel}31_1_Txh.root',\n    '${inputFilePrefixTTChannel}32_1_kqB.root',\n    '${inputFilePrefixTTChannel}33_1_93Z.root',\n    '${inputFilePrefixTTChannel}34_1_IAy.root',\n    '${inputFilePrefixTTChannel}35_1_dKK.root',\n    '${inputFilePrefixTTChannel}36_1_g85.root',\n    '${inputFilePrefixTTChannel}37_1_iKU.root',\n    '${inputFilePrefixTTChannel}38_1_GbY.root',\n    '${inputFilePrefixTTChannel}39_1_V4r.root',\n    '${inputFilePrefixTTChannel}3_1_yHQ.root',\n    '${inputFilePrefixTTChannel}40_1_vZJ.root',\n    '${inputFilePrefixTTChannel}41_1_qLi.root',\n    '${inputFilePrefixTTChannel}42_1_5ZF.root',\n    '${inputFilePrefixTTChannel}43_1_oiM.root',\n    '${inputFilePrefixTTChannel}44_1_2uY.root',\n    '${inputFilePrefixTTChannel}45_1_a9N.root',\n    '${inputFilePrefixTTChannel}46_1_644.root',\n    '${inputFilePrefixTTChannel}47_1_0d1.root',\n    '${inputFilePrefixTTChannel}48_1_Kg5.root',\n    '${inputFilePrefixTTChannel}49_1_WGm.root',\n    '${inputFilePrefixTTChannel}4_1_hVU.root',\n    '${inputFilePrefixTTChannel}50_1_nC2.root',\n    '${inputFilePrefixTTChannel}51_1_qFg.root',\n    '${inputFilePrefixTTChannel}52_1_PVM.root',\n    '${inputFilePrefixTTChannel}53_1_rFI.root',\n    '${inputFilePrefixTTChannel}54_1_wWT.root',\n    '${inputFilePrefixTTChannel}55_1_z0C.root',\n    '${inputFilePrefixTTChannel}56_1_Kqg.root',\n    '${inputFilePrefixTTChannel}57_1_2Cd.root',\n    '${inputFilePrefixTTChannel}58_1_cYy.root',\n    '${inputFilePrefixTTChannel}59_1_hM0.root',\n    '${inputFilePrefixTTChannel}5_1_pjG.root',\n    '${inputFilePrefixTTChannel}60_1_SfP.root',\n    '${inputFilePrefixTTChannel}61_1_PiG.root',\n    '${inputFilePrefixTTChannel}62_1_Nyj.root',\n    '${inputFilePrefixTTChannel}63_1_Zdk.root',\n    '${inputFilePrefixTTChannel}64_1_nO7.root',\n    '${inputFilePrefixTTChannel}65_1_InE.root',\n    '${inputFilePrefixTTChannel}66_1_9JV.root',\n    '${inputFilePrefixTTChannel}67_1_Bbv.root',\n    '${inputFilePrefixTTChannel}68_1_l6N.root',\n    '${inputFilePrefixTTChannel}69_1_Kn0.root',\n    '${inputFilePrefixTTChannel}6_1_wc2.root',\n    '${inputFilePrefixTTChannel}70_1_dT2.root',\n    '${inputFilePrefixTTChannel}71_1_oA5.root',\n    '${inputFilePrefixTTChannel}72_1_CMO.root',\n    '${inputFilePrefixTTChannel}73_1_QzH.root',\n    '${inputFilePrefixTTChannel}74_1_BTJ.root',\n    '${inputFilePrefixTTChannel}75_1_uoc.root',\n    '${inputFilePrefixTTChannel}76_1_eIb.root',\n    '${inputFilePrefixTTChannel}77_1_UVq.root',\n    '${inputFilePrefixTTChannel}7_1_tN2.root',\n    '${inputFilePrefixTTChannel}8_1_2nw.root',\n    '${inputFilePrefixTTChannel}9_1_VzE.root'\n    ])" "readFiles.extend([\n    '${inputFilePrefixTbarTChannel}10_1_zod.root',\n    '${inputFilePrefixTbarTChannel}11_1_F86.root',\n    '${inputFilePrefixTbarTChannel}12_1_qSU.root',\n    '${inputFilePrefixTbarTChannel}13_1_LrU.root',\n    '${inputFilePrefixTbarTChannel}14_1_3A6.root',\n    '${inputFilePrefixTbarTChannel}15_1_XeH.root',\n    '${inputFilePrefixTbarTChannel}16_1_jC2.root',\n    '${inputFilePrefixTbarTChannel}17_1_Il5.root',\n    '${inputFilePrefixTbarTChannel}18_1_psB.root',\n    '${inputFilePrefixTbarTChannel}19_1_Qd0.root',\n    '${inputFilePrefixTbarTChannel}1_1_ujA.root',\n    '${inputFilePrefixTbarTChannel}20_1_dS4.root',\n    '${inputFilePrefixTbarTChannel}21_1_MUT.root',\n    '${inputFilePrefixTbarTChannel}22_1_Eue.root',\n    '${inputFilePrefixTbarTChannel}23_1_nHB.root',\n    '${inputFilePrefixTbarTChannel}24_1_Aut.root',\n    '${inputFilePrefixTbarTChannel}25_1_i0x.root',\n    '${inputFilePrefixTbarTChannel}26_1_L0d.root',\n    '${inputFilePrefixTbarTChannel}27_1_IAp.root',\n    '${inputFilePrefixTbarTChannel}28_1_hI2.root',\n    '${inputFilePrefixTbarTChannel}29_1_Uau.root',\n    '${inputFilePrefixTbarTChannel}2_1_rEx.root',\n    '${inputFilePrefixTbarTChannel}30_1_FWY.root',\n    '${inputFilePrefixTbarTChannel}31_1_fF7.root',\n    '${inputFilePrefixTbarTChannel}32_1_KCb.root',\n    '${inputFilePrefixTbarTChannel}33_1_9k4.root',\n    '${inputFilePrefixTbarTChannel}34_1_ScX.root',\n    '${inputFilePrefixTbarTChannel}35_1_EPP.root',\n    '${inputFilePrefixTbarTChannel}36_1_Qfw.root',\n    '${inputFilePrefixTbarTChannel}37_1_RzM.root',\n    '${inputFilePrefixTbarTChannel}38_1_Yb0.root',\n    '${inputFilePrefixTbarTChannel}39_1_lZN.root',\n    '${inputFilePrefixTbarTChannel}3_1_mT6.root',\n    '${inputFilePrefixTbarTChannel}4_1_wp1.root',\n    '${inputFilePrefixTbarTChannel}5_1_clZ.root',\n    '${inputFilePrefixTbarTChannel}6_1_VeI.root',\n    '${inputFilePrefixTbarTChannel}7_1_zpQ.root',\n    '${inputFilePrefixTbarTChannel}8_1_jPO.root',\n    '${inputFilePrefixTbarTChannel}9_1_U0h.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_T_s-channel.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Tbar_s-channel.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_T_t-channel.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_Tbar_t-channel.root" )

#TauAnalyzer output files
highMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_T_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_Tbar_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_T_t-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_Tbar_t-channel_${version}.root" )
lowMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_T_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_Tbar_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_T_t-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_Tbar_t-channel_${version}.root" )
highMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_T_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_Tbar_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_T_t-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_Tbar_t-channel_${version}.root" )
lowMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_T_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_Tbar_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_T_t-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_Tbar_t-channel_${version}.root" )
highMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_T_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_Tbar_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_T_t-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_Tbar_t-channel_${version}.root" )
lowMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_T_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_Tbar_s-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_T_t-channel_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_Tbar_t-channel_${version}.root" )

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

#generate script that submits all jobs to LSF
cat <<EOF > submitSingleTopTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*T*channel*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitSingleTopTauAnalyzerJobs.sh

#generate script that copies all files locally from EOS
cat <<EOF > copySingleTopFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in "_s" "bar_s" "_t" "bar_t"
  do
  for cut in "Iso" "NonIso" ""
    do
    for MTBin in high low
      do
      cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_T\${sample}-channel_${version}.root /data1/`whoami`/SingleTop/analysis/
      cmsRm /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_T\${sample}-channel_${version}.root
    done
  done
done

exit 0
EOF
chmod a+x copySingleTopFromEOS.sh

exit 0
