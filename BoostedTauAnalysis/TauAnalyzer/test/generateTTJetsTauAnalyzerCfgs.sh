#!/bin/bash

if [ $# -gt 2 ]
    then
    echo "Usage: ./generateTTJetsTauAnalyzerCfgs.sh <version> <template cfg>"
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
inputFilePrefix="root://eoscms//eos/cms/store/user/friccita/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-Summer12_DR53X-PU_S10_START53_V7A-v2-AODSIM_skim_v3/data_no_selection_"

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
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}100_2_7xP.root',\n    '${inputFilePrefix}101_2_kLf.root',\n    '${inputFilePrefix}102_1_pcO.root',\n    '${inputFilePrefix}103_2_ejr.root',\n    '${inputFilePrefix}104_2_Ae8.root',\n    '${inputFilePrefix}105_2_agQ.root',\n    '${inputFilePrefix}106_2_I4m.root',\n    '${inputFilePrefix}107_2_mGQ.root',\n    '${inputFilePrefix}108_2_sWV.root',\n    '${inputFilePrefix}109_2_6Lz.root',\n    '${inputFilePrefix}10_2_VUX.root',\n    '${inputFilePrefix}110_2_lz7.root',\n    '${inputFilePrefix}111_2_x94.root',\n    '${inputFilePrefix}112_2_0va.root',\n    '${inputFilePrefix}113_2_aaN.root',\n    '${inputFilePrefix}114_2_aHH.root',\n    '${inputFilePrefix}115_2_M8D.root',\n    '${inputFilePrefix}116_2_KdI.root',\n    '${inputFilePrefix}117_2_4xP.root',\n    '${inputFilePrefix}118_2_iyJ.root',\n    '${inputFilePrefix}119_2_kQC.root',\n    '${inputFilePrefix}11_2_0el.root',\n    '${inputFilePrefix}120_2_pj5.root',\n    '${inputFilePrefix}121_2_p6r.root',\n    '${inputFilePrefix}122_2_ffd.root',\n    '${inputFilePrefix}123_2_rjC.root',\n    '${inputFilePrefix}124_2_UkE.root',\n    '${inputFilePrefix}125_2_1zn.root',\n    '${inputFilePrefix}126_2_tCU.root',\n    '${inputFilePrefix}127_2_yEW.root',\n    '${inputFilePrefix}128_2_EpA.root',\n    '${inputFilePrefix}129_2_4Hj.root',\n    '${inputFilePrefix}12_2_cNK.root',\n    '${inputFilePrefix}130_2_F99.root',\n    '${inputFilePrefix}131_2_EK7.root',\n    '${inputFilePrefix}132_2_jlu.root',\n    '${inputFilePrefix}133_2_GhP.root',\n    '${inputFilePrefix}134_2_QGs.root',\n    '${inputFilePrefix}135_2_M4o.root',\n    '${inputFilePrefix}136_1_RWk.root',\n    '${inputFilePrefix}137_2_wCj.root',\n    '${inputFilePrefix}13_2_ggm.root',\n    '${inputFilePrefix}14_1_qrz.root',\n    '${inputFilePrefix}15_2_Jkh.root',\n    '${inputFilePrefix}16_1_50J.root',\n    '${inputFilePrefix}17_2_Qmj.root',\n    '${inputFilePrefix}18_2_9FJ.root',\n    '${inputFilePrefix}19_2_WV8.root',\n    '${inputFilePrefix}1_1_x7M.root',\n    '${inputFilePrefix}20_2_cHG.root',\n    '${inputFilePrefix}21_2_218.root',\n    '${inputFilePrefix}22_1_Td1.root',\n    '${inputFilePrefix}23_2_riW.root',\n    '${inputFilePrefix}24_2_cZJ.root',\n    '${inputFilePrefix}25_2_uSY.root',\n    '${inputFilePrefix}26_2_u9s.root',\n    '${inputFilePrefix}27_2_kPQ.root',\n    '${inputFilePrefix}28_2_mpW.root',\n    '${inputFilePrefix}29_2_SvH.root',\n    '${inputFilePrefix}2_1_iqR.root',\n    '${inputFilePrefix}30_2_VtQ.root',\n    '${inputFilePrefix}31_2_XlB.root',\n    '${inputFilePrefix}32_2_Bvx.root',\n    '${inputFilePrefix}33_2_CT8.root',\n    '${inputFilePrefix}34_2_UUk.root',\n    '${inputFilePrefix}35_2_iJ8.root',\n    '${inputFilePrefix}36_2_FJH.root',\n    '${inputFilePrefix}37_2_V9d.root',\n    '${inputFilePrefix}38_2_nCT.root',\n    '${inputFilePrefix}39_2_6VN.root',\n    '${inputFilePrefix}3_1_aNo.root',\n    '${inputFilePrefix}40_2_FWT.root',\n    '${inputFilePrefix}41_2_CZn.root',\n    '${inputFilePrefix}42_2_8wh.root',\n    '${inputFilePrefix}43_2_21q.root',\n    '${inputFilePrefix}44_2_7A1.root',\n    '${inputFilePrefix}45_2_ra9.root',\n    '${inputFilePrefix}46_2_7hN.root',\n    '${inputFilePrefix}47_2_Uly.root',\n    '${inputFilePrefix}48_2_IJr.root',\n    '${inputFilePrefix}49_2_kGU.root',\n    '${inputFilePrefix}4_2_Jp2.root',\n    '${inputFilePrefix}50_2_YTN.root',\n    '${inputFilePrefix}51_2_jp7.root',\n    '${inputFilePrefix}52_2_6lj.root',\n    '${inputFilePrefix}53_2_Bn7.root',\n    '${inputFilePrefix}54_2_TJL.root',\n    '${inputFilePrefix}55_2_CcB.root',\n    '${inputFilePrefix}56_2_BVM.root',\n    '${inputFilePrefix}57_2_7Tv.root',\n    '${inputFilePrefix}58_2_Hfg.root',\n    '${inputFilePrefix}59_2_MJa.root',\n    '${inputFilePrefix}5_1_tGU.root',\n    '${inputFilePrefix}60_2_Odq.root',\n    '${inputFilePrefix}61_1_Cyv.root',\n    '${inputFilePrefix}62_2_7zu.root',\n    '${inputFilePrefix}63_1_yMI.root',\n    '${inputFilePrefix}64_1_ddF.root',\n    '${inputFilePrefix}65_2_be7.root',\n    '${inputFilePrefix}66_2_W4q.root',\n    '${inputFilePrefix}67_2_KUY.root',\n    '${inputFilePrefix}68_2_iR9.root',\n    '${inputFilePrefix}69_2_THx.root',\n    '${inputFilePrefix}6_1_hLg.root',\n    '${inputFilePrefix}70_2_yz8.root',\n    '${inputFilePrefix}71_2_DvQ.root',\n    '${inputFilePrefix}72_2_OV1.root',\n    '${inputFilePrefix}73_2_CdT.root',\n    '${inputFilePrefix}74_2_CdP.root',\n    '${inputFilePrefix}75_2_VSA.root',\n    '${inputFilePrefix}76_2_xVO.root',\n    '${inputFilePrefix}77_2_RJW.root',\n    '${inputFilePrefix}78_2_4tq.root',\n    '${inputFilePrefix}79_2_96p.root',\n    '${inputFilePrefix}7_2_eQk.root',\n    '${inputFilePrefix}80_2_7od.root',\n    '${inputFilePrefix}81_2_TV8.root',\n    '${inputFilePrefix}82_2_Y9U.root',\n    '${inputFilePrefix}83_2_IeC.root',\n    '${inputFilePrefix}84_2_nNp.root',\n    '${inputFilePrefix}85_2_3EB.root',\n    '${inputFilePrefix}86_2_qdt.root',\n    '${inputFilePrefix}87_2_AT9.root',\n    '${inputFilePrefix}88_2_fiV.root',\n    '${inputFilePrefix}89_2_Tr3.root',\n    '${inputFilePrefix}8_2_505.root',\n    '${inputFilePrefix}90_2_Tfp.root',\n    '${inputFilePrefix}91_2_vPB.root',\n    '${inputFilePrefix}92_2_O67.root',\n    '${inputFilePrefix}93_2_pju.root',\n    '${inputFilePrefix}94_2_qSI.root',\n    '${inputFilePrefix}95_2_t3V.root',\n    '${inputFilePrefix}96_2_njH.root',\n    '${inputFilePrefix}97_2_10I.root',\n    '${inputFilePrefix}98_2_YLm.root',\n    '${inputFilePrefix}99_2_KxG.root',\n    '${inputFilePrefix}9_1_bCS.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_TTJets.root" )

#TauAnalyzer output files
highMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_TTJets_${version}.root" )
lowMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_TTJets_${version}.root" )
highMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_TTJets_${version}.root" )
lowMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_TTJets_${version}.root" )
highMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_TTJets_${version}.root" )
lowMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_TTJets_${version}.root" )

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

#generate script that submits all jobs to LSF
cat <<EOF > submitTTJetsTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*TTJets*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitTTJetsTauAnalyzerJobs.sh

#generate script that copies all files locally from EOS
cat <<EOF > copyTTJetsFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  for cut in "Iso" "NonIso" ""
    do
    for MTBin in high low
      do
      cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_TTJets_${version}.root /data1/`whoami`/TTJets/analysis/
      cmsRm /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_TTJets_${version}.root
    done
  done
done

exit 0
EOF
chmod a+x copyTTJetsFromEOS.sh

exit 0
