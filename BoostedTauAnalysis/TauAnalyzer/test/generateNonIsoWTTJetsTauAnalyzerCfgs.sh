#!/bin/bash

if [ $# -gt 3 ]
    then
    echo "Usage: ./generateNonIsoWTTJetsTauAnalyzerCfgs.sh <version> <template cfg> [reweightOnly]"
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
inputFilePrefix="root://eoscms//eos/cms/store/user/friccita/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-Summer12_DR53X-PU_S10_START53_V7A-v2-AODSIM_skim_nonIsoW/data_no_selection"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/TTJets/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/nonIsoWTTJets/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}_100_1_DW7.root',\n    '${inputFilePrefix}_101_1_9h4.root',\n    '${inputFilePrefix}_102_1_U02.root',\n    '${inputFilePrefix}_103_1_lc5.root',\n    '${inputFilePrefix}_104_1_ce0.root',\n    '${inputFilePrefix}_105_1_c0U.root',\n    '${inputFilePrefix}_106_1_Ps2.root',\n    '${inputFilePrefix}_107_1_VjV.root',\n    '${inputFilePrefix}_108_1_i8j.root',\n    '${inputFilePrefix}_109_1_bKS.root',\n    '${inputFilePrefix}_10_1_Uk3.root',\n    '${inputFilePrefix}_110_1_Tm2.root',\n    '${inputFilePrefix}_111_1_eVU.root',\n    '${inputFilePrefix}_112_1_MDE.root',\n    '${inputFilePrefix}_113_1_P1z.root',\n    '${inputFilePrefix}_114_1_jN5.root',\n    '${inputFilePrefix}_115_1_TIB.root',\n    '${inputFilePrefix}_116_1_vZw.root',\n    '${inputFilePrefix}_117_1_FOh.root',\n    '${inputFilePrefix}_118_1_PKy.root',\n    '${inputFilePrefix}_119_1_gU0.root',\n    '${inputFilePrefix}_11_1_yTO.root',\n    '${inputFilePrefix}_120_1_wDZ.root',\n    '${inputFilePrefix}_121_1_Fc2.root',\n    '${inputFilePrefix}_122_1_xB5.root',\n    '${inputFilePrefix}_123_1_48j.root',\n    '${inputFilePrefix}_124_1_iT6.root',\n    '${inputFilePrefix}_125_1_zFs.root',\n    '${inputFilePrefix}_126_1_YlQ.root',\n    '${inputFilePrefix}_127_1_PKR.root',\n    '${inputFilePrefix}_128_1_UUO.root',\n    '${inputFilePrefix}_129_1_vuH.root',\n    '${inputFilePrefix}_12_1_CMS.root',\n    '${inputFilePrefix}_130_1_RTv.root',\n    '${inputFilePrefix}_131_1_uDD.root',\n    '${inputFilePrefix}_132_1_wSl.root',\n    '${inputFilePrefix}_133_1_WqA.root',\n    '${inputFilePrefix}_134_1_PaT.root',\n    '${inputFilePrefix}_135_1_BDz.root',\n    '${inputFilePrefix}_136_1_oVV.root',\n    '${inputFilePrefix}_137_1_qjK.root',\n    '${inputFilePrefix}_13_1_J2O.root',\n    '${inputFilePrefix}_14_1_uOq.root',\n    '${inputFilePrefix}_15_1_spS.root',\n    '${inputFilePrefix}_16_1_h9e.root',\n    '${inputFilePrefix}_17_1_H2F.root',\n    '${inputFilePrefix}_18_1_hLy.root',\n    '${inputFilePrefix}_19_1_3AA.root',\n    '${inputFilePrefix}_1_1_MLN.root',\n    '${inputFilePrefix}_20_1_NoE.root',\n    '${inputFilePrefix}_21_1_iM1.root',\n    '${inputFilePrefix}_22_1_dep.root',\n    '${inputFilePrefix}_23_1_wys.root',\n    '${inputFilePrefix}_24_1_oS1.root',\n    '${inputFilePrefix}_25_1_nya.root',\n    '${inputFilePrefix}_26_1_vCP.root',\n    '${inputFilePrefix}_27_1_ZXU.root',\n    '${inputFilePrefix}_28_1_tZJ.root',\n    '${inputFilePrefix}_29_1_o46.root',\n    '${inputFilePrefix}_2_1_ZFU.root',\n    '${inputFilePrefix}_30_1_eGv.root',\n    '${inputFilePrefix}_31_1_n67.root',\n    '${inputFilePrefix}_32_1_9fo.root',\n    '${inputFilePrefix}_33_1_Py3.root',\n    '${inputFilePrefix}_34_1_Rp0.root',\n    '${inputFilePrefix}_35_1_gnH.root',\n    '${inputFilePrefix}_36_1_tBM.root',\n    '${inputFilePrefix}_37_1_35i.root',\n    '${inputFilePrefix}_38_1_3UO.root',\n    '${inputFilePrefix}_39_1_JzZ.root',\n    '${inputFilePrefix}_3_1_Z6C.root',\n    '${inputFilePrefix}_40_1_Nva.root',\n    '${inputFilePrefix}_41_1_pP3.root',\n    '${inputFilePrefix}_42_1_mJ0.root',\n    '${inputFilePrefix}_43_1_AFe.root',\n    '${inputFilePrefix}_44_1_oRp.root',\n    '${inputFilePrefix}_45_1_cNR.root',\n    '${inputFilePrefix}_46_1_gF0.root',\n    '${inputFilePrefix}_47_1_clZ.root',\n    '${inputFilePrefix}_48_1_K9j.root',\n    '${inputFilePrefix}_49_1_uP2.root',\n    '${inputFilePrefix}_4_1_AbX.root',\n    '${inputFilePrefix}_50_1_t44.root',\n    '${inputFilePrefix}_51_1_VEe.root',\n    '${inputFilePrefix}_52_1_3Au.root',\n    '${inputFilePrefix}_53_1_hPm.root',\n    '${inputFilePrefix}_54_1_qpE.root',\n    '${inputFilePrefix}_55_1_Qs1.root',\n    '${inputFilePrefix}_56_1_FKN.root',\n    '${inputFilePrefix}_57_1_kum.root',\n    '${inputFilePrefix}_58_1_2Lm.root',\n    '${inputFilePrefix}_59_1_0n0.root',\n    '${inputFilePrefix}_5_1_4Gl.root',\n    '${inputFilePrefix}_60_1_7k4.root',\n    '${inputFilePrefix}_61_1_OWn.root',\n    '${inputFilePrefix}_62_1_v5N.root',\n    '${inputFilePrefix}_63_1_m5k.root',\n    '${inputFilePrefix}_64_1_Mdx.root',\n    '${inputFilePrefix}_65_1_cKL.root',\n    '${inputFilePrefix}_66_1_pzQ.root',\n    '${inputFilePrefix}_67_1_2gE.root',\n    '${inputFilePrefix}_68_1_GP5.root',\n    '${inputFilePrefix}_69_1_sHJ.root',\n    '${inputFilePrefix}_6_1_Jsm.root',\n    '${inputFilePrefix}_70_1_Cp9.root',\n    '${inputFilePrefix}_71_1_K57.root',\n    '${inputFilePrefix}_72_1_Skj.root',\n    '${inputFilePrefix}_73_1_KDN.root',\n    '${inputFilePrefix}_74_1_EhB.root',\n    '${inputFilePrefix}_75_1_dBM.root',\n    '${inputFilePrefix}_76_1_0Zp.root',\n    '${inputFilePrefix}_77_1_PtK.root',\n    '${inputFilePrefix}_78_1_W1E.root',\n    '${inputFilePrefix}_79_1_UxA.root',\n    '${inputFilePrefix}_7_1_SbH.root',\n    '${inputFilePrefix}_80_1_DuX.root',\n    '${inputFilePrefix}_81_1_AxF.root',\n    '${inputFilePrefix}_82_1_e5d.root',\n    '${inputFilePrefix}_83_1_jgs.root',\n    '${inputFilePrefix}_84_1_hTI.root',\n    '${inputFilePrefix}_85_1_KTt.root',\n    '${inputFilePrefix}_86_1_vtP.root',\n    '${inputFilePrefix}_87_1_mcX.root',\n    '${inputFilePrefix}_88_1_Af8.root',\n    '${inputFilePrefix}_89_1_Soe.root',\n    '${inputFilePrefix}_8_1_VHW.root',\n    '${inputFilePrefix}_90_1_pe1.root',\n    '${inputFilePrefix}_91_1_IFz.root',\n    '${inputFilePrefix}_92_1_1HO.root',\n    '${inputFilePrefix}_93_1_t6F.root',\n    '${inputFilePrefix}_94_1_zWR.root',\n    '${inputFilePrefix}_95_1_uEr.root',\n    '${inputFilePrefix}_96_1_fdK.root',\n    '${inputFilePrefix}_97_1_jsw.root',\n    '${inputFilePrefix}_98_1_ItX.root',\n    '${inputFilePrefix}_99_1_BMb.root',\n    '${inputFilePrefix}_9_1_403.root',\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_NonIsoWTTJets.root" )

#TauAnalyzer output files
highMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_NonIsoWTTJets_${version}.root" )
lowMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_NonIsoWTTJets_${version}.root" )
highMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_NonIsoWTTJets_${version}.root" )
lowMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_NonIsoWTTJets_${version}.root" )
nonIsoReweightTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoReweightAnalysis_NonIsoWTTJets_${version}.root" )
highMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_NonIsoWTTJets_${version}.root" )
lowMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_NonIsoWTTJets_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}NonIsoWTTJets${infoTag}_${version}.root" )

#samples
samples=( "NonIsoWTTJets" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%HIGHMTNONISOTAUANALYZEROUTFILE%${highMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTALLTAUANALYZEROUTFILE%${highMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTISOTAUANALYZEROUTFILE%${highMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTNONISOTAUANALYZEROUTFILE%${lowMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTALLTAUANALYZEROUTFILE%${lowMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTISOTAUANALYZEROUTFILE%${lowMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%HIGHMTSEQUENCE%process.highMTIsoTauAnalysisSequence%" -e "s%LOWMTSEQUENCE%process.lowMTIsoTauAnalysisSequence%" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_iso_cfg.py

  #generate cfg file for the non-isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%HIGHMTNONISOTAUANALYZEROUTFILE%${highMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTALLTAUANALYZEROUTFILE%${highMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTISOTAUANALYZEROUTFILE%${highMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTNONISOTAUANALYZEROUTFILE%${lowMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTALLTAUANALYZEROUTFILE%${lowMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTISOTAUANALYZEROUTFILE%${lowMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%HIGHMTSEQUENCE%process.highMTNonIsoTauAnalysisSequence%" -e "s%LOWMTSEQUENCE%process.lowMTNonIsoTauAnalysisSequence%" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_nonIso_cfg.py

  #generate cfg file for the non-isolated, reweighted sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%HIGHMTNONISOTAUANALYZEROUTFILE%${nonIsoReweightTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTALLTAUANALYZEROUTFILE%${highMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTISOTAUANALYZEROUTFILE%${highMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTNONISOTAUANALYZEROUTFILE%${lowMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTALLTAUANALYZEROUTFILE%${lowMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTISOTAUANALYZEROUTFILE%${lowMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%HIGHMTSEQUENCE%process.highMTNonIsoTauAnalysisSequence%" -e "s%LOWMTSEQUENCE%process.lowMTNonIsoTauAnalysisSequence%" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_nonIsoReweight_cfg.py

  #generate cfg file for the sample with no isolation cut
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%HIGHMTNONISOTAUANALYZEROUTFILE%${highMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTALLTAUANALYZEROUTFILE%${highMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTISOTAUANALYZEROUTFILE%${highMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTNONISOTAUANALYZEROUTFILE%${lowMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTALLTAUANALYZEROUTFILE%${lowMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTISOTAUANALYZEROUTFILE%${lowMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%HIGHMTSEQUENCE%process.highMTTauAnalysisSequence%" -e "s%LOWMTSEQUENCE%process.lowMTTauAnalysisSequence%" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_all_cfg.py

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
    cmsStage -f ${highMTNonIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
    cmsStage -f ${lowMTNonIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
    rm ${highMTNonIsoTauAnalyzerOutputFiles[${i}]} ${lowMTNonIsoTauAnalyzerOutputFiles[${i}]}
fi
#cmsRun \${fileNamePrefix}_nonIsoReweight_cfg.py
cmsStage -f ${highMTIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${lowMTIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
#cmsStage -f ${nonIsoReweightTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
rm ${highMTIsoTauAnalyzerOutputFiles[${i}]} ${lowMTIsoTauAnalyzerOutputFiles[${i}]}
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
cmsStage -f ${highMTAllTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${lowMTAllTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
rm ${highMTAllTauAnalyzerOutputFiles[${i}]} ${lowMTAllTauAnalyzerOutputFiles[${i}]}

exit 0
EOF
  chmod a+x tauanalyzer_${samples[${i}]}_all_cfg.sh
done

#generate run cfg that runs all files in the directory
cat <<EOF > runNonIsoWTTJetsTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *NonIsoWTTJets*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runNonIsoWTTJetsTauAnalyzerCfgs.sh

#generate script that submits all iso+nonIso+reweight jobs to LSF
cat <<EOF > submitNonIsoWTTJetsTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*NonIsoWTTJets*.sh | grep -v all | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitNonIsoWTTJetsTauAnalyzerJobs.sh

#generate script that submits all noIsoCut jobs to LSF
cat <<EOF > submitNonIsoWTTJetsAllTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*NonIsoWTTJets*all*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitNonIsoWTTJetsAllTauAnalyzerJobs.sh

#generate script that copies all iso+nonIso+reweight files locally from EOS
cat <<EOF > copyNonIsoWTTJetsFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  #for cut in Iso NonIso NonIsoReweight
  for cut in Iso NonIso
    do
    for MTBin in high low
      do
      if [ "\$cut" != "NonIso" ] || [ $reweightOnly -eq 0 ]
          then
          cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_NonIsoWTTJets_${version}.root /data1/`whoami`/nonIsoWTTJets/analysis/
          cmsRm /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_NonIsoWTTJets_${version}.root
      fi
    done
  done
done

exit 0
EOF
chmod a+x copyNonIsoWTTJetsFromEOS.sh

#generate script that copies all noIsoCut files locally from EOS
cat <<EOF > copyAllNonIsoWTTJetsFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in \`seq `expr $iBeg + 1` `expr $iEnd + 1`\`
  do
  for MTBin in high low
    do
    cmsStage -f /store/user/`whoami`/muHadAnalysis_\${MTBin}MT_NonIsoWTTJets_${version}.root /data1/`whoami`/nonIsoWTTJets/analysis/
    cmsRm /store/user/`whoami`/muHadAnalysis_\${MTBin}MT_NonIsoWTTJets_${version}.root
  done
done

exit 0
EOF
chmod a+x copyAllNonIsoWTTJetsFromEOS.sh

exit 0
