#!/bin/bash

####STUFF TO CONFIGURE####

#version
version="v10"
infoTag=""
dir=$version

#number of samples
nSamples=4
iBeg=0
iEnd=`expr $nSamples - 1`

#input file prefix
inputFilePrefix="root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/Wh1ToTauMuTauHadX/"

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
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_10_1_jaY.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_11_1_Pr8.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_12_1_C6u.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_13_1_HFc.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_14_1_43l.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_15_1_9jm.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_16_1_6cn.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_1_1_pH5.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_21_1_rQw.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_26_1_FFz.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_28_1_oWb.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_29_1_23b.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_2_1_QrZ.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_30_1_143.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_31_1_OWV.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_32_1_iQ6.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_3_1_6SS.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_4_1_mH5.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_5_1_Y1u.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_6_1_v0B.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_7_1_8Ou.root',\n    '${inputFilePrefix}SingleMu_Run2012A-22Jan2013-v1_AOD_skim_v2/data_no_selection_8_1_mzc.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_100_1_zKy.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_10_1_5pv.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_11_1_P8e.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_12_1_BIN.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_13_1_aKe.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_14_1_E4B.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_15_1_Xts.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_17_1_O1b.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_1_1_Yy4.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_20_1_zjw.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_2_1_mNS.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_3_1_eDd.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_46_1_ZP5.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_4_1_EQB.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_5_1_Unx.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_62_1_kWm.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_64_1_LNQ.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_65_1_whr.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_67_1_OyG.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_68_1_Erm.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_69_1_hFy.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_6_1_YX1.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_71_1_DlK.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_73_1_ToZ.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_74_1_KTR.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_77_1_rIW.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_78_1_TNS.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_79_1_bVM.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_7_1_Ewa.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_81_1_XSp.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_83_1_ic0.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_84_1_GPE.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_85_1_ggJ.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_86_1_d1h.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_87_1_qMM.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_8_1_OBW.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_91_1_ePD.root',\n    '${inputFilePrefix}SingleMu_Run2012B-22Jan2013-v1_AOD_skim_v2/data_no_selection_9_1_10u.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefix}SingleMu_Run2012C-22Jan2013-v1_AOD_skim_v2/data_no_selection_146_1_66O.root',\n    '${inputFilePrefix}SingleMu_Run2012C-22Jan2013-v1_AOD_skim_v2/data_no_selection_160_1_Jqv.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefix}SingleMu_Run2012D-22Jan2013-v1_AOD_skim_v2/data_no_selection_104_1_KhD.root',\n    '${inputFilePrefix}SingleMu_Run2012D-22Jan2013-v1_AOD_skim_v2/data_no_selection_105_1_A0A.root',\n    '${inputFilePrefix}SingleMu_Run2012D-22Jan2013-v1_AOD_skim_v2/data_no_selection_106_1_UNf.root',\n    '${inputFilePrefix}SingleMu_Run2012D-22Jan2013-v1_AOD_skim_v2/data_no_selection_109_1_DSr.root',\n    '${inputFilePrefix}SingleMu_Run2012D-22Jan2013-v1_AOD_skim_v2/data_no_selection_149_1_GnI.root',\n    '${inputFilePrefix}SingleMu_Run2012D-22Jan2013-v1_AOD_skim_v2/data_no_selection_35_1_Fld.root',\n    '${inputFilePrefix}SingleMu_Run2012D-22Jan2013-v1_AOD_skim_v2/data_no_selection_37_1_D4s.root',\n    '${inputFilePrefix}SingleMu_Run2012D-22Jan2013-v1_AOD_skim_v2/data_no_selection_38_1_ZRb.root',\n    '${inputFilePrefix}SingleMu_Run2012D-22Jan2013-v1_AOD_skim_v2/data_no_selection_39_1_t7V.root',\n    '${inputFilePrefix}SingleMu_Run2012D-22Jan2013-v1_AOD_skim_v2/data_no_selection_84_1_4vT.root',\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_SingleMu_Run2012A.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_SingleMu_Run2012B.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_SingleMu_Run2012C.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_SingleMu_Run2012D.root" )

#TauAnalyzer output files
isoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_SingleMu_Run2012A_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_SingleMu_Run2012B_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_SingleMu_Run2012C_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_SingleMu_Run2012D_${version}.root" )
nonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_SingleMu_Run2012A_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_SingleMu_Run2012B_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_SingleMu_Run2012C_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_SingleMu_Run2012D_${version}.root" )
allTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_SingleMu_Run2012A_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_SingleMu_Run2012B_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_SingleMu_Run2012C_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_SingleMu_Run2012D_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}SingleMu_Run2012A${infoTag}_${version}.root" "${EDMOutputFilePrefix}SingleMu_Run2012B${infoTag}_${version}.root" "${EDMOutputFilePrefix}SingleMu_Run2012C${infoTag}_${version}.root" "${EDMOutputFilePrefix}SingleMu_Run2012D${infoTag}_${version}.root" )

#samples
samples=( "SingleMu_Run2012A" "SingleMu_Run2012B" "SingleMu_Run2012C" "SingleMu_Run2012D" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file for the isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.isoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" ../tauanalyzer_data_template_cfg.py > tauanalyzer_${samples[${i}]}_iso_cfg.py

  #generate cfg file for the non-isolated sample
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%NONISOTAUANALYZEROUTFILE%${nonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%ALLTAUANALYZEROUTFILE%${allTauAnalyzerOutputFiles[${i}]}%" -e "s%ISOTAUANALYZEROUTFILE%${isoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%SEQUENCE%process.nonIsoTauAnalysisSequence%" -e "s%PUSCENARIO%S10%" ../tauanalyzer_data_template_cfg.py > tauanalyzer_${samples[${i}]}_nonIso_cfg.py

  #generate job submission script for LSF
  cat <<EOF > tauanalyzer_${samples[${i}]}_cfg.sh
#!/bin/bash

jobDir="`pwd`"
fileNamePrefix="tauanalyzer_${samples[${i}]}"

cd \$jobDir
eval \`scramv1 runtime -sh\`
cd -
cp \$jobDir/\${fileNamePrefix}_iso_cfg.py \$jobDir/\${fileNamePrefix}_nonIso_cfg.py .
#cmsRun \${fileNamePrefix}_iso_cfg.py #BLINDED!!!
cmsRun \${fileNamePrefix}_nonIso_cfg.py
#cmsStage -f ${isoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/ #BLINDED!!!
cmsStage -f ${nonIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
#rm ${isoTauAnalyzerOutputFiles[${i}]} #BLINDED!!!
rm ${nonIsoTauAnalyzerOutputFiles[${i}]} 

exit 0
EOF
  chmod a+x tauanalyzer_${samples[${i}]}_cfg.sh
done

#generate run cfg that runs all files in the directory
cat <<EOF > runDataTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *SingleMu*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  isIso=\`echo \$file | sed -e "s%.*\(iso\).*%\1%"\`
  if [ "\$isIso" != "iso" ] #BLINDED!!!
      then
      cmsRun \$file > \$outFile
  else
      echo "Not running script \$file due to blinding requirement"
  fi
done

exit 0
EOF
chmod a+x runDataTauAnalyzerCfgs.sh

#generate script that submits all jobs to LSF
cat <<EOF > submitDataTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*SingleMu*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 8nh -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitDataTauAnalyzerJobs.sh

#generate script that copies all files locally from EOS
cat <<EOF > copyDataFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in "A" "B" "C" "D"
  do
#  for cut in Iso NonIso #BLINDED!!!
  for cut in NonIso
    do
    cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_SingleMu_Run2012\${sample}_${version}.root /data1/`whoami`/data/analysis/
  done
done

exit 0
EOF
chmod a+x copyDataFromEOS.sh

exit 0
