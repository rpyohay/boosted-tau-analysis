#!/bin/bash

if [ $# -gt 4 ]
    then
    echo "Usage: ./generateJobFiles.sh <version> <MC template cfg> <data template cfg> <non-isolated data template cfg>"
    exit 0
fi

version=$1
MCTemplateCfg="tauanalyzer_WNJetsToLNu_Wh1_template_cfg.py"
if [ -n "$2" ]
    then
    MCTemplateCfg=$2
fi
dataTemplateCfg="tauanalyzer_data_template_cfg.py"
if [ -n "$3" ]
    then
    dataTemplateCfg=$3
fi
nonisodataTemplateCfg="tauanalyzer_nonisodata_template_cfg.py"
if [ -n "$4" ]
    then
    nonisodataTemplateCfg=$4
fi

./generateDYJetsToLLTauAnalyzerCfgs.sh $version $MCTemplateCfg
./generateDataTauAnalyzerCfgs.sh $version $dataTemplateCfg
./generateTTJetsTauAnalyzerCfgs.sh $version $MCTemplateCfg
#./generateWJetsToLNuTauAnalyzerCfgs.sh $version $MCTemplateCfg
./generateWNJetsToLNuTauAnalyzerCfgs.sh $version $MCTemplateCfg
#./generateWbbTauAnalyzerCfgs.sh $version $MCTemplateCfg
./generateSingleTopTauAnalyzerCfgs.sh $version $MCTemplateCfg
./generateWZTauAnalyzerCfgs.sh $version $MCTemplateCfg
./generateWWTauAnalyzerCfgs.sh $version $MCTemplateCfg
./generateZZTauAnalyzerCfgs.sh $version $MCTemplateCfg
#./generateQCDTauAnalyzerCfgs.sh $version $MCTemplateCfg
#./generateQCDBTauAnalyzerCfgs.sh $version $MCTemplateCfg
#./generateQCDBMuTauAnalyzerCfgs.sh $version $MCTemplateCfg
./generateWh1TauAnalyzerCfg.sh $version $MCTemplateCfg
./generateggTauAnalyzerCfg.sh $version $MCTemplateCfg
./generateZHTauAnalyzerCfg.sh $version $MCTemplateCfg
./generateVBFTauAnalyzerCfg.sh $version $MCTemplateCfg
./generateNonIsoWDataTauAnalyzerCfgs.sh $version $nonisodataTemplateCfg
#./generateSinglePhotonDataTauAnalyzerCfgs.sh $version $dataTemplateCfg
#./generateNonIsoWQCDTauAnalyzerCfgs.sh $version $MCTemplateCfg
#./generateNonIsoWQCDBTauAnalyzerCfgs.sh $version $MCTemplateCfg
#./generateNonIsoWQCDBMuTauAnalyzerCfgs.sh $version $MCTemplateCfg
#./generateNonIsoWDYJetsToLLTauAnalyzerCfgs.sh $version $MCTemplateCfg
#./generateNonIsoWTTJetsTauAnalyzerCfgs.sh $version $MCTemplateCfg
#./generateNonIsoWWNJetsToLNuTauAnalyzerCfgs.sh $version $MCTemplateCfg

#generate run cfg that runs all signals in the directory
cd $version
cat <<EOF > runSigTauAnalyzerCfgs.sh
#!/bin/bash

./runggTauAnalyzerCfgs.sh &
./runVBFTauAnalyzerCfgs.sh &
sleep 10m
./runWh1TauAnalyzerCfgs.sh &
./runZHTauAnalyzerCfgs.sh &

exit 0
EOF
chmod a+x runSigTauAnalyzerCfgs.sh

exit 0
