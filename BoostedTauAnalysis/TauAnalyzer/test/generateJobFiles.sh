#!/bin/bash

if [ $# -gt 3 ]
    then
    echo "Usage: ./generateJobFiles.sh <version> <MC template cfg> <data template cfg>"
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
./generateNonIsoWDataTauAnalyzerCfgs.sh $version $dataTemplateCfg
#./generateSinglePhotonDataTauAnalyzerCfgs.sh $version $dataTemplateCfg

exit 0
