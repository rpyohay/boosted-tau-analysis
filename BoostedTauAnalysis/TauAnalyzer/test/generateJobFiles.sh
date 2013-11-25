#!/bin/bash

if [ $# -ne 1 ]
    then
    echo "Usage: ./generateJobFiles.sh <version>"
    exit 0
fi

version=$1

./generateDYJetsToLLTauAnalyzerCfgs.sh $version
./generateDataTauAnalyzerCfgs.sh $version
./generateTTJetsTauAnalyzerCfgs.sh $version
./generateWJetsToLNuTauAnalyzerCfgs.sh $version
./generateWNJetsToLNuTauAnalyzerCfgs.sh $version
./generateSingleTopTauAnalyzerCfgs.sh $version
./generateWZTauAnalyzerCfgs.sh $version
./generateWWTauAnalyzerCfgs.sh $version
./generateZZTauAnalyzerCfgs.sh $version
./generateQCDTauAnalyzerCfgs.sh $version
./generateQCDBTauAnalyzerCfgs.sh $version
./generateQCDBMuTauAnalyzerCfgs.sh $version
./generateWh1TauAnalyzerCfg.sh $version

exit 0
