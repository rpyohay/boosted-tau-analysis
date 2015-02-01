#!/bin/bash

if [ $# -ne 1 ]
    then
    echo "Usage: ./submitAll.sh <version>"
    exit 0
fi

version=$1

if [ ! -d $version ]
    then
    echo "Run ./generateJobFiles.sh before this script."
    exit 0
fi

cd $version
./submitDYJetsToLLTauAnalyzerJobs.sh
./submitDataTauAnalyzerJobs.sh
./submitTTJetsTauAnalyzerJobs.sh
./submitWNJetsToLNuTauAnalyzerJobs.sh
#./submitWbbTauAnalyzerJobs.sh
#./submitWJetsToLNuTauAnalyzerJobs.sh
./submitSingleTopTauAnalyzerJobs.sh
./submitWZTauAnalyzerJobs.sh
./submitWWTauAnalyzerJobs.sh
./submitZZTauAnalyzerJobs.sh
#./submitQCDTauAnalyzerJobs.sh
#./submitQCDBTauAnalyzerJobs.sh
#./submitQCDBMuTauAnalyzerJobs.sh
./submitNonIsoWDataTauAnalyzerJobs.sh
#./submitSinglePhotonDataTauAnalyzerJobs.sh
#./submitNonIsoWQCDTauAnalyzerJobs.sh
#./submitNonIsoWQCDBTauAnalyzerJobs.sh
#./submitNonIsoWQCDBMuTauAnalyzerJobs.sh
#./submitNonIsoWDYJetsToLLTauAnalyzerJobs.sh
#./submitNonIsoWTTJetsTauAnalyzerJobs.sh
#./submitNonIsoWWNJetsToLNuTauAnalyzerJobs.sh
cd ..

exit 0
