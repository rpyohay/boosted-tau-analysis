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
./submitDYJetsToLLAllTauAnalyzerJobs.sh
./submitDataAllTauAnalyzerJobs.sh
./submitTTJetsAllTauAnalyzerJobs.sh
./submitWNJetsToLNuAllTauAnalyzerJobs.sh
#./submitWbbAllTauAnalyzerJobs.sh
#./submitWJetsToLNuAllTauAnalyzerJobs.sh
./submitSingleTopAllTauAnalyzerJobs.sh
./submitWZAllTauAnalyzerJobs.sh
./submitWWAllTauAnalyzerJobs.sh
./submitZZAllTauAnalyzerJobs.sh
#./submitQCDAllTauAnalyzerJobs.sh
#./submitQCDBAllTauAnalyzerJobs.sh
#./submitQCDBMuAllTauAnalyzerJobs.sh
./submitNonIsoWDataAllTauAnalyzerJobs.sh
#./submitSinglePhotonDataAllTauAnalyzerJobs.sh
./submitNonIsoWDYJetsToLLAllTauAnalyzerJobs.sh
./submitNonIsoWTTJetsAllTauAnalyzerJobs.sh
./submitNonIsoWWNJetsToLNuAllTauAnalyzerJobs.sh
cd ..

exit 0
