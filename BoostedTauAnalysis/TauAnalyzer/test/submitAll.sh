#!/bin/bash

cd v12
./submitDYJetsToLLTauAnalyzerJobs.sh
./submitDataTauAnalyzerJobs.sh
./submitTTJetsTauAnalyzerJobs.sh
./submitWNJetsToLNuTauAnalyzerJobs.sh
#./submitWJetsToLNuTauAnalyzerJobs.sh
./submitSingleTopTauAnalyzerJobs.sh
./submitWZTauAnalyzerJobs.sh
./submitZZTauAnalyzerJobs.sh
./submitQCDTauAnalyzerJobs.sh
./submitQCDBTauAnalyzerJobs.sh
./submitQCDBMuTauAnalyzerJobs.sh

cd ..

exit 0
