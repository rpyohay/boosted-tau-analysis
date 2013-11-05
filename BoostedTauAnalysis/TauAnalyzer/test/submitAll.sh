#!/bin/bash

cd v9
./submitDYJetsToLLTauAnalyzerJobs.sh
./submitDataTauAnalyzerJobs.sh
./submitTTJetsTauAnalyzerJobs.sh
./submitWNJetsToLNuTauAnalyzerJobs.sh
./submitSingleTopTauAnalyzerJobs.sh
./submitWZTauAnalyzerJobs.sh
./submitZZTauAnalyzerJobs.sh
./submitQCDTauAnalyzerJobs.sh
./submitQCDBTauAnalyzerJobs.sh
./submitQCDBMuTauAnalyzerJobs.sh

cd ..

exit 0
