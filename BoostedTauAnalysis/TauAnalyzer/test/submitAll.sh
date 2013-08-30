#!/bin/bash

cd v43
./submitDYJetsToLLTauAnalyzerJobs.sh
./submitDataTauAnalyzerJobs.sh
./submitTTJetsTauAnalyzerJobs.sh
./submitWNJetsToLNuTauAnalyzerJobs.sh
./submitWJetsToLNuTauAnalyzerJobs.sh
./submitSingleTopTauAnalyzerJobs.sh
./submitWZTauAnalyzerJobs.sh
./submitWWTauAnalyzerJobs.sh
./submitZZTauAnalyzerJobs.sh
cd ..

exit 0
