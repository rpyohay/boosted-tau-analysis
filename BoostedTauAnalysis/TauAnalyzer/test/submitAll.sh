#!/bin/bash

cd v51
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
