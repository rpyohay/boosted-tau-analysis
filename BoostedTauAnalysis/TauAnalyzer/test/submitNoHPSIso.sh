#!/bin/bash

cd v70
./submitDYJetsToLLAllTauAnalyzerJobs.sh
./submitDataAllTauAnalyzerJobs.sh
./submitTTJetsAllTauAnalyzerJobs.sh
./submitWNJetsToLNuAllTauAnalyzerJobs.sh
./submitWJetsToLNuAllTauAnalyzerJobs.sh
./submitSingleTopAllTauAnalyzerJobs.sh
./submitWZAllTauAnalyzerJobs.sh
./submitWWAllTauAnalyzerJobs.sh
./submitZZAllTauAnalyzerJobs.sh
cd ..

exit 0
