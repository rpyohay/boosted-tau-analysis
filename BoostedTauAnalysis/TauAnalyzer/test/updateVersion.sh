#!/bin/bash

sed -e "s%version=\"\(.*\)\"%version=\"v35\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWh1TauAnalyzerCfg.sh > generateWh1TauAnalyzerCfg.sh.tmp && mv generateWh1TauAnalyzerCfg.sh.tmp generateWh1TauAnalyzerCfg.sh
sed -e "s%version=\"\(.*\)\"%version=\"v35\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateDYJetsToLLTauAnalyzerCfgs.sh > generateDYJetsToLLTauAnalyzerCfgs.sh.tmp && mv generateDYJetsToLLTauAnalyzerCfgs.sh.tmp generateDYJetsToLLTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"v35\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateTTJetsTauAnalyzerCfgs.sh > generateTTJetsTauAnalyzerCfgs.sh.tmp && mv generateTTJetsTauAnalyzerCfgs.sh.tmp generateTTJetsTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"v35\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWNJetsToLNuTauAnalyzerCfgs.sh > generateWNJetsToLNuTauAnalyzerCfgs.sh.tmp && mv generateWNJetsToLNuTauAnalyzerCfgs.sh.tmp generateWNJetsToLNuTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"v35\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateDataTauAnalyzerCfgs.sh > generateDataTauAnalyzerCfgs.sh.tmp && mv generateDataTauAnalyzerCfgs.sh.tmp generateDataTauAnalyzerCfgs.sh

exit 0
