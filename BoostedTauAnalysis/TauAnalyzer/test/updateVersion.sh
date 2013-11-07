#!/bin/bash

version="v$1"

sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWh1TauAnalyzerCfg.sh > generateWh1TauAnalyzerCfg.sh.tmp && mv generateWh1TauAnalyzerCfg.sh.tmp generateWh1TauAnalyzerCfg.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateDYJetsToLLTauAnalyzerCfgs.sh > generateDYJetsToLLTauAnalyzerCfgs.sh.tmp && mv generateDYJetsToLLTauAnalyzerCfgs.sh.tmp generateDYJetsToLLTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateTTJetsTauAnalyzerCfgs.sh > generateTTJetsTauAnalyzerCfgs.sh.tmp && mv generateTTJetsTauAnalyzerCfgs.sh.tmp generateTTJetsTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWNJetsToLNuTauAnalyzerCfgs.sh > generateWNJetsToLNuTauAnalyzerCfgs.sh.tmp && mv generateWNJetsToLNuTauAnalyzerCfgs.sh.tmp generateWNJetsToLNuTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWJetsToLNuTauAnalyzerCfgs.sh > generateWJetsToLNuTauAnalyzerCfgs.sh.tmp && mv generateWJetsToLNuTauAnalyzerCfgs.sh.tmp generateWJetsToLNuTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateSingleTopTauAnalyzerCfgs.sh > generateSingleTopTauAnalyzerCfgs.sh.tmp && mv generateSingleTopTauAnalyzerCfgs.sh.tmp generateSingleTopTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWZTauAnalyzerCfgs.sh > generateWZTauAnalyzerCfgs.sh.tmp && mv generateWZTauAnalyzerCfgs.sh.tmp generateWZTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWWTauAnalyzerCfgs.sh > generateWWTauAnalyzerCfgs.sh.tmp && mv generateWWTauAnalyzerCfgs.sh.tmp generateWWTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateZZTauAnalyzerCfgs.sh > generateZZTauAnalyzerCfgs.sh.tmp && mv generateZZTauAnalyzerCfgs.sh.tmp generateZZTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateDataTauAnalyzerCfgs.sh > generateDataTauAnalyzerCfgs.sh.tmp && mv generateDataTauAnalyzerCfgs.sh.tmp generateDataTauAnalyzerCfgs.sh
sed -e "s%cd \(v.*\)%cd ${version}%" submitAll.sh > submitAll.sh.tmp && mv submitAll.sh.tmp submitAll.sh
sed -e "s%cd \(v.*\)%cd ${version}%" submitNoHPSIso.sh > submitNoHPSIso.sh.tmp && mv submitNoHPSIso.sh.tmp submitNoHPSIso.sh
sed -e "s%cd \(v.*\)%cd ${version}%" copyAll.sh > copyAll.sh.tmp && mv copyAll.sh.tmp copyAll.sh
sed -e "s%cd \(v.*\)%cd ${version}%" copyNoHPSIso.sh > copyNoHPSIso.sh.tmp && mv copyNoHPSIso.sh.tmp copyNoHPSIso.sh

chmod a+x generateWh1TauAnalyzerCfg.sh generateDYJetsToLLTauAnalyzerCfgs.sh generateTTJetsTauAnalyzerCfgs.sh generateWNJetsToLNuTauAnalyzerCfgs.sh generateWJetsToLNuTauAnalyzerCfgs.sh generateSingleTopTauAnalyzerCfgs.sh generateWZTauAnalyzerCfgs.sh generateWWTauAnalyzerCfgs.sh generateZZTauAnalyzerCfgs.sh generateDataTauAnalyzerCfgs.sh submitAll.sh submitNoHPSIso.sh copyAll.sh copyNoHPSIso.sh

exit 0
