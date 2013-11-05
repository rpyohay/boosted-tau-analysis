#!/bin/bash

version="v$1"

sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWh1TauAnalyzerCfg.sh > generateWh1TauAnalyzerCfg.sh.tmp && mv generateWh1TauAnalyzerCfg.sh.tmp generateWh1TauAnalyzerCfg.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateDYJetsToLLTauAnalyzerCfgs.sh > generateDYJetsToLLTauAnalyzerCfgs.sh.tmp && mv generateDYJetsToLLTauAnalyzerCfgs.sh.tmp generateDYJetsToLLTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateTTJetsTauAnalyzerCfgs.sh > generateTTJetsTauAnalyzerCfgs.sh.tmp && mv generateTTJetsTauAnalyzerCfgs.sh.tmp generateTTJetsTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWNJetsToLNuTauAnalyzerCfgs.sh > generateWNJetsToLNuTauAnalyzerCfgs.sh.tmp && mv generateWNJetsToLNuTauAnalyzerCfgs.sh.tmp generateWNJetsToLNuTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateSingleTopTauAnalyzerCfgs.sh > generateSingleTopTauAnalyzerCfgs.sh.tmp && mv generateSingleTopTauAnalyzerCfgs.sh.tmp generateSingleTopTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWZTauAnalyzerCfgs.sh > generateWZTauAnalyzerCfgs.sh.tmp && mv generateWZTauAnalyzerCfgs.sh.tmp generateWZTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateZZTauAnalyzerCfgs.sh > generateZZTauAnalyzerCfgs.sh.tmp && mv generateZZTauAnalyzerCfgs.sh.tmp generateZZTauAnalyzerCfgs.sh
<<<<<<< HEAD
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateQCDTauAnalyzerCfgs.sh > generateQCDTauAnalyzerCfgs.sh.tmp && mv generateQCDTauAnalyzerCfgs.sh.tmp generateQCDTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateQCDBTauAnalyzerCfgs.sh > generateQCDBTauAnalyzerCfgs.sh.tmp && mv generateQCDBTauAnalyzerCfgs.sh.tmp generateQCDBTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateQCDBMuTauAnalyzerCfgs.sh > generateQCDBMuTauAnalyzerCfgs.sh.tmp && mv generateQCDBMuTauAnalyzerCfgs.sh.tmp generateQCDBMuTauAnalyzerCfgs.sh
=======
>>>>>>> 46247d24eb20eea62a68284a461b2dc6bfa58c65
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateDataTauAnalyzerCfgs.sh > generateDataTauAnalyzerCfgs.sh.tmp && mv generateDataTauAnalyzerCfgs.sh.tmp generateDataTauAnalyzerCfgs.sh
sed -e "s%cd \(v.*\)%cd ${version}%" submitAll.sh > submitAll.sh.tmp && mv submitAll.sh.tmp submitAll.sh
sed -e "s%cd \(v.*\)%cd ${version}%" copyAll.sh > copyAll.sh.tmp && mv copyAll.sh.tmp copyAll.sh

<<<<<<< HEAD
chmod a+x generateWh1TauAnalyzerCfg.sh generateDYJetsToLLTauAnalyzerCfgs.sh generateTTJetsTauAnalyzerCfgs.sh generateWNJetsToLNuTauAnalyzerCfgs.sh generateSingleTopTauAnalyzerCfgs.sh generateWZTauAnalyzerCfgs.sh generateZZTauAnalyzerCfgs.sh generateQCDTauAnalyzerCfgs.sh generateDataTauAnalyzerCfgs.sh submitAll.sh copyAll.sh
=======
chmod a+x generateWh1TauAnalyzerCfg.sh generateDYJetsToLLTauAnalyzerCfgs.sh generateTTJetsTauAnalyzerCfgs.sh generateWNJetsToLNuTauAnalyzerCfgs.sh generateSingleTopTauAnalyzerCfgs.sh generateWZTauAnalyzerCfgs.sh generateZZTauAnalyzerCfgs.sh generateDataTauAnalyzerCfgs.sh submitAll.sh copyAll.sh
>>>>>>> 46247d24eb20eea62a68284a461b2dc6bfa58c65

exit 0
