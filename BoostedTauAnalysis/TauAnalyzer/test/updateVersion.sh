#!/bin/bash

version="v$1"
reweightOnly=0
if [ "$2" == "reweightOnly" ]
    then
    reweightOnly=1
fi

sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateWh1TauAnalyzerCfg.sh > generateWh1TauAnalyzerCfg.sh.tmp && mv generateWh1TauAnalyzerCfg.sh.tmp generateWh1TauAnalyzerCfg.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" -e "s%reweightOnly=.*%reweightOnly=${reweightOnly}%" generateDYJetsToLLTauAnalyzerCfgs.sh > generateDYJetsToLLTauAnalyzerCfgs.sh.tmp && mv generateDYJetsToLLTauAnalyzerCfgs.sh.tmp generateDYJetsToLLTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" -e "s%reweightOnly=.*%reweightOnly=${reweightOnly}%" generateTTJetsTauAnalyzerCfgs.sh > generateTTJetsTauAnalyzerCfgs.sh.tmp && mv generateTTJetsTauAnalyzerCfgs.sh.tmp generateTTJetsTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" -e "s%reweightOnly=.*%reweightOnly=${reweightOnly}%" generateWNJetsToLNuTauAnalyzerCfgs.sh > generateWNJetsToLNuTauAnalyzerCfgs.sh.tmp && mv generateWNJetsToLNuTauAnalyzerCfgs.sh.tmp generateWNJetsToLNuTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" -e "s%reweightOnly=.*%reweightOnly=${reweightOnly}%" generateWJetsToLNuTauAnalyzerCfgs.sh > generateWJetsToLNuTauAnalyzerCfgs.sh.tmp && mv generateWJetsToLNuTauAnalyzerCfgs.sh.tmp generateWJetsToLNuTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" -e "s%reweightOnly=.*%reweightOnly=${reweightOnly}%" generateSingleTopTauAnalyzerCfgs.sh > generateSingleTopTauAnalyzerCfgs.sh.tmp && mv generateSingleTopTauAnalyzerCfgs.sh.tmp generateSingleTopTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" -e "s%reweightOnly=.*%reweightOnly=${reweightOnly}%" generateWZTauAnalyzerCfgs.sh > generateWZTauAnalyzerCfgs.sh.tmp && mv generateWZTauAnalyzerCfgs.sh.tmp generateWZTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" -e "s%reweightOnly=.*%reweightOnly=${reweightOnly}%" generateWWTauAnalyzerCfgs.sh > generateWWTauAnalyzerCfgs.sh.tmp && mv generateWWTauAnalyzerCfgs.sh.tmp generateWWTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" -e "s%reweightOnly=.*%reweightOnly=${reweightOnly}%" generateZZTauAnalyzerCfgs.sh > generateZZTauAnalyzerCfgs.sh.tmp && mv generateZZTauAnalyzerCfgs.sh.tmp generateZZTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" -e "s%reweightOnly=.*%reweightOnly=${reweightOnly}%" generateDataTauAnalyzerCfgs.sh > generateDataTauAnalyzerCfgs.sh.tmp && mv generateDataTauAnalyzerCfgs.sh.tmp generateDataTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateQCDTauAnalyzerCfgs.sh > generateQCDTauAnalyzerCfgs.sh.tmp && mv generateQCDTauAnalyzerCfgs.sh.tmp generateQCDTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateQCDBTauAnalyzerCfgs.sh > generateQCDBTauAnalyzerCfgs.sh.tmp && mv generateQCDBTauAnalyzerCfgs.sh.tmp generateQCDBTauAnalyzerCfgs.sh
sed -e "s%version=\"\(.*\)\"%version=\"${version}\"%" -e "s%infoTag=\"\(.*\)\"%infoTag=\"\"%" generateQCDBMuTauAnalyzerCfgs.sh > generateQCDBMuTauAnalyzerCfgs.sh.tmp && mv generateQCDBMuTauAnalyzerCfgs.sh.tmp generateQCDBMuTauAnalyzerCfgs.sh
sed -e "s%cd \(v.*\)%cd ${version}%" submitAll.sh > submitAll.sh.tmp && mv submitAll.sh.tmp submitAll.sh
sed -e "s%cd \(v.*\)%cd ${version}%" submitNoHPSIso.sh > submitNoHPSIso.sh.tmp && mv submitNoHPSIso.sh.tmp submitNoHPSIso.sh
sed -e "s%cd \(v.*\)%cd ${version}%" copyAll.sh > copyAll.sh.tmp && mv copyAll.sh.tmp copyAll.sh
sed -e "s%cd \(v.*\)%cd ${version}%" copyNoHPSIso.sh > copyNoHPSIso.sh.tmp && mv copyNoHPSIso.sh.tmp copyNoHPSIso.sh

chmod a+x generateWh1TauAnalyzerCfg.sh generateDYJetsToLLTauAnalyzerCfgs.sh generateTTJetsTauAnalyzerCfgs.sh generateWNJetsToLNuTauAnalyzerCfgs.sh generateWJetsToLNuTauAnalyzerCfgs.sh generateSingleTopTauAnalyzerCfgs.sh generateWZTauAnalyzerCfgs.sh generateWWTauAnalyzerCfgs.sh generateZZTauAnalyzerCfgs.sh generateQCDTauAnalyzerCfgs.sh generateQCDBTauAnalyzerCfgs.sh generateQCDBMuTauAnalyzerCfgs.sh generateDataTauAnalyzerCfgs.sh submitAll.sh submitNoHPSIso.sh copyAll.sh copyNoHPSIso.sh

exit 0
