#!/bin/bash

if [ $# -ne 1 ]
    then
    echo "Usage: ./copyNoHPSIso.sh <version>"
    exit 0
fi

version=$1

if [ ! -d $version ]
    then
    echo "Run ./generateJobFiles.sh before this script."
    exit 0
fi

cd $version
./copyAllDYJetsToLLFromEOS.sh
./copyAllDataFromEOS.sh
./copyAllTTJetsFromEOS.sh
#./copyAllWJetsToLNuFromEOS.sh
./copyAllWNJetsToLNuFromEOS.sh
#./copyAllWbbFromEOS.sh
./copyAllSingleTopFromEOS.sh
./copyAllWZFromEOS.sh
./copyAllWWFromEOS.sh
./copyAllZZFromEOS.sh
#./copyAllQCDFromEOS.sh
#./copyAllQCDBFromEOS.sh
#./copyAllQCDBMuFromEOS.sh
./copyAllNonIsoWDataFromEOS.sh
#./copyAllSinglePhotonDataFromEOS.sh
./copyAllNonIsoWDYJetsToLLFromEOS.sh
./copyAllNonIsoWTTJetsFromEOS.sh
./copyAllNonIsoWWNJetsToLNuFromEOS.sh
cd ..

exit 0
