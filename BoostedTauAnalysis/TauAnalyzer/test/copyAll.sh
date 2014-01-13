#!/bin/bash

if [ $# -ne 1 ]
    then
    echo "Usage: ./copyAll.sh <version>"
    exit 0
fi

version=$1

if [ ! -d $version ]
    then
    echo "Run ./generateJobFiles.sh before this script."
    exit 0
fi

cd $version
./copyDYJetsToLLFromEOS.sh
./copyDataFromEOS.sh
./copyTTJetsFromEOS.sh
./copyWJetsToLNuFromEOS.sh
./copyWNJetsToLNuFromEOS.sh
./copySingleTopFromEOS.sh
./copyWZFromEOS.sh
./copyWWFromEOS.sh
./copyZZFromEOS.sh
./copyQCDFromEOS.sh
./copyQCDBFromEOS.sh
./copyQCDBMuFromEOS.sh
./copyNonIsoWDataFromEOS.sh
cd ..

exit 0
