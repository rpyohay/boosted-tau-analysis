#!/bin/bash

cd v12
./copyDYJetsToLLFromEOS.sh
./copyDataFromEOS.sh
./copyTTJetsFromEOS.sh
./copyWNJetsToLNuFromEOS.sh
#./copyWJetsToLNuFromEOS.sh
./copySingleTopFromEOS.sh
./copyWZFromEOS.sh
./copyZZFromEOS.sh
./copyQCDFromEOS.sh
./copyQCDBFromEOS.sh
./copyQCDBMuFromEOS.sh
cd ..

exit 0
