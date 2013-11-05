#!/bin/bash

cd v9
./copyDYJetsToLLFromEOS.sh
./copyDataFromEOS.sh
./copyTTJetsFromEOS.sh
./copyWNJetsToLNuFromEOS.sh
./copySingleTopFromEOS.sh
./copyWZFromEOS.sh
./copyZZFromEOS.sh
./copyQCDFromEOS.sh
./copyQCDBFromEOS.sh
./copyQCDBMuFromEOS.sh
cd ..

exit 0
