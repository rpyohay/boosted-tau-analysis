#!/bin/bash

cd v70
./copyAllDYJetsToLLFromEOS.sh
./copyAllDataFromEOS.sh
./copyAllTTJetsFromEOS.sh
./copyAllWJetsToLNuFromEOS.sh
./copyAllWNJetsToLNuFromEOS.sh
./copyAllSingleTopFromEOS.sh
./copyAllWZFromEOS.sh
./copyAllWWFromEOS.sh
./copyAllZZFromEOS.sh
cd ..

exit 0
