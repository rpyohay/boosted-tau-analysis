#!/bin/bash

cd v40
./copyDYJetsToLLFromEOS.sh
./copyDataFromEOS.sh
./copyTTJetsFromEOS.sh
./copyWNJetsToLNuFromEOS.sh
./copySingleTopFromEOS.sh
./copyWZFromEOS.sh
./copyZZFromEOS.sh
cd ..

exit 0
