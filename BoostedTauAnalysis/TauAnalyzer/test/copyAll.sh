#!/bin/bash

cd v51
./copyDYJetsToLLFromEOS.sh
./copyDataFromEOS.sh
./copyTTJetsFromEOS.sh
./copyWJetsToLNuFromEOS.sh
./copyWNJetsToLNuFromEOS.sh
./copySingleTopFromEOS.sh
./copyWZFromEOS.sh
./copyWWFromEOS.sh
./copyZZFromEOS.sh
cd ..

exit 0
