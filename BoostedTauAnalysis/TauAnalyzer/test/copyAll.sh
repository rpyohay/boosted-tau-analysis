#!/bin/bash

cd v35
./copyDYJetsToLLFromEOS.sh
./copyDataFromEOS.sh
./copyTTJetsFromEOS.sh
./copyWNJetsToLNuFromEOS.sh
cd ..

exit 0
