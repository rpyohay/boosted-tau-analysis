#!/bin/bash

cd v36
./copyDYJetsToLLFromEOS.sh
./copyDataFromEOS.sh
./copyTTJetsFromEOS.sh
./copyWNJetsToLNuFromEOS.sh
cd ..

exit 0
