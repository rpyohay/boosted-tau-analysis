#!/bin/bash

samples="WH gg"
for iSample in $samples
  do
  for iFirstFile in `seq 1 250 1000`
    do
    ./generateSkimCfgs.sh $iSample $iFirstFile `expr $iFirstFile + 249`
  done
done

exit 0
