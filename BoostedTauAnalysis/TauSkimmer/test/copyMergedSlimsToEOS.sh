#!/bin/bash

eval `scramv1 runtime -sh`
samples="WH gg"
for iSample in $samples
  do
  for iFirstFile in `seq 1 250 1000`
    do
    cmsStage -f /data1/yohay/NMSSMHiggs_${iSample}_files${iFirstFile}-`expr $iFirstFile + 249`_24Sep12.root /store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500_H2SMLike_${iSample}/
  done
done

exit 0
