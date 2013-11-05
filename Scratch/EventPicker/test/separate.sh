#!/bin/bash

num_lines=`cat /data1/yohay/Wh1_muon_jet_analysis.txt | wc -l`
for line in `seq 1 $num_lines`
  do
  string=`head -n $line /data1/yohay/Wh1_muon_jet_analysis.txt | tail -n 1`
  hist=`echo $string | sed -e "s%\([0-9]*\) [0-9]* [0-9]* [0-9]* [0-9]*%\1%"`
  x_bin=`echo $string | sed -e "s%[0-9]* \([0-9]*\) [0-9]* [0-9]* [0-9]*%\1%"`
  y_bin=`echo $string | sed -e "s%[0-9]* [0-9]* \([0-9]*\) [0-9]* [0-9]*%\1%"`
  run=`echo $string | sed -e "s%[0-9]* [0-9]* [0-9]* \([0-9]*\) [0-9]*%\1%"`
  evt=`echo $string | sed -e "s%[0-9]* [0-9]* [0-9]* [0-9]* \([0-9]*\)%\1%"`
  echo "'$run:$evt-$run:$evt'," | cat >> /data1/yohay/hist${hist}_xBin${x_bin}_yBin${y_bin}.txt
done

eval `scramv1 runtime -sh`

for iHist in `seq 0 3`
  do
  for iXBin in `seq 0 1`
    do
    for iYBin in `seq 0 1`
      do
      evts=`cat /data1/yohay/hist${iHist}_xBin${iXBin}_yBin${iYBin}_formatted.txt | tr '\n' ' '`
      sed -e "s%EVENTS%$evts%" -e "s%HIST%$iHist%" -e "s%XBIN%$iXBin%" -e "s%YBIN%$iYBin%" pickEvents.py > pickEvents_hist${iHist}_xBin${iXBin}_yBin${iYBin}.py
      cmsRun pickEvents_hist${iHist}_xBin${iXBin}_yBin${iYBin}.py
      #rm pickEvents_hist${iHist}_xBin${iXBin}_yBin${iYBin}.py /data1/yohay/hist${iHist}_xBin${iXBin}_yBin${iYBin}_formatted.txt /data1/yohay/hist${iHist}_xBin${iXBin}_yBin${iYBin}.txt
    done
  done
done

exit 0
