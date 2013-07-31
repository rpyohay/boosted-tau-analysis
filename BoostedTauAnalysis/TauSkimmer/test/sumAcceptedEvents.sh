#!/bin/bash

list=`grep -e "TrigReport Events" Summer12_DYToTauTau_skim/res/*.stdout | sed -e "s%.*passed = \([0-9]*\).*%\1%g"`
sum=0
for line in $list
  do
  sum=`expr $sum + $line`
done
echo "Sum: $sum"

exit 0
