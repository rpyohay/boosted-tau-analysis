#!/bin/bash

export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

outputFileName="kinematics_bkg_sel"
fDir="/afs/cern.ch/work/f/friccita"
finalFileName="merged.root"
inputFiles=`cmsLs /store/user/friccita/Testing/ | grep kinematics_bkg_sel | wc -l`
echo $inputFiles
cmsLs /store/user/friccita/Testing/ | grep kinematics_bkg_sel | awk '{ print $5 }' > kinFileList.txt

p=1
startp=1
endp=(${inputFiles} + 1)
echo $startp
echo $endp

for linenumber in `seq ${startp} ${endp}`
    do
    file=""
    file=`sed -n "${linenumber} p" kinFileList.txt `
    if [ $linenumber -le $endp ]; then
	echo $file
	cmsStage -f $file $fDir
    fi
done

hadd -f ${fDir}/${finalFileName} ${fDir}/${outputFileName}*.root
rm ${fDir}/${outputFileName}*.root

exit 0
