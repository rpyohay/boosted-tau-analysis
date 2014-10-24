 #!/bin/bash

export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

#outputFileName="kinematics_bkg_sel"
#WNJetsToLNu
inputFiles=`cmsLs /store/user/friccita/effTest | grep root | wc -l`
echo $inputFiles

cmsLs /store/user/friccita/effTest | grep root | awk '{ print $5 }' > fileList.txt

p=1
startp=1
endp=(${inputFiles} + 1)
echo $startp
echo $endp

for linenumber in `seq ${startp} ${endp}`
    do
    file=""
    file=`sed -n "${linenumber} p" fileList.txt `
    if [ $linenumber -le $endp ]; then
	echo $file
	cmsRm $file
    fi
done
