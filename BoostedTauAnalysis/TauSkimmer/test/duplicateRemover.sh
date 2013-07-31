#!/bin/bash

export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

eosDirectory="/store/user/friccita/Testing"
filenameFragment="kinematics_bkg_sel"

#numFiles=`cmsLs $eosDirectory | grep $filenameFragment | wc -l`
numFiles=10
echo "Number of files with desired string in filename: $numFiles" 
cmsLs $eosDirectory | grep $filenameFragment | awk '{ print $5 }' > outputFileList.txt

startp=1
endp=$(expr ${numFiles} + 1)
echo $startp $endp

filesToRemoveList=""
#for file1 in `cmsLs $eosDirectory | grep $filenameFragment | awk '{ print $5 }'`
# do
#    echo $file1
#     for file2 in `cmsLs $eosDirectory | grep $filenameFragment | awk '{ print $5 }'`
#      do
#	if [ $file1 != $file2 ]
#	 then
#	    DIFF=`diff "$file1" "$file2" -q`
#	    if [ "${DIFF%% *}" != "Files" ]
#	     then
#		echo "Same: $file1 $file2"
#		echo "Must remove: $file1"
#		fileName="${file1}\n"
#		filesToRemoveList="${filesToRemoveList}${fileName}"
#		break
#	    fi #ifDIFF
#	fi #if files are not the same
#      done
# done

for linenumber1 in `seq ${startp} ${endp}`
 do
    file1=`sed -n "${linenumber1} p" outputFileList.txt`
    startp2=$(expr ${startp} + 1)
     for linenumber2 in `seq ${startp2} ${endp}`
      do
       file2=`sed -n "${linenumber2} p" outputFileList.txt`
	if [ $file1 != $file2 ]
	 then
	    DIFF=`diff "$file1" "$file2" -q`
	    if [ "${DIFF%% *}" != "Files" ]
	     then
		echo "Same: $file1 $file2"
		echo "Must remove: $file1"
		fileName="${file1}\n"
		filesToRemoveList="${filesToRemoveList}${fileName}"
		break
	    fi #ifDIFF
	fi #if files are not the same
      done
 done


exit 0
