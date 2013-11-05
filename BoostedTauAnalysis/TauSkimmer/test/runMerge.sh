#!/bin/bash

#arguments
nFilesPerJob=$1
EOSDir=$2
queue=$3
startJob=$4

#set up CMSSW environment so cmsLs can be used
eval `scramv1 runtime -sh`

#macros
sample=`echo $EOSDir | sed -e "s%/store/user/yohay/\(.*\)%\1%" | sed -e "s%/%-%g"`
dir=$sample
prefix="'root://eoscms//eos/cms"
nFiles=`cmsLs $EOSDir | grep root | awk '{ print $5 }' | wc -l`
echo $nFiles

#make working directory
mkdir -p $dir
cd $dir

#loop over number of input files
fileCounter=1
inputFileList=""
nJobs=$startJob
for iFile in `cmsLs $EOSDir | grep root | awk '{ print $5 }'`
  do

  #format file name
  fileName="${prefix}${iFile}'"
  if [ `expr $fileCounter % $nFilesPerJob` != "0" ] && [ "$fileCounter" != "$nFiles" ]
      then

      #append comma
      fileName="${fileName},\n    "
      inputFileList="${inputFileList}${fileName}"

  #create cfg and sh files if you're at the max number of files per job
  else
      inputFileList="${inputFileList}${fileName}"
      outputFileName="${sample}_${nJobs}.root"
#      sed -e "s%FILES%$inputFileList%" -e "s%OUTPUT_FILE_NAME%$outputFileName%" ../merge.py > merge_${nJobs}.py
#      sed -e "s%DIR%$dir%" -e "s%CFG_FILE_NAME%merge_${nJobs}.py%" -e "s%OUTPUT_FILE_NAME%$outputFileName%" ../merge.sh > merge_${nJobs}.sh
#      chmod a+x merge_${nJobs}.sh
#      bsub -q $queue -J merge_${sample}_${nJobs} < merge_${nJobs}.sh
      inputFileList=""
      nJobs=`expr $nJobs + 1`
  fi
  fileCounter=`expr $fileCounter + 1`
done

exit 0
