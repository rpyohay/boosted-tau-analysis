#!/bin/bash

#arguments
sample=$1

#macros
dir="${sample}"

#make working directory
mkdir -p $dir
cd $dir

#for background case
if [ "$sample" = "WJets" ]
    then

    #set up CMSSW environment so cmsLs can be used
    eval `scramv1 runtime -sh`

    #macros
    #nFilesPerJob=88
    #first_file=1
    #last_file=1000
    #prefix="'root://eoscms//eos/cms"
    #EOSDir="/store/user/friccita/WToMuNu_skim/"
    #file_list=`cmsLs $EOSDir | grep root | awk '{ print $5 }'`
    #identifier="${EOSDir}Summer12_WJetsToLNu_WMuNuSkim"
    nFilesPerJob=$2
    first_file=$3
    last_file=$4
    prefix="'root://eoscms//eos/cms"
    EOSDir=$5
    filePrefix=$6
    file_list=`cmsLs $EOSDir | grep root | awk '{ print $5 }'`
    identifier="${EOSDir}${filePrefix}"
    queue=$7
    startJob=$8

    #loop over number of input files
    input_file_list=""
    nJobs=$startJob
    for iFile in `seq ${first_file} ${last_file}`
      do

      #if proposed file exists...
      file=`echo $file_list | sed -e "s%.*\(${identifier}_${iFile}_[0-9]*_[0-9A-Za-z]*\.root\).*%\1%g"`
      file_list_mod=`echo $file_list | sed -e "s%\(.*\)%\1%g"`
      if [ "$file" != "$file_list_mod" ]
	  then
	  file_name="${prefix}${file}'"
	  file_name=`echo $file_name | sed -e "s%/%\/%g" | sed -e "s%\.%\.%g"`

          #format file name
	  if [ `expr $iFile % $nFilesPerJob` != "0" ] && [ "$iFile" != "$last_file" ]
	      then

              #append comma
	      file_name="${file_name},\n    "
	      input_file_list="${input_file_list}${file_name}"

          #create cfg and sh files if you're at the max number of files per job
	  else
	      input_file_list="${input_file_list}${file_name}"
	      outputFilePrefix="${dir}_tau_analysis_${nJobs}"
	      sed -e "s%FILES%$input_file_list%" -e "s%JOB%$nJobs%g" -e "s%OUTPUTROOTFILE%${outputFilePrefix}.root%" -e "s%OUTPUTTEXTFILE%${outputFilePrefix}.txt%" ../tauanalyzer_template_cfg.py > tauanalyzer_template_cfg_${nJobs}.py
	      sed -e "s%DIR%$dir%" -e "s%JOB%$nJobs%g" ../tauanalyzer_template_cfg.sh > tauanalyzer_template_cfg_${nJobs}.sh
	      chmod a+x tauanalyzer_template_cfg_${nJobs}.sh
	      bsub -q $queue -J tauanalyzer_template_cfg_${nJobs} < tauanalyzer_template_cfg_${nJobs}.sh
	      input_file_list=""
	      nJobs=`expr $nJobs + 1`
	  fi
      fi
    done
fi

#for signal case
if [ "$sample" = "Wh1" ]
    then
    input_file_list="'file:/data1/yohay/NMSSMHiggs_WH_files1-250_24Sep12.root',\n    'file:/data1/yohay/NMSSMHiggs_WH_files251-500_24Sep12.root',\n    'file:/data1/yohay/NMSSMHiggs_WH_files501-750_24Sep12.root',\n    'file:/data1/yohay/NMSSMHiggs_WH_files751-1000_24Sep12.root'"
    input_file_list=`echo $input_file_list | sed -e "s%/%\/%g" | sed -e "s%\.%\.%g"`
    nJobs=0
    outputFilePrefix="/data1/yohay/${dir}/${dir}_tau_analysis_${nJobs}"
    sed -e "s%FILES%$input_file_list%" -e "s%JOB%$nJobs%g" -e "s%OUTPUTROOTFILE%${outputFilePrefix}.root%" -e "s%OUTPUTTEXTFILE%${outputFilePrefix}.txt%" ../tauanalyzer_template_cfg.py > tauanalyzer_template_cfg_${nJobs}.py
fi

exit 0
