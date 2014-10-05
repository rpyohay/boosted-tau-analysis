#!/bin/bash

#parse arguments
if [ $# -ne 7 ]
    then
    echo "Usage: ./generate.sh cfg_name script_name dir_name data_tier start_job num_jobs queue"
    exit 0
fi
cfg_name=$1
script_name=$2
dir_name=$3
data_tier=$4
start=$5
num_jobs=$6
queue=$7

#make directory on EOS
EOS_dir_query=`cmsLs /store/user/yohay/${dir_name}`
EOS_dir_query=`echo $EOS_dir_query | grep "No such file or directory"`
if [ "EOS_dir_query" != "" ]
    then
    cmsMkdir /store/user/yohay/${dir_name}
fi

#make local directory for holding all generated scripts and LSF output directories
mkdir -p $dir_name
cd $dir_name

#generate cfg and sh files for each job
for i in `seq $start $num_jobs`
  do
  sed -e "s%JOBNUM%$i%g" -e "s%DIRNAME%$dir_name%g" ../${cfg_name}.py > ${cfg_name}_${i}.py
  sed -e "s%JOBNUM%$i%g" -e "s%DIRNAME%$dir_name%g" -e "s%SCRIPTNAME%$cfg_name%g" -e "s%DATATIER%$data_tier%g" ../${script_name}.sh > ${cfg_name}_${i}.sh
  chmod a+x ${cfg_name}_${i}.sh
  bsub -q $queue -J ${cfg_name}_${i} < ${cfg_name}_${i}.sh
done

exit 0
