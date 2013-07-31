#!/bin/bash

#parse arguments
if [ $# -ne 3 ]
    then
    echo "Usage: ./generateSkimCfgs.sh sample first_file last_file"
    exit 0
fi
sample=$1
first_file=$2
last_file=$3

#set up the input file name prefix and suffix
prefix="root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500_H2SMLike_${sample}/Summer12_NMSSMHiggs_"
suffix=".root"

#set up the output file name
output_sample_name="$sample"
if [ "$sample" = "v2" ]
    then
    output_sample_name="ZH"
fi
output_file_name="/data1/yohay/NMSSMHiggs_${output_sample_name}_files${first_file}-${last_file}_24Sep12.root"

#set up headers and trailers for the file arrays

file_array_header="readFiles\.extend(\[\n    "
file_array_trailer="\n    \])\n"

#loop over number of input files
input_file_list=""
for iFile in `seq ${first_file} ${last_file}`
  do

  #format input file name
  file_name="${prefix}${iFile}${suffix}"
  file_name=`echo $file_name | sed -e "s%/%\/%g" | sed -e "s%\.%\.%g"`

  #prepend header if necessary
  if [ `expr $iFile % 256` = "1" ] || [ "$iFile" = "${first_file}" ]
      then
      file_name="${file_array_header}\"${file_name}\""

      #append comma if necessary
      if [ "$iFile" != "${last_file}" ]
	  then
	  file_name="${file_name},\n    "
      fi
  fi

  #append trailer if necessary
  if [ `expr $iFile % 256` = "0" ] || [ "$iFile" = "${last_file}" ]
      then
      file_name="\"${file_name}\"${file_array_trailer}"
  fi

  #append comma if necessary
  last_file_in_section=`echo $file_name | grep \)`
  has_comma=`echo $file_name | grep ,`
  if [ -z "$has_comma" ] && [ -z "$last_file_in_section" ]
      then
      file_name="\"${file_name}\",\n    "
  fi

  #append file name to list
  input_file_list="${input_file_list}${file_name}"
done

sed -e "s%FILES%$input_file_list%" -e "s%OUTPUT_FILE_NAME%$output_file_name%" merge.py > merge_${output_sample_name}_files${first_file}-${last_file}.py

exit 0
