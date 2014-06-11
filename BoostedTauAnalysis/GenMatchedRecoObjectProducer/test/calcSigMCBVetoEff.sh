#!/bin/bash

if [ $# -ne 3 ]
    then
    echo "Usage: ./calcSigMCBVetoEff.sh <in_version> <out_version> <template cfg>"
    exit 0
fi

####STUFF TO CONFIGURE####

#version
inVersion=$1
outVersion=$2
templateCfg=$3
dir=$outVersion

#a1 mass
a1Mass=( "5" "7" "9" "11" "13" "15" )
nJobs=${#a1Mass[@]}
iBeg=0
iEndJob=`expr $nJobs - 1`

#input file directories
dirs=( "Wh1_Medium" "gg" )
nSamples=${#dirs[@]}
iEndSample=`expr $nSamples - 1`

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for iSample in `seq $iBeg $iEndSample`
  do

  #loop over number of jobs
  for iJob in `seq $iBeg $iEndJob`
    do

    #generate cfg file
    sample=${dirs[${iSample}]}
    if [ ${sample} = "Wh1_Medium" ]
	then
	sample="Wh1"
    fi
    sample="${sample}_a${a1Mass[${iJob}]}"
    sed -e "s%DIR%${dirs[${iSample}]}%" -e "s%MASS%${a1Mass[${iJob}]}%" -e "s%VERSION%${inVersion}%" -e "s%SAMPLE%${sample}%" ../${templateCfg} > calcSigMCBVetoEff_${sample}.py
  done
done

#generate run cfg that runs all files in the directory
cat <<EOF > runSigMCBVetoEffCfgs.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for file in \`ls -alh *.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file >& \$outFile &
done

exit 0
EOF
chmod a+x runSigMCBVetoEffCfgs.sh

exit 0
