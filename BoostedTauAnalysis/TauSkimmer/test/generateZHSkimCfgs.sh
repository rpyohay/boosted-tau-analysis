#!/bin/bash

if [ $# -ne 2 ]
    then
    echo "Usage: ./generateZHSkimCfgs.sh <version> <template cfg>"
    exit 0
fi

####STUFF TO CONFIGURE####

#version
version=$1
templateCfg=$2
dir=$version

#number of samples
nSamples=6
iBeg=0
iEnd=`expr $nSamples - 1`

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#samples
samples=( "ZH_a5" "ZH_a7" "ZH_a9" "ZH_a11" "ZH_a13" "ZH_a15" )

#a1 masses
masses=( "a5" "a7" "a9" "a11" "a13" "a15" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

    #generate cfg file
  if [ ${samples[${i}]} = "ZH_a9" ]
      then
      sed -e "s%    'root.*_422\.root',%%" -e "s%    'root.*_7\.root',%%" -e "s%    'root.*_223\.root',%%" -e "s%    'root.*_54\.root',%%" -e "s%    'root.*_370\.root',%%" -e "s%    'root.*_523\.root',%%" -e "s%    'root.*_822\.root',%%" -e "s%    'root.*_117\.root',%%" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" -e "s%user\/yohay%group/phys_higgs/HiggsExo%" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  fi
done

#generate run cfg that runs all skim files in the directory
cat <<EOF > runZHSkimCfgs.sh
#!/bin/bash

for file in \`ls -alh *ZH*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file >& \$outFile &
done

exit 0
EOF
chmod a+x runZHSkimCfgs.sh

exit 0
