#!/bin/bash

if [ $# -ne 2 ]
    then
    echo "Usage: ./generateggSkimCfgs.sh <version> <template cfg>"
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
samples=( "gg_a5" "gg_a7" "gg_a9" "gg_a11" "gg_a13" "gg_a15" )

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
  if [ ${samples[${i}]} = "gg_a5" ]
      then
      sed -e "s%    'root.*_150\.root',%%" -e "s%    'root.*_368\.root',%%" -e "s%    'root.*_41\.root',%%" -e "s%    'root.*_487\.root',%%" -e "s%    'root.*_[5-9][0-9][1-9]\.root'.*%%" -e "s%    'root.*_[5-9][1-9][0-9]\.root'.*%%" -e "s%    'root.*_[6-9][0-9][0-9]\.root'.*%%" -e "s%    'root.*_1000\.root'.*%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "gg_a7" ]
      then
      sed -e "s%    'root.*_324\.root',%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" -e "s%user\/yohay%group/phys_higgs%" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "gg_a9" ]
      then
      sed -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "gg_a11" ]
      then
      sed -e "s%    'root.*_18\.root',%%" -e "s%    'root.*_125\.root',%%" -e "s%    'root.*_235\.root',%%" -e "s%    'root.*_359\.root',%%" -e "s%    'root.*_373\.root'.*%%" -e "s%    'root.*_436\.root'.*%%" -e "s%    'root.*_529\.root'.*%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" -e "s%user\/yohay%group/phys_higgs%" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "gg_a13" ]
      then
      sed -e "s%    'root.*_914\.root',%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "gg_a15" ]
      then
      sed -e "s%    'root.*_661\.root',%%" -e "s%    'root.*_761\.root',%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  else
      sed -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  fi
done

#generate run cfg that runs all skim files in the directory
cat <<EOF > runggSkimCfgs.sh
#!/bin/bash

for file in \`ls -alh *gg*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file >& \$outFile &
done

exit 0
EOF
chmod a+x runggSkimCfgs.sh

exit 0
