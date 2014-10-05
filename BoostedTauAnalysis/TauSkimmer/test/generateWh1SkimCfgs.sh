#!/bin/bash

if [ $# -ne 2 ]
    then
    echo "Usage: ./generateWh1SkimCfgs.sh <version> <template cfg>"
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
samples=( "Wh1_a5" "Wh1_a7" "Wh1_a9" "Wh1_a11" "Wh1_a13" "Wh1_a15" )

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
  if [ ${samples[${i}]} = "Wh1_a5" ]
      then
      sed -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "Wh1_a7" ]
      then
      sed -e "s%    'root.*_324\.root',%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" -e "s%user\/yohay%group/phys_higgs/HiggsExo%" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "Wh1_a9" ]
      then
      sed -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "Wh1_a11" ]
      then
      sed -e "s%    'root.*_18\.root',%%" -e "s%    'root.*_125\.root',%%" -e "s%    'root.*_235\.root',%%" -e "s%    'root.*_359\.root',%%" -e "s%    'root.*_373\.root'.*%%" -e "s%    'root.*_436\.root'.*%%" -e "s%    'root.*_529\.root'.*%%" -e "s%    'root.*_787\.root'.*%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" -e "s%user\/yohay%group/phys_higgs/HiggsExo%" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "Wh1_a13" ]
      then
      sed -e "s%    'root.*_914\.root',%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "Wh1_a15" ]
      then
      sed -e "s%    'root.*_661\.root',%%" -e "s%    'root.*_761\.root',%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  else
      sed -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  fi
done

#generate run cfg that runs all skim files in the directory
cat <<EOF > runWh1SkimCfgs.sh
#!/bin/bash

for file in \`ls -alh *Wh1*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file >& \$outFile &
done

exit 0
EOF
chmod a+x runWh1SkimCfgs.sh

exit 0
