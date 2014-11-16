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
      sed -e "s%    'root.*_119\.root',%%" -e "s%    'root.*_1055\.root',%%" -e "s%    'root.*_244\.root',%%" -e "s%    'root.*_245\.root',%%" -e "s%    'root.*_803\.root',%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "gg_a7" ]
      then
      sed -e "s%    'root.*_30\.root',%%" -e "s%    'root.*_74\.root',%%" -e "s%    'root.*_250\.root',%%" -e "s%    'root.*_385\.root',%%" -e "s%    'root.*_443\.root',%%" -e "s%    'root.*_524\.root',%%" -e "s%    'root.*_578\.root',%%" -e "s%    'root.*_651\.root',%%" -e "s%    'root.*_872\.root',%%" -e "s%    'root.*_946\.root',%%" -e "s%    'root.*_1076\.root',%%" -e "s%    'root.*_1181\.root',%%" -e "s%    'root.*_1275\.root',%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "gg_a9" ]
      then
      sed -e "s%    'root.*_36\.root',%%" -e "s%    'root.*_245\.root',%%" -e "s%    'root.*_714\.root',%%" -e "s%    'root.*_723\.root',%%" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "gg_a11" ]
      then
      sed -e "s%    'root.*_303\.root',%%" -e "s%    'root.*_380\.root',%%" -e "s%    'root.*_499\.root',%%" -e "s%    'root.*_48\.root',%%" -e "s%    'root.*_80\.root'.*%%" -e "s%    'root.*_286\.root'.*%%" -e "s%    'root.*_468\.root'.*%%" -e "s%    'root.*_497\.root'.*%%" -e "s%    'root.*_517\.root'.*%%" -e "s%    'root.*_591\.root'.*%%" -e "s%    'root.*_634\.root'.*%%" -e "s%    'root.*_748\.root'.*%%" -e "s%    'root.*_799\.root'.*%%" -e "s%    'root.*_863\.root'.*%%" -e "s%    'root.*_949\.root'.*%%" -e "s%    'root.*_1024\.root'.*%%" -e "s%    'root.*_1128\.root'.*%%" -e "s%    'root.*_1207\.root'.*%%" -e "s%    'root.*_1250\.root'.*%%" -e "s%    'root.*_1287\.root'.*%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" -e "s%user\/yohay%group/phys_higgs/HiggsExo%" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "gg_a13" ]
      then
      sed -e "s%    'root.*_973\.root',%%" -e "s%    'root.*_624\.root',%%" -e "s%    'root.*_336\.root',%%" -e "s%    'root.*_3\.root',%%" -e "s%    'root.*_127\.root',%%" -e "s%    'root.*_135\.root',%%" -e "s%    'root.*_332\.root',%%" -e "s%    'root.*_970\.root',%%" -e "s%    'root.*_1127\.root',%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" -e "s%user\/yohay%group/phys_higgs/HiggsExo%" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  elif [ ${samples[${i}]} = "gg_a15" ]
      then
      sed -e "s%    'root.*_46\.root',%%" -e "s%    'root.*_51\.root',%%" -e "s%    'root.*_136\.root',%%" -e "s%    'root.*_140\.root',%%" -e "s%    'root.*_259\.root',%%" -e "s%    'root.*_381\.root',%%" -e "s%    'root.*_489\.root',%%" -e "s%    'root.*_637\.root',%%" -e "s%    'root.*_668\.root',%%" -e "s%    'root.*_894\.root',%%" -e "s%    'root.*_1053\.root',%%" -e "s%a9%${masses[${i}]}%g" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" -e "s%user\/yohay%group/phys_higgs/HiggsExo%" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
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
