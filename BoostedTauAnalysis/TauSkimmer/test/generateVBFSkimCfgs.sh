#!/bin/bash

if [ $# -ne 2 ]
    then
    echo "Usage: ./generateVBFSkimCfgs.sh <version> <template cfg>"
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
samples=( "VBF_a5" "VBF_a7" "VBF_a9" "VBF_a11" "VBF_a13" "VBF_a15" )

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
  if [ ${samples[${i}]} = "VBF_a9" ]
      then
      sed -e "s%    'root.*_9\.root',%%" -e "s%    'root.*_28\.root',%%" -e "s%    'root.*_61\.root',%%" -e "s%    'root.*_87\.root',%%" -e "s%    'root.*_115\.root',%%" -e "s%    'root.*_118\.root',%%" -e "s%    'root.*_131\.root',%%" -e "s%    'root.*_142\.root',%%" -e "s%    'root.*_146\.root',%%" -e "s%    'root.*_147\.root',%%" -e "s%    'root.*_158\.root',%%" -e "s%    'root.*_160\.root',%%" -e "s%    'root.*_161\.root',%%" -e "s%    'root.*_170\.root',%%" -e "s%    'root.*_189\.root',%%" -e "s%    'root.*_202\.root',%%" -e "s%    'root.*_219\.root',%%" -e "s%    'root.*_242\.root',%%" -e "s%    'root.*_283\.root',%%" -e "s%    'root.*_301\.root',%%" -e "s%    'root.*_371\.root',%%" -e "s%    'root.*_398\.root',%%" -e "s%    'root.*_406\.root',%%" -e "s%    'root.*_414\.root',%%" -e "s%    'root.*_434\.root',%%" -e "s%    'root.*_442\.root',%%" -e "s%    'root.*_445\.root',%%" -e "s%    'root.*_462\.root'%%" -e "s%    'root.*_483\.root',%%" -e "s%    'root.*_497\.root',%%" -e "s%    'root.*_503\.root',%%" -e "s%    'root.*_514\.root',%%" -e "s%    'root.*_526\.root',%%" -e "s%    'root.*_543\.root',%%" -e "s%    'root.*_565\.root',%%" -e "s%    'root.*_572\.root',%%" -e "s%    'root.*_579\.root',%%" -e "s%    'root.*_603\.root',%%" -e "s%    'root.*_608\.root',%%" -e "s%    'root.*_780\.root',%%" -e "s%    'root.*_832\.root',%%" -e "s%    'root.*_852\.root',%%" -e "s%    'root.*_860\.root',%%" -e "s%    'root.*_876\.root',%%" -e "s%    'root.*_888\.root',%%" -e "s%    'root.*_958\.root',%%" -e "s%    'root.*_992\.root',%%" -e "s%    'root.*_1000\.root',%%" -e "s%    'root.*_1015\.root',%%" -e "s%    'root.*_1019\.root',%%" -e "s%    'root.*_1021\.root',%%" -e "s%    'root.*_1086\.root',%%" -e "s%    'root.*_1089\.root',%%" -e "s%    'root.*_1092\.root',%%" -e "s%    'root.*_1104\.root',%%" -e "s%    'root.*_1106\.root',%%" -e "s%    'root.*_1118\.root',%%" -e "s%    'root.*_1130\.root',%%" -e "s%    'root.*_1148\.root',%%" -e "s%    'root.*_1156\.root',%%" -e "s%    'root.*_1183\.root',%%" -e "s%    'root.*_1186\.root',%%" -e "s%    'root.*_1215\.root',%%" -e "s%    'root.*_1227\.root',%%" -e "s%    'root.*_1248\.root',%%" -e "s%    'root.*_1259\.root',%%" -e "s%    'root.*_1288\.root',%%" -e "s%    'root.*_1295\.root',%%" -e "s%    'root.*_1302\.root',%%" -e "s%    'root.*_1320\.root',%%" -e "s%    'root.*_1354\.root',%%" -e "s%    'root.*_617\.root',%%" -e "s%    'root.*_395\.root',%%" -e "s%SAMPLE%${samples[${i}]}%" -e "s%VERSION%${version}%g" -e "s%user\/yohay%group/phys_higgs/HiggsExo%" ../${templateCfg} > tauSelectionSkim_${samples[${i}]}.py
  fi
done

#generate run cfg that runs all skim files in the directory
cat <<EOF > runVBFSkimCfgs.sh
#!/bin/bash

for file in \`ls -alh *VBF*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file >& \$outFile &
done

exit 0
EOF
chmod a+x runVBFSkimCfgs.sh

exit 0
