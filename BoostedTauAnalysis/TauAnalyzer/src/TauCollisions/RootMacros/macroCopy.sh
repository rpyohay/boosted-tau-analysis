#
# macroCopy.sh 
# Daniel Coelho, Brenna Mockler, Michael Manighalam
# version 2.0 08/07/2013
# 
# This program is used by MultiEventGenerator.sh and singleEventGenerator.sh, but can be run independently of either of them
#
# This program uses templateRootMacro.C, root_macro_compiler.C, overlay_histograms_a01.C, overlay_histograms_h01.C, 2hist_overlay.C
#
# To run this script on it's own, type sh macroCopy.sh
#
# This Program is designed to copy the base code from the template file (templateRootMacro.C) to all the other root macros with different 
# combinations of h01 and a01 masses(it can either create the macros or overwrite existing ones, depending on whether they have been made 
# previously). 
# 
# It then compiles the macros using load_root_macro.C and overlays and saves the histograms using overlay_histograms_a01.C,
# overlay_histograms_h01.C or 2hist_overlay.C.
#
#! /bin/bash

var1=templateRootMacro #This is the template file you need to edit to edit all of them
var2='ma01("7")'
var3='mh01("115")'
varArray1=(TAU_h115_a7 TAU_h115_a9 TAU_h115_a11 TAU_h115_a13 TAU_h115_a15 TAU_h120_a7 TAU_h120_a9 TAU_h120_a11 TAU_h120_a13 TAU_h120_a15 TAU_h125_a7 TAU_h125_a9 TAU_h125_a11 TAU_h125_a13 TAU_h125_a15 ) # Names of files that will be made (or updated) 
varArray2=('ma01("7")' 'ma01("9")' 'ma01("11")' 'ma01("13")' 'ma01("15")' 'ma01("7")' 'ma01("9")' 'ma01("11")' 'ma01("13")' 'ma01("15")' 'ma01("7")' 'ma01("9")' 'ma01("11")' 'ma01("13")' 'ma01("15")' )
varArray3=('mh01("115")' 'mh01("115")' 'mh01("115")' 'mh01("115")' 'mh01("115")' 'mh01("120")' 'mh01("120")' 'mh01("120")' 'mh01("120")' 'mh01("120")' 'mh01("125")' 'mh01("125")' 'mh01("125")' 'mh01("125")' 'mh01("125")' )
# Make same length as varArray1

echo "---------------------------------------"
echo "Starting Script With Given Params"
echo "String to change: $var1"  # Variable you want to change
echo "String to change: $var2"  # Variable you want to change
echo "String to change: $var3"  # Variable you want to change
echo "---------------------------------------"


echo "Enter Loop"   #This loop updates the macros using the template (templateRootMacro.C) and then compiles them    

cd Output_MultiEvent/

for index in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 # Number of file copies you want to make.
do
    cp ../RootMacros/"templateRootMacro.C" "${varArray1[index]}.C"   #this line and the following three update the macros using the template file
    sed -i "s/$var1/${varArray1[index]}/" ${varArray1[index]}.C
    sed -i "s/$var2/${varArray2[index]}/" ${varArray1[index]}.C
    sed -i "s/$var3/${varArray3[index]}/" ${varArray1[index]}.C
    root -l -q -x ../RootMacros/"root_macro_compiler.C(\"${varArray1[index]}\")" #compiles the macros and gives us root files of the histograms
done

mv Histograms_h01_* ../HistsAndHistRootFiles_MultiEvent
cd ../HistsAndHistRootFiles_MultiEvent
	#only use one of the three following lines	
	#root -l -x ../OverlayHistMacros/"overlay_histograms_a01.C()" #overlays the histograms with the same a01 mass 
	root -l -x -q ../OverlayHistMacros/"overlay_histograms_h01.C()" #overlays the histograms with the same h01 mass
	#root -l -x ../OverlayHistMacros/"2hist_overlay.C()" #overlays two histograms specified in 2hist_overlay.C TO USE THIS COMMENT OUT THE PREVIOUS LINES UP 
	cd -

exit 0
