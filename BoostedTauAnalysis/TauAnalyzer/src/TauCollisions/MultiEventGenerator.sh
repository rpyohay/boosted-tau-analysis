#	
#	MultiEventGenerator.sh
#	
#	Michael Manighalam, Brenna Mockler, Daniel Coelho
#	August 7, 2013
#	Madgraph 5 version 1.5.11
#	
#	This Script uses macroCopy.sh, create_root_macro.sh, MadGraph5, and ExRootAnalysis.
#
#	Run using the following command in bash: sh multiEventGenerator.sh.
#	
#	What this script does:
#		1. Generates 15 .lhe files of the process (specified below) with h01 mass of 115, 120, and 125 GeV, and a01 mass of 7, 9, 11, 13, and 15 GeV.
#		2. Converts the .lhe files to .root files using ExRootAnalysis
#		3. Overlays different quantities of the processes using root and saves them as .png files (macroCopy.sh does this part)
#	
#	Process Generated: p p > h01 w+, w+ > e+ ve (h01 > a01 a01, a01 >tau+ tau-)
#	
#		***Specify paths to ExRootAnalysis and MadGraph in the code below***
#
#	Notes
#		-restrict_default.dat is what generates param_card.dat for each event, which is why we are editing it (we thought the param_card.dat in models/nmssm/ generated the param_card.dat for each event, but it's actually
#		 restrict_default.dat)
#		-if you want to use the cuts that we used, use our run_card.dat
#
#Don't delete the next line
#! /bin/bash


#Paths you may need to specify:
pathExRoot=~/ExRootAnalysis 	#Path to the ExRootAnalysis directory
pathMG5=~/MadGraph5_v1_5_11	#Path to the MadGraph directory

varArray1=('25 1.15e+02' '25 1.15e+02' '25 1.15e+02' '25 1.15e+02' '25 1.15e+02'
			  '25 1.20e+02' '25 1.20e+02' '25 1.20e+02' '25 1.20e+02' '25 1.20e+02'
			  '25 1.25e+02' '25 1.25e+02' '25 1.25e+02' '25 1.25e+02' '25 1.25e+02') 	 #Corresponds to the h01 mass defined in param_card.dat that will be updated each run

varArray2=('36 7.00e+00' '36 9.00e+00' '36 11.00e+00' '36 13.00e+00' '36 15.00e+00'
			  '36 7.00e+00' '36 9.00e+00' '36 11.00e+00' '36 13.00e+00' '36 15.00e+00'
			  '36 7.00e+00' '36 9.00e+00' '36 11.00e+00' '36 13.00e+00' '36 15.00e+00') #Corresponds to the a01 mass defined in param_card.dat that will be updated for each run

varArray3=('TAU_h115_a7' 'TAU_h115_a9' 'TAU_h115_a11' 'TAU_h115_a13' 'TAU_h115_a15'
			  'TAU_h120_a7' 'TAU_h120_a9' 'TAU_h120_a11' 'TAU_h120_a13' 'TAU_h120_a15' 
			  'TAU_h125_a7' 'TAU_h125_a9' 'TAU_h125_a11' 'TAU_h125_a13' 'TAU_h125_a15') #Names of the .lhe and .root files

MG5process='generate p p > w+ h01, w+ > e+ ve, (h01 > a01 a01, a01 > tau+ tau-)' #Define MG5 process

#Don't Change this, it is the process that is currently in process_card_template.
originalMG5process='generate p p > w+ h01, w+ > e+ ve, (h01 > a01 a01, a01 > tau+ tau-)'

#The following loop edits param_card.dat and the process_card.dat, generates 
for index in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 # Number of file copies you want to make.
do
	#Copies the template param card (restrict_default_template.dat) and template process card (process_card_template.dat) so they can be edited
	cp MadGraphCards/"restrict_default_template.dat" $pathMG5/models/nmssm/"restrict_default.dat"
	cp MadGraphCards/"process_card_template.dat" MadGraphCards/"process_card.dat"

		sed -i "s/${varArray1[0]}/${varArray1[$index]}/" $pathMG5/models/nmssm/"restrict_default.dat" #Updates the h01 mass
		sed -i "s/${varArray2[0]}/${varArray2[$index]}/" $pathMG5/models/nmssm/"restrict_default.dat" #Updates the a01 mass
		sed -i "s/${varArray3[0]}/${varArray3[$index]}/" MadGraphCards/"process_card.dat" #Updates the outputed event directory name
      	sed -i "s/$originalMG5process/$MG5process/" MadGraphCards/"process_card.dat"

	$pathMG5/bin/mg5 MadGraphCards/"process_card.dat" #Runs the updated process card

#Manual overwrite
	rm -r MadGraphEvents/"${varArray3[$index]}"/ 
	mv "${varArray3[$index]}" MadGraphEvents/

	gunzip -vf MadGraphEvents/"${varArray3[$index]}"/Events/run_01/"events.lhe.gz" #unzips the newly generated .lhe file
	cp MadGraphEvents/${varArray3[$index]}/Events/run_01/"events.lhe"	 Output_MultiEvent/"events_${varArray3[$index]}.lhe" #renames the .lhe file to correspond with the h01 mass and a01 mass for the event
	
	#Converts the LHE file to a root file
	$pathExRoot/ExRootLHEFConverter	Output_MultiEvent/"events_${varArray3[$index]}.lhe" Output_MultiEvent/"events_${varArray3[$index]}.root" 

	cd Output_MultiEvent/
	root -l -q -x ../RootMacros/"create_root_macro.C(\"${varArray3[$index]}\")" #Creates .C and .h files for each event (we really only do this for the .h files because we'll be generating our own .C files in Maccpy.sh). These .C and .h files will be located in Output_MultiEvent/ (because we cd into there on the previous line)
	cd -
done

sh RootMacros/macroCopy.sh

exit 0
