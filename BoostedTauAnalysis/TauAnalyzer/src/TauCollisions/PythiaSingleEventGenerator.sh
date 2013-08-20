#
#	PythiaSingleEventGenerator.sh
#
#	Michael Manighalam, Brenna Mockler, Daniel Coelho
#	August 7, 2013
#	Madgraph 5 version 1.5.11
#
#	Runs one madgraph process and one pythia process to produce showering
#
#	This Script uses create_root_macro.sh, MadGraph5, and ExRootAnalysis.
#
#	Run using the following command in bash: sh PythiaSingleEventGenerator.sh.
#	
#	What this script does:
#		1. Generates an .lhe file of the process specified in the variable MG5process with h01 mass and a01 mass specified in the variables mh01 and ma01. 
#		2. Converts the .lhe file to a .root file using ExRootAnalysis
#		3. Creates root macros (.C and .h files) to create histograms of some of the variables in the collision
#		4. Compiles the macros and creates the histograms in root
#	
#		***Specify paths to ExRootAnalysis and MadGraph in the code below***
#
#	Notes
#		-restrict_default.dat is what generates param_card.dat for each event, which is why we are editing it (we thought the param_card.dat in models/nmssm/ generated the param_card.dat for each event, but it's actually
#		 restrict_default.dat)
#		-if you want to use the cuts that we used, use our run_card.dat
#! /bin/bash


#Paths you need to specify:
pathExRoot=~/ExRootAnalysis 	#Path to the ExRootAnalysis directory
pathMG5=~/MadGraph5_v1_5_11	#Path to the MadGraph directory


mh01='115' #Define Higgs Mass for .C file
ma01='9' #Define a01 Mass for .C file
nevents='10000' #Define the number of events for MadGraph
fname="pythiaTau_Entries_$nevents" #Desired name of files 
MG5process="generate p p > w+ h01, w+ > e+ ve, (h01 > a01 a01, a01 > tau+ tau-)" #Define process for MadGraph
var1='1.15e+02' #Define Higgs Mass for param_card
var2='9.00e+00' #Define a01 Mass for param_card

#Don't change these. They are the Higgs mass, a01 mass, name of 
#files, MG5 process, and the number of events currently in restrict_default_template,
#run_card_template, process_card_template, and templateRootMacro.
hmass='25 1.15e+02'
amass='36 7.00e+00'
runEvents='10000 = nevents !'
originalFname='TAU_h115_a7'
originalMG5process='generate p p > w+ h01, w+ > e+ ve, (h01 > a01 a01, a01 > tau+ tau-)'
hmassC='115'
amassC='7'

	cp MadGraphCards/"restrict_default_template.dat" $pathMG5/models/nmssm/"restrict_default.dat"
	cp MadGraphCards/"run_card_template.dat" $pathMG5/Template/Cards/"run_card.dat"
	cp MadGraphCards/"process_card_template.dat" MadGraphCards/"process_card.dat"

     sed -i "s/$hmass/25 $var1/" $pathMG5/models/nmssm/"restrict_default.dat" #updates the h01 mass
     sed -i "s/$amass/36 $var2/" $pathMG5/models/nmssm/"restrict_default.dat" #updates the a01 mass
     sed -i "s/$runEvents/$nevents= nevents !/" $pathMG5/Template/Cards/"run_card.dat"
     sed -i "s/$originalFname/$fname/" MadGraphCards/"process_card.dat"
     sed -i "s/$originalMG5process/$MG5process/" MadGraphCards/"process_card.dat"
	$pathMG5/bin/mg5 MadGraphCards/"process_card.dat"
	
	
#Manual Overwrite	
	rm -r MadGraphEvents/$fname/ 
	mv $fname MadGraphEvents/
	cp madevent_interface_template.py MadGraphEvents/$fname/bin/internal/madevent_interface.py
	MadGraphEvents/$fname/bin/generate_events -f #generates pythia showering 	
cd ..
	gunzip -vf MadGraphEvents/$fname/Events/run_01/"tag_1_pythia_events.lhe.gz" #unzips the newly generated .lhe file
	cp MadGraphEvents/$fname/Events/run_01/"tag_1_pythia_events.lhe"	Output_singleEvent/"pythia_events_$fname.lhe" #renames the .lhe file to correspond with the h01 mass and a01 mass for the event
	
	#Converts the LHE file to a root file
	$pathExRoot/ExRootLHEFConverter	Output_singleEvent/"pythia_events_$fname.lhe" Output_singleEvent/"pythia_events_$fname.root"

	cd Output_singleEvent/
	root -l -q -x ../RootMacros/"create_root_macro.C(\"$fname\")" #Creates .C and .h files for each event (we really only do this for the .h files because we'll be generating our own .C files in Maccpy.sh). These .C and .h files will be located in Output_singleEvent/ (because we cd into there on the previous line)
	
#Now we generate the correct .C file from templateRootMacro and then run 
#root_macro_compiler to get a root file of histograms from the .C file
	 cp ../RootMacros/"pythiaTemplateRootMacro.C" $fname.C
    sed -i "s/pythiaTemplateRootMacro/$fname/" $fname.C
    sed -i "s/mh01(\"$hmassC\")/mh01(\"$mh01\")/" $fname.C
    sed -i "s/ma01(\"$amassC\")/ma01(\"$ma01\")/" $fname.C
    sed -i "s/Histograms_h01/Histograms_Pythia_Entries_""$nevents""_h01/" $fname.C

#We store all of the .C and .h files in C_h_files_singleEvent and we store all of the
#histograms and histogram root files in HistsAndHistRootFiles_singleEvent
	root -l -x -q ../RootMacros/"root_macro_compiler.C(\"$fname\")"
	cp "Histograms_Pythia_Entries_""$nevents""_h01_""$mh01""_a01_""$ma01"".root" ../HistsAndHistRootFiles_singleEvent
	rm "Histograms_Pythia_Entries_""$nevents""_h01_""$mh01""_a01_""$ma01"".root"
	cd ../

exit 0
