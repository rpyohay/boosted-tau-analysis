/*
	create_root_macro.C
	Daniel Coelho, Brenna Mockler, Michael Manighalam
	version 2.0 08/08/13
	code for root commands to make root macros from event files. Used by TauEventGenerator.sh
*/

#include <iostream>
#include <string>


void create_root_macro(string varArray3)

{    
	gROOT->ProcessLine(("TFile *f = new TFile(\"events_" + varArray3 + ".root\")").c_str()); 
	gROOT->ProcessLine(("LHEF->MakeClass(\"" + varArray3 + "\")").c_str());
}
