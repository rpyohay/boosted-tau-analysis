/*
	root_macro_compiler.C
	Daniel Coelho, Brenna Mockler, Michael Manighalam
	version 2.0 08/07/2013	
	
	Used by macroCopy.sh	
	
	Helper function to compile root macros and make root files of the histograms.

*/

#include <iostream>
#include <string>

void root_macro_compiler(string macroName)
{    
	gROOT->ProcessLine((".L " + macroName + ".C++").c_str()); 
	gROOT->ProcessLine((macroName + " f").c_str());
	gROOT->ProcessLine("f->Loop();"); 
}
