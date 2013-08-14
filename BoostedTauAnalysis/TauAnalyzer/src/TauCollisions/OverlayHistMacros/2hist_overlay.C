#include <iostream>
#include <string>
#include <stdio.h>
#include <TStyle.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <Rtypes.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TPave.h>

void simple_overlay()
{ //open function bracket




const int num_hists = 43;
TH1F *first_histogram[num_hists];
TH1F *second_histogram[num_hists]; 
string first_hname[num_hists]= {
"Energy_Positron","Energy_Neutrino", "Energy_First_a01", "Energy_Second_a01", "Energy_Tau_1", "Energy_Tau_2", "Energy_Tau_3", "Energy_Tau_4",
"PT_Positron","PT_Neutrino", "PT_First_a01", "PT_Second_a01", "PT_Tau_1", "PT_Tau_2", "PT_Tau_3", "PT_Tau_4",
"Eta_Tau_1","Eta_Tau_2", "Eta_Tau_3", "Eta_Tau_4",
"DeltaEta_Tau_12", "DeltaEta_Tau_13", "DeltaEta_Tau_14", "DeltaEta_Tau_23", "DeltaEta_Tau_24", "DeltaEta_Tau_34",
"DeltaPhi_Tau_12", "DeltaPhi_Tau_13", "DeltaPhi_Tau_14", "DeltaPhi_Tau_23", "DeltaPhi_Tau_24", "DeltaPhi_Tau_34",
"DeltaR_Tau_12", "DeltaR_Tau_13", "DeltaR_Tau_14", "DeltaR_Tau_23", "DeltaR_Tau_24", "DeltaR_Tau_34",
"Invariant_Mass_Pair12", "Invariant_Mass_Pair14", "Invariant_Mass_Pair32", "Invariant_Mass_Pair34", "Invariant_Mass_4_Taus"};          // enter info for histograms here 
string second_hname[num_hists]= {
"energy_ele_hist","energy_nu_hist","energy_a1_hist","energy_a2_hist","energy_tau1_hist","energy_tau2_hist","energy_tau3_hist", "energy_tau4_hist",
"pT_ele_hist","pT_nu_hist","pT_a1_hist","pT_a2_hist","pT_tau1_hist","pT_tau2_hist","pT_tau3_hist", "pT_tau4_hist",
"eta_tau1", "eta_tau2", "eta_tau3", "eta_tau4",
"deleta_tau12", "deleta_tau13", "deleta_tau14", "deleta_tau23", "deleta_tau24", "deleta_tau34",
"delphi_tau12", "delphi_tau13", "delphi_tau14", "delphi_tau23", "delphi_tau24", "delphi_tau34",
"delR_tau12", "delR_tau13", "delR_tau14", "delR_tau23", "delR_tau24", "delR_tau34",
"mass_12", "mass_14", "mass_32", "mass_34", "mass_1234"};

TFile *first_file = new TFile("Entries_3000_Histograms_h01_115_a01_9.root"); //put name of one file here
TFile *second_file = new TFile("KinematicsPlots_NewBinning.root"); //put name of other file here


for (int i=0; i < num_hists; i++)
	{
	TIter next(first_file->GetListOfKeys());
	TKey *key;
	string hName;

		while((key=(TKey*)next()))
			{	
			hName = key->GetName();
	
			if  (hName.find(first_hname[i]) != -1) 
				{
				first_histogram[i] = (TH1F*)(first_file->Get(hName.c_str()));
				cout << "first_hname[" << i << "]" << first_hname[i] << endl;
				break;
				}		
			}

	TIter next(second_file->GetListOfKeys());

		while((key=(TKey*)next()))
			{
			hName = key->GetName();
		
			if  (hName.find(second_hname[i]) != -1) 
				{
				second_histogram[i] = (TH1F*)(second_file->Get(hName.c_str()));
			cout << "second_hname[" << i << "]" << first_hname[i] << endl;
				break;				
				}					
			}

	TCanvas *c1 = new TCanvas();
	double statXstart[5] = {.8,.8,.8,.8,.8};  
	double statXend[5] = {1.,1.,1.,1.,1.};
	double statYstart[5] = {.80,.70,.60,.50,.40};
	double statYend[5] = {.70,.60,.50,.40,.30};
	Color_t color_array[15] = {kBlue+1,kRed+0,kViolet,kAzure+10,kGreen+2,kBlue+1,kRed+0,kViolet,kAzure+10,kGreen+2,kBlue+1,kRed+0,kViolet,kAzure+10,kGreen+2};
	size_t name_length;
	TPaveText *title = new TPaveText(0.306092,0.9390678,0.693908,0.995,"blNDC");
	gStyle->SetOptTitle(0);
	title->SetBorderSize(0);
	title->SetFillColor(kWhite);
	title->SetTextFont(42);
	TPaveStats *st;	
	
		string directory = "Francesca_Overlaid_Tau_Histograms";
		gSystem->MakeDirectory(directory.c_str());
		directory = directory + "/";

		c1->cd();
		gStyle->SetOptStat("nmre");
		gStyle->SetStatH(0.03);
		gStyle->SetStatW(0.20);

		//TH1F *first_histogram = new TH1F(first_histogram_title, first_histogram_title, 200, 0., 500);
		//TH1F *second_histogram = new TH1F(second_histogram_title, second_histogram_title, 200, 0., 500);	

		first_histogram[i]->SetLineColor(color_array[0]);
		first_histogram[i]->SetLineWidth(2);
	 	first_histogram[i]->SetTitle(first_hname[i].c_str());	// WHAT DOES THIS DO?	eg. does it just set the name for hist title, or also for legend title?
		first_histogram[i]->SetName(first_hname[i].c_str()); // WHAT DOES THIS DO? eg. does it just set name for legend title?
		first_histogram[i]->Draw();
	
		gPad->Update();	
	
		st = (TPaveStats*)first_histogram[i]->FindObject("stats");
		st->SetX1NDC(statXstart[0]); 
		st->SetX2NDC(statXend[0]);
		st->SetY1NDC(statYstart[0]); 
		st->SetY2NDC(statYend[0]);


		second_histogram[i]->SetLineColor(color_array[1]);
		second_histogram[i]->SetLineWidth(2);
		second_histogram[i]->SetTitle(second_hname[i].c_str()); 
		second_histogram[i]->SetName(second_hname[i].c_str());
		second_histogram[i]->Draw("sames");
	
		gPad->Update();		

		st = (TPaveStats*)second_histogram[i]->FindObject("stats");		
		st->SetX1NDC(statXstart[1]); 
		st->SetX2NDC(statXend[1]);
		st->SetY1NDC(statYstart[1]); 
		st->SetY2NDC(statYend[1]);				
		
		c1->BuildLegend(.85,.8,1.,1.);
		c1->SetTitle("Dwarves");
		title->Clear(); //this clears the title for the histogram
		title->AddText((first_hname[i]).c_str()); //this adds text to the TPaveText box that we created before the loop
		title->Draw(); //this draws the TPaveText box with the new text... e.g. it makes our new title
	
		c1->SaveAs((directory + first_hname[i] + ".png").c_str());
	} //close for loop

} //close function bracket
