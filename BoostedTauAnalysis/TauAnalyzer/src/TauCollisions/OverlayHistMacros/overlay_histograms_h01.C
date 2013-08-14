/*
	overlay_histograms_h01.C
	Daniel Coelho, Brenna Mockler, Michael Manighalam
	version 2.0 08/08/13
	Program to overlay histograms with same h01 mass and differing a01 masses
*/

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

void overlay_histograms_h01()
{ //open function bracket

int num_a01 = 7;
int num_h01 = 115;
int k = 0;
const int num_hist_files = 15;
const int num_variables = 8; //number of variables we want for the physical quantities 
int m = 0;


TFile *histogram_file[num_hist_files];
TH1F *histogram_variables[num_variables];//[physical_quantity][num_variables] --> we currently aren't using this variable because we don't actually need to save the original histograms again
string Histogram_Name[15][100];
stringstream ss;
stringstream tt;



for (int j = 0; j <= 2; j++)   
{		
	for (int i = 0; i <= 4; i++)   
	{
	ss << num_a01;
	tt << num_h01;
	histogram_file[k] = new TFile((("Histograms_h01_" + tt.str() + "_a01_" + ss.str() + ".root").c_str()));  // Insert Your Filename you want to put in
	ss.str("");
	tt.str("");
	num_a01 = num_a01 + 2;
	
	m = 0;
	TIter next(histogram_file[k]->GetListOfKeys());
	TKey *key;
	string hName;
	while((key=(TKey*)next()))
		{
		hName = key->GetName();
		if ( ((hName.find("Transverse_Energy") == -1) && (hName.find("Energy") != -1))
			|| (hName.find("Eta_Tau") != -1) 
			|| (hName.find("Invariant_Mass") != -1) 
			|| (hName.find("PT") != -1) 
			|| ((hName.find("a01s") == -1) && (hName.find("Delta") != -1))
 			|| (hName.find("Largest_Tau_Pt") != -1) )
			{
			Histogram_Name[k][m] = hName;
			//cout <<"Histogram_Name[" << k << "]["<< m <<"] = " << Histogram_Name[k][m] << endl;
			m++; //m goes from 0 to however many different variables we are plotting (eg. Transverse Energy, Delta_Eta, Delta_R)
			}
		}
	k++; //k goes from 0 to 15 and corresponds with each mass combination of the a01's and h01s (eg. "h01 115, a01 7")
	}
num_h01 = num_h01 + 5;
num_a01 = 7;
}


TCanvas *c1 = new TCanvas();
TH1F *arrayHis[num_hist_files];
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
	string directory = "Overlaid_h01_Tau_Histograms";
	gSystem->MakeDirectory(directory.c_str());
	directory = directory + "/";

for(int j = 0; j < m; j++)
	{
//115 GeV
	c1->cd();
	gStyle->SetOptStat("nmre");
	gStyle->SetStatH(0.03);
	gStyle->SetStatW(0.20);
	arrayHis[0] = (TH1F*)(histogram_file[0]->Get((Histogram_Name[0][j]).c_str()));
	arrayHis[0]->SetLineColor(kBlue+1);
	arrayHis[0]->SetLineWidth(2);
 	arrayHis[0]->SetTitle((Histogram_Name[0][j].substr(11,11)).c_str());		
	arrayHis[0]->SetName(Histogram_Name[0][j].substr(11,11).c_str()); // name of histogram (Histogram_Name[0][j].substr(0,10) + Histogram_Name[0][j].substr(21, Histogram_Name[0][j].length())).c_str()	
	arrayHis[0]->Draw();
	gPad->Update();	
	st = (TPaveStats*)arrayHis[0]->FindObject("stats");
	st->SetX1NDC(statXstart[0]); 
     st->SetX2NDC(statXend[0]);
	st->SetY1NDC(statYstart[0]); 
     st->SetY2NDC(statYend[0]);
	for(int i = 1; i < 5; i++)
		{
		arrayHis[i] = (TH1F*)(histogram_file[i]->Get((Histogram_Name[i][j]).c_str()));				
		arrayHis[i]->SetLineColor(color_array[i]);
		arrayHis[i]->SetLineWidth(2);
		arrayHis[i]->SetTitle((Histogram_Name[i][j].substr(11,11)).c_str()); //this sets labels for the legend and title for histogram as well, and because we are overlaying them it sets it to just one of them (the "h01 115, a01 7" one)
		arrayHis[i]->SetName((Histogram_Name[i][j].substr(11,11)).c_str());
		arrayHis[i]->Draw("sames");
		gPad->Update();		
		st = (TPaveStats*)arrayHis[i]->FindObject("stats");		
		st->SetX1NDC(statXstart[i]); 
		st->SetX2NDC(statXend[i]);
		st->SetY1NDC(statYstart[i]); 
		st->SetY2NDC(statYend[i]);				
		}
	c1->BuildLegend(.85,.8,1.,1.);
	c1->SetTitle("Elves");
	title->Clear(); //this clears the title for the histogram
	title->AddText(("Tau " + Histogram_Name[0][j].substr(0,10) + Histogram_Name[0][j].substr(21, Histogram_Name[0][j].length())).c_str()); //this adds text to the TPaveText box that we created before the loop
	title->Draw(); //this draws the TPaveText box with the new text... e.g. it makes our new title
	name_length = Histogram_Name[0][j].length();
	string newname = Histogram_Name[0][j].substr(0,10) + Histogram_Name[0][j].substr(21, Histogram_Name[0][j].length());	
	c1->SaveAs((directory + newname + ".png").c_str());

//120 GeV
	c1->cd();
	gStyle->SetOptStat("nmre");
	gStyle->SetStatH(0.03);
	gStyle->SetStatW(0.20);
	arrayHis[5] = (TH1F*)(histogram_file[5]->Get((Histogram_Name[5][j]).c_str()));
	arrayHis[5]->SetLineColor(kBlue+1);
	arrayHis[5]->SetLineWidth(2);
 	arrayHis[5]->SetTitle((Histogram_Name[5][j].substr(11,11)).c_str());
	arrayHis[5]->SetName((Histogram_Name[5][j].substr(11,11)).c_str());
	arrayHis[5]->Draw();
	gPad->Update();	
	st = (TPaveStats*)arrayHis[5]->FindObject("stats");	
	st->SetX1NDC(statXstart[0]); 
	st->SetX2NDC(statXend[0]);
	st->SetY1NDC(statYstart[0]); 
	st->SetY2NDC(statYend[0]);	
	
	
	for(int i = 6; i < 10; i++)
		{
		arrayHis[i] = (TH1F*)(histogram_file[i]->Get((Histogram_Name[i][j]).c_str()));				
		arrayHis[i]->SetLineColor(color_array[i]);
		arrayHis[i]->SetLineWidth(2);
		arrayHis[i]->SetTitle((Histogram_Name[i][j].substr(11,11)).c_str()); //this sets title for histogram as well, and because we are overlaying them it sets it to just one of them (the "h01 115, a01 7" one)
		arrayHis[i]->SetName((Histogram_Name[i][j].substr(11,11)).c_str());
		arrayHis[i]->Draw("sames");
		gPad->Update();		
		st = (TPaveStats*)arrayHis[i]->FindObject("stats");		
		st->SetX1NDC(statXstart[i-5]); 
		st->SetX2NDC(statXend[i-5]);
		st->SetY1NDC(statYstart[i-5]); 
		st->SetY2NDC(statYend[i-5]);		
		}
	c1->BuildLegend(.85,.8,1.,1.);
	title->Clear();
	title->AddText(("Tau " + Histogram_Name[5][j].substr(0,10) + Histogram_Name[0][j].substr(21, Histogram_Name[0][j].length())).c_str());
	title->Draw();
	name_length = Histogram_Name[5][j].length();
	string newname = Histogram_Name[5][j].substr(0,10) + Histogram_Name[0][j].substr(21, Histogram_Name[0][j].length());	
	c1->SaveAs((directory + newname + ".png").c_str());

//125 GeV
	c1->cd();
	gStyle->SetOptStat("nmre");
	gStyle->SetStatH(0.03);
	gStyle->SetStatW(0.20);
	arrayHis[10] = (TH1F*)(histogram_file[10]->Get((Histogram_Name[10][j]).c_str()));
	arrayHis[10]->SetLineColor(kBlue+1);
	arrayHis[10]->SetLineWidth(2);
 	arrayHis[10]->SetTitle((Histogram_Name[10][j].substr(11,11)).c_str());
	arrayHis[10]->SetName((Histogram_Name[10][j].substr(11,11)).c_str());
	arrayHis[10]->Draw();
	gPad->Update();
	st = (TPaveStats*)arrayHis[10]->FindObject("stats");		
	st->SetX1NDC(statXstart[0]); 
	st->SetX2NDC(statXend[0]);
	st->SetY1NDC(statYstart[0]); 
	st->SetY2NDC(statYend[0]);	
	
	for(int i = 11; i < num_hist_files; i++)
		{
		arrayHis[i] = (TH1F*)(histogram_file[i]->Get((Histogram_Name[i][j]).c_str()));				
		arrayHis[i]->SetLineColor(color_array[i]);
		arrayHis[i]->SetLineWidth(2);
		arrayHis[i]->SetTitle((Histogram_Name[i][j].substr(11,11)).c_str()); //this sets title for histogram as well, and because we are overlaying them it sets it to just one of them (the "h01 115, a01 7" one)
		arrayHis[i]->SetName((Histogram_Name[i][j].substr(11,11)).c_str());
		arrayHis[i]->Draw("sames");
		gPad->Update();		
		st = (TPaveStats*)arrayHis[i]->FindObject("stats");		
		st->SetX1NDC(statXstart[i-10]); 
		st->SetX2NDC(statXend[i-10]);
		st->SetY1NDC(statYstart[i-10]); 
		st->SetY2NDC(statYend[i-10]);		
		}
	c1->BuildLegend(.85,.8,1.,1.);
	title->Clear();
	title->AddText(("Tau " + Histogram_Name[10][j].substr(0,10) + Histogram_Name[0][j].substr(21, Histogram_Name[0][j].length())).c_str());
	title->Draw();
	name_length = Histogram_Name[10][j].length();
	string newname = Histogram_Name[10][j].substr(0,10) + Histogram_Name[0][j].substr(21, Histogram_Name[0][j].length());
	c1->SaveAs((directory + newname + ".png").c_str());
	}

} //close function bracket

