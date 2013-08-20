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

/*
Program to overlay histograms with same a01 mass and differing h01 masses
*/

void overlay_histograms_a01()
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
	while((key=(TKey*)next())) //goes through keys and puts the ones whose names correspond to the histograms we want to plot into the array Histogram_Name
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
Color_t color_array[15] = {kBlue+1,kAzure+10,kRed+0,kGreen+2,kViolet,kBlue+1,kRed+0,kViolet,kAzure+10,kGreen+2,kBlue+1,kRed+0,kViolet,kAzure+10,kGreen+2};
size_t name_length;
TPaveText *title = new TPaveText(0.306092,0.9390678,0.693908,0.995,"blNDC");
gStyle->SetOptTitle(0);
title->SetBorderSize(0);
title->SetFillColor(kWhite);
title->SetTextFont(42);
TPaveStats *st;	
	string directory = "Overlaid_a01_Tau_Histograms";
	gSystem->MakeDirectory(directory.c_str());
	directory = directory + "/";

for(int j = 0; j < m; j++)

	{
//a01 7 GeV
	c1->cd();
	gStyle->SetOptStat("nmre"); //sets what we want in the stats boxes
	gStyle->SetStatH(0.03);
	gStyle->SetStatW(0.20);
	arrayHis[0] = (TH1F*)(histogram_file[0]->Get((Histogram_Name[0][j]).c_str()));
	arrayHis[0]->SetLineColor(kBlue+1);
	arrayHis[0]->SetLineWidth(2);
 	arrayHis[0]->SetTitle((Histogram_Name[0][j].substr(0,11)).c_str());		
	arrayHis[0]->SetName(Histogram_Name[0][j].substr(0,11).c_str()); 
	arrayHis[0]->Draw();
	gPad->Update();	
	st = (TPaveStats*)arrayHis[0]->FindObject("stats");
	st->SetX1NDC(statXstart[0]); 
     st->SetX2NDC(statXend[0]);
	st->SetY1NDC(statYstart[0]); 
     st->SetY2NDC(statYend[0]);
		
	for(int i = 5; i < 15; i = i+5)
		{
		arrayHis[i] = (TH1F*)(histogram_file[i]->Get((Histogram_Name[i][j]).c_str()));				
		arrayHis[i]->SetLineColor(color_array[(i/5)]);
		arrayHis[i]->SetLineWidth(2);
		arrayHis[i]->SetTitle((Histogram_Name[i][j].substr(0,11)).c_str()); //this sets labels for the legend and title for histogram as well, and because we are overlaying them it sets it to just one of them 
		arrayHis[i]->SetName((Histogram_Name[i][j].substr(0,11)).c_str());
		arrayHis[i]->Draw("sames");
		gPad->Update();		
		st = (TPaveStats*)arrayHis[i]->FindObject("stats");		
		st->SetX1NDC(statXstart[(i/5)]); 
		st->SetX2NDC(statXend[(i/5)]);
		st->SetY1NDC(statYstart[(i/5)]); 
		st->SetY2NDC(statYend[(i/5)]);				
		}
	c1->BuildLegend(.85,.8,1.,1.);
	c1->SetTitle("Elves");
	title->Clear(); //this clears the title for the histogram
	title->AddText(("Tau " + Histogram_Name[0][j].substr(11, Histogram_Name[0][j].length())).c_str()); //this adds text to the TPaveText box that we created before the loop
	title->Draw(); //this draws the TPaveText box with the new text... e.g. it makes our new title
	name_length = Histogram_Name[0][j].length();
	string newname = Histogram_Name[0][j].substr(11, Histogram_Name[0][j].length());  //file name
	c1->SaveAs((directory + newname + ".png").c_str());

//a01 9 GeV
	c1->cd();
	gStyle->SetOptStat("nmre");
	gStyle->SetStatH(0.03);
	gStyle->SetStatW(0.20);
	arrayHis[1] = (TH1F*)(histogram_file[1]->Get((Histogram_Name[1][j]).c_str()));
	arrayHis[1]->SetLineColor(kBlue+1);
	arrayHis[1]->SetLineWidth(2);
 	arrayHis[1]->SetTitle((Histogram_Name[1][j].substr(0,11)).c_str());
	arrayHis[1]->SetName((Histogram_Name[1][j].substr(0,11)).c_str());
	arrayHis[1]->Draw();
	gPad->Update();	
	st = (TPaveStats*)arrayHis[1]->FindObject("stats");	
	st->SetX1NDC(statXstart[0]); 
	st->SetX2NDC(statXend[0]);
	st->SetY1NDC(statYstart[0]); 
	st->SetY2NDC(statYend[0]);	
	
	
	for(int i = 6; i < 15; i=i+5)
		{
		arrayHis[i] = (TH1F*)(histogram_file[i]->Get((Histogram_Name[i][j]).c_str()));				
		arrayHis[i]->SetLineColor(color_array[(i-1)/5]);
		arrayHis[i]->SetLineWidth(2);
		arrayHis[i]->SetTitle((Histogram_Name[i][j].substr(0,11)).c_str()); //this sets title for histogram as well, and because we are overlaying them it sets it to just one of them 
		arrayHis[i]->SetName((Histogram_Name[i][j].substr(0,11)).c_str());
		arrayHis[i]->Draw("sames");
		gPad->Update();		
		st = (TPaveStats*)arrayHis[i]->FindObject("stats");		
		st->SetX1NDC(statXstart[(i-1)/5]); 
		st->SetX2NDC(statXend[(i-1)/5]);
		st->SetY1NDC(statYstart[(i-1)/5]); 
		st->SetY2NDC(statYend[(i-1)/5]);		
		}
	c1->BuildLegend(.85,.8,1.,1.);
	title->Clear();
	title->AddText(("Tau " + Histogram_Name[1][j].substr(11, Histogram_Name[1][j].length())).c_str());
	title->Draw();
	name_length = Histogram_Name[1][j].length();
	string newname = Histogram_Name[1][j].substr(11, Histogram_Name[1][j].length());	
	c1->SaveAs((directory + newname + ".png").c_str());

//a01 11
	c1->cd();
	gStyle->SetOptStat("nmre");
	gStyle->SetStatH(0.03);
	gStyle->SetStatW(0.20);
	arrayHis[2] = (TH1F*)(histogram_file[2]->Get((Histogram_Name[2][j]).c_str()));
	arrayHis[2]->SetLineColor(kBlue+1);
	arrayHis[2]->SetLineWidth(2);
 	arrayHis[2]->SetTitle((Histogram_Name[2][j].substr(0,11)).c_str());
	arrayHis[2]->SetName((Histogram_Name[2][j].substr(0,11)).c_str());
	arrayHis[2]->Draw();
	gPad->Update();
	st = (TPaveStats*)arrayHis[2]->FindObject("stats");		
	st->SetX1NDC(statXstart[0]); 
	st->SetX2NDC(statXend[0]);
	st->SetY1NDC(statYstart[0]); 
	st->SetY2NDC(statYend[0]);	
	
	for(int i = 7; i < num_hist_files; i=i+5)
		{
		arrayHis[i] = (TH1F*)(histogram_file[i]->Get((Histogram_Name[i][j]).c_str()));				
		arrayHis[i]->SetLineColor(color_array[(i-2)/5]);
		arrayHis[i]->SetLineWidth(2);
		arrayHis[i]->SetTitle((Histogram_Name[i][j].substr(0,11)).c_str()); //this sets title for histogram as well, and because we are overlaying them it sets it to just one of them 
		arrayHis[i]->SetName((Histogram_Name[i][j].substr(0,11)).c_str());
		arrayHis[i]->Draw("sames");
		gPad->Update();		
		st = (TPaveStats*)arrayHis[i]->FindObject("stats");		
		st->SetX1NDC(statXstart[(i-2)/5]); 
		st->SetX2NDC(statXend[(i-2)/5]);
		st->SetY1NDC(statYstart[(i-2)/5]); 
		st->SetY2NDC(statYend[(i-2)/5]);		
		}
	c1->BuildLegend(.85,.8,1.,1.);
	title->Clear();
	title->AddText(("Tau " + Histogram_Name[2][j].substr(11, Histogram_Name[0][j].length())).c_str());
	title->Draw();
	name_length = Histogram_Name[2][j].length();
	string newname = Histogram_Name[2][j].substr(11, Histogram_Name[2][j].length());
	c1->SaveAs((directory + newname + ".png").c_str());

//a01 13	
	c1->cd();
	gStyle->SetOptStat("nmre");
	gStyle->SetStatH(0.03);
	gStyle->SetStatW(0.20);
	arrayHis[3] = (TH1F*)(histogram_file[3]->Get((Histogram_Name[3][j]).c_str()));
	arrayHis[3]->SetLineColor(kBlue+1);
	arrayHis[3]->SetLineWidth(2);
 	arrayHis[3]->SetTitle((Histogram_Name[3][j].substr(0,11)).c_str());		
	arrayHis[3]->SetName(Histogram_Name[3][j].substr(0,11).c_str()); 
	arrayHis[3]->Draw();
	gPad->Update();	
	st = (TPaveStats*)arrayHis[3]->FindObject("stats");
	st->SetX1NDC(statXstart[0]); 
     st->SetX2NDC(statXend[0]);
	st->SetY1NDC(statYstart[0]); 
     st->SetY2NDC(statYend[0]);
	for(int i = 8; i < 15; i=i+5)
		{
		arrayHis[i] = (TH1F*)(histogram_file[i]->Get((Histogram_Name[i][j]).c_str()));				
		arrayHis[i]->SetLineColor(color_array[(i-3)/5]);
		arrayHis[i]->SetLineWidth(2);
		arrayHis[i]->SetTitle((Histogram_Name[i][j].substr(0,11)).c_str()); //this sets labels for the legend and title for histogram as well, and because we are overlaying them it sets it to just one of them 
		arrayHis[i]->SetName((Histogram_Name[i][j].substr(0,11)).c_str());
		arrayHis[i]->Draw("sames");
		gPad->Update();		
		st = (TPaveStats*)arrayHis[i]->FindObject("stats");		
		st->SetX1NDC(statXstart[(i-3)/5]); 
		st->SetX2NDC(statXend[(i-3)/5]);
		st->SetY1NDC(statYstart[(i-3)/5]); 
		st->SetY2NDC(statYend[(i-3)/5]);				
		}
	c1->BuildLegend(.85,.8,1.,1.);
	c1->SetTitle("Elves");
	title->Clear(); //this clears the title for the histogram
	title->AddText(("Tau " + Histogram_Name[3][j].substr(11, Histogram_Name[3][j].length())).c_str()); //this adds text to the TPaveText box that we created before the loop
	title->Draw(); //this draws the TPaveText box with the new text... e.g. it makes our new title
	name_length = Histogram_Name[3][j].length();
	string newname = Histogram_Name[3][j].substr(11, Histogram_Name[3][j].length());	
	c1->SaveAs((directory + newname + ".png").c_str());

//a01 15

c1->cd();
	gStyle->SetOptStat("nmre");
	gStyle->SetStatH(0.03);
	gStyle->SetStatW(0.20);
	arrayHis[4] = (TH1F*)(histogram_file[4]->Get((Histogram_Name[4][j]).c_str()));
	arrayHis[4]->SetLineColor(kBlue+1);
	arrayHis[4]->SetLineWidth(2);
 	arrayHis[4]->SetTitle((Histogram_Name[4][j].substr(0,11)).c_str());		
	arrayHis[4]->SetName(Histogram_Name[4][j].substr(0,11).c_str());
	arrayHis[4]->Draw();
	gPad->Update();	
	st = (TPaveStats*)arrayHis[4]->FindObject("stats");
	st->SetX1NDC(statXstart[0]); 
     st->SetX2NDC(statXend[0]);
	st->SetY1NDC(statYstart[0]); 
     st->SetY2NDC(statYend[0]);
	for(int i = 9; i < 15; i=i+5)
		{
		arrayHis[i] = (TH1F*)(histogram_file[i]->Get((Histogram_Name[i][j]).c_str()));				
		arrayHis[i]->SetLineColor(color_array[(i-4)/5]);
		arrayHis[i]->SetLineWidth(2);
		arrayHis[i]->SetTitle((Histogram_Name[i][j].substr(0,11)).c_str()); //this sets labels for the legend and title for histogram as well, and because we are overlaying them it sets it to just one of them 
		arrayHis[i]->SetName((Histogram_Name[i][j].substr(0,11)).c_str());
		arrayHis[i]->Draw("sames");
		gPad->Update();		
		st = (TPaveStats*)arrayHis[i]->FindObject("stats");		
		st->SetX1NDC(statXstart[(i-4)/5]); 
		st->SetX2NDC(statXend[(i-4)/5]);
		st->SetY1NDC(statYstart[(i-4)/5]); 
		st->SetY2NDC(statYend[(i-4)/5]);				
		}
	c1->BuildLegend(.85,.8,1.,1.);
	c1->SetTitle("Elves");
	title->Clear(); //this clears the title for the histogram
	title->AddText(("Tau " + Histogram_Name[4][j].substr(11, Histogram_Name[4][j].length())).c_str()); //this adds text to the TPaveText box that we created before the loop
	title->Draw(); //this draws the TPaveText box with the new text... e.g. it makes our new title
	name_length = Histogram_Name[4][j].length();
	string newname = Histogram_Name[4][j].substr(11, Histogram_Name[4][j].length());	//name of file
	c1->SaveAs((directory + newname + ".png").c_str());

	

}


} //close function bracket

