#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TProfile.h"
#include "TF1.h"
#include "TKey.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooBinning.h"
#include "Error.C"

//default drawing options
const Double_t defaultXAxisLabelSize = 0.05;
const Int_t defaultCanvasWidth = 600;
const Int_t defaultCanvasHeight = 600;

//error codes
enum ErrorCode {SUCCESS = 0, CANNOT_RETRIEVE_OBJECT};

//error: file cannot be opened
string errorCannotOpenFile(const string& fnName, const string& fileName)
{
  return ("Error in " + fnName + ":\nCannot open file \"" + fileName + "\".\n");
}

//error: object could not be retrieved from file
string errorCannotRetrieveObject(const string& fnName, const string& fileName, 
				 const string& objName)
{
  return ("Error in " + fnName + ":\nCannot retrieve object \"" + objName + "\" from file \"" + 
	  fileName + "\".\n");
}

//add trailing / to save path if missing
void formatSavePath(string& savePath)
{
  if (savePath.find('/') != (savePath.length() - 1)) savePath+="/";
}

//set canvas drawing options
void setCanvasOptions(TVirtualPad& canvas, const Int_t grid, const Int_t logY, const Int_t logZ)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
}

//set canvas drawing options
void setCanvasOptions(TCanvas& canvas, const Int_t grid, const Int_t logY, const Int_t logZ)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
}

//set canvas margins
void setCanvasMargins(TVirtualPad& canvas, const float left, const float top, const float right, 
		      const float bottom)
{
  canvas.cd()->SetLeftMargin(left);
  canvas.cd()->SetTopMargin(top);
  canvas.cd()->SetRightMargin(right);
  canvas.cd()->SetBottomMargin(bottom);
}

//set canvas margins
void setCanvasMargins(TCanvas& canvas, const float left, const float top, const float right, 
		      const float bottom)
{
  canvas.cd()->SetLeftMargin(left);
  canvas.cd()->SetTopMargin(top);
  canvas.cd()->SetRightMargin(right);
  canvas.cd()->SetBottomMargin(bottom);
}

//set axis drawing options
void setAxisOptions(TAxis* axis, const Float_t labelSize, const Float_t titleOffset, 
		    const char* title)
{
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetTitleSize(0.06);
  axis->SetTitleOffset(titleOffset);
  axis->SetTitle(title);
}

//set axis drawing options
void setAxisOptions(TAxis* axis, const Float_t labelSize, const Float_t titleSize, 
		    const Float_t titleOffset, const char* title)
{
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetTitleSize(titleSize);
  axis->SetTitleOffset(titleOffset);
  axis->SetTitle(title);
}

//set axis drawing options
void setAxisOptions(TAxis* axis, const Float_t labelSize, const Float_t titleSize, 
		    const Float_t titleOffset)
{
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetTitleSize(titleSize);
  axis->SetTitleOffset(titleOffset);
}

//set axis labels
void setAxisLabels(TAxis* axis, const vector<string>& binLabels)
{
  for (Int_t iBin = 1; iBin <= axis->GetNbins(); ++iBin) {
    if (iBin <= (int)binLabels.size()) axis->SetBinLabel(iBin, binLabels[iBin - 1].c_str());
  }
}

//set graph drawing options
void setGraphOptions(TGraphAsymmErrors& graph, const Color_t color, const Size_t size, 
		     const Style_t style, const char* xAxisTitle, const char* yAxisTitle)
{
  graph.SetMarkerColor(color);
  graph.SetMarkerSize(size);
  graph.SetMarkerStyle(style);
  graph.SetLineColor(color);
  graph.SetTitle("");
  setAxisOptions(graph.GetXaxis(), defaultXAxisLabelSize, 0.9, xAxisTitle);
  setAxisOptions(graph.GetYaxis(), defaultXAxisLabelSize, 1.05, yAxisTitle);
}

//set 1D histogram drawing options
void setHistogramOptions(TH1F* histogram, const Color_t color, const Size_t size, 
			 const Style_t style, const Double_t scale, const char* xAxisTitle, 
			 const char* yAxisTitle, 
			 const Double_t xAxisLabelSize = defaultXAxisLabelSize)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(2);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), xAxisLabelSize, 0.9, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), defaultXAxisLabelSize, 1.05, yAxisTitle);
  histogram->Scale(scale);
}

//set 2D histogram options
void setHistogramOptions(TH2F* histogram, const Color_t color, const Size_t size, 
			 const Style_t style, const Float_t yAxisTitleOffset, 
			 const Double_t scale, const char* xAxisTitle, const char* yAxisTitle)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), 0.04, 0.04, 1.1, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), 0.04, 0.04, yAxisTitleOffset, yAxisTitle);
  setAxisOptions(histogram->GetZaxis(), 0.04, 0.04, histogram->GetZaxis()->GetTitleOffset(), "");
  histogram->Scale(scale);
}

//set 2D histogram options
void setHistogramOptions(TH2F* histogram, const Color_t color, const Size_t size, 
			 const Style_t style, const Float_t yAxisTitleOffset, const Double_t scale)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), 0.04, 0.04, 1.1);
  setAxisOptions(histogram->GetYaxis(), 0.04, 0.04, yAxisTitleOffset);
  setAxisOptions(histogram->GetZaxis(), 0.04, 0.04, histogram->GetZaxis()->GetTitleOffset(), "");
  histogram->Scale(scale);
}

//set legend drawing options
void setLegendOptions(TLegend& legend, const char* header)
{
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextFont(42);
  legend.SetHeader(header);
}


//for agreement plotting/assessment:
//set bin content equal to sqrt(bin^2 + sigma_bin^2)
TH1F* denomErrorScale(TH1F* histInput)
{
  TH1F* hist = (TH1F*)histInput->Clone();
  for (int i = 1; i < (hist->GetNbinsX() - 2); ++i)
    {
      double bin = hist->GetBinContent(i);
      double binErr = hist->GetBinError(i);
      double newContent = sqrt((bin*bin) + (binErr*binErr));
      hist->SetBinContent(i, newContent);
    }
  return hist;
}

//make efficiency plots from input TH1s and save them to the output file
ErrorCode plotEfficiency(TFile& in, TFile& out, const map<string, pair<string, string> >& 
			 effHistMap, const map<string, vector<string> >& binLabelMap, 
			 string savePath)
{
  string fnName("ErrorCode plotEfficiency(const TH1F& in, TH1F& out, const map<string, ");
  fnName+="pair<string, string> >& effHistMap, const map<string, vector<string> >& binLabelMap, ";
  fnName+="string savePath)";

  //loop over efficiency histogram map
  for (map<string, pair<string, string> >::const_iterator iEffHist = effHistMap.begin(); 
       iEffHist != effHistMap.end(); ++iEffHist) {

    //get histograms from file
    TH1F* numeratorHist = NULL;
    TH1F* denominatorHist = NULL;
    in.GetObject(iEffHist->first.c_str(), numeratorHist);
    in.GetObject(iEffHist->second.first.c_str(), denominatorHist);
    if ((numeratorHist == NULL) || (denominatorHist == NULL)) {
      cerr << errorCannotRetrieveObject(fnName, in.GetName(), iEffHist->first + "\" or \"" + 
					iEffHist->second.first);
      return CANNOT_RETRIEVE_OBJECT;
    }

    //set axis labels and titles
    Double_t xAxisLabelSize = defaultXAxisLabelSize;
    map<string, vector<string> >::const_iterator iHistLabels = 
      binLabelMap.find(iEffHist->second.first);
    if (iHistLabels != binLabelMap.end()) {
      setAxisLabels(numeratorHist->GetXaxis(), iHistLabels->second);
      setAxisLabels(denominatorHist->GetXaxis(), iHistLabels->second);
      xAxisLabelSize = 0.04;
    }
    setHistogramOptions(numeratorHist, kBlack, 0.7, 20, 1.0, iEffHist->second.second.c_str(), 
			"#epsilon", xAxisLabelSize);
    setHistogramOptions(denominatorHist, kBlack, 0.7, 20, 1.0, iEffHist->second.second.c_str(), 
			"#epsilon", xAxisLabelSize);
    numeratorHist->GetYaxis()->SetRangeUser(0.0, 1.1);
    denominatorHist->GetYaxis()->SetRangeUser(0.0, 1.1);

//     //rebin
//     cerr << "Rebinning\n";
//     numeratorHist->Rebin(2);
//     denominatorHist->Rebin(2);

    //make efficiency histogram
    TGraphAsymmErrors effGraph(numeratorHist, denominatorHist);
    setGraphOptions(effGraph, kBlack, 0.7, 20, iEffHist->second.second.c_str(), "#epsilon");
    effGraph.GetYaxis()->SetRangeUser(0.0, 1.1);

    //draw efficiency histogram
    out.cd();
    string effCanvasName("eff_" + iEffHist->first + "_over_" + iEffHist->second.first);
    Int_t canvasWidth = defaultCanvasWidth;
    Int_t canvasHeight = defaultCanvasHeight;
    if (iHistLabels != binLabelMap.end()) {
      canvasWidth = 1200;
      canvasHeight = 600;
    }
    TCanvas effCanvas(effCanvasName.c_str(), "", canvasWidth, canvasHeight);
    setCanvasOptions(effCanvas, 1, 0, 0);
    if (iHistLabels != binLabelMap.end()) effCanvas.cd()->SetRightMargin(0.1);
    numeratorHist->Draw("AXIS");
    numeratorHist->Draw("AXIGSAME");
    effGraph.Draw("PSAME");
    effCanvas.Write();

    //save PDF of efficiency plot
    if (savePath != "noPDF") {
      formatSavePath(savePath);
      effCanvas.SaveAs((savePath + effCanvasName + ".pdf").c_str());
    }
  }

  //success
  return SUCCESS;
}

//make efficiency plots from input TH2s and save them to the output file
ErrorCode plot2DEfficiency(TFile& in, TFile& out, 
			   const map<pair<string, string>, pair<string, string> >& effHistMap, 
			   string savePath)
{
  string fnName("ErrorCode plot2DEfficiency(const TH1F& in, TH1F& out, const map<pair<string, ");
  fnName+="string, pair<string, string> >& effHistMap, string savePath)";

  //loop over efficiency histogram map
  for (map<pair<string, string>, pair<string, string> >::const_iterator iEffHist = 
	 effHistMap.begin(); iEffHist != effHistMap.end(); ++iEffHist) {

    //get histograms from file
    TH2F* numeratorHist = NULL;
    TH2F* denominatorHist = NULL;
    in.GetObject(iEffHist->first.first.c_str(), numeratorHist);
    in.GetObject(iEffHist->first.second.c_str(), denominatorHist);
    if ((numeratorHist == NULL) || (denominatorHist == NULL)) {
      cerr << errorCannotRetrieveObject(fnName, in.GetName(), iEffHist->first.first + "\" or \"" + 
					iEffHist->first.second);
      return CANNOT_RETRIEVE_OBJECT;
    }

    //set axis labels and titles
    setHistogramOptions(numeratorHist, kBlack, 0.7, 20, 1.35, 1.0, 
			iEffHist->second.first.c_str(), iEffHist->second.second.c_str());
    setHistogramOptions(denominatorHist, kBlack, 0.7, 20, 1.35, 1.0, 
			iEffHist->second.first.c_str(), iEffHist->second.second.c_str());
    numeratorHist->GetZaxis()->SetRangeUser(0.0, 1.1);
    denominatorHist->GetZaxis()->SetRangeUser(0.0, 1.1);

    //make efficiency histogram
    TH2F* effHist = (TH2F*)numeratorHist->Clone();
    effHist->Divide(denominatorHist);

    //print efficiencies
    for (Int_t iPTBin = 1; iPTBin <= effHist->GetNbinsX(); ++iPTBin) {
      for (Int_t iAbsEtaBin = 1; iAbsEtaBin <= effHist->GetNbinsY(); ++iAbsEtaBin) {
	cout << "    mistagEffVsPTAndEta_->SetBinContent(" << iPTBin << ", " << iAbsEtaBin << ", ";
	cout << effHist->GetBinContent(iPTBin, iAbsEtaBin) << ");\n";
      }
    }

    //draw efficiency histogram
    out.cd();
    string effCanvasName("eff_" + iEffHist->first.first + "_over_" + iEffHist->first.second);
    Int_t canvasWidth = defaultCanvasWidth;
    Int_t canvasHeight = defaultCanvasHeight;
    TCanvas effCanvas(effCanvasName.c_str(), "", canvasWidth, canvasHeight);
    setCanvasOptions(effCanvas, 1, 0, 0);
    setCanvasMargins(effCanvas, 0.125, 0.05, 0.125, 0.1);
    effHist->Draw("COLZ");
    effCanvas.Write();

    //save PDF of efficiency plot
    if (savePath != "noPDF") {
      formatSavePath(savePath);
      effCanvas.SaveAs((savePath + effCanvasName + ".pdf").c_str());
    }
  }

  //success
  return SUCCESS;
}

//make 1D histogram plots from input TH1s and save them to the output file
ErrorCode plot1DHistograms(TFile& in, TFile& out, const map<string, string>& hist1DMap, 
			   const map<string, vector<string> >& binLabelMap, string savePath)
{
  string fnName("ErrorCode plot1DHistograms(TFile& in, TFile& out, const map<string, string>& ");
  fnName+="hist1DMap, const map<string, vector<string> >& binLabelMap, string savePath)";

  //loop over 1D histogram map
  for (map<string, string>::const_iterator iHist1D = hist1DMap.begin(); 
       iHist1D != hist1DMap.end(); ++iHist1D) {

    //get histograms from file
    TH1F* hist = NULL;
    in.GetObject(iHist1D->first.c_str(), hist);
    if (hist == NULL) {
      cerr << errorCannotRetrieveObject(fnName, in.GetName(), iHist1D->first);
      return CANNOT_RETRIEVE_OBJECT;
    }

    //make 1D histogram
    Double_t xAxisLabelSize = defaultXAxisLabelSize;
    map<string, vector<string> >::const_iterator iHistLabels = 
      binLabelMap.find(iHist1D->first);
    if (iHistLabels != binLabelMap.end()) {
      setAxisLabels(hist->GetXaxis(), iHistLabels->second);
      xAxisLabelSize = 0.04;
    }
    setHistogramOptions(hist, kBlack, 0.7, 20, 1.0, iHist1D->second.c_str(), "", xAxisLabelSize);

    //draw 1D histogram
    out.cd();
    string canvas1DName(iHist1D->first);
    Int_t canvasWidth = defaultCanvasWidth;
    Int_t canvasHeight = defaultCanvasHeight;
    if (iHistLabels != binLabelMap.end()) {
      canvasWidth = 1200;
      canvasHeight = 600;
    }
    TCanvas canvas1D(canvas1DName.c_str(), "", canvasWidth, canvasHeight);
    setCanvasOptions(canvas1D, 1, 0, 0);
    if (iHistLabels != binLabelMap.end()) canvas1D.cd()->SetRightMargin(0.1);
    hist->Draw();
    canvas1D.Write();

    //save PDF of efficiency plot
    if (savePath != "noPDF") {
      formatSavePath(savePath);
      canvas1D.SaveAs((savePath + canvas1DName + ".pdf").c_str());
    }
  }

  //success
  return SUCCESS;
}

//make plots with efficiency and histogram overlaid
ErrorCode 
plot1DHistogramAndEfficiencyOverlaid(TFile& in, TFile& out, 
				     const map<string, pair<string, string> >& effHistMap, 
				     const map<string, vector<string> >& binLabelMap, 
				     string savePath, const float weight)
{
  string fnName("ErrorCode plot1DHistogramAndEfficiencyOverlaid(/*TFile& in, */TFile& out, ");
  fnName+="const map<string, pair<string, string> >& effHistMap, ";
  fnName+="const map<string, vector<string> >& binLabelMap, string savePath)";

  //loop over efficiency histogram map
  for (map<string, pair<string, string> >::const_iterator iEffHist = effHistMap.begin(); 
       iEffHist != effHistMap.end(); ++iEffHist) {

    //create canvas
    out.cd();
    string effAndHistCanvasName("eff_" + iEffHist->first + "_over_" + iEffHist->second.first + 
				"_and_hist_" + iEffHist->second.first);
    Int_t canvasWidth = defaultCanvasWidth;
    Int_t canvasHeight = defaultCanvasHeight;
    map<string, vector<string> >::const_iterator iHistLabels = 
      binLabelMap.find(iEffHist->second.first);
    if (iHistLabels != binLabelMap.end()) {
      canvasWidth = 1200;
      canvasHeight = 600;
    }
    TCanvas effAndHistCanvas(effAndHistCanvasName.c_str(), "", canvasWidth, canvasHeight);
    setCanvasOptions(effAndHistCanvas, 1, 0, 0);
    if (iHistLabels != binLabelMap.end()) effAndHistCanvas.cd()->SetRightMargin(0.1);

    //create pads, one transparent, for the two objects and axes
    string effCanvasName("eff_" + iEffHist->first + "_over_" + iEffHist->second.first);
    string effGraphName("divide_" + iEffHist->first + "_by_" + iEffHist->second.first);
    TPad graphPad(effCanvasName.c_str(), "", 0, 0, 1, 1);
    graphPad.SetRightMargin(0.2);
    graphPad.SetLeftMargin(0.2);
    graphPad.SetTopMargin(0.2);
    graphPad.SetBottomMargin(0.2);
    TPad histPad(effGraphName.c_str(), "", 0, 0, 1, 1);
    histPad.SetFillStyle(4000);
    histPad.SetRightMargin(0.2);
    histPad.SetLeftMargin(0.2);
    histPad.SetTopMargin(0.2);
    histPad.SetBottomMargin(0.2);
    graphPad.Draw();
    graphPad.cd();

    //get efficiency graph from file
    TCanvas* effCanvas = NULL;
    TGraphAsymmErrors* effGraph = NULL;
    out.GetObject(effCanvasName.c_str(), effCanvas);
    if (effCanvas != NULL) {
      effGraph = (TGraphAsymmErrors*)effCanvas->GetPrimitive(effGraphName.c_str());
      if (effGraph == NULL) {
	cerr << errorCannotRetrieveObject(fnName, out.GetName(), effGraphName);
	return CANNOT_RETRIEVE_OBJECT;
      }
    }
    else {
      cerr << errorCannotRetrieveObject(fnName, out.GetName(), effCanvasName);
      return CANNOT_RETRIEVE_OBJECT;
    }

    //get histogram
    TH1F* denominatorHist = NULL;
    in.cd();
    in.GetObject(iEffHist->second.first.c_str(), denominatorHist);

    //draw efficiency graph on canvas with left side axis
    out.cd();
    effAndHistCanvas.cd();
    graphPad.cd();
    graphPad.SetTicks(0, 0);
    effGraph->GetXaxis()->SetRangeUser(denominatorHist->GetXaxis()->GetXmin(), 
				       denominatorHist->GetXaxis()->GetXmax());
    effGraph->Draw("AP");
    graphPad.SetTicks(0, 0);
    graphPad.Update();
    effAndHistCanvas.cd();

    //draw histogram on canvas with right side axis
    out.cd();
    effAndHistCanvas.cd();
    histPad.Draw();
    histPad.cd();
    histPad.SetTicks(0, 0);
    denominatorHist->Scale(weight);
    denominatorHist->GetXaxis()->SetLabelColor(kWhite);
    denominatorHist->GetXaxis()->SetTitleColor(kWhite);
    denominatorHist->GetYaxis()->SetTitle("Entries / 5 GeV");
    denominatorHist->GetYaxis()->SetTitleOffset(1.6);
    denominatorHist->Draw("X+Y+");
    histPad.SetTicks(0, 0);
    effAndHistCanvas.Write();

    //save PDF of efficiency plot
    if (savePath != "noPDF") {
      formatSavePath(savePath);
      effAndHistCanvas.SaveAs((savePath + effAndHistCanvasName + ".pdf").c_str());
    }
  }

  //success
  return SUCCESS;
}

/*make canvases with multiple 1D histogram plots from different input files and save them to the 
  output file*/
ErrorCode plotMultiple1DHistogramsDifferentFiles(TFile& out, map<pair<string, string>, 
						 vector<pair<pair<TFile*, Option_t*>, 
						 pair<Color_t, Style_t> > > >& canvasMap, 
						 const map<string, vector<string> >& binLabelMap, 
						 string savePath)
{
  string fnName("ErrorCode plotMultiple1DHistogramsDifferentFiles(TFile& out, map<pair<string, ");
  fnName+="string>, vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > > >& ";
  fnName+="canvasMap, const map<string, vector<string> >& binLabelMap, string savePath)";

  //loop over 1D histogram map
  for (map<pair<string, string>, vector<pair<pair<TFile*, Option_t*>, 
	 pair<Color_t, Style_t> > > >::iterator iCanvas = canvasMap.begin(); 
       iCanvas != canvasMap.end(); ++iCanvas) {

    //unpack this map
    string canvas1DName(iCanvas->first.first);
    string unit(iCanvas->first.second);

    //get a pointer to the custom bin labels if applicable
    map<string, vector<string> >::const_iterator iHistLabels = 
      binLabelMap.find(canvas1DName);

    //make canvas
    out.cd();
    Int_t canvasWidth = defaultCanvasWidth;
    Int_t canvasHeight = defaultCanvasHeight;
    if (iHistLabels != binLabelMap.end()) {
      canvasWidth = 1200;
      canvasHeight = 600;
    }
    TCanvas canvas1D(canvas1DName.c_str(), "", canvasWidth, canvasHeight);
    setCanvasOptions(canvas1D, 1, 0, 0);
    if (iHistLabels != binLabelMap.end()) canvas1D.cd()->SetRightMargin(0.1);

    //loop over vector of histograms and their options to put on one canvas
    for (vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > >::iterator iHist = 
	   iCanvas->second.begin(); iHist != iCanvas->second.end(); ++iHist) {

      //unpack this vector
      TFile* in = iHist->first.first;
      Option_t* drawOption = iHist->first.second;
      Color_t markerColor = iHist->second.first;
      Style_t markerStyle = iHist->second.second;

      //get histogram from file
      TH1F* hist = NULL;
      in->GetObject(canvas1DName.c_str(), hist);
      if (hist == NULL) {
	cerr << errorCannotRetrieveObject(fnName, in->GetName(), canvas1DName);
	return CANNOT_RETRIEVE_OBJECT;
      }

      //make 1D histogram
      Double_t xAxisLabelSize = defaultXAxisLabelSize;
      if (iHistLabels != binLabelMap.end()) {
	setAxisLabels(hist->GetXaxis(), iHistLabels->second);
	xAxisLabelSize = 0.04;
      }
      setHistogramOptions(hist, markerColor, 0.7, markerStyle, 1.0/hist->GetEntries(), 
			  unit.c_str(), "", xAxisLabelSize);

      //draw 1D histogram
      out.cd();
      hist->Draw(drawOption);
    }

    //write canvas
    canvas1D.Write();

    //save PDF of efficiency plot
    if (savePath != "noPDF") {
      formatSavePath(savePath);
      canvas1D.SaveAs((savePath + canvas1DName + ".pdf").c_str());
    }
  }

  //success
  return SUCCESS;
}

//nicely format plots and save them to an output file, optionally printing them as PDFs
void plotNice(const string& inputFileName, const map<string, pair<string, string> >& effHistMap1D, 
	      const map<pair<string, string>, pair<string, string> >& effHistMap2D, 
	      const map<string, vector<string> >& binLabelMap, 
	      const map<string, string>& hist1DMap, const string& outputFileName, 
	      const string& savePath, const float weight)
{
  string fnName("const string& inputFileName, const map<string, pair<string, string> >& ");
  fnName+="effHistMap1D, const map<pair<string, string>, pair<string, string>& effHistMap2D, ";
  fnName+="const map<string, vector<string> >& binLabelMap, ";
  fnName+="const map<string, string>& hist1DMap, const string& outputFileName, ";
  fnName+="const string& savePath)";

  //open input file
  TFile in(inputFileName.c_str());
  if (!in.IsOpen()) {
    cerr << errorCannotOpenFile(fnName, inputFileName);
    cerr << inputFileName << ".\n";
    return;
  }

  //open output file
  TFile out(outputFileName.c_str(), "RECREATE");
  if (!out.IsOpen()) {
    cerr << errorCannotOpenFile(fnName, outputFileName);
    in.Close();
    return;
  }
  cout << "Processing file: " << inputFileName.c_str() << endl;

  //make efficiency plots
  plotEfficiency(in, out, effHistMap1D, binLabelMap, savePath);
  plot2DEfficiency(in, out, effHistMap2D, savePath);

  //make 1D histogram plots
  plot1DHistograms(in, out, hist1DMap, binLabelMap, savePath);

  //make plots with efficiency and histogram overlaid
  plot1DHistogramAndEfficiencyOverlaid(in, out, effHistMap1D, binLabelMap, savePath, weight);

  //write output file
  out.cd();
  out.Write();

  //close files
  out.Close();
  in.Close();
}

/*nicely format plots from different input files and save them to an output file, optionally 
  printing them as PDFs*/
void plotNiceDifferentFiles(map<pair<string, string>, vector<pair<pair<TFile*, Option_t*>, 
			    pair<Color_t, Style_t> > > >& canvasMap, 
			    const map<string, vector<string> >& binLabelMap, 
			    const string& outputFileName, const string& savePath)
{
  string fnName("void plotNiceDifferentFiles(const map<pair<string, string>, ");
  fnName+="vector<pair<pair<string, Option_t*>, pair<Color_t, Style_t> > > >& canvasMap, ";
  fnName+="const map<string, vector<string> >& binLabelMap, const string& outputFileName, ";
  fnName+="const string& savePath)";

  //open output file
  TFile out(outputFileName.c_str(), "RECREATE");
  if (!out.IsOpen()) {
    cerr << errorCannotOpenFile(fnName, outputFileName);
    return;
  }

  //make 1D histogram plots
  plotMultiple1DHistogramsDifferentFiles(out, canvasMap, binLabelMap, savePath);

  //write output file
  out.cd();
  out.Write();

  //close files
  out.Close();
  for (map<pair<string, string>, vector<pair<pair<TFile*, Option_t*>, 
	 pair<Color_t, Style_t> > > >::iterator iCanvas = canvasMap.begin(); 
       iCanvas != canvasMap.end(); ++iCanvas) {
    for (vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > >::iterator iHist = 
	   iCanvas->second.begin(); iHist != iCanvas->second.end(); ++iHist) {
      iHist->first.first->Close();
      delete iHist->first.first;
    }
  }
}

//return the histogram with the largest bin content
TH1F* histWithLargestBinContent(const vector<TH1F*>& hists)
{
  TH1F* pHistWithMaxMaxBin = NULL;
  Double_t maxBinContent = 0.0;
  for (vector<TH1F*>::const_iterator iHist = hists.begin(); iHist != hists.end(); ++iHist) {
    Int_t histMaxBin = (*iHist)->GetMaximumBin();
    Double_t histMaxBinContent = (*iHist)->GetBinContent(histMaxBin);
    if (histMaxBinContent >= maxBinContent) {
      maxBinContent = histMaxBinContent;
      pHistWithMaxMaxBin = *iHist;
    }
  }
  return pHistWithMaxMaxBin;
}

template<typename T>
void setup(const vector<string>& canvasNames, vector<TCanvas*>& outputCanvases, 
	   const bool setLogY, vector<TLegend*>& legends, vector<THStack*>& stacks, 
	   const vector<string>& legendHeaders, vector<vector<T*> >& hists, 
	   const unsigned int size, const bool dataMC, const bool twoDim)
{
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    if (dataMC) {
      outputCanvases.push_back(new TCanvas(iCanvasName->c_str(), "", 600, 900));
      outputCanvases[outputCanvases.size() - 1]->Divide(1, 2);
      outputCanvases[outputCanvases.size() - 1]->cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
      setCanvasOptions(*outputCanvases[outputCanvases.size() - 1]->cd(1), 1, setLogY, 0);
      if (twoDim) {
	setCanvasMargins(*outputCanvases[outputCanvases.size() - 1]->cd(1), 0.2, 0.2, 0.2, 0.2);
      }
      outputCanvases[outputCanvases.size() - 1]->cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
      setCanvasOptions(*outputCanvases[outputCanvases.size() - 1]->cd(2), 1, 0, 0);
      if (twoDim) {
	setCanvasMargins(*outputCanvases[outputCanvases.size() - 1]->cd(2), 0.2, 0.2, 0.2, 0.2);
      }
    }
    else {
      outputCanvases.push_back(new TCanvas(iCanvasName->c_str(), "", 600, 600));
      setCanvasOptions(*outputCanvases[outputCanvases.size() - 1], 1, setLogY, 0);
      if (twoDim) setCanvasMargins(*outputCanvases[outputCanvases.size() - 1], 0.2, 0.2, 0.2, 0.2);
    }
    legends.push_back(new TLegend(0.4, 0.5, 0.8, 0.9));
    string stackName(*iCanvasName);
    if (stackName.find("Canvas") != string::npos) {
      stackName.replace(stackName.find("Canvas"), 6, "Stack");
    }
    else stackName+="Stack";
    stacks.push_back(new THStack(stackName.c_str(), ""));
    if (legendHeaders.size() > canvasIndex) setLegendOptions(*legends[legends.size() - 1], 
							     legendHeaders[canvasIndex].c_str());
    hists.push_back(vector<T*>(size, NULL));
  }
}

template<typename T>
void setup(const vector<string>& canvasNames, vector<TCanvas*>& outputCanvases, 
	   vector<T*>& hists, const bool twoDim = false)
{
  vector<TLegend*> legends;
  vector<THStack*> stacks;
  vector<string> legendHeaders;
  vector<vector<TH1F*> > dummyHists;
  setup(canvasNames, outputCanvases, false, legends, stacks, legendHeaders, dummyHists, 1, false, 
	twoDim);
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) { hists.push_back(NULL); }
}

template<typename T>
void scaleAndAdd(const vector<string>& canvasNames, TFile* in, const vector<string>& graphNames, 
		 const float weight, vector<T*>& hists, const unsigned int fileIndex, 
		 const vector<Int_t>& blindLow, const vector<Int_t>& blindHigh)
{
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    TCanvas* pCanvas;
    in->GetObject(iCanvasName->c_str(), pCanvas);
    T* pHist = (T*)pCanvas->GetPrimitive(graphNames[canvasIndex].c_str());

    //include overflows in last bin
    Int_t lastBin = pHist->GetNbinsX();
    Double_t lastBinErr = pHist->GetBinError(lastBin);
    Double_t overflowBinErr = pHist->GetBinError(lastBin + 1);
    pHist->
      SetBinContent(lastBin, pHist->GetBinContent(lastBin) + pHist->GetBinContent(lastBin + 1));
    pHist->SetBinError(lastBin, sqrt(lastBinErr*lastBinErr + overflowBinErr*overflowBinErr));
    pHist->SetBinContent(lastBin + 1, 0.0);
    pHist->SetBinError(lastBin + 1, 0.0);
    pHist->Scale(weight);
    if (fileIndex == 0) hists[canvasIndex] = pHist;
    else hists[canvasIndex]->Add(pHist);
    Int_t iBin = blindLow[canvasIndex];
    Int_t endBin = blindHigh[canvasIndex];
    if ((blindLow[canvasIndex] < 0) || (blindHigh[canvasIndex] < -2)) {
      cerr << "Error: blindLow should be in [0, inf) and blindHigh should be in [-2, inf).\n";
      cerr << "Blinding entire histogram.\n";
      cerr << blindLow[canvasIndex] << " " << blindHigh[canvasIndex] << endl;
      iBin = 0;
      endBin = hists[canvasIndex]->GetNbinsX() + 1;
    }
    else {
      if (blindHigh[canvasIndex] == -1) endBin = hists[canvasIndex]->GetNbinsX() + 1;
      if (blindHigh[canvasIndex] == -2) endBin = blindLow[canvasIndex] - 1;
    }
    while (iBin <= endBin) {
      hists[canvasIndex]->SetBinContent(iBin, 0.0);
      hists[canvasIndex]->SetBinError(iBin, 0.0);
      ++iBin;
    }
  }
}

void draw(const vector<string>& canvasNames, TFile& outStream, vector<TCanvas*>& outputCanvases, 
	  vector<TH1F*>& hists)
{
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    outStream.cd();
    outputCanvases[canvasIndex]->cd();
    if (hists[canvasIndex] != NULL) hists[canvasIndex]->Draw();
  }
}

void draw(const vector<string>& canvasNames, TFile& outStream, vector<TCanvas*>& outputCanvases, 
	  vector<TH2F*>& hists)
{
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    outStream.cd();
    outputCanvases[canvasIndex]->cd();
    outputCanvases[canvasIndex]->SetLogz();
    if (hists[canvasIndex] != NULL) {
      setHistogramOptions(hists[canvasIndex], kBlack, 0.7/*4.2*/, 20, 1.6, 1.0);
//       hists[canvasIndex]->GetZaxis()->SetRangeUser(0.00001, 1.0);
      hists[canvasIndex]->GetZaxis()->SetRangeUser(0.1, 10000.0);
      hists[canvasIndex]->Draw("COLZ");
//       TProfile* profileX = 
// 	hists[canvasIndex]->ProfileX((string(hists[canvasIndex]->GetName()) + "_pfx").c_str());
//       profileX->SetLineWidth(3);
//       profileX->Draw("HISTSAME");
//       profileX->Draw("ESAME");
    }
  }
}

void write(vector<TCanvas*>& outputCanvases)
{
  for (vector<TCanvas*>::iterator iOutputCanvas = outputCanvases.begin(); 
       iOutputCanvas != outputCanvases.end(); ++iOutputCanvas) {
    (*iOutputCanvas)->Write();
//     stringstream fileName;
//     fileName << "/afs/cern.ch/user/y/yohay/AN-13-254/notes/AN-13-254/trunk/regB-data-MC-v149-";
//     fileName << (*iOutputCanvas)->GetName() << ".pdf";
//     (*iOutputCanvas)->SaveAs(fileName.str().c_str());
  }
}

template<typename T>
void deleteObjects(vector<T*>& objs)
{
  for (typename vector<T*>::iterator iObj = objs.begin(); iObj != objs.end(); 
       ++iObj) { delete *iObj; }
}

void deleteStreams(vector<TFile*>& streams)
{
  for (vector<TFile*>::const_iterator iStream = streams.begin(); 
       iStream != streams.end(); ++iStream) {
    (*iStream)->Close();
    delete *iStream;
  }
}

/*create a histogram of (data - bkg.)/sqrt(jetFakeStatErr^2 + resBkgStatErr^2)
  for all cases, use the nominal error as the denominator*/
TH1F* makeDataBkgAgreementHist(const TH1F* data, const TH1F* bkg)
{
  if ((data->GetNbinsX() + 2) != (bkg->GetNbinsX() + 2)) {
    cerr << "Error: size mismatch.\n";
    return NULL;
  }
  TH1F* ratioHist = (TH1F*)data->Clone();
  ratioHist->Add(bkg, -1.0);
  for (Int_t iBin = 0; iBin <= (bkg->GetNbinsX() + 1); ++iBin) {
    Double_t dataBkgStatErr = ratioHist->GetBinError(iBin);
    Double_t bkgStatErr = bkg->GetBinError(iBin);
    Double_t binContent = bkgStatErr == 0.0 ? 0.0 : (ratioHist->GetBinContent(iBin)/bkgStatErr);
    ratioHist->SetBinContent(iBin, binContent);
    ratioHist->SetBinError(iBin, dataBkgStatErr/bkgStatErr);
  }
  ratioHist->GetYaxis()->SetTitle("#frac{Data - Bkg.}{Bkg. stat. error}");
  ratioHist->GetYaxis()->SetRangeUser(-3.0, 3.0);
  return ratioHist;
}

/*create a histogram of (data - bkg.)/
  sqrt(jetFakeStatErr^2 + resBkgStatErr^2 + jetFakeSystErr^2)*/
TH1F* makeDataBkgAgreementHist(const TH1F* data, const TH1F* nomBkg, const TH1F* allQCDBkg, 
			       const TH1F* allEWBkg)
{
  if (((data->GetNbinsX()) != (nomBkg->GetNbinsX())) || 
      ((nomBkg->GetNbinsX()) != (allQCDBkg->GetNbinsX())) || 
      ((allQCDBkg->GetNbinsX()) != (allEWBkg->GetNbinsX()))) {
    cerr << "Error: size mismatch.\n";
    return NULL;
  }
  TH1F* ratioHist = (TH1F*)data->Clone();
  ratioHist->Add(nomBkg, -1.0);
  for (Int_t iBin = 0; iBin <= (nomBkg->GetNbinsX() + 1); ++iBin) {
    Double_t dataNomBkgStatErr = ratioHist->GetBinError(iBin);
    Double_t nomBkgStatErr = nomBkg->GetBinError(iBin);
    Double_t nomBkgVal = nomBkg->GetBinContent(iBin);
    Double_t allQCDBkgSystErr2 = allQCDBkg->GetBinContent(iBin) - nomBkgVal;
    allQCDBkgSystErr2*=allQCDBkgSystErr2;
    Double_t allEWBkgSystErr2 = allEWBkg->GetBinContent(iBin) - nomBkgVal;
    allEWBkgSystErr2*=allEWBkgSystErr2;
    Double_t totBkgErr = 0.0;
    if (allQCDBkgSystErr2 > allEWBkgSystErr2) {
      totBkgErr = sqrt(nomBkgStatErr*nomBkgStatErr + allQCDBkgSystErr2);
    }
    else totBkgErr = sqrt(nomBkgStatErr*nomBkgStatErr + allEWBkgSystErr2);
    Double_t binContent = totBkgErr == 0.0 ? 0.0 : (ratioHist->GetBinContent(iBin)/totBkgErr);
    ratioHist->SetBinContent(iBin, binContent);
    ratioHist->SetBinError(iBin, dataNomBkgStatErr/nomBkgStatErr); //meaningless
  }
  ratioHist->GetYaxis()->SetTitle("#frac{Data - Bkg.}{Bkg. error (excl. resonances)}");
  ratioHist->GetYaxis()->SetRangeUser(-3.0, 3.0);
  return ratioHist;
}

//merge efficiency graphs from different files into one canvas
void drawMultipleEfficiencyGraphsOn1Canvas(const string& outputFileName, 
					   const vector<string>& inputFiles, 
					   const vector<string>& canvasNames, 
					   const vector<string>& graphNames, 
					   const vector<string>& legendHeaders, 
					   const vector<Color_t>& colors, 
					   const vector<Style_t>& styles, 
					   const vector<string>& legendEntries, 
					   const vector<float>& weights, const bool setLogY, 
					   const bool drawStack, const bool dataMC)
{
  if ((inputFiles.size() > colors.size()) || (inputFiles.size() > styles.size()) || 
      (inputFiles.size() > legendEntries.size()) || (inputFiles.size() > weights.size()) || 
      (canvasNames.size() > graphNames.size()) || (canvasNames.size() > legendHeaders.size())) {
    cerr << "Error: vector size mismatch.\n";
    cerr << inputFiles.size() << endl;
    cerr << colors.size() << endl;
    cerr << styles.size() << endl;
    cerr << legendEntries.size() << endl;
    for (vector<string>::const_iterator i = legendEntries.begin(); i != legendEntries.end(); ++i) {
      cerr << *i << endl;
    }
    cerr << weights.size() << endl;
    cerr << legendHeaders.size() << endl;
    cerr << canvasNames.size() << endl;
    cerr << graphNames.size() << endl;
    return;
  }
  TFile outStream(outputFileName.c_str(), "RECREATE");
  vector<TFile*> inputStreams;
  vector<TCanvas*> outputCanvases;
  vector<TLegend*> legends;
  vector<THStack*> stacks;
  vector<vector<TH1F*> > hists;
  setup(canvasNames, outputCanvases, setLogY, legends, stacks, legendHeaders, hists, 
	inputFiles.size(), dataMC, false);
  bool data = true;
  bool gg = false;
  bool ZH = false;
  bool VBF = false;
  for (vector<string>::const_iterator iInputFile = inputFiles.begin(); 
       iInputFile != inputFiles.end(); ++iInputFile) {
    const unsigned int fileIndex = iInputFile - inputFiles.begin();
    if ((fileIndex == 0) && (iInputFile->find("Wh1") != string::npos)) data = false;
    if ((fileIndex == 1) && (iInputFile->find("gg") != string::npos)) gg = true;
    if ((fileIndex == 2) && (iInputFile->find("ZH") != string::npos)) ZH = true;
    if ((fileIndex == 3) && (iInputFile->find("VBF") != string::npos)) VBF = true;
    inputStreams.push_back(new TFile(iInputFile->c_str()));
    for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
	 iCanvasName != canvasNames.end(); ++iCanvasName) {
      const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
      TCanvas* pCanvas;
      inputStreams[inputStreams.size() - 1]->GetObject(iCanvasName->c_str(), pCanvas);
      TGraphAsymmErrors* pGraph = NULL;
      TH1F* pHist = NULL;
      if (graphNames[canvasIndex].find("divide") != string::npos) {
	pGraph = (TGraphAsymmErrors*)pCanvas->GetPrimitive(graphNames[canvasIndex].c_str());
	setGraphOptions(*pGraph, colors[fileIndex], 0.7, styles[fileIndex], 
			pGraph->GetXaxis()->GetTitle(), pGraph->GetYaxis()->GetTitle());
	legends[canvasIndex]->AddEntry(pGraph, legendEntries[fileIndex].c_str(), "l");
      }
      else {
	pHist = (TH1F*)pCanvas->GetPrimitive(graphNames[canvasIndex].c_str());
	float weight = weights[fileIndex];
	if (weights[fileIndex] == 0.0) weight = 1.0/pHist->Integral(0, -1);
	if (weights[fileIndex] == -1.0) {
	  if (fileIndex == 0) weight = 1.0;
	  else weight = hists[canvasIndex][0]->Integral(0, -1)/pHist->Integral(0, -1);
	}
	setHistogramOptions(pHist, colors[fileIndex], 0.7, styles[fileIndex], 
			    weight, 
			    string(pHist->GetName()) == "muHadPTOverMuHadMass" ? 
			    "p_{T}/m" : pHist->GetXaxis()->GetTitle(), 
			    pHist->GetYaxis()->GetTitle());
	if (string(pHist->GetName()) == "muHadMass"/*"tauHadIso"*/) {
	  cout << "Processing file " << *iInputFile << endl;
	  Double_t mGeq4GeVErr = 0.0;
	  Double_t mGeq4GeVInt = pHist->IntegralAndError(5/*17*/, -1, mGeq4GeVErr);
	  Double_t mGeq4GeVErrTerm = 
	    mGeq4GeVInt == 0.0 ? 0.0 : (mGeq4GeVErr*mGeq4GeVErr)/(mGeq4GeVInt*mGeq4GeVInt);
	  cout << "m >= 4 GeV: " << mGeq4GeVInt << " +/- " << mGeq4GeVErr << endl;
	  Double_t mLe2GeVErr = 0.0;
	  Double_t mLe2GeVInt = pHist->IntegralAndError(1, 2/*8*/, mLe2GeVErr);
	  Double_t mLe2GeVErrTerm = 
	    mLe2GeVInt == 0.0 ? 0.0 : (mLe2GeVErr*mLe2GeVErr)/(mLe2GeVInt*mLe2GeVInt);
	  cout << "m < 2 GeV: " << mLe2GeVInt << " +/- " << mLe2GeVErr << endl;
	  Double_t ratio = mGeq4GeVInt == 0.0 ? 0.0 : (mLe2GeVInt/mGeq4GeVInt);
	  double errorRatio = ratio*sqrt(mGeq4GeVErrTerm + mLe2GeVErrTerm);
	  cout << "(m < 2 GeV)/(m >= 4 GeV): " << setprecision(3) << ratio << " +/- ";
	  cout << errorRatio << endl;
	  if (!data)
	    cout << "Total m integral for dz: " << pHist->Integral() << endl;
	}
	string histName(pHist->GetName());
	if (histName == "jet_pt_etacut") pHist->GetXaxis()->SetTitle("p_{T} (GeV)");
	if (histName == "jet_eta") pHist->GetXaxis()->SetTitle("#eta");
	if (histName == "jet_phi") pHist->GetXaxis()->SetTitle("#phi");
	if (histName == "jet_mass_etacut") pHist->GetXaxis()->SetTitle("m (GeV)");
	if (histName == "jet_ptmj_etacut") {
	  pHist->GetXaxis()->SetTitle("#frac{p_{T}}{m}");
	}
// 	if ((histName == "muHadMass") && 
// 	    (outputFileName.find("dataVsMC_muHadNonIsoAnalysis") != string::npos))
// 	  {
// 	    cout << "For " << iInputFile->c_str() << ", " << histName << " total integral = ";
// 	    cout << pHist->Integral() << endl;
// 	  }
	if (setLogY) pHist->GetYaxis()->SetRangeUser(0.1, 10000.0);
	string legendStyle("l");
	if (drawStack) {
// 	  pHist->SetFillStyle(0);
// 	  pHist->SetFillColor(0);
	  pHist->SetFillStyle(1001);
	  if (((fileIndex == 0) || (gg && (fileIndex == 1)) || (ZH && (fileIndex == 2)) || 
	       (VBF && (fileIndex == 3))) && !data) pHist->SetFillColor(0);
	  else pHist->SetFillColor(colors[fileIndex]);
	  if (fileIndex == 0) {
	    if(!data)
	      {
		pHist->SetLineWidth(4);
		pHist->SetLineStyle(2);
	      }
	    if (data) legendStyle = "lp";
	    else legendStyle = "l";
	  }
	  else if (((gg && (fileIndex == 1)) || (ZH && (fileIndex == 2)) || 
		    (VBF && (fileIndex == 3))) && !data)
	    {
	      legendStyle = "l";
	      pHist->SetLineWidth(4);
	      pHist->SetLineStyle(2);
	    }
	  else legendStyle = "f";
	}
	hists[canvasIndex][fileIndex] = pHist;
	legends[canvasIndex]->
	  AddEntry(pHist, legendEntries[fileIndex].c_str(), legendStyle.c_str());
	if ((data && (fileIndex != 0)) || 
	    (!data && (fileIndex != 0) && ((gg && (fileIndex != 1)) || !gg) && 
	     ((ZH && (fileIndex != 2)) || !ZH) && ((VBF && (fileIndex != 3)) || !VBF))) {
	  stacks[canvasIndex]->Add(pHist, "HIST");
	}
      }
      outStream.cd();
      outputCanvases[canvasIndex]->cd(dataMC ? 1 : 0);
      if (fileIndex == 0) { if (pGraph != NULL) pGraph->Draw("AP"); }
      else if (pGraph != NULL) pGraph->Draw("PSAME");
    }
  }
  if (!drawStack) {
    for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
	 iCanvasName != canvasNames.end(); ++iCanvasName) {
      const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
      outStream.cd();
      outputCanvases[canvasIndex]->cd(dataMC ? 1 : 0);
      TH1F* pHistWithMaxMaxBin = histWithLargestBinContent(hists[canvasIndex]);
      if (pHistWithMaxMaxBin != NULL) pHistWithMaxMaxBin->Draw("HIST");
    }
  }
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    outStream.cd();
    outputCanvases[canvasIndex]->cd(dataMC ? 1 : 0);
    if (!drawStack) {
      TH1F* ratioHist = NULL;
      for (vector<string>::const_iterator iInputFile = inputFiles.begin(); 
	   iInputFile != inputFiles.end(); ++iInputFile) {
	const unsigned int fileIndex = iInputFile - inputFiles.begin();
	if (hists[canvasIndex][fileIndex] != NULL) {
	  hists[canvasIndex][fileIndex]->Draw("HISTESAME");
	  if (dataMC) {
	    if (fileIndex == 0) ratioHist = (TH1F*)hists[canvasIndex][fileIndex]->Clone();
	    else ratioHist->Divide(hists[canvasIndex][fileIndex]);
	  }
	}
      }
      if (ratioHist != NULL) {
	ratioHist->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratioHist->GetYaxis()->SetTitle("#frac{Non-isolated}{Isolated}");
	outputCanvases[canvasIndex]->cd(dataMC ? 2 : 0);
	ratioHist->Draw();
      }
    }
    else {
      if (setLogY) {
	stacks[canvasIndex]->SetMinimum(1.0);
	stacks[canvasIndex]->SetMaximum(10000000.0);
      }
      //stacks[canvasIndex]->Draw("e");
      stacks[canvasIndex]->Draw();
      TList* stackedHists = stacks[canvasIndex]->GetHists();
      TH1F* hist = (TH1F*)stackedHists->First();
      stacks[canvasIndex]->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
      stacks[canvasIndex]->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
      string drawOpt("SAME");
      if (!data) drawOpt = "HISTSAME";
      hists[canvasIndex][0]->Draw(drawOpt.c_str());
      if (!data && gg) hists[canvasIndex][1]->Draw(drawOpt.c_str());
      if (!data && ZH) hists[canvasIndex][2]->Draw(drawOpt.c_str());
      if (!data && VBF) hists[canvasIndex][3]->Draw(drawOpt.c_str());
      TH1F* stackSumHist = NULL;
      for (Int_t i = 0; i < stackedHists->GetEntries(); ++i) {
	TH1F* stackHist = (TH1F*)stackedHists->At(i)->Clone();
	stackHist->SetTitle("");
	if (i == 0) stackSumHist = stackHist;
	else stackSumHist->Add(stackHist);
      }
      if (*iCanvasName == "muHadMassCanvas") {
	Double_t MCQCDErrSigReg, MCQCDErrNormReg;
	Double_t MCQCDSigReg = stackSumHist->IntegralAndError(5/*17*/, -1, MCQCDErrSigReg);
	Double_t MCQCDNormReg = stackSumHist->IntegralAndError(1, 2/*8*/, MCQCDErrNormReg);
	Double_t MCQCDNormRegFracErr = 
	  ((MCQCDNormReg > 0.0) ? (MCQCDErrNormReg/MCQCDNormReg) : 0.0);
	Double_t MCQCDSigRegFracErr = ((MCQCDSigReg > 0.0) ? (MCQCDErrSigReg/MCQCDSigReg) : 0.0);
	Double_t MCQCDNormToSigRatio = ((MCQCDSigReg > 0.0) ? (MCQCDNormReg/MCQCDSigReg) : 0.0);
	Double_t MCQCDNormToSigRatioErr2 = MCQCDNormToSigRatio*MCQCDNormToSigRatio*
	  (MCQCDNormRegFracErr*MCQCDNormRegFracErr + MCQCDSigRegFracErr*MCQCDSigRegFracErr);
	string reg = (dataMC ? "B" : "A");
	cout << endl << "MC (optionally + QCD) prediction in region ";
	cout << reg << ", m >= 4 GeV: ";
	cout << MCQCDSigReg << " +/- " << MCQCDErrSigReg << endl;
	cout << "MC (optionally + QCD) prediction in region ";
	cout << reg << ", m < 2 GeV: ";
	cout << setprecision(3) << MCQCDNormReg << " +/- " << MCQCDErrNormReg << endl;
	cout << "(m < 2 GeV)/(m >= 4 GeV): " << setprecision(3) << MCQCDNormToSigRatio << " +/- ";
	cout << sqrt(MCQCDNormToSigRatioErr2) << endl;
      }
      TH1F* dataHist = (TH1F*)hists[canvasIndex][0]->Clone();

      // New: calculate pull //
      TH1F* pullHist = makeDataBkgAgreementHist(dataHist, stackSumHist);

      // Old data-MC agreement histogram //
     //  TH1F* denomHist = denomErrorScale((TH1F*)stackSumHist->Clone()); 
     //  stackSumHist->Add(dataHist, -1.0);
     // //stackSumHist->Divide(dataHist);
     //  stackSumHist->Divide(denomHist);
     //  setHistogramOptions(stackSumHist, kBlack, 0.7, 20, 1.0, 
     // 			  stackSumHist->GetXaxis()->GetTitle(), 
     // 			  //"#frac{N_{MC} - N_{data}}{N_{data}}");
     // 			  "#frac{N_{MC} - N_{data}}{sqrt(N_{MC}^2 + #sigma_{MC}^2}");
     //  stackSumHist->GetYaxis()->SetRangeUser(-1.0, 1.0);

      if (dataMC) {
	outputCanvases[canvasIndex]->cd(dataMC ? 2 : 0);

	// New: draw pull //
	pullHist->GetYaxis()->SetTitle("#frac{Data - MC}{MC stat. error}");
	pullHist->SetLineColor(kBlack);
	pullHist->SetFillColor(kAzure);
	pullHist->SetFillStyle(1001);
	pullHist->Draw("HIST");

	// Old data-MC agreement histogram //
	// stackSumHist->Draw("e");
      }
    }
    outputCanvases[canvasIndex]->cd(dataMC ? 1 : 0);
    legends[canvasIndex]->Draw();
  }
  outStream.cd();
  write(outputCanvases);
  outStream.Write();
  outStream.Close();
  deleteObjects(legends);
  deleteObjects(stacks);
  deleteObjects(outputCanvases);
  deleteStreams(inputStreams);
}

//predicates needed for find_if
bool lessThan0(int i) { return (i < 0); }
bool lessThanNeg2(int i) { return (i < -2); }

//hadd histograms drawn on canvases
//blindLow = -1 ==> overflow bin
//blindHigh = -2 ==> don't blind anything
void haddCanvases(const string& outputFileName, const vector<string>& inputFiles, 
		  const vector<float>& weights, const vector<string>& canvasNames1D, 
		  const vector<string>& graphNames1D, const vector<string>& canvasNames2D, 
		  const vector<string>& graphNames2D, const vector<Int_t>& blindLow, 
		  const vector<Int_t>& blindHigh)
{
  if ((inputFiles.size() > weights.size()) || (canvasNames1D.size() != graphNames1D.size()) || 
      (graphNames1D.size() != blindLow.size()) || (blindLow.size() != blindHigh.size()) || 
      (canvasNames2D.size() != graphNames2D.size())) {
    cerr << "Error: vector size mismatch.\n";
    return;
  }
  if ((find_if(blindLow.begin(), blindLow.end(), lessThan0) != blindLow.end()) || 
      (find_if(blindHigh.begin(), blindHigh.end(), lessThanNeg2) != blindHigh.end())) {
    cerr << "Error: blindLow should be in [0, inf) and blindHigh should be in [-2, inf).\n";
    return;
  }
  TFile outStream(outputFileName.c_str(), "RECREATE");
  vector<TFile*> inputStreams;
  vector<TCanvas*> outputCanvases1D;
  vector<TCanvas*> outputCanvases2D;
  vector<TH1F*> hists1D;
  vector<TH2F*> hists2D;
  setup(canvasNames1D, outputCanvases1D, hists1D);
  setup(canvasNames2D, outputCanvases2D, hists2D, true);
  vector<Int_t> nullBlindLow(canvasNames2D.size(), 0);
  vector<Int_t> nullBlindHigh(canvasNames2D.size(), -2);
  for (vector<string>::const_iterator iInputFile = inputFiles.begin(); 
       iInputFile != inputFiles.end(); ++iInputFile) {
    const unsigned int fileIndex = iInputFile - inputFiles.begin();
    inputStreams.push_back(new TFile(iInputFile->c_str()));
    scaleAndAdd(canvasNames1D, inputStreams[inputStreams.size() - 1], graphNames1D, 
		weights[fileIndex], hists1D, fileIndex, blindLow, blindHigh);
    scaleAndAdd(canvasNames2D, inputStreams[inputStreams.size() - 1], graphNames2D, 
		weights[fileIndex], hists2D, fileIndex, nullBlindLow, nullBlindHigh);
  }
  draw(canvasNames1D, outStream, outputCanvases1D, hists1D);
  draw(canvasNames2D, outStream, outputCanvases2D, hists2D);
  outStream.cd();
  write(outputCanvases1D);
  write(outputCanvases2D);
  outStream.Write();
  outStream.Close();
  deleteObjects(outputCanvases1D);
  deleteObjects(outputCanvases2D);
  deleteStreams(inputStreams);
}

//plot histograms with different names in the same file on the same canvas
void mergePlotsIn1File(TFile& in, TFile& out, const vector<string>& canvasNames, 
		       const vector<string>& histNames, const float* weights, 
		       const Color_t* colors, const Style_t* styles, const char* canvasName, 
		       const char* legendHeader, const char** legendEntries, const bool setLogY)
{
  TCanvas canvas(canvasName, "", 600, 600);
  TLegend legend(0.5, 0.6, 0.9, 0.8);
  vector<TH1F*> hists(canvasNames.size(), NULL);
  setCanvasOptions(canvas, 1, setLogY, 0);
  setLegendOptions(legend, legendHeader);
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    TCanvas* pCanvas = NULL;
    in.GetObject(iCanvasName->c_str(), pCanvas);
    TH1F* pHist = NULL;
    if (pCanvas != NULL) {
      pHist = (TH1F*)pCanvas->GetPrimitive(histNames[canvasIndex].c_str());
      float weight = 
	weights[canvasIndex] == 0.0 ? 1.0/pHist->Integral(0, -1) : weights[canvasIndex];
      setHistogramOptions(pHist, colors[canvasIndex], 0.7, styles[canvasIndex], 
			  weight, pHist->GetXaxis()->GetTitle(), 
			  pHist->GetYaxis()->GetTitle());
      if (setLogY) pHist->GetYaxis()->SetRangeUser(0.1, 50000.0);
      hists[canvasIndex] = pHist;
      legend.AddEntry(pHist, legendEntries[canvasIndex], "l");
    }
  }
  out.cd();
  canvas.cd();
  TH1F* pHistWithMaxMaxBin = histWithLargestBinContent(hists);
  if (pHistWithMaxMaxBin != NULL) pHistWithMaxMaxBin->Draw();
  for (vector<TH1F*>::const_iterator iHist = hists.begin(); 
       iHist != hists.end(); ++iHist) {
    if (*iHist != NULL) (*iHist)->Draw("SAME");
  }
  legend.Draw();
  canvas.Write();
}

//plot histograms with different names in the same file on the same canvas
void mergePlotsIn1File(const char* inputFileName, const char* outputFileName)
{
  TFile in(inputFileName);
  TFile out(outputFileName, "RECREATE");
  vector<string> canvasNames;
  canvasNames.push_back("dRWMuSoftMuCanvas");
  canvasNames.push_back("dRWMuSoftMuMuHadMassGe2Canvas");
  vector<string> histNames;
  histNames.push_back("dRWMuSoftMu");
  histNames.push_back("dRWMuSoftMuMuHadMassGe2");
  const float weights[2] = {0.0, 0.0};
  const Color_t colors[2] = {kBlack, kRed};
  const Style_t styles[2] = {20, 21};
  const char* legendEntries[2] = {"m_{#mu+had} > 0 GeV", "m_{#mu+had} > 2 GeV"};
  mergePlotsIn1File(in, out, canvasNames, histNames, weights, colors, styles, "dRWMuSoftMu", 
		    "Normalized to 1", legendEntries, false);
  in.Close();
  out.Close();
}

//draw difference between data and MC
void drawDifferenceGraphsOn1Canvas(const string& outputFileName, 
					   const vector<string>& inputFiles, 
					   const vector<string>& canvasNames, 
					   const vector<string>& graphNames, 
					   const vector<string>& legendHeaders, 
					   const vector<Color_t>& colors, 
					   const vector<Style_t>& styles, 
					   const vector<string>& legendEntries, 
					   const vector<float>& weights, const bool setLogY, 
					   const bool dataMC)
{
  if ((inputFiles.size() > colors.size()) || (inputFiles.size() > styles.size()) || 
      (inputFiles.size() > legendEntries.size()) || (inputFiles.size() > weights.size()) || 
      (canvasNames.size() > graphNames.size()) || (canvasNames.size() > legendHeaders.size())) {
    cerr << "Error: vector size mismatch.\n";
    cerr << inputFiles.size() << endl;
    cerr << colors.size() << endl;
    cerr << styles.size() << endl;
    cerr << legendEntries.size() << endl;
    for (vector<string>::const_iterator i = legendEntries.begin(); i != legendEntries.end(); ++i) {
      cerr << *i << endl;
    }
    cerr << weights.size() << endl;
    return;
  }
  TFile outStream(outputFileName.c_str(), "RECREATE");
  vector<TFile*> inputStreams;
  vector<TCanvas*> outputCanvases;
  vector<TLegend*> legends;
  vector<THStack*> stacks;
  vector<vector<TH1F*> > hists;
  vector<TH1F*> histDiff(canvasNames.size());
  setup(canvasNames, outputCanvases, setLogY, legends, stacks, legendHeaders, hists, 
	inputFiles.size(), dataMC, false);
  for (vector<string>::const_iterator iInputFile = inputFiles.begin(); 
       iInputFile != inputFiles.end(); ++iInputFile) {
    const unsigned int fileIndex = iInputFile - inputFiles.begin();
    inputStreams.push_back(new TFile(iInputFile->c_str()));
    for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
	 iCanvasName != canvasNames.end(); ++iCanvasName) {
      const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
      TCanvas* pCanvas;
      inputStreams[inputStreams.size() - 1]->GetObject(iCanvasName->c_str(), pCanvas);
      TGraphAsymmErrors* pGraph = NULL;
      TH1F* pHist = NULL;
      if (graphNames[canvasIndex].find("divide") != string::npos) {
	pGraph = (TGraphAsymmErrors*)pCanvas->GetPrimitive(graphNames[canvasIndex].c_str());
	setGraphOptions(*pGraph, colors[fileIndex], 0.7, styles[fileIndex], 
			pGraph->GetXaxis()->GetTitle(), pGraph->GetYaxis()->GetTitle());
	legends[canvasIndex]->AddEntry(pGraph, legendEntries[fileIndex].c_str(), "l");
      }
      else {
	pHist = (TH1F*)pCanvas->GetPrimitive(graphNames[canvasIndex].c_str());
	float weight = weights[fileIndex] == 0.0 ? 1.0/pHist->Integral(0, -1) : weights[fileIndex];
	setHistogramOptions(pHist, colors[fileIndex], 0.7, styles[fileIndex], 
			    weight, 
			    string(pHist->GetName()) == "muHadPTOverMuHadMass" ? 
			    "p_{T}/m" : pHist->GetXaxis()->GetTitle(), 
			    pHist->GetYaxis()->GetTitle());
	string histName(pHist->GetName());
	if (histName == "jet_pt_etacut") pHist->GetXaxis()->SetTitle("p_{T} (GeV)");
	if (histName == "jet_eta") pHist->GetXaxis()->SetTitle("#eta");
	if (histName == "jet_phi") pHist->GetXaxis()->SetTitle("#phi");
	if (histName == "jet_mass_etacut") pHist->GetXaxis()->SetTitle("m (GeV)");
	if (histName == "jet_ptmj_etacut") {
	  pHist->GetXaxis()->SetTitle("#frac{p_{T}}{m}");
	}
	if (setLogY) pHist->GetYaxis()->SetRangeUser(0.1, 10000.0);
	string legendStyle("l");
	hists[canvasIndex][fileIndex] = pHist;
	legends[canvasIndex]->
	  AddEntry(pHist, legendEntries[fileIndex].c_str(), legendStyle.c_str());
// 	if (fileIndex == (inputFiles.size() - 1)) stacks[canvasIndex]->Add(pHist, "HIST");
// 	else stacks[canvasIndex]->Add(pHist, "HISTE");
	if (fileIndex != 0) stacks[canvasIndex]->Add(pHist, "HIST");
      }
      outStream.cd();
      outputCanvases[canvasIndex]->cd(dataMC ? 1 : 0);
      //outputCanvases[canvasIndex]->cd(0);
      if (fileIndex == 0) { if (pGraph != NULL) pGraph->Draw("AP"); }
      else if (pGraph != NULL) pGraph->Draw("PSAME");
    }
  }
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    outStream.cd();
    outputCanvases[canvasIndex]->cd(dataMC ? 1 : 0);
    TH1F* pHistWithMaxMaxBin = histWithLargestBinContent(hists[canvasIndex]);
    if (pHistWithMaxMaxBin != NULL) pHistWithMaxMaxBin->Draw();
  }
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    outStream.cd();
    outputCanvases[canvasIndex]->cd(dataMC ? 1 : 0);
    //outputCanvases[canvasIndex]->cd(0);
    for (vector<string>::const_iterator iInputFile = inputFiles.begin(); 
	 iInputFile != inputFiles.end(); ++iInputFile) { // loop over input files
      const unsigned int fileIndex = iInputFile - inputFiles.begin();
      if (hists[canvasIndex][fileIndex] != NULL)
	{ // if there is a histogram for this filename and canvas
	  if (fileIndex == 0)
	    { // if data
	      histDiff[canvasIndex] = (TH1F*)hists[canvasIndex][fileIndex]->Clone();
	    } // if data
	  else
	    { // if MC
	      if (inputFiles[fileIndex].find("QCD") != string::npos)
		{ // do nothing
		  if (graphNames[canvasIndex].find("muHadMass") != string::npos)
		    {
		    }
		  
		} // do nothing
	      else
		{
		  if (graphNames[canvasIndex].find("muHadMass") != string::npos)
		    {
		    }
		  histDiff[canvasIndex]->Add(hists[canvasIndex][fileIndex], -1.); /* subtract from 
										     data*/
		}
	    } // if MC
	} // if there is a histogram for this filename and canvas
    } // loop over input files
  }
    
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    outputCanvases[canvasIndex]->cd(dataMC ? 1 : 0);
    //outputCanvases[canvasIndex]->cd(0);
    histDiff[canvasIndex]->Draw();
    //legends[canvasIndex]->Draw("same");
  } 

  outStream.cd();
  write(outputCanvases);
  outStream.Write();
  outStream.Close();
  deleteObjects(legends);
  deleteObjects(stacks);
  deleteObjects(outputCanvases);
  deleteStreams(inputStreams);
}

//get object from canvas
template<typename T>
T* getObjectFromCanvas(TFile& file, const string& objName, const string& canvasName, 
		       const unsigned int pad)
{
  TCanvas* canvas = NULL;
  file.GetObject(canvasName.c_str(), canvas);
  T* obj = NULL;
  if (canvas != NULL) {
    canvas->Draw();
    obj = (T*)canvas->cd(pad)->GetPrimitive(objName.c_str())/*->Clone()*/;
  }
  return obj;
}

//get integral and error from a histogram
pair<Double_t, Double_t> getIntegralAndError(const string& fileName, 
					     const pair<Int_t, Int_t>& bins, const string& var, 
					     const unsigned int pad)
{
  Double_t integral = 0.0;
  Double_t err = 0.0;
  TFile file(fileName.c_str());
  if (file.IsOpen()) {
    TH1F* hist = getObjectFromCanvas<TH1F>(file, var, var + "Canvas", pad);
    if (hist != NULL) {
      integral = hist->IntegralAndError(bins.first, bins.second, err);
    }
    else {
      cerr << "Error getting histogram " << var << " from canvas " << var << "Canvas from file ";
      cerr << fileName << ".\n";
    }
//     delete hist;
//     hist = NULL;
  }
  else cerr << "Error opening file " << fileName << ".\n";
  file.Close();
  return pair<Double_t, Double_t>(integral, err);
}

//compute normalization factor and error from two histograms and a normalization region
pair<Double_t, Double_t> normFactorAndError(const pair<string, string>& numeratorSample, 
					    const pair<string, string>& denominatorSample, 
					    const pair<Int_t, Int_t>& normReg)
{
  //get no. events and stat. error in numerator normalization region
  pair<Double_t, Double_t> numeratorNormRegIntAndErr = 
    getIntegralAndError(numeratorSample.first, normReg, numeratorSample.second, 0);
  Double_t numeratorNormReg = numeratorNormRegIntAndErr.first;
  Double_t numeratorNormRegErr = numeratorNormRegIntAndErr.second;

  //get no. events and stat. error in denominator normalization region
  pair<Double_t, Double_t> denominatorNormRegIntAndErr = 
    getIntegralAndError(denominatorSample.first, normReg, denominatorSample.second, 0);
  Double_t denominatorNormReg = denominatorNormRegIntAndErr.first;
  Double_t denominatorNormRegErr = denominatorNormRegIntAndErr.second;

  //compute normalization factor
  const Double_t normFactor = numeratorNormReg/denominatorNormReg;

  //compute normalization factor error
  Double_t numeratorFracErr2 = 
    numeratorNormReg == 0.0 ? 0.0 : (numeratorNormRegErr/numeratorNormReg);
  numeratorFracErr2*=numeratorFracErr2;
  Double_t denominatorFracErr2 = 
    denominatorNormReg == 0.0 ? 0.0 : (denominatorNormRegErr/denominatorNormReg);
  denominatorFracErr2*=denominatorFracErr2;
  const Double_t normErr = normFactor*sqrt(numeratorFracErr2 + denominatorFracErr2);

  //return
  return pair<Double_t, Double_t>(normFactor, normErr);
}

//get Region A QCD histograms
void drawQCDRegionAHistograms(const string& outputFileA, 
			      const string& inputFileNameB, 
			      const vector<string>& canvasNames, 
			      const vector<string>& graphNames, 
			      const vector<string>& legendHeaders, 
			      const vector<Color_t>& colors, 
			      const vector<Style_t>& styles, 
			      const vector<string>& legendEntries, 
			      const vector<float>& weights, const bool setLogY, 
			      const bool dataMC, const pair<Double_t, Double_t>& SFAndErr)
{ // start routine
  
  if ((canvasNames.size() > graphNames.size()) || (canvasNames.size() > legendHeaders.size())) {
    cerr << "Error: vector size mismatch.\n";
    cerr << colors.size() << endl;
    cerr << styles.size() << endl;
    cerr << legendEntries.size() << endl;
    for (vector<string>::const_iterator i = legendEntries.begin(); i != legendEntries.end(); ++i) {
      cerr << *i << endl;
    }
    cerr << weights.size() << endl;
    return;
  }

  TFile outStream(outputFileA.c_str(), "RECREATE");
  TFile inputFileB(inputFileNameB.c_str(), "read");
  vector<TCanvas*> outputCanvases;
  vector<TLegend*> legends;
  vector<THStack*> stacks;
  vector<vector<TH1F*> > hists;
  setup(canvasNames, outputCanvases, setLogY, legends, stacks, legendHeaders, hists, 3, dataMC, 
	false);

  //compute normalization factor and error from muHadMass histogram
  const Double_t SF = SFAndErr.first;
  const Double_t SFErr = SFAndErr.second;
  Double_t SFFracErr2 = SF == 0.0 ? 0.0 : (SFErr/SF);
  SFFracErr2*=SFFracErr2;

//   //debug
//   cout << "Normalization factor: " << setprecision(3) << SF << " +/- " << SFErr << endl;

  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) { // loop over canvases
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    TCanvas* pCanvasB;
    inputFileB.GetObject(iCanvasName->c_str(), pCanvasB);
    TH1F* pHistB = NULL;
    pHistB = (TH1F*)pCanvasB->GetPrimitive(graphNames[canvasIndex].c_str())->Clone();

    //scale histogram and set total statistical error
    for (Int_t iBin = 0; iBin <= (pHistB->GetNbinsX() + 1); ++iBin) {
      Double_t nBErr = pHistB->GetBinError(iBin);
      Double_t nB = pHistB->GetBinContent(iBin);
      Double_t nBFracStatErr2 = nB == 0.0 ? 0.0 : ((nBErr*nBErr)/(nB*nB));
      Double_t totErr = nB*SF*sqrt(nBFracStatErr2 + SFFracErr2);
      pHistB->SetBinContent(iBin, nB*SF);
      pHistB->SetBinError(iBin, totErr);

//       //debug
//       if (string(pHistB->GetName()) == "muHadMass") {
// 	Double_t statErr = nB*SF*sqrt(nBFracStatErr2);
// 	Double_t normErr = nB*SF*sqrt(SFFracErr2);
// 	Double_t fullFracErr = (totErr/(nB*SF));
// 	Double_t statFracErr = (statErr/(nB*SF));
// 	cout << "Bin " << iBin << endl << setprecision(3) << pHistB->GetBinContent(iBin);
// 	cout << " +/- " << statErr << "(stat.) +/- " << normErr << "(norm.)" << endl;
// 	cout << "Full fractional error: " << fullFracErr << endl;
// 	cout << "Fractional error neglecting normalization term: " << statFracErr << endl;
//       }
    }

    float weight = 1.0;
    setHistogramOptions(pHistB, colors[0], 0.7, styles[0], 
			weight, 
			string(pHistB->GetName()) == "muHadPTOverMuHadMass" ? 
			"p_{T}/m" : pHistB->GetXaxis()->GetTitle(), 
			pHistB->GetYaxis()->GetTitle());
    string histName(pHistB->GetName());
    if (histName == "jet_pt_etacut") pHistB->GetXaxis()->SetTitle("p_{T} (GeV)");
    if (histName == "jet_eta") pHistB->GetXaxis()->SetTitle("#eta");
    if (histName == "jet_phi") pHistB->GetXaxis()->SetTitle("#phi");
    if (histName == "jet_mass_etacut") pHistB->GetXaxis()->SetTitle("m (GeV)");
    if (histName == "jet_ptmj_etacut") {
      pHistB->GetXaxis()->SetTitle("#frac{p_{T}}{m}");
    }
    if (setLogY) pHistB->GetYaxis()->SetRangeUser(0.1, 10000.0);
    string legendStyle("l");
    legends[canvasIndex]->
      AddEntry(pHistB, legendEntries[0].c_str(), legendStyle.c_str());
    outStream.cd();
    outputCanvases.at(canvasIndex)->cd(dataMC ? 1 : 0);
    pHistB->Draw();
  } // loop over canvases

  outStream.cd();
  write(outputCanvases);
  outStream.Write();
  outStream.Close();
  inputFileB.Close();
  deleteObjects(legends);
  deleteObjects(stacks);
  deleteObjects(outputCanvases);
} // end routine

//replace muHadMass histogram in region A QCD-only file with the resonance-subtracted estimate
void setRegAQCDMuHadMassEstToResSubtrRegC(const string& regAQCDFileName, 
					  const string& resBkgFileName)
{
  //clone resonance-subtracted histogram out of resonance background file
  TH1F* nonResRegCHist = NULL;
  TFile resBkgFile(resBkgFileName.c_str());
  if (resBkgFile.IsOpen()) {
    TIter iKey(resBkgFile.GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)iKey()) && (nonResRegCHist == NULL)) {
      string objName(key->GetName());
      if (objName.find("weightedAvgkFixed_nonResBkg") != string::npos) {
	TH1F* obj = NULL;
	resBkgFile.GetObject(objName.c_str(), obj);
	if (obj != NULL) {
	  nonResRegCHist = (TH1F*)obj->Clone();
	  nonResRegCHist->SetName("muHadMass");
	  nonResRegCHist->SetDirectory(0);
	}
      }
    }
    if (nonResRegCHist == NULL) {
      cerr << "Error finding histogram with name containing the phrase ";
      cerr << "\"weightedAvgkFixed_nonResBkg\" from file " << resBkgFileName << ".\n";
      resBkgFile.ls();
    }
  }
  else cerr << "Error opening file " << resBkgFileName << ".\n";
  resBkgFile.Close();

  //write resonance-subtracted histogram over existing histogram on muHadMass canvas in region A 
  //QCD file
  TFile regAQCDFile(regAQCDFileName.c_str(), "UPDATE");
  if (regAQCDFile.IsOpen()) {
    regAQCDFile.cd();
    TCanvas resBkgCanvas("muHadMassCanvas", "muHadMassCanvas", 600, 600);
    setCanvasOptions(resBkgCanvas, 0, 0, 0);
    resBkgCanvas.cd();
    nonResRegCHist->Draw();
    resBkgCanvas.Write("muHadMassCanvas", TObject::kOverwrite);
    regAQCDFile.Write();
  }
  else cerr << "Error opening file " << regAQCDFileName << ".\n";
  regAQCDFile.Close();
}

//test data-driven background estimation method with MC on a given variable
void addClosurePlot(TFile& sigVsBkgIsoStream, const string& var, const string& unit, 
                    TFile& bkgNonIsoStream, const float nonIsoScale, 
                    const int normRegionLowerBin, const int normRegionUpperBin, TFile& outStream)
{
  //top level declarations
  string canvasName(var + "Canvas");
  string stackName(var + "Stack");

  //get plots of signal and background MC, isolated tau sample
  TCanvas* canvasIso = NULL;
  sigVsBkgIsoStream.GetObject(canvasName.c_str(), canvasIso);

  //get the signal histograms
  TH1F* histSig1 = NULL;
  THStack* stackBkgIso = NULL;
  if (canvasIso != NULL) {
    TList* sigs = canvasIso->GetListOfPrimitives();
    histSig1 = (TH1F*)sigs->At(2)->Clone();
    histSig1->GetYaxis()->SetRangeUser(0.1, 10000.0);

    /*get the background stack histogram in the signal region (isolated taus) and convert it to a 
      TH1*/
    stackBkgIso = (THStack*)canvasIso->GetPrimitive(stackName.c_str());
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << sigVsBkgIsoStream.GetName() << ".\n";
    return;
  }
  TList* stackedHistsIso = NULL;
  if (stackBkgIso != NULL) stackedHistsIso = stackBkgIso->GetHists();
  else {
    cerr << "Error opening stack " << stackName << " from canvas " << canvasName << " from file ";
    cerr << sigVsBkgIsoStream.GetName() << ".\n";
    return;
  }
  double sum2 = 0.;
  double err2 = 0.;
  double sum4 = 0.;
  double err4 = 0.;
  TH1F* histBkgIso = NULL;
  bool isItMass = false;
  if (stackedHistsIso != NULL) {
    for (Int_t i = 0; i < stackedHistsIso->GetEntries(); ++i) {
      TH1F* stackHist = (TH1F*)stackedHistsIso->At(i)->Clone();
      if (canvasName.find("muHadMassCanvas") != string::npos)
	{
	  isItMass = true;
	  Double_t sampleErr2 = 0.0;
	  Double_t sampleErr4 = 0.0;
	  sum2 += stackHist->IntegralAndError(1, 2/*8*/, sampleErr2);
	  sum4 += stackHist->IntegralAndError(5/*17*/, -1, sampleErr4);
	  err2+=(sampleErr2*sampleErr2);
	  err4+=(sampleErr4*sampleErr4);
	}
      if (i == 0) histBkgIso = stackHist;
      else histBkgIso->Add(stackHist);
    }
    if (isItMass)
      {
	cout << "Total MC in Region A (m < 2 GeV): " << setprecision(3) << sum2 << " +/- ";
	cout << sqrt(err2) << endl;
	cout << "Total MC in Region A (m >= 4 GeV): " << setprecision(3) << sum4 << " +/- ";
	cout << sqrt(err4) << endl;
	Double_t ratioRegA = sum4 == 0.0 ? 0.0 : sum2/sum4;
	Double_t mLe2GeVFracErr2RegA = sum2 == 0.0 ? 0.0 : err2/(sum2*sum2);
	Double_t mGeq4GeVFracErr2RegA = sum4 == 0.0 ? 0.0 : err4/(sum4*sum4);
	cout << "Total MC in Region A (m < 2 GeV)/(m >= 4 GeV): " << setprecision(3) << ratioRegA;
	cout << " +/- " << ratioRegA*sqrt(mLe2GeVFracErr2RegA + mGeq4GeVFracErr2RegA) << endl;
      }
    setHistogramOptions(histBkgIso, kBlue, 0.7, 21, 1.0, unit.c_str(), "");
    histBkgIso->GetYaxis()->SetRangeUser(0.1, 10000.0);
  }
  else {
    cerr << "Error opening histogram list from stack " << stackName << " from canvas ";
    cerr << canvasName << " from file " << sigVsBkgIsoStream.GetName() << ".\n";
    return;
  }

  //get plots of background MC, non-isolated tau sample
  TCanvas* canvasNonIso = NULL;
  bkgNonIsoStream.GetObject(canvasName.c_str(), canvasNonIso);
  THStack* stackBkgNonIso = NULL;
  if (canvasNonIso != NULL) {
    canvasNonIso->Draw();
    stackBkgNonIso = (THStack*)canvasNonIso->cd(1)->GetPrimitive(stackName.c_str());
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file " << bkgNonIsoStream.GetName();
    cerr << ".\n";
    return;
  }
  TList* stackedHistsNonIso = NULL;
  if (stackBkgNonIso != NULL) stackedHistsNonIso = stackBkgNonIso->GetHists();
  else {
    cerr << "Error opening stack " << stackName << " from canvas " << canvasName << " from file ";
    cerr << bkgNonIsoStream.GetName() << ".\n";
    return;
  }
  TH1F* histBkgNonIso = NULL;
  double sum2B = 0.;
  double err2B = 0.;
  double sum4B = 0.;
  double err4B = 0.;
  bool isItMassB = false;
  if (stackedHistsNonIso != NULL) {
    for (Int_t i = 0; i < stackedHistsNonIso->GetEntries(); ++i) {
      TH1F* stackHist = (TH1F*)stackedHistsNonIso->At(i)->Clone();
      if (canvasName.find("muHadMassCanvas") != string::npos)
	{
	  isItMassB = true;
	  Double_t sampleErr2 = 0.0;
	  Double_t sampleErr4 = 0.0;
	  sum2B += stackHist->IntegralAndError(1, 2/*8*/, sampleErr2);
	  sum4B += stackHist->IntegralAndError(5/*17*/, -1, sampleErr4);
	  err2B+=(sampleErr2*sampleErr2);
	  err4B+=(sampleErr4*sampleErr4);
	}
      if (i == 0) histBkgNonIso = stackHist;
      else histBkgNonIso->Add(stackHist);
    }
    if (isItMassB)
      {
	cout << "Total MC in Region B (m < 2 GeV): " << setprecision(3) << sum2B << " +/- ";
	cout << sqrt(err2B) << endl;
	cout << "Total MC in Region B (m >= 4 GeV): " << setprecision(3) << sum4B << " +/- ";
	cout << sqrt(err4B) << endl;
	Double_t ratioRegB = sum4B == 0.0 ? 0.0 : sum2B/sum4B;
	Double_t mLe2GeVFracErr2RegB = sum2B == 0.0 ? 0.0 : err2B/(sum2B*sum2B);
	Double_t mGeq4GeVFracErr2RegB = sum4B == 0.0 ? 0.0 : err4B/(sum4B*sum4B);
	cout << "Total MC in Region B (m < 2 GeV)/(m >= 4 GeV): " << setprecision(3) << ratioRegB;
	cout << " +/- " << ratioRegB*sqrt(mLe2GeVFracErr2RegB + mGeq4GeVFracErr2RegB) << endl;
      }
    histBkgNonIso->Scale(nonIsoScale);
  }
  else {
    cerr << "Error opening histogram list from stack " << stackName << " from canvas ";
    cerr << canvasName << " from file " << bkgNonIsoStream.GetName() << ".\n";
    return;
  }

  /*normalize non-isolated background MC histogram to isolated background MC in signal-depleted 
    region*/
  histBkgNonIso->Scale(histBkgIso->Integral(normRegionLowerBin, normRegionUpperBin)/
                       histBkgNonIso->Integral(normRegionLowerBin, normRegionUpperBin));
  setHistogramOptions(histBkgNonIso, kRed, 0.7, 20, 1.0, unit.c_str(), "");
  histBkgNonIso->GetYaxis()->SetRangeUser(0.1, 10000.0);

  //write to file
  outStream.cd();
  TCanvas outCanvas(canvasName.c_str(), "", 600, 900);
  outCanvas.Divide(1, 2);
  outCanvas.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  setCanvasOptions(*outCanvas.cd(1), 1, 1, 0);
  outCanvas.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  setCanvasOptions(*outCanvas.cd(2), 1, 0, 0);
  outCanvas.cd(1);
  histBkgIso->Draw("HISTE");
  histBkgNonIso->Draw("HISTESAME");
  outCanvas.cd(2);
  TH1F* isoOverNonIso = (TH1F*)histBkgIso->Clone();
  isoOverNonIso->Divide(histBkgNonIso);
  isoOverNonIso->GetYaxis()->SetTitle("#frac{MC (A)}{MC (B)}");
  isoOverNonIso->GetYaxis()->SetRangeUser(0.0, 2.0);
  isoOverNonIso->Draw();
  outCanvas.Write();
}

//test data-driven background estimation method with MC
void makeMCClosurePlots(const string& sigVsBkgIsoFileName, const vector<string>& vars, 
			const vector<string>& units, const string& bkgNonIsoFileName, 
			const float nonIsoScale, const vector<int>& normRegionLowerBins, 
			const vector<int>& normRegionUpperBins, const string& outputFileName)
{
  //open files
  TFile sigVsBkgIsoStream(sigVsBkgIsoFileName.c_str());
  TFile bkgNonIsoStream(bkgNonIsoFileName.c_str());
  TFile outStream(outputFileName.c_str(), "RECREATE");
  if (sigVsBkgIsoStream.IsOpen() && bkgNonIsoStream.IsOpen() && outStream.IsOpen()) {
    for (vector<string>::const_iterator iVar = vars.begin(); iVar != vars.end(); ++iVar) {
      const unsigned int varIndex = iVar - vars.begin();
      addClosurePlot(sigVsBkgIsoStream, *iVar, units[varIndex], bkgNonIsoStream, 
		     nonIsoScale, normRegionLowerBins[varIndex], normRegionUpperBins[varIndex], 
		     outStream);
    }
  }
  else {
    cerr << "Error opening files " << sigVsBkgIsoFileName << " or " << bkgNonIsoFileName;
    cerr << " or " << outputFileName << ".\n";
    return;
  }

  //write to file
  outStream.Write();
  outStream.Close();
  sigVsBkgIsoStream.Close();
  bkgNonIsoStream.Close();
}

void QCDVsMCClosurePlots(const vector<string>& QCDVsMCInputFileNames, const string& var, 
			 const string& units, const pair<string, string>& legend, 
			 const int normRegionLowerBin, const int normRegionUpperBin, 
			 const Double_t xMin, const Double_t xMax, 
			 const string& outputFileName, const string& outputFileMode)
{

  TFile outStream(outputFileName.c_str(), outputFileMode.c_str());
  vector<TFile*> inputStreams;
  vector<TCanvas*> outputCanvases;
  vector<TH1F*> hists(2);
  string canvasName(var + "Canvas");

  for (vector<string>::const_iterator iInputFile = QCDVsMCInputFileNames.begin(); 
       iInputFile != QCDVsMCInputFileNames.end(); ++iInputFile) { // loop over input files
    const unsigned int fileIndex = iInputFile - QCDVsMCInputFileNames.begin();
    inputStreams.push_back(new TFile(iInputFile->c_str()));
    TCanvas* pCanvas;
    inputStreams[inputStreams.size() - 1]->GetObject(canvasName.c_str(), pCanvas);
    TH1F* pHist = (TH1F*)pCanvas->GetPrimitive(var.c_str());
    if (fileIndex == 0)
      { // if this is the QCD file
	hists[0] = (TH1F*)pHist->Clone();
      } // if this is the QCD file
    else if (fileIndex == 1)
      { // if this is the first MC bkg file
	hists[1] = (TH1F*)pHist->Clone();
      } // if this is the first MC bkg file
    else
      { // if this is another MC bkg file
	hists[1]->Add(pHist,1.);
      } // if this is another MC bkg file
    
  } // loop over input files

  hists[0]->SetLineColor(4); // blue for QCD
  hists[1]->SetLineColor(2); // red for bkg
  hists[0]->SetMarkerColor(4); // blue for QCD
  hists[1]->SetMarkerColor(2); // red for bkg
  hists[0]->Scale(hists[1]->Integral(normRegionLowerBin,normRegionUpperBin)/
		  hists[0]->Integral(normRegionLowerBin,normRegionUpperBin));
  //hists[0]->Scale(1./hists[0]->Integral());
  //hists[1]->Scale(1./hists[1]->Integral());
  hists[1]->GetXaxis()->SetTitle(units.c_str());

  TH1F* ratioHist = (TH1F*)hists[0]->Clone();
//   ratioHist->Add(hists[1],-1.);
  ratioHist->Divide(hists[1]);
  ratioHist->GetYaxis()->SetTitle(/*"Fractional difference"*/"Ratio");
  ratioHist->GetYaxis()->SetRangeUser(/*-2.0,2.0*/0.0, 2.0);
  ratioHist->SetLineColor(1);
  ratioHist->SetMarkerColor(1);

  TLegend *leg = new TLegend(0.55, 0.65, 0.95, 0.85, "");
  setLegendOptions(*leg, "");
  leg->AddEntry(hists[0], legend.first.c_str(), "lp");
  leg->AddEntry(hists[1], legend.second.c_str(), "lp");

  //find out which histogram has the maximum bin content so that it can be drawn first
  TH1F* pHistWithMaxMaxBin = histWithLargestBinContent(hists);

  //write to file
  outStream.cd();
  TCanvas outCanvas(canvasName.c_str(), "", 600, 600);
  outCanvas.Divide(1,2);
  outCanvas.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  outCanvas.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  outCanvas.cd(1);
  pHistWithMaxMaxBin->Draw("HISTE");
  hists[0]->Draw("HISTESAME");
  hists[1]->Draw("HISTESAME");
  pHistWithMaxMaxBin->GetXaxis()->SetRangeUser(xMin, xMax);
  hists[0]->GetXaxis()->SetRangeUser(xMin, xMax);
  hists[1]->GetXaxis()->SetRangeUser(xMin, xMax);
  leg->Draw();
  outCanvas.cd(2);
  ratioHist->Draw();
  ratioHist->GetXaxis()->SetRangeUser(xMin, xMax);
  outCanvas.Write();
  outStream.Write();
  outStream.Close();
}

void compareTotalMCBToA(const vector<string>& QCDVsMCInputFileNames1,
			    const vector<string>& QCDVsMCInputFileNames2,
			    const string& var, 
			    const string& units, 
			    const int normRegionLowerBin, 
			    const int normRegionUpperBin, const string& outputFileName)
{

  TFile outStream(outputFileName.c_str(), "RECREATE");
  vector<TFile*> inputStreams;
  vector<TCanvas*> outputCanvases;
  vector<TH1F*> hists(2);
  string canvasName(var + "Canvas");

  // loop over Region B files
  for (vector<string>::const_iterator iInputFile = QCDVsMCInputFileNames1.begin(); 
       iInputFile != QCDVsMCInputFileNames1.end(); ++iInputFile) { // loop over input files
    const unsigned int fileIndex = iInputFile - QCDVsMCInputFileNames1.begin();
    inputStreams.push_back(new TFile(iInputFile->c_str()));
    TCanvas* pCanvas;
    inputStreams[inputStreams.size() - 1]->GetObject(canvasName.c_str(), pCanvas);
    TH1F* pHist = (TH1F*)pCanvas->GetPrimitive(var.c_str());
    if (fileIndex == 0)
      { // if this is the first file
	hists[0] = (TH1F*)pHist->Clone();
      } // if this is the first file
    else
      { // if this is not the first file
	hists[0]->Add(pHist,1.);
      } // if this is not the first file
    
  } // loop over input files

  // loop over Region A files
  for (vector<string>::const_iterator iInputFile = QCDVsMCInputFileNames2.begin(); 
       iInputFile != QCDVsMCInputFileNames2.end(); ++iInputFile) { // loop over input files
    const unsigned int fileIndex = iInputFile - QCDVsMCInputFileNames2.begin();
    inputStreams.push_back(new TFile(iInputFile->c_str()));
    TCanvas* pCanvas;
    inputStreams[inputStreams.size() - 1]->GetObject(canvasName.c_str(), pCanvas);
    TH1F* pHist = (TH1F*)pCanvas->GetPrimitive(var.c_str());
    if (fileIndex == 0)
      { // if this is the first file
	hists[1] = (TH1F*)pHist->Clone();
      } // if this is the first file
    else
      { // if this is not the first file
	hists[1]->Add(pHist,1.);
      } // if this is not the first file
    
  } // loop over input files


  hists[0]->SetLineColor(4); // blue for Region B
  hists[1]->SetLineColor(2); // red for Region A
  hists[0]->Scale(hists[1]->Integral(normRegionLowerBin,normRegionUpperBin)/
		  hists[0]->Integral(normRegionLowerBin,normRegionUpperBin)); // scale B to A
  //hists[0]->Scale(1./hists[0]->Integral());
  //hists[1]->Scale(1./hists[1]->Integral());
  hists[1]->GetXaxis()->SetTitle(units.c_str());

  TH1F* ratioHist = (TH1F*)hists[0]->Clone();
//   ratioHist->Add(hists[1],-1.);
  ratioHist->Divide(hists[1]);
  ratioHist->GetYaxis()->SetTitle(/*"#frac{B-A}{A}"*/"Reg. B / Reg. A");
  ratioHist->GetYaxis()->SetRangeUser(/*-2.0,2.0*/0.0, 2.0);
  ratioHist->SetLineColor(1);

  TLegend *leg = new TLegend(0.35, 0.55, 0.75, 0.75, "");
  leg->AddEntry(hists[0], "Total MC + data-driven QCD, Region B", "lp");
  leg->AddEntry(hists[1], "Total MC + data-driven QCD, Region A", "lp");

  //write to file
  outStream.cd();
  TCanvas outCanvas(canvasName.c_str(), "", 600, 600);
  outCanvas.Divide(1,2);
  outCanvas.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  outCanvas.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  outCanvas.cd(1);
  hists[0]->Draw("HISTE");
  hists[1]->Draw("HISTESAME");
  leg->Draw();
  outCanvas.cd(2);
  ratioHist->Draw();

  outCanvas.Write();
  outStream.Write();
  outStream.Close();
}
  
/*make final plot showing:
  - Wh1 expected signal
  - gg expected signal
  - full background estimate from region B data
  - background estimate from MC
  - QCD estimate from region B/C/D data*/
/*options for displaying MC background:
  - each sample separately
  - WW/WZ/ZZ combined into diboson, tt and single top combined into top
  - all samples combined
 */
void addFinalPlot(pair<TFile*, float>& isoSigBkgFile, TFile& isoDataFile, 
		  pair<TFile*, float>& nonIsoDataFile, const string& var, const string& unit, 
		  const int normRegionLowerBin, const int normRegionUpperBin, 
		  const string& option, TFile& outStream, const bool ma9GeV)
{
  //top level declarations
  string canvasName(var + "Canvas");
  string stackName(var + "Stack");
  string canvasReweightErrSqName(var + "ReweightErrSqCanvas");
  string varReweightErrSq(var + "ReweightErrSq");

  //get plots of signal and background MC and data-driven QCD estimate, isolated tau sample
  TCanvas* canvasIsoSigBkg = NULL;
  isoSigBkgFile.first->GetObject(canvasName.c_str(), canvasIsoSigBkg);

  //get the signal histograms
  //for ma = 9 GeV, there are 4 signals
  //for all other pseudoscalar masses, there are 2 signals
  vector<TH1F*> isoSig(4, NULL);
  if (!ma9GeV) isoSig.erase(isoSig.begin() + 2, isoSig.begin() + 4);
  THStack* isoBkg = NULL;
  TLegend legendBkgSep(0.35, 0.55, 0.75, 0.95);
  TLegend legendBkgMain5(0.35, 0.55, 0.75, 0.95);
  TLegend legendBkgAll(0.35, 0.55, 0.75, 0.95);
  setLegendOptions(legendBkgSep, "CMS 19.7 fb^{-1}");
  setLegendOptions(legendBkgMain5, "CMS 19.7 fb^{-1}");
  setLegendOptions(legendBkgAll, "CMS 19.7 fb^{-1}");
  if (canvasIsoSigBkg != NULL) {
    TList* sigs = canvasIsoSigBkg->GetListOfPrimitives();
    for (vector<TH1F*>::iterator iIsoSig = isoSig.begin(); iIsoSig != isoSig.end(); 
	 ++iIsoSig) {
      const unsigned int i = iIsoSig - isoSig.begin();
      *iIsoSig = (TH1F*)sigs->At(i + 2)->Clone();
      (*iIsoSig)->GetYaxis()->SetRangeUser(0.01, 10000.0);
      for (int b = 5/*17*/; b < (*iIsoSig)->GetNbinsX(); ++b)
	{
	  cout << "stat error on signal sample " << i << " = " << (*iIsoSig)->GetBinError(b);
	  cout << endl;
	}
    }
    legendBkgSep.AddEntry(isoSig[0], "Wh_{1}", "l");
    legendBkgSep.AddEntry(isoSig[1], "gg fusion", "l");
    legendBkgMain5.AddEntry(isoSig[0], "Wh_{1}", "l");
    legendBkgMain5.AddEntry(isoSig[1], "gg fusion", "l");
    legendBkgAll.AddEntry(isoSig[0], "Wh_{1}", "l");
    legendBkgAll.AddEntry(isoSig[1], "gg fusion", "l");
    isoSig[0]->SetName((var + "Wh1Search").c_str());
    isoSig[1]->SetName((var + "ggSearch").c_str());
    setHistogramOptions(isoSig[0], kSpring - 1, 0.7, 20, isoSigBkgFile.second, unit.c_str(), "");
    setHistogramOptions(isoSig[1], kBlue, 0.7, 20, isoSigBkgFile.second, unit.c_str(), "");
    isoSig[0]->SetLineWidth(4);
    isoSig[0]->SetLineStyle(2);
    isoSig[1]->SetLineWidth(4);
    isoSig[1]->SetLineStyle(2);
    if (ma9GeV) {
      legendBkgSep.AddEntry(isoSig[2], "ZH", "l");
      legendBkgSep.AddEntry(isoSig[3], "VBF", "l");
      legendBkgMain5.AddEntry(isoSig[2], "ZH", "l");
      legendBkgMain5.AddEntry(isoSig[3], "VBF", "l");
      legendBkgAll.AddEntry(isoSig[2], "ZH", "l");
      legendBkgAll.AddEntry(isoSig[3], "VBF", "l");
      isoSig[2]->SetName((var + "ZHSearch").c_str());
      isoSig[3]->SetName((var + "VBFSearch").c_str());
      setHistogramOptions(isoSig[2], kRed, 0.7, 20, isoSigBkgFile.second, unit.c_str(), "");
      setHistogramOptions(isoSig[3], kMagenta, 0.7, 20, isoSigBkgFile.second, unit.c_str(), "");
      isoSig[2]->SetLineWidth(4);
      isoSig[2]->SetLineStyle(2);
      isoSig[3]->SetLineWidth(4);
      isoSig[3]->SetLineStyle(2);
    }

    /*get the MC + data-driven QCD background stack histogram and make some new stacks with 
      different combinations of the backgrounds*/
    isoBkg = (THStack*)canvasIsoSigBkg->GetPrimitive(stackName.c_str());
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << isoSigBkgFile.first->GetName() << ".\n";
    return;
  }
  TList* isoBkgHists = NULL;
  if (isoBkg != NULL) isoBkgHists = isoBkg->GetHists();
  else {
    cerr << "Error opening stack " << stackName << " from canvas " << canvasName << " from file ";
    cerr << isoSigBkgFile.first->GetName() << ".\n";
    return;
  }
  TH1F* isoBkgAllHist = NULL;
  TH1F* isoBkgDibosonHist = NULL;
  TH1F* isoBkgWNJetsHist = NULL;
  TH1F* isoBkgTopHist = NULL;
  TH1F* isoBkgDrellYanHist = NULL;
  TH1F* isoBkgQCDHist = NULL;
  THStack isoBkgSep("isoBkgSep", "");
  THStack isoBkgAll("isoBkgAll", "");
  THStack isoBkgMain5("isoBkgMain5", "");
  if (isoBkgHists != NULL) {
    for (Int_t i = 0; i < isoBkgHists->GetEntries(); ++i) {
      TH1F* stackHist = (TH1F*)isoBkgHists->At(i)->Clone();
      stackHist->Scale(isoSigBkgFile.second);
      stackHist->GetYaxis()->SetRangeUser(0.01, 10000.0);
      isoBkgSep.Add(stackHist, "HIST");
      Double_t err = 0.0;
      Double_t val = stackHist->IntegralAndError(5/*17*/, -1, err);
      if (i == 0) {
	isoBkgAllHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgAllHist->Scale(isoSigBkgFile.second);
	isoBkgDibosonHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgDibosonHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "WW (MC)", "f");
	cout << "prediction from WW: " << val << " +/- " << err << endl;

      }
      else isoBkgAllHist->Add(stackHist);
      if (i == 1)
	{
	  legendBkgSep.AddEntry(stackHist, "ZZ (MC)", "f");
	  cout << "prediction from ZZ: " << val << " +/- " << err << endl;
	}
      if (i == 2)
	{
	  legendBkgSep.AddEntry(stackHist, "WZ (MC)", "f");
	  cout << "prediction from WZ: " << val << " +/- " << err << endl;
	}
      if ((i == 1) || (i == 2)) isoBkgDibosonHist->Add(stackHist);
      if (i == 3) {
	isoBkgWNJetsHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgWNJetsHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "W + #geq1 jet (MC)", "f");
	cout << "prediction from WNJets: " << val << " +/- " << err << endl;
      }
      if (i == 4) {
	isoBkgTopHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgTopHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "t/#bar{t} (MC)", "f");
	cout << "prediction from t/tbar: " << val << " +/- " << err << endl;
      }
      if (i == 5) {
	isoBkgTopHist->Add(stackHist);
	legendBkgSep.AddEntry(stackHist, "t#bar{t} + jets (MC)", "f");
	cout << "prediction from TTJets: " << val << " +/- " << err << endl;
      }
      if (i == 6) {
	isoBkgDrellYanHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgDrellYanHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "Drell-Yan + jets (MC)", "f");
	cout << "prediction from DY: " << val << " +/- " << err << endl;
      }
      if (i == 7) {
	isoBkgQCDHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgQCDHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "QCD (data)", "f");
	cout << "prediction from QCD (data): " << val << " +/- " << err << endl;
      }
    }
    isoBkgAllHist->GetYaxis()->SetRangeUser(0.01, 10000.0);
    isoBkgAll.Add(isoBkgAllHist, "HISTE C");
    isoBkgMain5.Add(isoBkgDibosonHist, "HIST");
    isoBkgMain5.Add(isoBkgWNJetsHist, "HIST");
    isoBkgMain5.Add(isoBkgTopHist, "HIST");
    isoBkgMain5.Add(isoBkgDrellYanHist, "HIST");
    isoBkgMain5.Add(isoBkgQCDHist, "HIST");
    legendBkgAll.AddEntry(isoBkgAllHist, "MC EW + data QCD", "f");
    legendBkgMain5.AddEntry(isoBkgDibosonHist, "MC diboson", "f");
    legendBkgMain5.AddEntry(isoBkgWNJetsHist, "MC W + #geq1 jet", "f");
    legendBkgMain5.AddEntry(isoBkgTopHist, "MC top + jets", "f");
    legendBkgMain5.AddEntry(isoBkgDrellYanHist, "MC Drell-Yan + jets", "f");
    legendBkgMain5.AddEntry(isoBkgQCDHist, "Data QCD", "f");
  }
  else {
    cerr << "Error opening histogram list from stack " << stackName << " from canvas ";
    cerr << canvasName << " from file " << isoSigBkgFile.first->GetName() << ".\n";
    return;
  }

  //get the data histogram, non-isolated tau sample
  TCanvas* canvasNonIsoData = NULL;
  nonIsoDataFile.first->GetObject(canvasName.c_str(), canvasNonIsoData);
  TH1F* nonIsoData = NULL;
  if (canvasNonIsoData != NULL) {
    canvasNonIsoData->Draw();
    nonIsoData = (TH1F*)canvasNonIsoData->cd(1)->GetPrimitive(var.c_str())->Clone();
    nonIsoData->SetName((var + "DataControl").c_str());
    setHistogramOptions(nonIsoData, kRed, 0.7, 20, nonIsoDataFile.second, unit.c_str(), "");
    nonIsoData->GetYaxis()->SetRangeUser(0.01, 10000.0);
    legendBkgSep.AddEntry(nonIsoData, "Background (from data)", "lp");
    legendBkgMain5.AddEntry(nonIsoData, "Background (from data)", "lp");
    legendBkgAll.AddEntry(nonIsoData, "Background (from data)", "lp");
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << nonIsoDataFile.first->GetName() << ".\n";
    return;
  }

  //get the reweighting error histogram, non-isolated tau sample
  TCanvas* canvasNonIsoDataReweightErrSq = NULL;
  nonIsoDataFile.first->GetObject(canvasReweightErrSqName.c_str(), canvasNonIsoDataReweightErrSq);
  TH1F* nonIsoDataReweightErrSq = NULL;
  if (canvasNonIsoDataReweightErrSq != NULL) {
    canvasNonIsoDataReweightErrSq->Draw();
    nonIsoDataReweightErrSq = 
      (TH1F*)canvasNonIsoDataReweightErrSq->cd(1)->GetPrimitive(varReweightErrSq.c_str())->Clone();
  }
  else {
    cerr << "Error opening canvas " << canvasReweightErrSqName << " from file ";
    cerr << nonIsoDataFile.first->GetName() << ".\n";
    return;
  }

  //get the data histogram, isolated tau sample
  TCanvas* canvasIsoData = NULL;
  isoDataFile.GetObject(canvasName.c_str(), canvasIsoData);
  TH1F* isoData = NULL;
  if (canvasIsoData != NULL) {
    isoData = (TH1F*)canvasIsoData->GetPrimitive(var.c_str())->Clone();
    isoData->SetName((var + "DataSearch").c_str());
    legendBkgSep.AddEntry(isoData, "Data", "p");
    legendBkgMain5.AddEntry(isoData, "Data", "p");
    legendBkgAll.AddEntry(isoData, "Data", "p");
    isoData->GetYaxis()->SetRangeUser(0.01, 10000.0);
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file " << isoDataFile.GetName();
    cerr << ".\n";
    return;
  }

  //calculate the normalization factor
  Double_t norm = isoData->Integral(normRegionLowerBin, normRegionUpperBin)/
    nonIsoData->Integral(normRegionLowerBin, normRegionUpperBin);
  cout << "The normalization constant is: " << norm << endl;

  /*calculate the statistical error on the background prediction from the non-isolated data, 
    including the term from the error on the normalization factor*/
  const TH1* nonIsoDataPtrCast = dynamic_cast<const TH1*>(nonIsoData);
  const TH1* nonIsoDataReweightErrSqPtrCast = dynamic_cast<const TH1*>(nonIsoDataReweightErrSq);
  vector<Double_t> nonIsoDataStatErrSq;
  for (Int_t iBin = 1; iBin <= nonIsoData->GetNbinsX(); ++iBin) {
    nonIsoDataStatErrSq.
      push_back(bkgErrSq(nonIsoDataPtrCast, nonIsoDataReweightErrSqPtrCast, iBin, norm, 
			 normErrSq(dynamic_cast<const TH1*>(isoData), nonIsoDataPtrCast, 
				   nonIsoDataReweightErrSqPtrCast, normRegionLowerBin, 
				   normRegionUpperBin, norm)));
  }

  //normalize non-isolated data histogram to isolated data in signal-depleted region
  cout << "Region B m > 4 integral before normalization = " << nonIsoData->Integral(5/*17*/,-1);
  cout << endl;
  nonIsoData->Scale(norm);

  //set statistical error in each bin of the non-isolated data histogram
  for (Int_t iBin = 1; iBin <= nonIsoData->GetNbinsX(); ++iBin) {
    nonIsoData->SetBinError(iBin, sqrt(nonIsoDataStatErrSq[iBin - 1]));
  }

  //write to file
  outStream.cd();
  TCanvas outCanvas(canvasName.c_str(), "", 600, 900);
  outCanvas.Divide(1, 2);
  outCanvas.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  setCanvasOptions(*outCanvas.cd(1), 1, 1, 0);
  outCanvas.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  setCanvasOptions(*outCanvas.cd(2), 1, 0, 0);
  outCanvas.cd(1);
  Double_t B = 0.;
  if (option == "separate") {
    isoBkgSep.Draw();
    isoBkgSep.SetMinimum(0.01);
    isoBkgSep.SetMaximum(10000.0);
    isoBkgSep.GetHistogram()->GetXaxis()->SetTitle(unit.c_str());
  }
  else if (option == "main 5") {
    isoBkgMain5.Draw();
    Float_t sum = 0.0;
    Float_t err = 0.0;
    for (Int_t iHist = 0; iHist < isoBkgMain5.GetHists()->GetEntries(); ++iHist) {
      TH1F* hist = (TH1F*)isoBkgMain5.GetHists()->At(iHist);
      for (Int_t iBin = 5/*17*/; iBin <= (hist->GetNbinsX() + 1); ++iBin) {
	sum+=hist->GetBinContent(iBin);
	err+=(hist->GetBinError(iBin)*hist->GetBinError(iBin));
      }
    }
    B = sum;
    
    cout << "Region A MC + data-driven QCD, m > 4: " << setprecision(3) << sum << " +/- ";
    cout << setprecision(3) << sqrt(err) <<  endl;
    isoBkgMain5.SetMinimum(0.01);
    isoBkgMain5.SetMaximum(10000.0);
    isoBkgMain5.GetHistogram()->GetXaxis()->SetTitle(unit.c_str());
  }
  else if (option == "combined") {
    isoBkgAll.Draw();
    isoBkgAll.SetMinimum(0.01);
    isoBkgAll.SetMaximum(10000.0);
    isoBkgAll.GetHistogram()->GetXaxis()->SetTitle(unit.c_str());
  }
  nonIsoData->Draw("HISTESAME");
  Double_t statErrNonIsoData;
  nonIsoData->IntegralAndError(5/*17*/, -1, statErrNonIsoData);
  cout << "Region B data, m > 4: " << setprecision(3) << nonIsoData->Integral(5/*17*/, -1);
  cout << " +/- " << setprecision(3) << statErrNonIsoData << endl;
  for (vector<TH1F*>::iterator iIsoSig = isoSig.begin(); iIsoSig != isoSig.end(); 
       ++iIsoSig) {
    (*iIsoSig)->Draw("HISTSAME");
    Double_t statErrIsoSig;
    Double_t S = (*iIsoSig)->IntegralAndError(5/*17*/, -1, statErrIsoSig);
    cout << "Region A signal " << iIsoSig - isoSig.begin() << ", m > 4: ";
    cout << setprecision(3) << S << " +/- " << setprecision(3) << statErrIsoSig;
    cout << " (S/sqrt(S+B) = " << S/sqrt(S + B) << ")" << endl;
  }
  isoData->Draw("ESAME");
  Double_t regANormErr;
  Double_t regANorm = 
    isoData->IntegralAndError(normRegionLowerBin, normRegionUpperBin, regANormErr);
  cout << "Region A data, m < 2: " << setprecision(3) << regANorm << " +/- " << regANormErr;
  cout << endl;
  if (option == "separate") legendBkgSep.Draw();
  else if (option == "main 5") legendBkgMain5.Draw();
  else if (option == "combined") legendBkgAll.Draw();
  outCanvas.cd(2);
  TH1F* nonIsoDataMinusIsoBkgAll = (TH1F*)nonIsoData->Clone();
  nonIsoDataMinusIsoBkgAll->Add(isoBkgAllHist, -1.0);
  TH1F* nonIsoDataMinusIsoBkgAllDenom = denomErrorScale((TH1F*)nonIsoData->Clone());
  //  nonIsoDataMinusIsoBkgAll->Divide(nonIsoData);
  nonIsoDataMinusIsoBkgAll->Divide(nonIsoDataMinusIsoBkgAllDenom);
  //  nonIsoDataMinusIsoBkgAll->GetYaxis()->SetTitle("#frac{Data (B) - MC (A)}{Data (B)}");
  nonIsoDataMinusIsoBkgAll->GetYaxis()->SetTitle("#frac{Data (B) - Data (A)}{sqrt(Data(B)^2 + #sigma_{DataB}^2)}");
  nonIsoDataMinusIsoBkgAll->GetYaxis()->SetRangeUser(-1.0, 1.0);
  nonIsoDataMinusIsoBkgAll->Draw();
  Double_t nonIsoDataMinusIsoBkgAllErr = 0.0;
  Double_t nonIsoDataMinusIsoBkgAllVal = 
    nonIsoDataMinusIsoBkgAll->IntegralAndError(5/*17*/, -1, nonIsoDataMinusIsoBkgAllErr);
  cout << "percent deviation in final bin: " << setprecision(4) << nonIsoDataMinusIsoBkgAllVal;
  cout << " +/- " << nonIsoDataMinusIsoBkgAllErr << endl;
  TH1F* nonIsoDataMinusIsoData = (TH1F*)nonIsoData->Clone();
  nonIsoDataMinusIsoData->Add(isoData, -1.0);
  TH1F* nonIsoDataMinusIsoDataDenom = denomErrorScale((TH1F*)nonIsoData->Clone());
  //  nonIsoDataMinusIsoData->Divide(nonIsoData);
  nonIsoDataMinusIsoData->Divide(nonIsoDataMinusIsoDataDenom);
  //  nonIsoDataMinusIsoData->GetYaxis()->SetTitle("#frac{Data (B) - Data (A)}{Data (B)}");
  nonIsoDataMinusIsoData->GetYaxis()->SetTitle("#frac{Data (B) - Data (A)}{sqrt(Data(B)^2 + #sigma_{DataB}^2)}");
  nonIsoDataMinusIsoData->GetYaxis()->SetRangeUser(-1.0, 1.0);
  nonIsoDataMinusIsoData->SetMarkerColor(kBlack);
  nonIsoDataMinusIsoData->SetLineColor(kBlack);
//   nonIsoDataMinusIsoData->Draw("SAME");
  outCanvas.Write();
}

/*make plot showing:
  - jet fake background estimate from region B data
  - jet fake background estimate from MC
  - jet fake QCD estimate from region B/C/D data
  options for displaying MC background:
  - each sample separately
  - WW/WZ/ZZ combined into diboson, tt and single top combined into top
  - all samples combined*/
void addJetFakeBkgFinalPlot(pair<TFile*, float>& isoSigBkgFile, TFile& isoDataFile, 
			    pair<TFile*, float>& nonIsoDataFile, const string& var, 
			    const string& unit, const int normRegionLowerBin, 
			    const int normRegionUpperBin, const string& option, TFile& outStream)
{
  //top level declarations
  string canvasName(var + "Canvas");
  string stackName(var + "Stack");

  //get plots of background MC, isolated tau sample
  TLegend legendBkgSep(0.35, 0.55, 0.75, 0.95);
  TLegend legendBkgMain5(0.35, 0.55, 0.75, 0.95);
  TLegend legendBkgAll(0.35, 0.55, 0.75, 0.95);
  setLegendOptions(legendBkgSep, "CMS 19.7 fb^{-1}");
  setLegendOptions(legendBkgMain5, "CMS 19.7 fb^{-1}");
  setLegendOptions(legendBkgAll, "CMS 19.7 fb^{-1}");
  TCanvas* canvasIsoSigBkg = NULL;
  THStack* isoBkg = NULL;
  isoSigBkgFile.first->GetObject(canvasName.c_str(), canvasIsoSigBkg);
  if (canvasIsoSigBkg != NULL) isoBkg = (THStack*)canvasIsoSigBkg->GetPrimitive(stackName.c_str());
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << isoSigBkgFile.first->GetName() << ".\n";
    return;
  }
  TList* isoBkgHists = NULL;
  if (isoBkg != NULL) isoBkgHists = isoBkg->GetHists();
  else {
    cerr << "Error opening stack " << stackName << " from canvas " << canvasName << " from file ";
    cerr << isoSigBkgFile.first->GetName() << ".\n";
    return;
  }
  TH1F* isoBkgAllHist = NULL;
  TH1F* isoBkgDibosonHist = NULL;
  TH1F* isoBkgWNJetsHist = NULL;
  TH1F* isoBkgTopHist = NULL;
  TH1F* isoBkgDrellYanHist = NULL;
  TH1F* isoBkgQCDHist = NULL;
  THStack isoBkgSep("isoBkgSep", "");
  THStack isoBkgAll("isoBkgAll", "");
  THStack isoBkgMain5("isoBkgMain5", "");
  if (isoBkgHists != NULL) {
    for (Int_t i = 0; i < isoBkgHists->GetEntries(); ++i) {
      TH1F* stackHist = (TH1F*)isoBkgHists->At(i)->Clone();
      stackHist->Scale(isoSigBkgFile.second);
      stackHist->GetYaxis()->SetRangeUser(0.01, 10000.0);
      isoBkgSep.Add(stackHist, "HIST");
      Double_t err = 0.0;
      Double_t val = stackHist->IntegralAndError(5/*17*/, -1, err);
      if (i == 0) {
	isoBkgAllHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgAllHist->Scale(isoSigBkgFile.second);
	isoBkgDibosonHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgDibosonHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "WW (MC)", "f");
	cout << "prediction from WW: " << val << " +/- " << err << endl;
      }
      else isoBkgAllHist->Add(stackHist);
      if (i == 1)
	{
	  legendBkgSep.AddEntry(stackHist, "ZZ (MC)", "f");
	  cout << "prediction from ZZ: " << val << " +/- " << err << endl;
	}
      if (i == 2)
	{
	  legendBkgSep.AddEntry(stackHist, "WZ (MC)", "f");
	  cout << "prediction from WZ: " << val << " +/- " << err << endl;
	}
      if ((i == 1) || (i == 2)) isoBkgDibosonHist->Add(stackHist);
      if (i == 3) {
	isoBkgWNJetsHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgWNJetsHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "W + #geq1 jet (MC)", "f");
	cout << "prediction from WNJets: " << val << " +/- " << err << endl;
      }
      if (i == 4) {
	isoBkgTopHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgTopHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "t/#bar{t} (MC)", "f");
	cout << "prediction from t/tbar: " << val << " +/- " << err << endl;
      }
      if (i == 5) {
	isoBkgTopHist->Add(stackHist);
	legendBkgSep.AddEntry(stackHist, "t#bar{t} + jets (MC)", "f");
	cout << "prediction from TTJets: " << val << " +/- " << err << endl;
      }
      if (i == 6) {
	isoBkgDrellYanHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgDrellYanHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "Drell-Yan + jets (MC)", "f");
	cout << "prediction from DY: " << val << " +/- " << err << endl;
      }
      if (i == 7) {
	isoBkgQCDHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgQCDHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "QCD (data)", "f");
	cout << "prediction from non-resonant QCD (data): " << val << " +/- " << err << endl;
      }
    }
    isoBkgAllHist->GetYaxis()->SetRangeUser(0.01, 10000.0);
    isoBkgAll.Add(isoBkgAllHist, "HISTE C");
    isoBkgMain5.Add(isoBkgDibosonHist, "HIST");
    isoBkgMain5.Add(isoBkgWNJetsHist, "HIST");
    isoBkgMain5.Add(isoBkgTopHist, "HIST");
    isoBkgMain5.Add(isoBkgDrellYanHist, "HIST");
    isoBkgMain5.Add(isoBkgQCDHist, "HIST");
    legendBkgAll.AddEntry(isoBkgAllHist, "MC EW + data QCD", "f");
    legendBkgMain5.AddEntry(isoBkgDibosonHist, "MC diboson", "f");
    legendBkgMain5.AddEntry(isoBkgWNJetsHist, "MC W + #geq1 jet", "f");
    legendBkgMain5.AddEntry(isoBkgTopHist, "MC top + jets", "f");
    legendBkgMain5.AddEntry(isoBkgDrellYanHist, "MC Drell-Yan + jets", "f");
    legendBkgMain5.AddEntry(isoBkgQCDHist, "Data non-resonant QCD", "f");
  }
  else {
    cerr << "Error opening histogram list from stack " << stackName << " from canvas ";
    cerr << canvasName << " from file " << isoSigBkgFile.first->GetName() << ".\n";
    return;
  }

  //get the data histogram, non-isolated tau sample
  TCanvas* canvasNonIsoData = NULL;
  nonIsoDataFile.first->GetObject(canvasName.c_str(), canvasNonIsoData);
  TH1F* nonIsoData = NULL;
  if (canvasNonIsoData != NULL) {
    canvasNonIsoData->Draw();
    nonIsoData = (TH1F*)canvasNonIsoData->cd(1)->GetPrimitive(var.c_str())->Clone();
    nonIsoData->SetName((var + "DataControl").c_str());
    setHistogramOptions(nonIsoData, kRed, 0.7, 20, nonIsoDataFile.second, unit.c_str(), "");
    nonIsoData->GetYaxis()->SetRangeUser(0.01, 10000.0);
    legendBkgSep.AddEntry(nonIsoData, "Jet fake background (from data)", "lp");
    legendBkgMain5.AddEntry(nonIsoData, "Jet fake background (from data)", "lp");
    legendBkgAll.AddEntry(nonIsoData, "Jet fake background (from data)", "lp");
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << nonIsoDataFile.first->GetName() << ".\n";
    return;
  }

  //get the data histogram, isolated tau sample
  TCanvas* canvasIsoData = NULL;
  isoDataFile.GetObject(canvasName.c_str(), canvasIsoData);
  TH1F* isoData = NULL;
  if (canvasIsoData != NULL) {
    isoData = (TH1F*)canvasIsoData->GetPrimitive(var.c_str())->Clone();
    isoData->SetName((var + "DataSearch").c_str());
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file " << isoDataFile.GetName();
    cerr << ".\n";
    return;
  }

  //calculate the normalization factor
  Double_t norm = isoData->Integral(normRegionLowerBin, normRegionUpperBin)/
    nonIsoData->Integral(normRegionLowerBin, normRegionUpperBin);

  /*calculate the statistical error on the background prediction from the non-isolated data, 
    including the term from the error on the normalization factor*/
  pair<Double_t, Double_t> regANormIntAndErr = 
    getIntegralAndError(string(isoDataFile.GetName()), 
			pair<Int_t, Int_t>(normRegionLowerBin, normRegionUpperBin), "muHadMass", 
			0);
  Double_t normRegA = regANormIntAndErr.first;
  Double_t normRegAErr = regANormIntAndErr.second;
  pair<Double_t, Double_t> regBNormIntAndErr = 
    getIntegralAndError(string(nonIsoDataFile.first->GetName()), 
			pair<Int_t, Int_t>(normRegionLowerBin, normRegionUpperBin), "muHadMass", 
			1);
  Double_t normRegB = regBNormIntAndErr.first;
  Double_t normRegBErr = regBNormIntAndErr.second;
  vector<Double_t> nonIsoDataStatErrSq;
  for (Int_t iBin = 0; iBin <= (nonIsoData->GetNbinsX() + 1); ++iBin) {
    Double_t regBBinContentRaw = nonIsoData->GetBinContent(iBin);
    Double_t regBBinErrRaw = nonIsoData->GetBinError(iBin);
    Double_t regBBinErrSqRaw = regBBinContentRaw == 0.0 ? 
      0.0 : (regBBinErrRaw*regBBinErrRaw)/(regBBinContentRaw*regBBinContentRaw);
    Double_t normRegBErrSq = normRegB == 0.0 ? 0.0 : (normRegBErr*normRegBErr)/(normRegB*normRegB);
    Double_t normRegAErrSq = normRegA == 0.0 ? 0.0 : (normRegAErr*normRegAErr)/(normRegA*normRegA);
    Double_t regBBinContentSqFinal = 
      nonIsoData->GetBinContent(iBin)*nonIsoData->GetBinContent(iBin)*norm*norm;
    Double_t regBBinErrSqFinal = 
      regBBinContentSqFinal*(regBBinErrSqRaw + normRegBErrSq + normRegAErrSq);
    nonIsoDataStatErrSq.push_back(regBBinErrSqFinal);
  }

  //normalize non-isolated data histogram to isolated data in signal-depleted region
  nonIsoData->Scale(norm);

  //set statistical error in each bin of the non-isolated data histogram
  for (Int_t iBin = 0; iBin <= (nonIsoData->GetNbinsX() + 1); ++iBin) {
    nonIsoData->SetBinError(iBin, sqrt(nonIsoDataStatErrSq[iBin]));
  }

  //write to file
  outStream.cd();
  TCanvas outCanvas(canvasName.c_str(), "", 600, 900);
  outCanvas.Divide(1, 2);
  outCanvas.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  setCanvasOptions(*outCanvas.cd(1), 1, 1, 0);
  outCanvas.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  setCanvasOptions(*outCanvas.cd(2), 1, 0, 0);
  outCanvas.cd(1);
  if (option == "separate") {
    isoBkgSep.Draw();
    isoBkgSep.SetMinimum(0.01);
    isoBkgSep.SetMaximum(10000.0);
    isoBkgSep.GetHistogram()->GetXaxis()->SetTitle(unit.c_str());
  }
  else if (option == "main 5") {
    isoBkgMain5.Draw();
    Float_t sum = 0.0;
    Float_t err = 0.0;
    for (Int_t iHist = 0; iHist < isoBkgMain5.GetHists()->GetEntries(); ++iHist) {
      TH1F* hist = (TH1F*)isoBkgMain5.GetHists()->At(iHist);
      for (Int_t iBin = 5/*17*/; iBin <= (hist->GetNbinsX() + 1); ++iBin) {
	sum+=hist->GetBinContent(iBin);
	err+=(hist->GetBinError(iBin)*hist->GetBinError(iBin));
      }
    }
    cout << "Region A MC + data-driven QCD, m > 4: " << setprecision(3) << sum << " +/- ";
    cout << sqrt(err) <<  endl;
    isoBkgMain5.SetMinimum(0.01);
    isoBkgMain5.SetMaximum(10000.0);
    isoBkgMain5.GetHistogram()->GetXaxis()->SetTitle(unit.c_str());
  }
  else if (option == "combined") {
    isoBkgAll.Draw();
    isoBkgAll.SetMinimum(0.01);
    isoBkgAll.SetMaximum(10000.0);
    isoBkgAll.GetHistogram()->GetXaxis()->SetTitle(unit.c_str());
  }
  nonIsoData->Draw("HISTESAME");
  Double_t statErrNonIsoData;
  nonIsoData->IntegralAndError(5/*17*/, -1, statErrNonIsoData);
  cout << "Region B data, m > 4: " << setprecision(3) << nonIsoData->Integral(5/*17*/, -1);
  cout << " +/- " << statErrNonIsoData << endl;
  Double_t regANormErr;
  Double_t regANorm = 
    isoData->IntegralAndError(normRegionLowerBin, normRegionUpperBin, regANormErr);
  cout << "Region A data, m < 2: " << setprecision(3) << regANorm << " +/- " << regANormErr;
  cout << endl;
  if (option == "separate") legendBkgSep.Draw();
  else if (option == "main 5") legendBkgMain5.Draw();
  else if (option == "combined") legendBkgAll.Draw();
  outCanvas.cd(2);
  TH1F* nonIsoDataOverIsoBkgAll = (TH1F*)nonIsoData->Clone();
  nonIsoDataOverIsoBkgAll->Divide(isoBkgAllHist);
  nonIsoDataOverIsoBkgAll->GetYaxis()->SetTitle("#frac{Data (B)}{A}");
  nonIsoDataOverIsoBkgAll->GetYaxis()->SetRangeUser(0.0, 2.0);
  nonIsoDataOverIsoBkgAll->Draw();
  outCanvas.Write();
}

void arcQuest(const vector<string>& QCDVsMCInputFileNames, const string& isoDataFileName, 
		  const string& nonIsoDataFileName, const string& var, const string& unit, 
		  const int normRegionLowerBin, const int normRegionUpperBin, 
		  const string& outputFileName)
{ // beginning
  //top level declarations
  TFile outStream(outputFileName.c_str(), "RECREATE");
  TFile isoDataFile(isoDataFileName.c_str(), "read");
  TFile nonIsoDataFile(nonIsoDataFileName.c_str(), "read");
  string canvasName(var + "Canvas");
  string canvasReweightErrSqName(var + "ReweightErrSqCanvas");
  string varReweightErrSq(var + "ReweightErrSq");
  vector<TFile*> inputStreams;
  vector<TH1F*> hists(2);

  TLegend legendBkgQCDEWK(0.35, 0.55, 0.75, 0.95);
  setLegendOptions(legendBkgQCDEWK, "CMS 19.7 fb^{-1}");

  cout << "arcQuest: looping over QCDMC input files" << endl;
  for (vector<string>::const_iterator iInputFile = QCDVsMCInputFileNames.begin(); 
       iInputFile != QCDVsMCInputFileNames.end(); ++iInputFile) { // loop over input files
    const unsigned int fileIndex = iInputFile - QCDVsMCInputFileNames.begin();
    inputStreams.push_back(new TFile(iInputFile->c_str()));
    cout << "arcQuest: filename = " << iInputFile->c_str() << endl;
    TCanvas* pCanvas;
    inputStreams[inputStreams.size() - 1]->GetObject(canvasName.c_str(), pCanvas);
    TH1F* pHist = (TH1F*)pCanvas->GetPrimitive(var.c_str());
    if (fileIndex == 0)
      { // if this is the QCD file
	hists[0] = (TH1F*)pHist->Clone();
	hists[0]->SetLineColor(3); // green QCD
	hists[0]->SetMarkerColor(3); // green QCD
      } // if this is the QCD file
    else if (fileIndex == 1)
      { // if this is the first MC bkg file
	hists[1] = (TH1F*)pHist->Clone();
	hists[1]->SetLineColor(4); // blue EWK
	hists[1]->SetMarkerColor(4); // blue EWK
      } // if this is the first MC bkg file
    else
      { // if this is another MC bkg file
	hists[1]->Add(pHist,1.);
      } // if this is another MC bkg file
    
  } // loop over input files

  if (hists[0] != NULL)
    legendBkgQCDEWK.AddEntry(hists[0], "QCD (data-driven, from Region D)", "lp");
  if (hists[1] != NULL)
    legendBkgQCDEWK.AddEntry(hists[1], "EWK (MC, from Region B)", "lp");

  //get the data histogram, non-isolated tau sample
  TCanvas* canvasNonIsoData = NULL;
  nonIsoDataFile.GetObject(canvasName.c_str(), canvasNonIsoData);
  TH1F* nonIsoData = NULL;
  if (canvasNonIsoData != NULL) {
    canvasNonIsoData->Draw();
    nonIsoData = (TH1F*)canvasNonIsoData->cd(1)->GetPrimitive(var.c_str())->Clone();
    nonIsoData->SetName((var + "DataControl").c_str());
    setHistogramOptions(nonIsoData, kRed, 0.7, 20, 1.0, unit.c_str(), "");
    nonIsoData->GetYaxis()->SetRangeUser(0.01, 10000.0);
    legendBkgQCDEWK.AddEntry(nonIsoData, "Background (from data)", "lp");
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << nonIsoDataFile.GetName() << ".\n";
    return;
  }

  //get the reweighting error histogram, non-isolated tau sample
  TCanvas* canvasNonIsoDataReweightErrSq = NULL;
  nonIsoDataFile.GetObject(canvasReweightErrSqName.c_str(), canvasNonIsoDataReweightErrSq);
  TH1F* nonIsoDataReweightErrSq = NULL;
  if (canvasNonIsoDataReweightErrSq != NULL) {
    canvasNonIsoDataReweightErrSq->Draw();
    nonIsoDataReweightErrSq = 
      (TH1F*)canvasNonIsoDataReweightErrSq->cd(1)->GetPrimitive(varReweightErrSq.c_str())->Clone();
  }
  else {
    cerr << "Error opening canvas " << canvasReweightErrSqName << " from file ";
    cerr << nonIsoDataFile.GetName() << ".\n";
    return;
  }

  //get the data histogram, isolated tau sample
  TCanvas* canvasIsoData = NULL;
  isoDataFile.GetObject(canvasName.c_str(), canvasIsoData);
  TH1F* isoData = NULL;
  if (canvasIsoData != NULL) {
    isoData = (TH1F*)canvasIsoData->GetPrimitive(var.c_str())->Clone();
    isoData->SetName((var + "DataSearch").c_str());
    //    legendBkgQCDEWK.AddEntry(isoData, "Data", "p");
//     setHistogramOptions(isoData, kBlack, 0.7, 20, 1.0, unit.c_str(), "");
    isoData->GetYaxis()->SetRangeUser(0.01, 10000.0);
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file " << isoDataFile.GetName();
    cerr << ".\n";
    return;
  }

  //calculate the normalization factor
  Double_t norm = isoData->Integral(normRegionLowerBin, normRegionUpperBin)/
    nonIsoData->Integral(normRegionLowerBin, normRegionUpperBin);
  cout << "The normalization constant is: " << norm << endl;
  /*calculate the statistical error on the background prediction from the non-isolated data, 
    including the term from the error on the normalization factor*/
  const TH1* nonIsoDataPtrCast = dynamic_cast<const TH1*>(nonIsoData);
  const TH1* nonIsoDataReweightErrSqPtrCast = dynamic_cast<const TH1*>(nonIsoDataReweightErrSq);
  vector<Double_t> nonIsoDataStatErrSq;
  for (Int_t iBin = 1; iBin <= nonIsoData->GetNbinsX(); ++iBin) {
    nonIsoDataStatErrSq.
      push_back(bkgErrSq(nonIsoDataPtrCast, nonIsoDataReweightErrSqPtrCast, iBin, norm, 
			 normErrSq(dynamic_cast<const TH1*>(isoData), nonIsoDataPtrCast, 
				   nonIsoDataReweightErrSqPtrCast, normRegionLowerBin, 
				   normRegionUpperBin, norm)));
  }

  //normalize non-isolated data histogram to isolated data in signal-depleted region
  cout << "Region B m > 4 integral before normalization = " << nonIsoData->Integral(5,-1) << endl;
  nonIsoData->Scale(norm);
  cout << "Normalization constant for region B = " << norm << endl;

  //set statistical error in each bin of the non-isolated data histogram
  for (Int_t iBin = 1; iBin <= nonIsoData->GetNbinsX(); ++iBin) {
    nonIsoData->SetBinError(iBin, sqrt(nonIsoDataStatErrSq[iBin - 1]));
  }

  //normalize region B QCD histogram to isolated data in signal-depleted region
  Double_t normQCD = isoData->Integral(normRegionLowerBin, normRegionUpperBin)/
    hists[0]->Integral(normRegionLowerBin, normRegionUpperBin);
  hists[0]->Scale(normQCD);

  //normalize region B total EWK histogram to isolated data in signal-depleted region
  //normalize region B QCD histogram to isolated data in signal-depleted region
  Double_t normEWK = isoData->Integral(normRegionLowerBin, normRegionUpperBin)/
    hists[1]->Integral(normRegionLowerBin, normRegionUpperBin);
  hists[1]->Scale(normEWK);

  //errors on m > 4 integral for QCD and EWK histograms
  Float_t sum0 = 0.0;
  Float_t err0 = 0.0;
  Float_t sum1 = 0.0;
  Float_t err1 = 0.0;
  for (Int_t iBin = 5; iBin <= (hists[0]->GetNbinsX() + 1); ++iBin) {
    sum0+=hists[0]->GetBinContent(iBin);
    err0+=(hists[0]->GetBinError(iBin)*hists[0]->GetBinError(iBin));
  }
  for (Int_t iBin = 5; iBin <= (hists[1]->GetNbinsX() + 1); ++iBin) {
    sum1+=hists[1]->GetBinContent(iBin);
    err1+=(hists[1]->GetBinError(iBin)*hists[1]->GetBinError(iBin));
  }

  //write to file
  outStream.cd();
  TCanvas outCanvasAllQCDEWK("Region B data vs all QCD/EWK bkg", "", 600, 900);
  outCanvasAllQCDEWK.Divide(1, 2);
  outCanvasAllQCDEWK.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  setCanvasOptions(*outCanvasAllQCDEWK.cd(1), 1, 1, 0);
  outCanvasAllQCDEWK.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  setCanvasOptions(*outCanvasAllQCDEWK.cd(2), 1, 0, 0);

  outCanvasAllQCDEWK.cd(1);
  nonIsoData->Draw("HISTE");
  hists[0]->Draw("HISTESAME");
  hists[1]->Draw("HISTESAME");
  legendBkgQCDEWK.Draw();
  Double_t statErrNonIsoData;
  nonIsoData->IntegralAndError(5, -1, statErrNonIsoData);
  cout << "Region B data, m > 4: " << setprecision(3) << nonIsoData->Integral(5, -1) << " +/- ";
  cout << setprecision(3) << statErrNonIsoData << endl;
  cout << "Region B QCD normalized, m > 4: " << setprecision(3) << sum0 << " +/- " << setprecision(3) << sqrt(err0) << endl;
  cout << "Region B EWK normalized, m > 4: " << setprecision(3) << sum1 << " +/- " << setprecision(3) << sqrt(err1) << endl;
  outCanvasAllQCDEWK.cd(2);
  TH1F* nonIsoDataMinusRegBQCD = (TH1F*)nonIsoData->Clone();
  nonIsoDataMinusRegBQCD->Add(hists[0], -1.0);
//   TH1F* nonIsoDataMinusRegBQCDDenom = denomErrorScale((TH1F*)nonIsoData->Clone());
  nonIsoDataMinusRegBQCD->Divide(/*nonIsoDataMinusRegBQCDDenom*/nonIsoData);
  nonIsoDataMinusRegBQCD->GetYaxis()->
    SetTitle(/*"#frac{Data(B) - Exp(B)}{sqrt(Data(B)^2 + #sigma_{DataB}^2)}"*/"Syst. error envelope");
  nonIsoDataMinusRegBQCD->GetYaxis()->SetRangeUser(-1.0, 1.0);
  nonIsoDataMinusRegBQCD->SetLineColor(3);
  nonIsoDataMinusRegBQCD->SetMarkerColor(3);
  nonIsoDataMinusRegBQCD->Draw("pe");
  Double_t regBQCDerr = 0.0;
  Double_t regBQCDintegral = nonIsoDataMinusRegBQCD->IntegralAndError(5, -1, regBQCDerr);
  cout << "percent deviation in final bin (QCD): " << regBQCDintegral << " +/- " << regBQCDerr << endl;
  TH1F* nonIsoDataMinusRegBEWK = (TH1F*)nonIsoData->Clone();
  nonIsoDataMinusRegBEWK->Add(hists[1], -1.0);
//   TH1F* nonIsoDataMinusRegBEWKDenom = denomErrorScale((TH1F*)nonIsoData->Clone());
  nonIsoDataMinusRegBEWK->Divide(/*nonIsoDataMinusRegBEWKDenom*/nonIsoData);
  nonIsoDataMinusRegBEWK->GetYaxis()->
    SetTitle(/*"#frac{Data(B) - Exp(B)}{sqrt(Data(B)^2 + #sigma_{DataB}^2)}"*/"Syst. error envelope");
  nonIsoDataMinusRegBEWK->GetYaxis()->SetRangeUser(-1.0, 1.0);
  nonIsoDataMinusRegBEWK->SetLineColor(4);
  nonIsoDataMinusRegBEWK->SetMarkerColor(4);
  nonIsoDataMinusRegBEWK->Draw("samepe");
  Double_t regBEWKerr = 0.0;
  Double_t regBEWKintegral = nonIsoDataMinusRegBEWK->IntegralAndError(5, -1, regBEWKerr);
  cout << "percent deviation in final bin (EWK): " << regBEWKintegral << " +/- " << regBEWKerr << endl;

  outCanvasAllQCDEWK.Write();

} // end

void GetVBFZHEstimate(const string& vbfOutputFileName, const string& isoggHFileName, double ggHWeight, double VBFWeight)
{ // start VBFEstimateFromGGH
  // N.B. this routine can also get ZH from WH
  TFile outStream(vbfOutputFileName.c_str(), "RECREATE");
  TFile isoggHFile(isoggHFileName.c_str(), "read"); // make sure this is the ggH iso file after hadding!

  string var = "muHadMass";
  string canvasName(var + "Canvas");
  string canvasReweightErrSqName(var + "ReweightErrSqCanvas");
  string varReweightErrSq(var + "ReweightErrSq");

  //get the ggH m_(mu+had) histogram, isolated tau sample
  TCanvas* canvasIsoGGH = NULL;
  isoggHFile.GetObject(canvasName.c_str(), canvasIsoGGH);
  TH1F* isoGGH = NULL;
  if (canvasIsoGGH != NULL) {
    isoGGH = (TH1F*)canvasIsoGGH->GetPrimitive(var.c_str())->Clone();
    isoGGH->SetName(var.c_str());
    isoGGH->GetYaxis()->SetRangeUser(0.01, 10000.0);
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file " << isoggHFile.GetName();
    cerr << ".\n";
    return;
  }

  //calculate normalization
  double norm = VBFWeight/ggHWeight;

  isoGGH->Scale(norm);
  cout << "Sample used for estimate: " << isoggHFileName.c_str() << endl;
  cout << "Estimate: m > 4 = " << isoGGH->Integral(5,-1) << " +/- " << isoGGH->GetBinError(5) << endl;

  TCanvas outCanvas(canvasName.c_str(), "", 600, 600);
  outCanvas.cd();
  isoGGH->Draw();

  outStream.cd();
  outCanvas.Write();
  outStream.Write();
  outStream.Close();

} // end VBFEstimateFromGGH

//normalize plot and propagate statistical error from normalization term
//don't use getIntegralAndError because they don't work for MC stacks
void normalizeHistogram(const TH1F* histToNormalizeTo, TH1F* histToNormalize, 
			const int normRegionLowerBin, const int normRegionUpperBin)
{
  /*calculate the statistical error on the data to be normalized, including the term from the 
    error on the normalization factor*/
  Double_t normToErr, toBeNormErr;
  const Double_t normTo = 
    histToNormalizeTo->IntegralAndError(normRegionLowerBin, normRegionUpperBin, normToErr);
  const Double_t toBeNorm = 
    histToNormalize->IntegralAndError(normRegionLowerBin, normRegionUpperBin, toBeNormErr);
  const Double_t norm = normTo/toBeNorm;
  vector<Double_t> histToNormalizeStatErrSq;
  for (Int_t iBin = 0; iBin <= (histToNormalize->GetNbinsX() + 1); ++iBin) {
    const Double_t histToNormalizeBinContentRaw = histToNormalize->GetBinContent(iBin);
    const Double_t histToNormalizeBinErrRaw = histToNormalize->GetBinError(iBin);
    const Double_t histToNormalizeBinErrSqRaw = histToNormalizeBinContentRaw == 0.0 ? 
      0.0 : (histToNormalizeBinErrRaw*histToNormalizeBinErrRaw)/
      (histToNormalizeBinContentRaw*histToNormalizeBinContentRaw);
    const Double_t toBeNormErrSq = 
      toBeNorm == 0.0 ? 0.0 : (toBeNormErr*toBeNormErr)/(toBeNorm*toBeNorm);
    const Double_t normToErrSq = normTo == 0.0 ? 0.0 : (normToErr*normToErr)/(normTo*normTo);
    const Double_t histToNormalizeBinContentSqFinal = 
      histToNormalize->GetBinContent(iBin)*histToNormalize->GetBinContent(iBin)*norm*norm;
    const Double_t histToNormalizeBinErrSqFinal = 
      histToNormalizeBinContentSqFinal*(histToNormalizeBinErrSqRaw + toBeNormErrSq + normToErrSq);
    histToNormalizeStatErrSq.push_back(histToNormalizeBinErrSqFinal);

//     //debug
//     const Double_t sigmaBNorm = 
//       histToNormalize->GetBinContent(iBin)*norm*sqrt(toBeNormErrSq + normToErrSq);
//     const Double_t sigmaBStat = 
//       histToNormalize->GetBinContent(iBin)*norm*sqrt(histToNormalizeBinErrSqRaw);
//     const Double_t sigmaB = sqrt(histToNormalizeBinErrSqFinal);
//     cout << "Bin " << iBin << endl;
//     cout << "sigmaBNorm = " << setprecision(3) << sigmaBNorm << endl;
//     cout << "sigmaBStat = " << setprecision(3) << sigmaBStat << endl;
//     cout << "sigmaB = " << setprecision(3) << sigmaB << endl;
  }

  //normalize histogram
  histToNormalize->Scale(norm);

  //set statistical error
  for (Int_t iBin = 0; iBin <= (histToNormalize->GetNbinsX() + 1); ++iBin) {
    histToNormalize->SetBinError(iBin, sqrt(histToNormalizeStatErrSq[iBin]));

    //debug
//     cout << "Bin error: " << histToNormalize->GetBinError(iBin) << endl;
//     cout << "Bin " << iBin << " content: " << histToNormalize->GetBinContent(iBin) << endl;
  }
}

//make and format pull histogram for bkg. stat. error only
TH1F* makeAndFormatPullHistogram(TH1F* data, TH1F* resBkgTemplate, TH1F* jetFakeBkgTemplate, 
				 const Color_t fillColor, Option_t* drawOpt)
{
  TH1F* totBkgHist = (TH1F*)jetFakeBkgTemplate->Clone();
  Double_t jetFakeBkgSigRegErr;
  Double_t jetFakeBkgSigReg = totBkgHist->IntegralAndError(5/*17*/, -1, jetFakeBkgSigRegErr);
  cout << "Jet fake bkg.: " << setprecision(3) << jetFakeBkgSigReg << " +/- ";
  cout << jetFakeBkgSigRegErr << endl;
  totBkgHist->Add((TH1F*)resBkgTemplate->Clone());
  TH1F* pull = makeDataBkgAgreementHist(data, totBkgHist);

//   //debug
//   for (Int_t iBin = 9; iBin <= pull->GetNbinsX(); ++iBin) { pull->SetBinContent(iBin, 0.0); }

  pull->SetLineColor(kBlack);
  pull->SetFillColor(fillColor);
  pull->SetFillStyle(1001);
  pull->Draw(drawOpt);
  return pull;
}

//make and format pull histogram for bkg. stat. + syst. error
TH1F* makeAndFormatPullHistogram(TH1F* data, TH1F* resBkgTemplate, TH1F* jetFakeBkgTemplate, 
				 const TH1F* allQCDBkgTemplate, const TH1F* allEWBkgTemplate, 
				 const Color_t fillColor, Option_t* drawOpt)
{
  TH1F* totBkgHist = (TH1F*)jetFakeBkgTemplate->Clone();
  Double_t jetFakeBkgSigRegErr;
  Double_t jetFakeBkgSigReg = totBkgHist->IntegralAndError(5/*17*/, -1, jetFakeBkgSigRegErr);
  cout << "Jet fake bkg.: " << setprecision(3) << jetFakeBkgSigReg << " +/- ";
  cout << jetFakeBkgSigRegErr << endl;
  totBkgHist->Add((TH1F*)resBkgTemplate->Clone());
  TH1F* pull = makeDataBkgAgreementHist(data, totBkgHist, allQCDBkgTemplate, allEWBkgTemplate);

//   //debug
//   for (Int_t iBin = 9; iBin <= pull->GetNbinsX(); ++iBin) { pull->SetBinContent(iBin, 0.0); }

  pull->SetLineColor(kBlack);
  pull->SetFillColor(fillColor);
  pull->SetFillStyle(1001);
  pull->Draw(drawOpt);
  return pull;
}

/*make final plot showing:
  - expected signal
  - full background estimate from regions B and C data
  - background systematic error*/
void addFinalPlot2(pair<TFile*, float>& isoSigBkgFile, TFile& isoDataFile, 
		   pair<TFile*, float>& nonIsoDataFile, pair<TFile*, float>& resBkgFile, 
		   pair<TFile*, float>& nonIsoWNonIsoDataFile, 
		   const int normRegionLowerBin, const int normRegionUpperBin, TFile& outStream, 
		   const bool ma9GeV)
{
  //top level declarations
  string canvasName("muHadMassCanvas");

  //get plots of signal MC, isolated tau sample
  TCanvas* canvasIsoSigBkg = NULL;
  isoSigBkgFile.first->GetObject(canvasName.c_str(), canvasIsoSigBkg);

  //get the signal histograms
  //for ma = 9 GeV, there are 4 signals
  //for all other pseudoscalar masses, there are 2 signals
  vector<TH1F*> isoSig(4, NULL);
  if (!ma9GeV) isoSig.erase(isoSig.begin() + 2, isoSig.end());
  vector<string> sigSamples;
  sigSamples.push_back("WH");
  sigSamples.push_back("ggH");
  sigSamples.push_back("ZH");
  sigSamples.push_back("VBF");
  vector<Color_t> sigColors;
  sigColors.push_back(kSpring - 1);
  sigColors.push_back(kBlue);
  sigColors.push_back(kSpring - 7);
  sigColors.push_back(kMagenta);
  TLegend legendBkgMain5(0.6, 0.55, 0.9, 0.95);
  setLegendOptions(legendBkgMain5, "CMS 19.7 fb^{-1}");
  if (canvasIsoSigBkg != NULL) {
    TList* sigs = canvasIsoSigBkg->GetListOfPrimitives();
    for (vector<TH1F*>::iterator iIsoSig = isoSig.begin(); iIsoSig != isoSig.end(); 
	 ++iIsoSig) {
      const unsigned int i = iIsoSig - isoSig.begin();
      *iIsoSig = (TH1F*)sigs->At(i + 2)->Clone();
      legendBkgMain5.AddEntry(*iIsoSig, sigSamples[i].c_str(), "l");
      (*iIsoSig)->SetName(("muHadMass" + sigSamples[i] + "Search").c_str());
      setHistogramOptions(*iIsoSig, sigColors[i], 0.7, 20, isoSigBkgFile.second, 
			  "m_{#mu+had} (GeV)", "");
      (*iIsoSig)->SetLineWidth(4);
      (*iIsoSig)->SetLineStyle(2);
      (*iIsoSig)->GetYaxis()->SetRangeUser(0.01, 100000.0);
    }
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << isoSigBkgFile.first->GetName() << ".\n";
    return;
  }

  //get the data histogram, non-isolated tau sample
  TCanvas* canvasNonIsoData = NULL;
  nonIsoDataFile.first->GetObject(canvasName.c_str(), canvasNonIsoData);
  TH1F* nonIsoData = NULL;
  if (canvasNonIsoData != NULL) {
    canvasNonIsoData->Draw();
    nonIsoData = (TH1F*)canvasNonIsoData->cd(1)->GetPrimitive("muHadMass")->Clone();
    nonIsoData->SetName("muHadMassDataControl");
    setHistogramOptions(nonIsoData, kRed, 0.7, 20, nonIsoDataFile.second, 
			"m_{#mu+had} (GeV)", "");
    nonIsoData->SetFillColor(kRed);
    nonIsoData->SetFillStyle(3002);
    nonIsoData->GetYaxis()->SetRangeUser(0.01, 100000.0);
    legendBkgMain5.AddEntry(nonIsoData, "Jet fake background", "f");
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << nonIsoDataFile.first->GetName() << ".\n";
    return;
  }

  //get the MC stack, non-isolated tau sample
  THStack* nonIsoMCStack = NULL;
  if (canvasNonIsoData != NULL) {
    canvasNonIsoData->Draw();
    nonIsoMCStack = (THStack*)canvasNonIsoData->cd(1)->GetPrimitive("muHadMassStack")->Clone();
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << nonIsoDataFile.first->GetName() << ".\n";
    return;
  }

  //sum the region B MC into a single histogram
  TH1F* nonIsoMCHist = NULL;
  if (nonIsoMCStack != NULL) {
    TList* nonIsoMCHists = nonIsoMCStack->GetHists();
    for (Int_t iHist = 0; iHist < nonIsoMCHists->GetEntries(); ++iHist) {
      TH1F* hist = (TH1F*)nonIsoMCHists->At(iHist)->Clone();
      if (iHist == 0) nonIsoMCHist = hist;
      else nonIsoMCHist->Add(hist);
    }
  }
  else {
    cerr << "Error opening stack muHadMassStack from canvas " << canvasName;
    cerr << " from file " << nonIsoDataFile.first->GetName() << ".\n";
  }
  if (nonIsoMCHist != NULL) nonIsoMCHist->SetName("muHadMassRegBMCSyst");

  //get the data histogram, non-isolated W + non-isolated tau sample
  TCanvas* canvasNonIsoWNonIsoData = NULL;
  nonIsoWNonIsoDataFile.first->GetObject(canvasName.c_str(), canvasNonIsoWNonIsoData);
  TH1F* nonIsoWNonIsoData = NULL;
  if (canvasNonIsoWNonIsoData != NULL) {
    canvasNonIsoWNonIsoData->Draw();
    nonIsoWNonIsoData = (TH1F*)canvasNonIsoWNonIsoData->cd(1)->GetPrimitive("muHadMass")->Clone();
    nonIsoWNonIsoData->SetName("muHadMassRegDDataSyst");
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << nonIsoWNonIsoDataFile.first->GetName() << ".\n";
    return;
  }

  //get the data histogram, isolated tau sample
  TCanvas* canvasIsoData = NULL;
  isoDataFile.GetObject(canvasName.c_str(), canvasIsoData);
  TH1F* isoData = NULL;
  if (canvasIsoData != NULL) {
    isoData = (TH1F*)canvasIsoData->GetPrimitive("muHadMass")->Clone();
    isoData->SetName("muHadMassDataSearch");
    isoData->GetYaxis()->SetRangeUser(0.01, 100000.0);
    legendBkgMain5.AddEntry(isoData, "Data", "p");
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file " << isoDataFile.GetName();
    cerr << ".\n";
    return;
  }

  //normalize region B MC to region A data
  normalizeHistogram(const_cast<const TH1F*>(isoData), nonIsoMCHist, 
		     normRegionLowerBin, normRegionUpperBin);

  //normalize region D data to region A data
  normalizeHistogram(const_cast<const TH1F*>(isoData), nonIsoWNonIsoData, 
		     normRegionLowerBin, normRegionUpperBin);

  //get the resonance background histograms: nominal, 1st fit type, and 2nd fit type
  TH1F* resBkg = NULL;
  TH1F* resBkgRegCkFixed = NULL;
  TH1F* resBkgRegDkFixed = NULL;
  TIter iKey(resBkgFile.first->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)iKey()) && 
	 ((resBkg == NULL) || (resBkgRegCkFixed == NULL) || (resBkgRegDkFixed == NULL))) {
    string objName(key->GetName());
    if (objName.find("weightedAvgkFixed_resBkg") != string::npos) {
      TH1F* obj = NULL;
      resBkgFile.first->GetObject(objName.c_str(), obj);
      if (obj != NULL) {
	resBkg = (TH1F*)obj->Clone();
	resBkg->SetName("muHadMassResBkg");
	setHistogramOptions(resBkg, kCyan + 2, 0.7, 20, resBkgFile.second, 
			    "m_{#mu+had} (GeV)", "");
	resBkg->SetFillColor(kCyan + 2);
	resBkg->SetFillStyle(1001);
	resBkg->GetYaxis()->SetRangeUser(0.01, 100000.0);
	legendBkgMain5.AddEntry(resBkg, "Resonance background", "f");
      }
    }
    else if (objName.find("regCkFixed_resBkg") != string::npos) {
      TH1F* obj = NULL;
      resBkgFile.first->GetObject(objName.c_str(), obj);
      if (obj != NULL) {
	resBkgRegCkFixed = (TH1F*)obj->Clone();
	resBkgRegCkFixed->SetName("muHadMassResBkgRegCkFixedSyst");
      }
    }
    else if (objName.find("regDkFixed_resBkg") != string::npos) {
      TH1F* obj = NULL;
      resBkgFile.first->GetObject(objName.c_str(), obj);
      if (obj != NULL) {
	resBkgRegDkFixed = (TH1F*)obj->Clone();
	resBkgRegDkFixed->SetName("muHadMassResBkgRegDkFixedSyst");
      }
    }
  }

  //normalize region B data to region A data
  normalizeHistogram(const_cast<const TH1F*>(isoData), nonIsoData, 
		     normRegionLowerBin, normRegionUpperBin);

  //stacked histogram for background components
  THStack* totBkg = new THStack("muHadMassStack", "");
  totBkg->Add(resBkg, "HIST");
  totBkg->Add(nonIsoData, "HIST");

  //write to file
  outStream.cd();
  TCanvas outCanvas(canvasName.c_str(), "", 600, 900);
  outCanvas.Divide(1, 2);
  outCanvas.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  setCanvasOptions(*outCanvas.cd(1), 1, 1, 0);
  outCanvas.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  setCanvasOptions(*outCanvas.cd(2), 1, 0, 0);
  outCanvas.cd(1);
  totBkg->Draw();
  totBkg->SetMinimum(0.01);
  totBkg->SetMaximum(100000.0);
  Double_t B = 0.0;
  Double_t err = 0.0;
  Double_t upsilonBkg = 0.0;
  Double_t upsilonBkgErr = 0.0;
  Double_t JPsiBkg = 0.0;
  Double_t JPsiBkgErr = 0.0;
  Double_t jetFakeBkg = 0.0;
  Double_t jetFakeBkgErr = 0.0;
  Double_t jetFakeBkgAroundJPsi = 0.0;
  Double_t jetFakeBkgAroundJPsiErr = 0.0;
  for (Int_t iHist = 0; iHist < totBkg->GetHists()->GetEntries(); ++iHist) {
    TH1F* hist = (TH1F*)totBkg->GetHists()->At(iHist);
    for (Int_t iBin = 5/*17*/; iBin <= (hist->GetNbinsX() + 1); ++iBin) {
      B+=hist->GetBinContent(iBin);
      if (hist->GetBinContent(iBin) != 0.0) err+=(hist->GetBinError(iBin)*hist->GetBinError(iBin));
      if (string(hist->GetName()) == "muHadMassResBkg") {
	upsilonBkg+=hist->GetBinContent(iBin);
	if (hist->GetBinContent(iBin) != 0.0) {
	  upsilonBkgErr+=(hist->GetBinError(iBin)*hist->GetBinError(iBin));
	}
      }
      if (string(hist->GetName()) == "muHadMassDataControl") {
	jetFakeBkg+=hist->GetBinContent(iBin);
	if (hist->GetBinContent(iBin) != 0.0) {
	  jetFakeBkgErr+=(hist->GetBinError(iBin)*hist->GetBinError(iBin));
	}
      }
    }
    for (Int_t iBin = 1; iBin <= 4/*16*/; ++iBin) {
      if (string(hist->GetName()) == "muHadMassResBkg") {
	JPsiBkg+=hist->GetBinContent(iBin);
	if (hist->GetBinContent(iBin) != 0.0) {
	  JPsiBkgErr+=(hist->GetBinError(iBin)*hist->GetBinError(iBin));
	}
      }
      if ((string(hist->GetName()) == "muHadMassDataControl") && (iBin > 2/*8*/)) {
	jetFakeBkgAroundJPsi+=hist->GetBinContent(iBin);
	if (hist->GetBinContent(iBin) != 0.0) {
	  jetFakeBkgAroundJPsiErr+=(hist->GetBinError(iBin)*hist->GetBinError(iBin));
	}
      }
    }
  }
  Double_t errSqrt = sqrt(err);
  Double_t upsilonBkgErrSqrt = sqrt(upsilonBkgErr);
  Double_t JPsiBkgErrSqrt = sqrt(JPsiBkgErr);
  Double_t jetFakeBkgErrSqrt = sqrt(jetFakeBkgErr);
  Double_t jetFakeBkgAroundJPsiErrSqrt = sqrt(jetFakeBkgAroundJPsiErr);
  cout << "Tot. bkg.: " << setprecision(3) << B << " +/- " << errSqrt << endl;
  cout << "Upsilon bkg.: " << setprecision(3) << upsilonBkg << " +/- " << upsilonBkgErrSqrt;
  cout << endl;
  cout << "J/psi bkg.: " << setprecision(3) << JPsiBkg << " +/- " << JPsiBkgErrSqrt;
  cout << endl;
  cout << "Jet fake bkg.: " << setprecision(3) << jetFakeBkg << " +/- " << jetFakeBkgErrSqrt;
  cout << endl;
  cout << "Jet fake bkg. around J/psi: " << setprecision(3) << jetFakeBkgAroundJPsi << " +/- ";
  cout << jetFakeBkgAroundJPsiErrSqrt << endl;
  Double_t ratioUpsilonToJetFakeBkg = jetFakeBkg == 0.0 ? 0.0 : (upsilonBkg/jetFakeBkg);
  Double_t upsilonBkgFracErr2 = upsilonBkg == 0.0 ? 0.0 : (upsilonBkgErr/(upsilonBkg*upsilonBkg));
  Double_t jetFakeBkgFracErr2 = jetFakeBkg == 0.0 ? 0.0 : (jetFakeBkgErr/(jetFakeBkg*jetFakeBkg));
  Double_t ratioUpsilonToJetFakeBkgErr = 
    ratioUpsilonToJetFakeBkg*sqrt(upsilonBkgFracErr2 + jetFakeBkgFracErr2);
  cout << "Ratio upsilon bkg. to jet fake bkg.: " << setprecision(3) << ratioUpsilonToJetFakeBkg;
  cout << " +/- " << ratioUpsilonToJetFakeBkgErr << endl;
  Double_t ratioJPsiToJetFakeBkg = 
    jetFakeBkgAroundJPsi == 0.0 ? 0.0 : (JPsiBkg/jetFakeBkgAroundJPsi);
  Double_t JPsiBkgFracErr2 = JPsiBkg == 0.0 ? 0.0 : (JPsiBkgErr/(JPsiBkg*JPsiBkg));
  Double_t jetFakeBkgAroundJPsiFracErr2 = jetFakeBkgAroundJPsi == 0.0 ? 
    0.0 : (jetFakeBkgAroundJPsiErr/(jetFakeBkgAroundJPsi*jetFakeBkgAroundJPsi));
  Double_t ratioJPsiToJetFakeBkgErr = 
    ratioJPsiToJetFakeBkg*sqrt(JPsiBkgFracErr2 + jetFakeBkgAroundJPsiFracErr2);
  cout << "Ratio JPsi bkg. to jet fake bkg.: " << setprecision(3) << ratioJPsiToJetFakeBkg;
  cout << " +/- " << ratioJPsiToJetFakeBkgErr << endl;
  for (vector<TH1F*>::iterator iIsoSig = isoSig.begin(); iIsoSig != isoSig.end(); 
       ++iIsoSig) {
    (*iIsoSig)->Draw("HISTSAME");
    Double_t statErrIsoSig;
    Double_t S = (*iIsoSig)->IntegralAndError(5/*17*/, -1, statErrIsoSig);
    cout << "Region A signal " << iIsoSig - isoSig.begin() << ", m > 4: ";
    cout << setprecision(3) << S << " +/- " << statErrIsoSig;
    cout << " (S/sqrt(S+B) = " << S/sqrt(S + B) << ")" << endl;
  }
  isoData->Draw("ESAME");
  isoData->GetYaxis()->SetRangeUser(0.01, 100000.0);
  legendBkgMain5.Draw();

  //make pull plots for nominal and systematic variations of jet fake and resonance backgrounds
  outCanvas.cd(2);

  /*- there are 3 jet fake background shape templates: region B MC, nominal, and region D data
    - there are 3 resonance background shape templates: exponent fixed from region C mass 
    sidebands, nominal, and exponent fixed from region D
    - define pull as pull(resonance bkg. template, jet fake bkg. template) ==> 9 possible pulls
    -                  reg. B MC | nom. | reg. D data
    -                  ------------------------------
    - reg. C sideband |__________|______|____________|
    - nom.            |__________|______|____________|
    - reg. D          |          |      |            |
    -                  ------------------------------
    - option 1 is to show all 9 separately (current choice)
    - option 2 is to show the nominal, the min, and the max
    - this applies to the high-MT bin only*/

  /*- pull(nom., nom.)
    - note that there is a small correlation between the normalization errors of the jet fake and 
    resonance backgrounds, but they are far outweighed by the statistical errors so the 
    correlation is ignored*/
  cout << "Nominal -- ";
//   makeAndFormatPullHistogram(isoData, resBkg, nonIsoData, kAzure, "HIST");
  makeAndFormatPullHistogram(isoData, resBkg, nonIsoData, nonIsoWNonIsoData, nonIsoMCHist, kAzure, 
			     "HIST");

//   //high-MT bin only
//   if (resBkgRegCkFixed != NULL) {

//     //pull(reg. C sideband, reg. B MC)
//     cout << "J/psi bkg. shape from reg. C, jet fake bkg. shape from reg. B MC -- ";
//     makeAndFormatPullHistogram(isoData, resBkgRegCkFixed, nonIsoMCHist, kAzure - 1, 
// 			       "HISTSAME");

//     //pull(reg. C sideband, nom.)
//     cout << "J/psi bkg. shape from reg. C, nominal jet fake bkg. shape -- ";
//     makeAndFormatPullHistogram(isoData, resBkgRegCkFixed, nonIsoData, kAzure - 2, 
// 			       "HISTSAME");

//     //pull(reg. C sideband, reg. D data)
//     cout << "J/psi bkg. shape from reg. C, jet fake bkg. shape from reg. D MC -- ";
//     makeAndFormatPullHistogram(isoData, resBkgRegCkFixed, nonIsoWNonIsoData, kAzure - 3, 
// 			       "HISTSAME");
//   }

//   //pull(nom, reg. B MC)
//   cout << "Nominal J/psi bkg. shape, jet fake bkg. shape from reg. B MC -- ";
//   makeAndFormatPullHistogram(isoData, resBkg, nonIsoMCHist, kAzure - 4, "HISTSAME");

//   //pull(nom, reg. D data)
//   cout << "Nominal J/psi bkg. shape, jet fake bkg. shape from reg. D data -- ";
//   makeAndFormatPullHistogram(isoData, resBkg, nonIsoWNonIsoData, kAzure - 5, "HISTSAME");

//   //high-MT bin only
//   if (resBkgRegDkFixed != NULL) {

//     //pull(reg. D, reg. B MC)
//     cout << "J/psi bkg. shape from reg. D, jet fake bkg. shape from reg. B MC -- ";
//     makeAndFormatPullHistogram(isoData, resBkgRegDkFixed, nonIsoMCHist, kAzure - 6, 
// 			       "HISTSAME");

//     //pull(reg. D, nom.)
//     cout << "J/psi bkg. shape from reg. D, nominal jet fake bkg. shape -- ";
//     makeAndFormatPullHistogram(isoData, resBkgRegDkFixed, nonIsoData, kAzure - 7, 
// 			       "HISTSAME");

//     //pull(reg. D, reg. D data)
//     cout << "J/psi bkg. shape from reg. D, jet fake bkg. shape from reg. D data -- ";
//     makeAndFormatPullHistogram(isoData, resBkgRegDkFixed, nonIsoWNonIsoData, kAzure - 8, 
// 			       "HISTSAME");
//   }

  //write the canvas to file
  outCanvas.Write();

  //plot the nominal and alternative shapes for the jet fake background on the same axes
  TCanvas jetFakeBkgSystCanvas("jetFakeBkgSystCanvas", "", 600, 600);
  setCanvasOptions(jetFakeBkgSystCanvas, 0, 1, 0);
  setHistogramOptions(nonIsoMCHist, kRed, 0.7, 20, 1.0, "m_{#mu+had} (GeV)", "");
  setHistogramOptions(nonIsoWNonIsoData, kBlue, 0.7, 20, 1.0, "m_{#mu+had} (GeV)", "");
  setHistogramOptions(nonIsoData, kBlack, 0.7, 20, 1.0, "m_{#mu+had} (GeV)", "");
  nonIsoData->SetFillColor(kBlack);
  nonIsoMCHist->SetFillColor(kRed);
  nonIsoWNonIsoData->SetFillColor(kBlue);
  nonIsoData->SetFillStyle(3001);
  nonIsoMCHist->SetFillStyle(3001);
  nonIsoWNonIsoData->SetFillStyle(3001);
  nonIsoData->Draw("E");
  nonIsoMCHist->Draw("E2SAME");
  nonIsoWNonIsoData->Draw("E2SAME");
  TLegend legend(0.6, 0.6, 0.9, 0.9);
  setLegendOptions(legend, "#splitline{CMS 19.7 fb^{-1}}{Jet fake bkg. estimate}");
  legend.AddEntry(nonIsoMCHist, "All-EW (region B MC)", "pf");
  legend.AddEntry(nonIsoWNonIsoData, "All-QCD (region D data)", "pf");
  legend.AddEntry(nonIsoData, "Nominal", "lp");
  legend.Draw();
  jetFakeBkgSystCanvas.Write();
}

//create a file of properly formatted final plots
void makeFinalPlot(const pair<string, float>& isoMC, const string& isoDataFileName, 
		   const pair<string, float>& nonIsoData, const pair<string, float>& resBkg, 
		   const pair<string, float>& nonIsoWNonIsoData, const vector<string>& vars, 
		   const vector<string>& units, const vector<int>& normRegionLowerBins, 
		   const vector<int>& normRegionUpperBins, const string& outputFileName, 
		   const string& option, const bool ma9GeV)
{
  //open files
  pair<TFile*, float> isoSigBkgFile(new TFile(isoMC.first.c_str()), isoMC.second);
  TFile isoDataFile(isoDataFileName.c_str());
  pair<TFile*, float> nonIsoDataFile(new TFile(nonIsoData.first.c_str()), nonIsoData.second);
  pair<TFile*, float> resBkgFile(new TFile(resBkg.first.c_str()), resBkg.second);
  pair<TFile*, float> 
    nonIsoWNonIsoDataFile(new TFile(nonIsoWNonIsoData.first.c_str()), nonIsoWNonIsoData.second);
  TFile outStream(outputFileName.c_str(), "RECREATE");
  if (isoSigBkgFile.first->IsOpen() && isoDataFile.IsOpen() && nonIsoDataFile.first->IsOpen() && 
      outStream.IsOpen()) {

    //add plots
    for (vector<string>::const_iterator iVar = vars.begin(); iVar != vars.end(); ++iVar) {
      const unsigned int varIndex = iVar - vars.begin();
      addFinalPlot2(isoSigBkgFile, isoDataFile, nonIsoDataFile, resBkgFile, nonIsoWNonIsoDataFile, 
		    normRegionLowerBins[varIndex], normRegionUpperBins[varIndex], outStream, 
		    ma9GeV);
      addJetFakeBkgFinalPlot(isoSigBkgFile, isoDataFile, nonIsoDataFile, *iVar, units[varIndex], 
			     normRegionLowerBins[varIndex], normRegionUpperBins[varIndex], option, 
			     outStream);
      // addFinalPlot(isoSigBkgFile, isoDataFile, nonIsoDataFile, *iVar, units[varIndex], 
      // 		   normRegionLowerBins[varIndex], normRegionUpperBins[varIndex], option, 
      // 		   outStream, ma9GeV);
    }
  }
  else {
    cerr << "Error opening files " << isoMC.first << " or " << isoDataFileName << " or ";
    cerr << nonIsoData.first << " or " << outputFileName << ".\n";
    return;
  }

  //write to file
  outStream.Write();
  outStream.Close();
  isoSigBkgFile.first->Close();
  isoDataFile.Close();
  nonIsoDataFile.first->Close();
  delete isoSigBkgFile.first;
  delete nonIsoDataFile.first;
}

TH2F* twoDimTotalBackgroundHistogram(const vector<string>& files, const string& histName)
{
  //set canvas name
  string canvasName(histName + "Canvas");

  //set total background histogram pointer
  TH2F* bkgTotHist = NULL;

  //loop over files
  for (vector<string>::const_iterator iFile = files.begin(); iFile != files.end(); ++iFile) {
    const unsigned int fileIndex = iFile - files.begin();
    TFile* file = new TFile(iFile->c_str());
    TCanvas* canvas = NULL;
    TH2F* hist = NULL;
    if (file->IsOpen()) {

      //add  histogram to running tally
      file->GetObject(canvasName.c_str(), canvas);
      if (canvas != NULL) {
	hist = (TH2F*)canvas->GetPrimitive(histName.c_str());
	if (hist != NULL) {
	  if (fileIndex == 0) bkgTotHist = (TH2F*)hist->Clone();
	  else bkgTotHist->Add(hist);
	}
	else cerr << "Error: could not get " << histName << " from " << canvasName << endl;
      }
      else cerr << "Error: could not get " << canvasName << " from " << *iFile << endl;
    }
    else cerr << "Error: could not open " << *iFile << endl;

//     //close  file
//     file->Close();
//     delete file;
  }

  //return
  return bkgTotHist;
}

void printWeightsAndErrors(const string& numeratorFile, const string& denominatorFile)
{
  TFile numeratorStream(numeratorFile.c_str());
  TFile denominatorStream(denominatorFile.c_str());
  TCanvas* numeratorCanvas = NULL;
  TCanvas* denominatorCanvas = NULL;
  if (numeratorStream.IsOpen() && denominatorStream.IsOpen()) {
    numeratorStream.GetObject("tauHadPTCanvas", numeratorCanvas);
    denominatorStream.GetObject("tauHadPTCanvas", denominatorCanvas);
    if ((numeratorCanvas != NULL) && (denominatorCanvas != NULL)) {
      TH1F* numeratorHist = (TH1F*)numeratorCanvas->GetPrimitive("tauHadPT");
      TH1F* denominatorHist = (TH1F*)denominatorCanvas->GetPrimitive("tauHadPT");
      const Int_t nBinsNumerator = numeratorHist->GetNbinsX();
      const Int_t nBinsDenominator = denominatorHist->GetNbinsX();
      if (nBinsNumerator == nBinsDenominator) {
	const Double_t nEvtsNumerator = numeratorHist->Integral(0, -1);
	const Double_t nEvtsDenominator = denominatorHist->Integral(0, -1);
	if ((nEvtsNumerator != 0.0) && (nEvtsDenominator != 0.0)) {
	  const Double_t coefficient = nEvtsDenominator/nEvtsNumerator;
	  for (Int_t iBin = 1; iBin <= nBinsNumerator; ++iBin) {
	    const Double_t binContentNumerator = numeratorHist->GetBinContent(iBin);
	    const Double_t binContentDenominator = denominatorHist->GetBinContent(iBin);
	    if ((binContentNumerator != 0.0) && (binContentDenominator != 0.0)) {
	      const Double_t weight = coefficient*binContentNumerator/binContentDenominator;
	      const Double_t error = 
		weight*sqrt(1.0/binContentNumerator + 1.0/binContentDenominator + 
			    1.0/nEvtsDenominator + 1.0/nEvtsNumerator);
	      cout << iBin << " " << weight << " " << error << endl;
	    }
	    else cerr << "Error: numerator or denominator bin " << iBin << " has no entries.\n";
	  }
	}
	else cerr << "Error: numerator or denominator histogram has no entries.\n";
      }
      else {
	cerr << "Error: no. of numerator bins " << nBinsNumerator;
	cerr << " unequal to no. of denominator bins " << nBinsDenominator << ".\n";
      }
    }
    else {
      cerr << "Error: could not get canvas tauHadPTCanvas from " << numeratorFile << " or ";
      cerr << denominatorFile << ".\n";
    }
  }
  else {
    cerr << "Error: could not open files " << numeratorFile << " or " << denominatorFile;
    cerr << ".\n";
  }
  numeratorStream.Close();
  denominatorStream.Close();
}

void print2DWeights(const vector<string>& isoFiles, const vector<string>& nonIsoFiles, 
		    const string& hist)
{
  //sanity check
  const unsigned int isoFilesSize = isoFiles.size();
  const unsigned int nonIsoFilesSize = nonIsoFiles.size();
  if (isoFilesSize > nonIsoFilesSize) {
    cerr << "Error: size of isoFiles (" << isoFilesSize << ") > size of nonIsoFiles (";
    cerr << nonIsoFilesSize << ")\n";
    return;
  }

  //fill total background histograms
  TH2F* isoBkgTot = twoDimTotalBackgroundHistogram(isoFiles, hist);
  TH2F* nonIsoBkgTot = twoDimTotalBackgroundHistogram(nonIsoFiles, hist);

  /*normalize both histograms, then divide iso histogram by non-iso histogram and print result bin 
    by bin*/
  isoBkgTot->Scale(1.0/isoBkgTot->Integral(0, -1, 0, -1));
  nonIsoBkgTot->Scale(1.0/nonIsoBkgTot->Integral(0, -1, 0, -1));
  isoBkgTot->Divide(nonIsoBkgTot);
  TAxis* xAxis = isoBkgTot->GetXaxis();
  TAxis* yAxis = isoBkgTot->GetYaxis();
  cout.width(15);
  cout << left << "Tau pT (GeV)";
  cout.width(15);
  cout << left << "Muon pT (GeV)";
  cout << endl;
  for (Int_t iBinX = 1; iBinX <= isoBkgTot->GetNbinsX(); ++iBinX) {
    Double_t xBinLowEdge = xAxis->GetBinLowEdge(iBinX);
    stringstream xStream;
    xStream << xBinLowEdge << "-" << (xBinLowEdge + xAxis->GetBinWidth(iBinX));
    for (Int_t iBinY = 1; iBinY <= isoBkgTot->GetNbinsY(); ++iBinY) {
      Double_t yBinLowEdge = yAxis->GetBinLowEdge(iBinY);
      stringstream yStream;
      yStream << yBinLowEdge << "-" << (yBinLowEdge + yAxis->GetBinWidth(iBinY));
      cout.width(12);
      cout << right << xStream.str();
      cout.width(15);
      cout << right << yStream.str();
      cout.precision(3);
      cout << scientific << left << ": " << isoBkgTot->GetBinContent(iBinX, iBinY);
      cout << endl;
    }
  }
}

//make plots of hadronic tau pT to support reweighting
void plotTauHadPT(const vector<string>& fileNames, const vector<pair<Color_t, Color_t> >& colors, 
		  const vector<pair<Style_t, Style_t> >& styles, const string& outputFileName)
{
  //setup
  string thisFunction("const vector<string>& fileNames, ");
  thisFunction+="const vector<pair<Color_t, Color_t> >& colors, ";
  thisFunction+="const vector<pair<Style_t, Style_t> >& styles, const string&";
  const string histogramName("tauHadPT");
  const string canvasName(histogramName + "Canvas");
  const string stackName(histogramName + "Stack");
  const string unit("p_{T} (GeV)");
  vector<TFile*> files(fileNames.size(), NULL);
  vector<TCanvas*> canvases(fileNames.size(), NULL);
  vector<pair<TH1F*, TH1F*> > histograms(fileNames.size(), pair<TH1F*, TH1F*>(NULL, NULL));

  //loop over file names
  for (vector<string>::const_iterator iFileName = fileNames.begin(); iFileName != fileNames.end(); 
       ++iFileName) {
    const unsigned int i = iFileName - fileNames.begin();

    //open files
    files[i] = new TFile(iFileName->c_str());

    //get canvases
    if (files[i]->IsOpen()) files[i]->GetObject(canvasName.c_str(), canvases[i]);
    else cerr << errorCannotOpenFile(thisFunction, *iFileName);

    //get TH1Fs (data for regions B and C, Wh1 for region A)
    THStack* stack = NULL;
    if (canvases[i] != NULL) {
      if (i == 1) {
	canvases[i]->Draw();
	histograms[i].first = (TH1F*)canvases[i]->cd(1)->GetPrimitive(histogramName.c_str());
      }
      else histograms[i].first = (TH1F*)canvases[i]->GetPrimitive(histogramName.c_str());
      setHistogramOptions(histograms[i].first, colors[i].first, 0.7, styles[i].first, 
			  1.0/histograms[i].first->Integral(0, -1), unit.c_str(), "");

      //get stacks (MC for regions A and B, nonexistent for region C)
      TList* primitives = NULL;
      if (i == 1) primitives = canvases[i]->cd(1)->GetListOfPrimitives();
      else primitives = canvases[i]->GetListOfPrimitives();
      Int_t j = 0;
      while ((j < primitives->GetEntries()) && (stack == NULL)) {
	if (string(primitives->At(j)->ClassName()) == "THStack") {
	  if (i == 1) {
	    canvases[i]->Draw();
	    stack = (THStack*)canvases[i]->cd(1)->GetPrimitive(stackName.c_str());
	  }
	  else stack = (THStack*)canvases[i]->GetPrimitive(stackName.c_str());
	}
	++j;
      }
    }
    else cerr << errorCannotRetrieveObject(thisFunction, *iFileName, canvasName);
    TList* stackedHistograms = NULL;
    if (stack != NULL) {
      stackedHistograms = stack->GetHists();
      for (Int_t j = 0; j < stackedHistograms->GetEntries(); ++j) {
	if (j == 0) histograms[i].second = (TH1F*)stackedHistograms->At(j);
	else histograms[i].second->Add((TH1F*)stackedHistograms->At(j));
      }
      setHistogramOptions(histograms[i].second, colors[i].second, 0.7, styles[i].second, 
			  1.0/histograms[i].second->Integral(0, -1), unit.c_str(), "");
    }
  }

  //open output file
  TFile outputFile(outputFileName.c_str(), "RECREATE");
  if (outputFile.IsOpen()) {

    //plot region A background MC vs. region B background MC
    TCanvas AMCVsBMCCanvas("AMCVsBMCCanvas", "", 600, 600);
    setCanvasOptions(AMCVsBMCCanvas, 1, 0, 0);
    TLegend AMCVsBMCLegend(0.35, 0.55, 0.75, 0.75);
    setLegendOptions(AMCVsBMCLegend, "MC (excluding QCD) normalized to 1");
    AMCVsBMCLegend.AddEntry(histograms[0].second, "Signal region", "lp");
    AMCVsBMCLegend.AddEntry(histograms[1].second, "Control region", "lp");
    AMCVsBMCCanvas.cd();
    histograms[0].second->Draw("E");
    histograms[1].second->Draw("ESAME");
    histograms[0].second->GetXaxis()->SetRange(1, histograms[0].second->GetNbinsX() - 1);
    histograms[1].second->GetXaxis()->SetRange(1, histograms[1].second->GetNbinsX() - 1);
    AMCVsBMCLegend.Draw();
    outputFile.cd();
    AMCVsBMCCanvas.Write();

    //plot region A background MC vs. region C data
    TCanvas AMCVsCDataCanvas("AMCVsCDataCanvas", "", 600, 600);
    setCanvasOptions(AMCVsCDataCanvas, 1, 0, 0);
    TLegend AMCVsCDataLegend(0.35, 0.55, 0.75, 0.75);
    setLegendOptions(AMCVsCDataLegend, "Normalized to 1");
    AMCVsCDataLegend.AddEntry(histograms[0].second, "Signal region", "lp");
    AMCVsCDataLegend.AddEntry(histograms[2].first, "Isolated W muon QCD control region", "lp");
    AMCVsCDataCanvas.cd();
    histograms[0].second->Draw("E");
    histograms[2].first->Draw("ESAME");
    histograms[0].second->GetXaxis()->SetRange(1, histograms[0].second->GetNbinsX() - 1);
    histograms[2].first->GetXaxis()->SetRange(1, histograms[2].first->GetNbinsX() - 1);
    AMCVsCDataLegend.Draw();
    outputFile.cd();
    AMCVsCDataCanvas.Write();

    //plot region A background MC vs. region B background MC vs. region B data
    TCanvas AMCVsBMCVsBDataCanvas("AMCVsBMCVsBDataCanvas", "", 600, 600);
    setCanvasOptions(AMCVsBMCVsBDataCanvas, 1, 0, 0);
    TLegend AMCVsBMCVsBDataLegend(0.35, 0.55, 0.75, 0.75);
    setLegendOptions(AMCVsBMCVsBDataLegend, "Normalized to 1");
    AMCVsBMCVsBDataLegend.AddEntry(histograms[1].first, "Data control region", "lp");
    AMCVsBMCVsBDataLegend.AddEntry(histograms[0].second, "MC signal region", "lp");
    AMCVsBMCVsBDataLegend.AddEntry(histograms[1].second, "MC control region", "lp");
    AMCVsBMCVsBDataCanvas.cd();
    histograms[1].first->Draw("E");
    histograms[0].second->Draw("ESAME");
    histograms[1].second->Draw("ESAME");
//     histograms[1].first->GetXaxis()->SetRange(1, histograms[1].first->GetNbinsX() - 1);
//     histograms[0].second->GetXaxis()->SetRange(1, histograms[0].second->GetNbinsX() - 1);
//     histograms[1].second->GetXaxis()->SetRange(1, histograms[1].second->GetNbinsX() - 1);
    AMCVsBMCVsBDataLegend.Draw();
    outputFile.cd();
    AMCVsBMCVsBDataCanvas.Write();
  }

  //close and delete
  outputFile.Write();
  outputFile.Close();
  deleteStreams(files);
}

//get stack from canvas and create a histogram of the sum of the stacks
template<typename T>
T* getStackSumHist(TFile& file, const string& stackName, const string& canvasName, 
		   const unsigned int pad)
{
  T* stackSumHist = NULL;
  THStack* stack = getObjectFromCanvas<THStack>(file, stackName, canvasName, pad);
  TList* hists = NULL;
  if (stack != NULL) {
    hists = stack->GetHists();
    if (hists != NULL) {
      for (Int_t iHist = 0; iHist < hists->GetEntries(); ++iHist) {
	T* hist = (T*)hists->At(iHist)->Clone();
	if (iHist == 0) stackSumHist = hist;
	else stackSumHist->Add(hist);
      }
    }
  }
  return stackSumHist;
}

//divide and scale histograms, rebinning optionally
TH1* divideAndScale(TH1* hist1, TH1* hist2, const double scale, const vector<double>& bins)
{
  const Double_t hist1Integral = hist1->Integral(0, -1);
  TH1* rebinnedHist1 = NULL;
  TH1* rebinnedHist2 = NULL;
  if (bins.size() > 0) {
    rebinnedHist1 = hist1->Rebin(bins.size() - 1, hist1->GetName(), &bins[0]);
    rebinnedHist2 = hist2->Rebin(bins.size() - 1, hist2->GetName(), &bins[0]);
  }
  else {
    rebinnedHist1 = hist1;
    rebinnedHist2 = hist2;
  }
  if (rebinnedHist1 != NULL) {
    rebinnedHist1->Divide(rebinnedHist2);
    rebinnedHist1->Scale(scale == 0.0 ? rebinnedHist2->Integral(0, -1)/hist1Integral : scale);
  }
  return rebinnedHist1;
}

//plot the fraction of events in one part of histogram to the total vs. x
TH1* integralRatioX(TH2* hist, const double scale, const unsigned int factor, const int yBinLow, 
		    const int yBinHigh)
{
  TH1D* integralRatioHist = NULL;
  TH2* rebinnedHist = NULL;
  if (factor > 1) rebinnedHist = hist->RebinX(factor, hist->GetName());
  else rebinnedHist = hist;
  if (rebinnedHist != NULL) {
    rebinnedHist->Scale(scale == 0.0 ? 1.0/rebinnedHist->Integral(0, -1, 0, -1, "") : scale);
    integralRatioHist = 
      rebinnedHist->ProjectionX((string(rebinnedHist->GetName()) + "IntegralRatio").c_str(), 
				yBinLow, yBinHigh, "e");
    for (Int_t iXBin = 1; iXBin <= rebinnedHist->GetNbinsX(); ++iXBin) {
      Double_t partIntegralErr = integralRatioHist->GetBinError(iXBin);
      Double_t partIntegral = integralRatioHist->GetBinContent(iXBin);
      Double_t fullIntegralErr = 0.0;
      Double_t fullIntegral = 
	rebinnedHist->IntegralAndError(iXBin, iXBin, 0, -1, fullIntegralErr, "");
      integralRatioHist->SetBinContent(iXBin, partIntegral/fullIntegral);
      integralRatioHist->
	SetBinError(iXBin, (partIntegral/fullIntegral)*
		    sqrt((partIntegralErr*partIntegralErr)/(partIntegral*partIntegral) + 
			 (fullIntegralErr*fullIntegralErr)/(fullIntegral*fullIntegral)));
    }
  }
  return integralRatioHist;
}

//draw histogram
void draw(TFile& file, TCanvas& canvas, TH1* hist, const string& fitFunction, 
	  const string& drawOpt)
{
  file.cd();
  canvas.cd();
//   TF1 fit("fit", "[0]*exp([1]*x) + [2]", 0.0, 600.0);
//   fit.SetParameter(0, 3.0);
//   fit.SetParameter(1, -0.05);
//   fit.SetParameter(2, 3.0);
  if (fitFunction != "") hist->Fit(fitFunction.c_str());
  hist->Draw(drawOpt.c_str());
}

//format histogram
void format(TH1* hist, const float yAxisLow, const float yAxisHigh, const Color_t color, 
	    const double scale)
{
  hist->GetYaxis()->SetRangeUser(yAxisLow, yAxisHigh);
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
  hist->Scale(scale == 0.0 ? 1.0/hist->Integral(0, -1) : scale);
}

//get version from file name
string version(const string& fileName)
{
  string retVal;
  const size_t begPos = fileName.find("_v");
  const size_t endPos1 = fileName.find("_", begPos + 1);
  const size_t endPos2 = fileName.find(".root", begPos + 1);
  if (begPos != string::npos) {
    size_t endPos = string::npos;
    if (endPos1 != string::npos) endPos = endPos1;
    else if (endPos2 != string::npos) endPos = endPos2;
    if (endPos != string::npos) retVal = fileName.substr(begPos + 1, endPos - begPos - 1);
  }
  return retVal;
}

//plot ratio of 2 histograms and fit
void divideAndFit(const string& fileName1, const string& fileName2, const string& outputFileName, 
		  const vector<string>& outputCanvasName, const bool stack, 
		  const vector<string>& histName1, const vector<string>& histName2, 
		  const vector<string>& canvasName1, const vector<string>& canvasName2, 
		  const unsigned int pad, const double scale, const vector<vector<double> >& bins, 
		  const string& fitFunction, const float yAxisLow, const float yAxisHigh)
{
  TFile file1(fileName1.c_str());
  TFile file2(fileName2.c_str());
  TFile outputFile(outputFileName.c_str(), "RECREATE");
  if (file1.IsOpen() && file2.IsOpen()) {
    for (vector<string>::const_iterator iHist = histName1.begin(); iHist != histName1.end(); 
	 ++iHist) {
      const unsigned int i = iHist - histName1.begin();
      TH1F* hist1 = NULL;
      TH1F* hist2 = NULL;
      if (stack) {
	hist1 = getStackSumHist<TH1F>(file1, histName1[i].c_str(), canvasName1[i].c_str(), pad);
	hist2 = getStackSumHist<TH1F>(file2, histName2[i].c_str(), canvasName2[i].c_str(), pad);
      }
      else {
	hist1 = 
	  getObjectFromCanvas<TH1F>(file1, histName1[i].c_str(), canvasName1[i].c_str(), pad);
	hist2 = 
	  getObjectFromCanvas<TH1F>(file2, histName2[i].c_str(), canvasName2[i].c_str(), pad);
      }
      if ((hist1 != NULL) && (hist2 != NULL)) {
	TH1* outputHist = divideAndScale(hist1, hist2, scale, bins[i]);
	TCanvas outputCanvas(outputCanvasName[i].c_str(), "", 600, 600);
	setCanvasOptions(outputCanvas, 1, 0, 0);
	draw(outputFile, outputCanvas, outputHist, fitFunction.c_str(), "");
	format(outputHist, yAxisLow, yAxisHigh, kBlack, 1.0);
	outputCanvas.Write();
// 	cerr << "-------\n";
// 	for (Int_t iBin = 1; iBin <= outputHist->GetNbinsX(); ++iBin) {
// 	  cerr << iBin << " " << outputHist->GetBinContent(iBin) << " +/- ";
// 	  cerr << outputHist->GetBinError(iBin) << endl;
// 	}
// 	cerr << "-------\n";
	outputFile.Write();
      }
    }
  }
  else cerr << "Error: could not open files " << fileName1 << " or " << fileName2 << ".\n";
  file1.Close();
  file2.Close();
  outputFile.Close();
}

//plot 1D histogram from 2D histogram
void plot1DFrom2D(const vector<string>& fileNames, const string& outputFileName, 
		  const vector<string>& outputCanvasTags, const vector<string>& histNames, 
		  const vector<string>& canvasNames, const double scale, 
		  const unsigned int factor, const int yBinLow, const int yBinHigh, 
		  const float yAxisLow, const float yAxisHigh)
{
  if ((fileNames.size() != outputCanvasTags.size()) || (histNames.size() != canvasNames.size())) {
    cerr << "Error: fileNames.size() = " << fileNames.size() << " but outputCanvasTags.size() = ";
    cerr << outputCanvasTags.size() << "; histNames.size() = " << histNames.size() << " but ";
    cerr << "canvasNames.size() = " << canvasNames.size() << ".\n";
    return;
  }
  TFile outputFile(outputFileName.c_str(), "RECREATE");
  for (vector<string>::const_iterator iFileName = fileNames.begin(); iFileName != fileNames.end(); 
       ++iFileName) {
    const unsigned int i = iFileName - fileNames.begin();
    TFile file(iFileName->c_str());
    if (file.IsOpen()) {
      for (vector<string>::const_iterator iHist = histNames.begin(); iHist != histNames.end(); 
	   ++iHist) {
	const unsigned int j = iHist - histNames.begin();
	TH2F* hist = NULL;
	hist = getObjectFromCanvas<TH2F>(file, histNames[j].c_str(), canvasNames[j].c_str(), 0);
	if (hist != NULL) {
	  TH1* outputHist = integralRatioX(hist, scale, factor, yBinLow, yBinHigh);
	  TCanvas outputCanvas((string(outputHist->GetName()) + outputCanvasTags[i]).c_str(), "", 
			       600, 600);
	  setCanvasOptions(outputCanvas, 1, 0, 0);
	  draw(outputFile, outputCanvas, outputHist, "", "");
	  format(outputHist, yAxisLow, yAxisHigh, kBlack, 1.0);
	  outputCanvas.Write();
	}
      }
    }
    else cerr << "Error: could not open file " << *iFileName << ".\n";
    file.Close();
  }
  outputFile.Close();
}

//plot 2 histograms on one canvas
void plot2Histograms(const vector<string>& fileName1, const vector<string>& fileName2, 
		     const string& outputFileName, const vector<string>& outputCanvasTags, 
		     const vector<bool>& stack, const vector<string>& histName1, 
		     const vector<string>& histName2, const vector<string>& canvasName1, 
		     const vector<string>& canvasName2, const vector<unsigned int>& pad, 
		     const double scale, const vector<float>& yAxisLow, 
		     const vector<float>& yAxisHigh)
{
  if ((fileName1.size() != fileName2.size()) || (fileName2.size() != outputCanvasTags.size()) || 
      (outputCanvasTags.size() != stack.size()) || (stack.size() != pad.size()) || 
      (histName1.size() != histName2.size()) || (histName2.size() != yAxisLow.size()) || 
      (yAxisLow.size() != yAxisHigh.size())) {
    cerr << "Error:\nfileName1.size() = " << fileName1.size() << endl;
    cerr << "fileName2.size() = " << fileName2.size() << endl;
    cerr << "outputCanvasTags.size() = " << outputCanvasTags.size() << endl;
    cerr << "stack.size() = " << stack.size() << endl;
    cerr << "histName1.size() = " << histName1.size() << endl;
    cerr << "histName2.size() = " << histName2.size() << endl;
    cerr << "pad.size() = " << pad.size() << endl;
    cerr << "yAxisLow.size() = " << yAxisLow.size() << endl;
    cerr << "yAxisHigh.size() = " << yAxisHigh.size() << endl;
    return;
  }
  TFile outputFile(outputFileName.c_str(), "RECREATE");
  for (vector<string>::const_iterator iFileName1 = fileName1.begin(); 
       iFileName1 != fileName1.end(); ++iFileName1) {
    const unsigned int i = iFileName1 - fileName1.begin();
    TFile file1(iFileName1->c_str());
    TFile file2(fileName2[i].c_str());
    if (file1.IsOpen() && file2.IsOpen()) {
      for (vector<string>::const_iterator iHist = histName1.begin(); iHist != histName1.end(); 
	   ++iHist) {
	const unsigned int j = iHist - histName1.begin();
	TH1F* hist1 = NULL;
	TH1F* hist2 = NULL;
	if (stack[i]) {
	  hist1 = getStackSumHist<TH1F>(file1, (histName1[j] + "Stack").c_str(), 
					canvasName1[j].c_str(), pad[i]);
	  hist2 = getStackSumHist<TH1F>(file2, (histName2[j] + "Stack").c_str(), 
					canvasName2[j].c_str(), pad[i]);
	}
	else {
	  hist1 = 
	    getObjectFromCanvas<TH1F>(file1, histName1[j].c_str(), canvasName1[j].c_str(), pad[i]);
	  hist2 = 
	    getObjectFromCanvas<TH1F>(file2, histName2[j].c_str(), canvasName2[j].c_str(), pad[i]);
	}
	if ((hist1 != NULL) && (hist2 != NULL)) {
	  TCanvas 
	    outputCanvas((string(hist1->GetName()) + outputCanvasTags[i]).c_str(), "", 600, 600);
	  setCanvasOptions(outputCanvas, 1, 0, 0);
	  TLegend outputLegend(0.4, 0.55, 0.8, 0.75);
	  setLegendOptions(outputLegend, (string(hist1->GetName()) + outputCanvasTags[i]).c_str());
	  outputLegend.AddEntry(hist1, version(*iFileName1).c_str(), "lp");
	  outputLegend.AddEntry(hist2, version(fileName2[i]).c_str(), "lp");
	  draw(outputFile, outputCanvas, hist1, "", "E");
	  draw(outputFile, outputCanvas, hist2, "", "ESAME");
	  format(hist1, yAxisLow[j], yAxisHigh[j], kBlack, scale);
	  format(hist2, yAxisLow[j], yAxisHigh[j], kRed, 
		 hist1->GetBinContent(1)/hist2->GetBinContent(1));
	  outputLegend.Draw();
	  outputCanvas.Write();
	}
      }
    }
    else cerr << "Error: could not open files " << *iFileName1 << " or " << fileName2[i] << ".\n";
    file1.Close();
    file2.Close();
  }
  outputFile.Close();
}

//plot ratio of region C and B data tau pT spectra and fit
void calculateTauPTWeightsFromFit(const string& regionCDataFileName, 
				  const string& regionBDataFileName, 
				  const string& tauPTWeightsFileName)
{
  const string histName("tauHadPT");
  const string canvasName(histName + "Canvas");
  const vector<string> histNames(1, histName);
  const vector<string> canvasNames(1, canvasName);
  divideAndFit(regionCDataFileName, regionBDataFileName, tauPTWeightsFileName, 
	       vector<string>(1, "tauPTWeightsCanvas"), false, histNames, histNames, canvasNames, 
	       canvasNames, 1, 0.0, vector<vector<double> >(1, vector<double>()), "expo", 0.0, 
	       2.0);
}

//plot fake rate (ratio of isolated to non-isolated tau events)
void plotFakeRate(const string& isoFileName, const string& nonIsoFileName, 
		  const string& fakeRateFileName, const bool stack = false)
{
  vector<string> histNames;
  histNames.push_back("tauHadPT");
  histNames.push_back("tauHadEta");
  histNames.push_back("muHadMass");
  histNames.push_back("tauHadDecayMode");
  vector<string> canvasNames;
  vector<string> outputCanvasNames;
  for (vector<string>::const_iterator iHist = histNames.begin(); iHist != histNames.end(); 
       ++iHist) {
    const unsigned int i = iHist - histNames.begin();
    canvasNames.push_back(histNames[i] + "Canvas");
    outputCanvasNames.push_back(histNames[i] + "FakeRateCanvas");
    if (stack) {
      histNames[i]+="Stack";
    }
  }
  vector<double> tauHadEtaBins;
  tauHadEtaBins.push_back(-2.3);
  tauHadEtaBins.push_back(-1.9);
  tauHadEtaBins.push_back(-1.5);
  tauHadEtaBins.push_back(-1.1);
  tauHadEtaBins.push_back(-0.7);
  tauHadEtaBins.push_back(-0.3);
  tauHadEtaBins.push_back(0.1);
  tauHadEtaBins.push_back(0.5);
  tauHadEtaBins.push_back(0.9);
  tauHadEtaBins.push_back(1.3);
  tauHadEtaBins.push_back(1.7);
  tauHadEtaBins.push_back(2.1);
  tauHadEtaBins.push_back(2.3);
  vector<vector<double> > bins;
  bins.push_back(vector<double>());
  bins.push_back(tauHadEtaBins);
  bins.push_back(vector<double>());
  bins.push_back(vector<double>());
  divideAndFit(isoFileName, nonIsoFileName, fakeRateFileName, outputCanvasNames, stack, histNames, 
	       histNames, canvasNames, canvasNames, 1, 1.0, bins, "", 0.0001, 0.3);
}

//plot ratio between two fake rates
void plotFakeRateRatio(const string& dataFileName, const string& MCFileName, 
		       const string& fakeRateRatioFileName)
{
  vector<string> histNames;
  histNames.push_back("tauHadPT");
  histNames.push_back("tauHadEta");
  histNames.push_back("muHadMass");
  histNames.push_back("tauHadDecayMode");
  vector<string> canvasNames;
  vector<string> outputCanvasNames;
  for (vector<string>::const_iterator iHist = histNames.begin(); iHist != histNames.end(); 
       ++iHist) {
    const unsigned int i = iHist - histNames.begin();
    canvasNames.push_back(histNames[i] + "FakeRateCanvas");
    outputCanvasNames.push_back(histNames[i] + "FakeRateRatioCanvas");
  }
  divideAndFit(dataFileName, MCFileName, fakeRateRatioFileName, outputCanvasNames, false, 
	       histNames, histNames, canvasNames, canvasNames, 1, 1.0, 
	       vector<vector<double> >(histNames.size(), vector<double>()), "", 0.1, 10.0);
}

//plot ratio of tail of mu+had mass histogram to total vs. CSV score
void plotMuHadMassTailToFullRatioVsCSVScore(const vector<string>& fileNames, 
					    const string& outputFileName, 
					    const vector<string>& outputCanvasTags)
{
  plot1DFrom2D(fileNames, outputFileName, outputCanvasTags, 
	       vector<string>(1, "muHadMassVsCSVScore"), 
	       vector<string>(1, "muHadMassVsCSVScoreCanvas"), 1.0, 5, 5, -1, 0.0, 0.5);
}

//compare the same plot from 2 versions of the analysis
void compare2Versions(const vector<string>& fileName1, const vector<string>& fileName2, 
		      const string& outputFileName, const vector<string>& outputCanvasTags, 
		      const vector<bool>& stack, const vector<unsigned int>& pad)
{
  vector<string> histNames;
  histNames.push_back("muHadMass");
  histNames.push_back("muHadMass1Prong");
  histNames.push_back("muHadMass1Prong1Pi0");
  histNames.push_back("muHadMass1Prong2Pi0");
  histNames.push_back("muHadMass3Prong");
  vector<string> canvasNames;
  for (vector<string>::const_iterator iHist = histNames.begin(); iHist != histNames.end(); 
       ++iHist) {
    const unsigned int i = iHist - histNames.begin();
    canvasNames.push_back(histNames[i] + "Canvas");
  }
  plot2Histograms(fileName1, fileName2, outputFileName, outputCanvasTags, stack, 
		  histNames, histNames, canvasNames, canvasNames, pad, 1.0, 
		  vector<float>(5, 0.001), vector<float>(5, 1000.0));
}

//plot fit information for different selections
void plotFitInfo()
{
  Double_t sels[8] = {0, 1, 2, 3, 4, 5, 6, 7};
//   string selNames[8] = {"Baseline", "+ jet", "+#slash{E}_{T}", "+ jet + #slash{E}_{T}", "+M_{T}", 
// 			"+M_{T} + jet", "+ jet (40 GeV), "+ jet (20 GeV)"};
  Double_t selErrs[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Double_t par0[8] = {1.2, 1.4, 1.4, 1.1, 1.2, 1.2, 1.4, 1.3};
  Double_t par1[8] = {-0.059, -0.067, -0.07, -0.04, -0.05, -0.05, -0.06, -0.061};
  Double_t par0Err[8] = {0.1, 0.2, 0.2, 0.4, 0.2, 0.5, 0.2, 0.1};
  Double_t par1Err[8] = {0.005, 0.009, 0.01, 0.02, 0.01, 0.02, 0.01, 0.007};
  Double_t chi2PerDOF[8] = {1.18, 1.17, 1.09, 0.72, 0.39, 0.93, 1.02, 1.07};
  TGraphErrors par0Graph(8, sels, par0, selErrs, par0Err);
  TGraphErrors par1Graph(8, sels, par1, selErrs, par1Err);
  TGraph chi2PerDOFGraph(8, sels, chi2PerDOF);
  par0Graph.SetMarkerColor(kRed);
  par1Graph.SetMarkerColor(kRed);
  chi2PerDOFGraph.SetMarkerColor(kRed);
  par0Graph.SetLineColor(kRed);
  par1Graph.SetLineColor(kRed);
  chi2PerDOFGraph.SetLineColor(kRed);
  TCanvas par0Canvas("par0Canvas", "", 600, 600);
  TCanvas par1Canvas("par1Canvas", "", 600, 600);
  TCanvas chi2PerDOFCanvas("chi2PerDOFCanvas", "", 600, 600);
  setCanvasOptions(par0Canvas, 1, 0, 0);
  setCanvasOptions(par1Canvas, 1, 0, 0);
  setCanvasOptions(chi2PerDOFCanvas, 1, 0, 0);
  TFile out("fit_info.root", "RECREATE");
  out.cd();
  par0Canvas.cd();
  par0Graph.Draw("AP");
//   par0Graph.GetHistogram()->GetXaxis()->SetBinLabel(1, selNames[0].c_str());
//   par0Graph.GetHistogram()->GetXaxis()->SetBinLabel(25, selNames[1].c_str());
//   par0Graph.GetHistogram()->GetXaxis()->SetBinLabel(75, selNames[2].c_str());
//   par0Graph.GetHistogram()->GetXaxis()->SetBinLabel(100, selNames[3].c_str());
  par0Canvas.Write();
  par1Canvas.cd();
  par1Graph.Draw("AP");
//   par1Graph.GetHistogram()->GetXaxis()->SetBinLabel(1, selNames[0].c_str());
//   par1Graph.GetHistogram()->GetXaxis()->SetBinLabel(25, selNames[1].c_str());
//   par1Graph.GetHistogram()->GetXaxis()->SetBinLabel(75, selNames[2].c_str());
//   par1Graph.GetHistogram()->GetXaxis()->SetBinLabel(100, selNames[3].c_str());
  par1Canvas.Write();
  chi2PerDOFCanvas.cd();
  chi2PerDOFGraph.Draw("AP");
//   chi2PerDOFGraph.GetHistogram()->GetXaxis()->SetBinLabel(1, selNames[0].c_str());
//   chi2PerDOFGraph.GetHistogram()->GetXaxis()->SetBinLabel(25, selNames[1].c_str());
//   chi2PerDOFGraph.GetHistogram()->GetXaxis()->SetBinLabel(75, selNames[2].c_str());
//   chi2PerDOFGraph.GetHistogram()->GetXaxis()->SetBinLabel(100, selNames[3].c_str());
  chi2PerDOFCanvas.Write();
  out.Write();
  out.Close();
}

//conservative estimate of the J/psi-->mumu or upsilon-->mumu/tautau background from region C
void estimatePeakingBackground(const string& regBQCDFileName, const string& regCDataFileName, 
			       const string& regDDataFileName, const string& var, 
			       const pair<Int_t, Int_t>& normBins, 
			       const pair<Int_t, Int_t>& peakBins)
{
  //get no. QCD events and stat. error in region B normalization region
  Double_t normRegB = 0.0;
  Double_t normRegBErr = 0.0;
  TFile regBQCDFile(regBQCDFileName.c_str());
  if (regBQCDFile.IsOpen()) {
    TH1F* hist = getObjectFromCanvas<TH1F>(regBQCDFile, var, var + "Canvas", 0);
    if (hist != NULL) {
      normRegB = hist->IntegralAndError(normBins.first, normBins.second, normRegBErr);
    }
    else {
      cerr << "Error getting histogram " << var << " from canvas " << var << "Canvas from file ";
      cerr << regBQCDFileName << ".\n";
    }
  }
  else cerr << "Error opening file " << regBQCDFileName << ".\n";
  regBQCDFile.Close();

  //get no. events and stat. error in region C peaking region
  Double_t peakRegC = 0.0;
  Double_t peakRegCErr = 0.0;
  TFile regCDataFile(regCDataFileName.c_str());
  if (regCDataFile.IsOpen()) {
    TH1F* hist = getObjectFromCanvas<TH1F>(regCDataFile, var, var + "Canvas", 0);
    if (hist != NULL) {
      peakRegC = hist->IntegralAndError(peakBins.first, peakBins.second, peakRegCErr);
    }
    else {
      cerr << "Error getting histogram " << var << " from canvas " << var << "Canvas from file ";
      cerr << regCDataFileName << ".\n";
    }
  }
  else cerr << "Error opening file " << regCDataFileName << ".\n";
  regCDataFile.Close();

  //get no. events and stat. error in region D normalization region
  Double_t normRegD = 0.0;
  Double_t normRegDErr = 0.0;
  TFile regDDataFile(regDDataFileName.c_str());
  if (regDDataFile.IsOpen()) {
    TH1F* hist = getObjectFromCanvas<TH1F>(regDDataFile, var, var + "Canvas", 0);
    if (hist != NULL) {
      normRegD = hist->IntegralAndError(normBins.first, normBins.second, normRegDErr);
    }
    else {
      cerr << "Error getting histogram " << var << " from canvas " << var << "Canvas from file ";
      cerr << regDDataFileName << ".\n";
    }
  }
  else cerr << "Error opening file " << regDDataFileName << ".\n";
  regDDataFile.Close();

  //calculate peaking background estimate and stat. error
  Double_t bkg = peakRegC*(normRegB/normRegD);
  Double_t bkgErr = bkg*sqrt((peakRegCErr*peakRegCErr)/(peakRegC*peakRegC) + 
			     (normRegBErr*normRegBErr)/(normRegB*normRegB) + 
			     (normRegDErr*normRegDErr)/(normRegD*normRegD));
  cout << "Peaking background: " << bkg << " +/- " << bkgErr << endl;
}

//make a string with decimal point in number replaced with 'p' and negative sign 
//replaced with 'm' for ROOT object name
string replacePeskyCharacters(Double_t val)
{
  stringstream valStream;
  valStream << val;
  string valString(valStream.str());
  size_t dotPos = valString.find('.');
  if (dotPos != string::npos) valString.replace(dotPos, 1, "p");
  size_t negSgnPos = valString.find('-');
  if (negSgnPos != string::npos) valString.replace(negSgnPos, 1, "m");
  return valString;
}

//save canvas as PDF
void saveCanvasAsPDF(const string& canvasName, const string& saveName, const string& dir, 
		     TFile& file)
{
  TCanvas* canvas = NULL;
  file.GetObject(canvasName.c_str(), canvas);
  if (canvas != NULL) {
    canvas->SaveAs((dir + saveName + ".pdf").c_str());
  }
  else {
    cerr << "Error getting canvas with name " << canvasName << " from file ";
    cerr << file.GetName() << ".\n";
  }
}

//get object from RooPlot
template <typename T>
T* getObjectFromRooPlot(const string& plotName, const string& objName, TFile& file)
{
  RooPlot* plot = NULL;
  file.GetObject(plotName.c_str(), plot);
  T* obj = NULL;
  if (plot != NULL) {
    obj = (T*)plot->findObject(objName.c_str());
  }
  return obj;
}

//fit region D for QCD jet fake background shape
RooFitResult* fitRegionD(const string& regDDataFileName, TFile& outFile)
{
  //return value
  RooFitResult* fitRes = NULL;

  //get region D histogram
  TH1F* regDHistRaw = NULL;
  TFile regDDataFile(regDDataFileName.c_str());
  if (regDDataFile.IsOpen()) {
    TH1F* hist = getObjectFromCanvas<TH1F>(regDDataFile, "muHadMass", "muHadMassCanvas", 0);
    if (hist != NULL) {
      regDHistRaw = (TH1F*)hist->Clone();

      //observable
      RooRealVar muHadMass("muHadMass", "m_{#mu+had} (GeV)", 0.0, 4.0);

      //build the PDF
      RooRealVar decayConst("decayConst", "Exponential decay constant", -10.0, -20.0, 0.0);
      RooExponential PDF("PDF", "Background PDF", muHadMass, decayConst);
 
      //construct a RooDataHist from the region D histogram
      RooDataHist regDRooDataHistRaw(regDHistRaw->GetName(), regDHistRaw->GetTitle(), 
				     RooArgList(muHadMass), regDHistRaw);

      //fit the PDF to the region D data
      fitRes = 
	PDF.fitTo(regDRooDataHistRaw, RooFit::SumW2Error(kFALSE), RooFit::Range(1.5, 4.0/*16.0*/), 
		  RooFit::Save(kTRUE));

      //plot region D data and PDF overlaid
      RooPlot* muHadMassFrame = muHadMass.frame();
      muHadMassFrame->SetName("frame_muHadMass_regD");
      muHadMassFrame->SetTitle("");
      muHadMassFrame->GetYaxis()->SetTitle("");
      regDRooDataHistRaw.plotOn(muHadMassFrame, RooFit::Name(regDRooDataHistRaw.GetName()));
      PDF.plotOn(muHadMassFrame, RooFit::Name(PDF.GetName()));
      muHadMassFrame->GetYaxis()->SetTitle("");

      //add legend
      TLegend legend(0.6, 0.7, 0.9, 0.9);
      setLegendOptions(legend, "CMS 19.7 fb^{-1}");
      legend.
	AddEntry(muHadMassFrame->findObject(regDRooDataHistRaw.GetName()), "Region D data", "lep");
      legend.AddEntry(muHadMassFrame->findObject(PDF.GetName()), "Exponential fit", "l");

      //draw RooPlot and legend on canvas
      TCanvas canvas(muHadMassFrame->GetName(), muHadMassFrame->GetName(), 600, 600);
      setCanvasOptions(canvas, 0, 0, 0);
      muHadMassFrame->Draw();
      legend.Draw();

      //write to file
      if (outFile.IsOpen()) {
	outFile.cd();
	muHadMassFrame->Write();
	canvas.Write();
      }
      else cerr << "Error: file " << outFile.GetName() << " is not open.\n";
    }
    else {
      cerr << "Error getting histogram muHadMass from canvas muHadMassCanvas from file ";
      cerr << regDDataFileName << ".\n";
    }
    regDDataFile.Close();
  }
  else cerr << "Error opening file " << regDDataFileName << ".\n";

  //return the fit result
  return fitRes;
}

//fit J/psi peak for configurable background shapes and return fit result
RooFitResult* fitJPsiPeakAndBackground(TFile& outputFile, const Double_t bkgDecayConstVal, 
				       const Double_t bkgDecayConstMin, 
				       const Double_t bkgDecayConstMax, const TH1* hist, 
				       const string& label)
{
  //observable
  RooRealVar muHadMass("muHadMass", "m_{#mu+had} (GeV)", /*0.0, 4.0*/1.5, 3.75);
 
  //build the signal PDF
  RooRealVar sigMean("sigMean", "m_{J/#psi} (GeV)", 3.1, 3.0, 3.2);
  RooRealVar sigWidth("sigWidth", "#sigma_{experimental}", 0.5, 0.1, 1.0);
  RooGaussian sigPDF("sigPDF", "Signal PDF", muHadMass, sigMean, sigWidth); 

  //build the background PDF
  RooRealVar bkgDecayConst("bkgDecayConst", "Exponential decay constant", bkgDecayConstVal, 
			   bkgDecayConstMin, bkgDecayConstMax);
  RooExponential bkgPDF("bkgPDF", "Background PDF", muHadMass, bkgDecayConst);

  //build signal+background PDF
  RooRealVar nSig("nSig", "No. J/#psi events", 5.0, 0.0, 50.0);
  RooRealVar nBkg("nBkg", "No. non-resonant QCD events", 60.0, 40.0, 80.0/*55.0, 30.0, 500.0*/);
  RooAddPdf sigBkgModel("sigBkgModel", "Gaussian + exponential model", RooArgList(sigPDF, bkgPDF), 
			RooArgList(nSig, nBkg));

  //construct a RooDataHist from the histogram assuming unweighted data
  RooDataHist RooDataHistRaw(hist->GetName(), hist->GetTitle(), RooArgList(muHadMass), hist);

  //fit the signal+background PDF to the data
  RooFitResult* fitRes = sigBkgModel.fitTo(RooDataHistRaw, RooFit::SumW2Error(kFALSE), 
					   RooFit::Range(1.5, 3.75), RooFit::Extended(kTRUE), 
					   RooFit::SumCoefRange("fit"), RooFit::Save(kTRUE));

  //create a histogram of the signal PDF, scaled by nSig
  //bin errors are from fit error on nSig, error due to sampling the PDF at the bin center (e.g. 
  //due to bin width) neglected because bin width is small enough that it's negligible
  TH1F* sigHist = (TH1F*)sigPDF.createHistogram(sigPDF.GetName(), muHadMass, RooFit::Binning(90));
  sigHist->Rebin(10);
  sigHist->SetName((label + "_sigPDF").c_str());
  RooRealVar* nSigPostFitRooRealVar = (RooRealVar*)fitRes->floatParsFinal().find("nSig");
  Double_t nSigPostFit = nSigPostFitRooRealVar->getVal();
  sigHist->Scale(nSigPostFit/sigHist->Integral());
  for (Int_t iBin = 0; iBin <= (sigHist->GetNbinsX() + 1); ++iBin) {
    sigHist->SetBinError(iBin, sigHist->GetBinContent(iBin)*
			 (nSigPostFitRooRealVar->getError()/nSigPostFit));
  }

  //plot data and signal+background PDF overlaid
  RooPlot* muHadMassFrame = muHadMass.frame();
  muHadMassFrame->SetName(label.c_str());
  muHadMassFrame->SetTitle("");
  RooDataHistRaw.plotOn(muHadMassFrame, RooFit::Name(RooDataHistRaw.GetName()));
  sigBkgModel.plotOn(muHadMassFrame, RooFit::Name(sigBkgModel.GetName()));
  sigBkgModel.plotOn(muHadMassFrame, RooFit::Components(bkgPDF), RooFit::LineStyle(kDashed), 
		     RooFit::Name(bkgPDF.GetName()));
  sigBkgModel.plotOn(muHadMassFrame, RooFit::Components(sigPDF), RooFit::LineStyle(kDashed), 
		     RooFit::LineColor(kRed), RooFit::Name(sigPDF.GetName()));
  muHadMassFrame->GetYaxis()->SetTitle("");
  setAxisOptions(muHadMassFrame->GetXaxis(), 0.04, 0.9, muHadMassFrame->GetXaxis()->GetTitle());

  //add legend
  TLegend legend(0.6, 0.7, 0.9, 0.9);
  setLegendOptions(legend, "CMS 19.7 fb^{-1}");
  legend.AddEntry(muHadMassFrame->findObject(RooDataHistRaw.GetName()), "Region C data", "lep");
  legend.AddEntry(muHadMassFrame->findObject(sigBkgModel.GetName()), "Total fit", "l");
  legend.AddEntry(muHadMassFrame->findObject(bkgPDF.GetName()), "Background fit", "l");
  legend.AddEntry(muHadMassFrame->findObject(sigPDF.GetName()), "Signal fit", "l");

  //draw RooPlot and legend on canvas
  TCanvas canvas(muHadMassFrame->GetName(), muHadMassFrame->GetName(), 600, 600);
  setCanvasOptions(canvas, 0, 0, 0);
  muHadMassFrame->Draw();
  legend.Draw();

  //write the RooPlot to a file
  if (outputFile.IsOpen()) {
    outputFile.cd();
    muHadMassFrame->Write();
    sigHist->Write();
    canvas.Write();
  }
  else cerr << "File " << outputFile.GetName() << " is not open.\n";

  //return the fit result
  return fitRes;
}

/*create a histogram with only the upsilon background component, write to file, and return 
  histogram*/
TH1F* createUpsilonBackgroundHistogram(const TH1F* inputHist, TFile& file, const string& label, 
				       const Int_t nNewBins = -1, 
				       const vector<Double_t>& newBinEdges = vector<Double_t>())
{
  //clone the input histogram, assumed to be the region C histogram scaled properly for the 
  //upsilon bins
  TH1F* resBkg = (TH1F*)inputHist->Clone();
  resBkg->SetName((label + "_resBkg").c_str());

  //zero out the non-upsilon bins
  for (Int_t iBin = 0; iBin < 17; ++iBin) {
    resBkg->SetBinContent(iBin, 0.0);
    resBkg->SetBinError(iBin, 0.0);
  }

  //optionally rebin
  TH1F* resBkgRebinned = NULL;
  if ((nNewBins > 0) && (newBinEdges.size() == (nNewBins + 1))) {
    resBkgRebinned = (TH1F*)resBkg->Rebin(nNewBins, resBkg->GetName(), &newBinEdges[0]);
  }

  //write to file and return
  file.cd();
  if (resBkgRebinned != NULL) resBkgRebinned->Write();
  else resBkg->Write();
  return ((resBkgRebinned != NULL) ? resBkgRebinned : resBkg);
}

/*- create a histogram with only the resonant background components, write to file, and return 
  histogram*/
TH1F* createResonanceBackgroundHistogram(const TH1F* inputHist, TFile& file, const string& label, 
					 const Double_t ABCDNorm, const Double_t ABCDNormErr, 
					 const Int_t nNewBins = -1, 
					 const vector<Double_t>& newBinEdges = vector<Double_t>())
{
  //calculate the ABCD normalization relative error squared
  Double_t ABCDNormRelErr2 = ABCDNorm == 0.0 ? 0.0 : (ABCDNormErr/ABCDNorm);
  ABCDNormRelErr2*=ABCDNormRelErr2;

  //clone the input histogram, assumed to be the region C histogram scaled properly for the 
  //upsilon bins
  TH1F* resBkg = (TH1F*)inputHist->Clone();
  resBkg->SetName((label + "_resBkg").c_str());

  //zero out the non-upsilon bins
  for (Int_t iBin = 0; iBin < 17; ++iBin) {
    resBkg->SetBinContent(iBin, 0.0);
    resBkg->SetBinError(iBin, 0.0);
  }

  //get the J/psi shape and add to the upsilon shape
  TH1F* JPsiHist = NULL;
  file.GetObject((label + "_sigPDF").c_str(), JPsiHist);
  if (JPsiHist != NULL) {
    for (Int_t iJPsiBin = 1; iJPsiBin <= JPsiHist->GetNbinsX(); ++iJPsiBin) {
      bool foundBin = false;
      Int_t iUpsilonBin = 1;
      while ((iUpsilonBin <= resBkg->GetNbinsX()) && !foundBin) {
	if (JPsiHist->GetBinLowEdge(iJPsiBin) == resBkg->GetBinLowEdge(iUpsilonBin)) {
	  Double_t nJPsiPreABCDNorm = JPsiHist->GetBinContent(iJPsiBin);
	  Double_t nJPsiErrPreABCDNorm = JPsiHist->GetBinError(iJPsiBin);
	  Double_t nJPsiRelErr2PreABCDNorm = 
	    nJPsiPreABCDNorm == 0.0 ? 0.0 : (nJPsiErrPreABCDNorm/nJPsiPreABCDNorm);
	  nJPsiRelErr2PreABCDNorm*=nJPsiRelErr2PreABCDNorm;
	  Double_t nJPsiPostABCDNorm = ABCDNorm*nJPsiPreABCDNorm;
	  resBkg->SetBinContent(iUpsilonBin, nJPsiPostABCDNorm);
	  resBkg->SetBinError(iUpsilonBin, 
			      nJPsiPostABCDNorm*sqrt(nJPsiRelErr2PreABCDNorm + ABCDNormRelErr2));
	  foundBin = true;
	}
	++iUpsilonBin;
      }
    }
  }
  else {
    cerr << "Error getting histogram " << label << "_sigPDF from file " << file.GetName() << ".\n";
  }

  //optionally rebin
  TH1F* resBkgRebinned = NULL;
  if ((nNewBins > 0) && (newBinEdges.size() == (nNewBins + 1))) {
    resBkgRebinned = (TH1F*)resBkg->Rebin(nNewBins, resBkg->GetName(), &newBinEdges[0]);
  }

  //print J/psi yield
  TH1F* nonNullHist = resBkgRebinned != NULL ? resBkgRebinned : resBkg;
  Int_t upperLim = 0.0;
  Int_t iBin = 1;
  while ((iBin <= nonNullHist->GetNbinsX()) && (upperLim == 0.0)) {
    if (nonNullHist->GetBinLowEdge(iBin) == 4.0/*GeV*/) upperLim = iBin - 1;
    ++iBin;
  }
  Double_t nJPsiErr = 0.0;
  Double_t nJPsi = nonNullHist->IntegralAndError(1, upperLim, nJPsiErr);
  cout << "No. J/psi events: " << setprecision(3) << nJPsi << " +/- " << nJPsiErr << endl;
  cout << "-------------------------------\n";

  //write to file and return
  file.cd();
  if (resBkgRebinned != NULL) resBkgRebinned->Write();
  else resBkg->Write();
  return ((resBkgRebinned != NULL) ? resBkgRebinned : resBkg);
}

/*- create a histogram with only the non-resonant background component, write to file, and return 
  histogram*/
TH1F* createNonResonantBackgroundHistogram(const TH1F* inputHist, TFile& file, 
					   const string& label, const Int_t nNewBins = -1, 
					   const vector<Double_t>& newBinEdges = 
					   vector<Double_t>())
{
  //clone the input histogram, assumed to be the region C histogram scaled properly for the 
  //upsilon bins
  TH1F* nonResBkg = (TH1F*)inputHist->Clone();
  nonResBkg->SetName((label + "_nonResBkg").c_str());

  //get the resonance background histogram and subtract it from the full region C histogram
  TH1F* resBkg = NULL;
  file.GetObject((label + "_resBkg").c_str(), resBkg);
  if (resBkg != NULL) {
    nonResBkg->Add(resBkg, -1.0);
    for (Int_t iBin = 1; iBin <= nonResBkg->GetNbinsX(); ++iBin) {
      Double_t nonResBkgErr2 = nonResBkg->GetBinError(iBin);
      nonResBkgErr2*=nonResBkgErr2;
      Double_t resBkgErr2 = resBkg->GetBinError(iBin);
      resBkgErr2*=resBkgErr2;
      nonResBkg->SetBinError(iBin, sqrt(nonResBkgErr2 + resBkgErr2));
    }
  }
  else {
    cerr << "Error getting histogram " << label << "_resBkg from file " << file.GetName() << ".\n";
  }

  //optionally rebin
  TH1F* nonResBkgRebinned = NULL;
  if ((nNewBins > 0) && (newBinEdges.size() == (nNewBins + 1))) {
    nonResBkgRebinned = (TH1F*)nonResBkg->Rebin(nNewBins, nonResBkg->GetName(), &newBinEdges[0]);
  }

  //write to file and return
  file.cd();
  if (nonResBkgRebinned != NULL) nonResBkgRebinned->Write();
  else nonResBkg->Write();
  return ((nonResBkgRebinned != NULL) ? nonResBkgRebinned : nonResBkg);
}

//fit region C for J/psi and upsilon signal and convert to a histogram
void estimateBoostedResonanceBackground(const string& outFileName, const string& regBQCDFileName, 
					const string& regCDataFileName, 
					const string& regDDataFileName, 
					const pair<Int_t, Int_t>& normBins, const string& HLTPath, 
					const string& MTBin, const Int_t nNewBins = -1, 
					const vector<Double_t>& newBinEdges = vector<Double_t>())
{
  //get no. QCD events and stat. error in region B normalization region
  pair<Double_t, Double_t> normRegBIntAndErr = 
    getIntegralAndError(regBQCDFileName, normBins, "muHadMass", 0);
  Double_t normRegB = normRegBIntAndErr.first;
  Double_t normRegBErr = normRegBIntAndErr.second;
  Double_t regBNormFracErr2 = normRegB == 0.0 ? 
    0.0 : (normRegBErr*normRegBErr)/(normRegB*normRegB);

  //get no. events and stat. error in region D normalization region
  pair<Double_t, Double_t> normRegDIntAndErr = 
    getIntegralAndError(regDDataFileName, normBins, "muHadMass", 0);
  Double_t normRegD = normRegDIntAndErr.first;
  Double_t normRegDErr = normRegDIntAndErr.second;
  Double_t regDNormFracErr2 = normRegD == 0.0 ? 
    0.0 : (normRegDErr*normRegDErr)/(normRegD*normRegD);

  //compute ABCD normalization factor and error
  Double_t ABCDNorm = normRegB/normRegD;
  Double_t ABCDNormErr = ABCDNorm*sqrt(regBNormFracErr2 + regDNormFracErr2);

  //get region C histogram
  TH1F* regCHistFinal = NULL;
  TH1F* regCHistRaw = NULL;
  TFile regCDataFile(regCDataFileName.c_str());
  if (regCDataFile.IsOpen()) {
    TH1F* hist = getObjectFromCanvas<TH1F>(regCDataFile, "muHadMass", "muHadMassCanvas", 0);
    if (hist != NULL) {
      regCHistFinal = (TH1F*)hist->Clone();
      regCHistRaw = (TH1F*)hist->Clone();
    }
    else {
      cerr << "Error getting histogram muHadMass from canvas muHadMassCanvas from file ";
      cerr << regCDataFileName << ".\n";
    }
  }
  else cerr << "Error opening file " << regCDataFileName << ".\n";

  //scale region C histogram by B/D
  regCHistFinal->Scale(ABCDNorm);

  //print information for gmN modeling in datacard
  cout << "---------------------------------------------------------\n";
  cout << "No. reg. C events with m >= 4 GeV before normalization: ";
  cout << regCHistRaw->Integral(17, -1) << endl;
  cout << "Normalization constant: " << ABCDNorm << endl;
  cout << "---------------------------------------------------------\n";

  //propagate stat. error on regions B and D
  for (Int_t iBin = 0; iBin <= (regCHistFinal->GetNbinsX() + 1); ++iBin) {
    Double_t regCBinContentRaw = regCHistRaw->GetBinContent(iBin);
    Double_t regCBinErrRaw = regCHistRaw->GetBinError(iBin);
    Double_t regCBinFracErr2 = regCBinContentRaw == 0.0 ? 
      0.0 : (regCBinErrRaw*regCBinErrRaw)/(regCBinContentRaw*regCBinContentRaw);
    Double_t regCBinErrFinal = regCHistFinal->GetBinContent(iBin)*
      sqrt(regCBinFracErr2 + regBNormFracErr2 + regDNormFracErr2);
    regCHistFinal->SetBinError(iBin, regCBinErrFinal);

//     //debug
//     Double_t statErr = regCHistFinal->GetBinContent(iBin)*sqrt(regCBinFracErr2);
//     Double_t normErr = 
//       regCHistFinal->GetBinContent(iBin)*sqrt(regBNormFracErr2 + regDNormFracErr2);
//     Double_t fullFracErr = regCHistFinal->GetBinContent(iBin) == 0.0 ? 0.0 : 
//       (regCBinErrFinal/regCHistFinal->GetBinContent(iBin));
//     Double_t statFracErr = sqrt(regCBinFracErr2);
//     cout << "Bin " << iBin << endl << setprecision(3) << regCHistFinal->GetBinContent(iBin);
//     cout << " +/- " << statErr << "(stat.) +/- " << normErr << "(norm.)\nFull fractional error: ";
//     cout << fullFracErr << endl << "Fractional error neglecting normalization term: ";
//     cout << statFracErr << endl;
  }

  //create output file
  TFile outFile(outFileName.c_str(), "RECREATE");

  /*in the high-MT bin, fit the J/psi peak and subtract it from the non-resonant background in 
    region C before preparing resonant and non-resonant background histograms*/
  stringstream label2;
  stringstream label4;
  stringstream label5;
  if (MTBin == "_highMT") {

    //step 1: fit with all parameters floating
    stringstream label1;
    label1 << "frame_muHadMass_-10_-20_0_regCAllParsFloat";
    RooFitResult* floatingFitRes = 
      fitJPsiPeakAndBackground(outFile, -10.0, -20.0, 0.0, regCHistRaw, label1.str());

    //step 2: fit with background shape fixed to that found in step 1
    RooRealVar* bkgDecayConstFloating = 
      (RooRealVar*)floatingFitRes->floatParsFinal().find("bkgDecayConst");
    Double_t k2 = bkgDecayConstFloating->getVal();
    label2 << "frame_muHadMass_" << replacePeskyCharacters(k2) << "_0_0_regCkFixed";
    if (bkgDecayConstFloating != NULL) {
      fitJPsiPeakAndBackground(outFile, k2, 0.0, 0.0, regCHistRaw, label2.str());
    }
    else {
      cerr << "Error: parameter bkgDecayConst not found in fit result " << floatingFitRes << ".\n";
    }

    //step 3: fit region D for the background shape only
    RooFitResult* regDBkgFitRes = fitRegionD(regDDataFileName, outFile);

    //step 4: fit with background shape fixed to that found in step 3
    RooRealVar* bkgDecayConstFixedToRegD = 
      (RooRealVar*)regDBkgFitRes->floatParsFinal().find("decayConst");
    Double_t k4 = bkgDecayConstFixedToRegD->getVal();
    label4 << "frame_muHadMass_" << replacePeskyCharacters(k4) << "_0_0_regDkFixed";
    if (bkgDecayConstFixedToRegD != NULL) {
      fitJPsiPeakAndBackground(outFile, k4, 0.0, 0.0, regCHistRaw, label4.str());
    }
    else cerr << "Error: parameter decayConst not found in fit result " << regDBkgFitRes << ".\n";

    //step 5: fit with background shape fixed to weighted average of those found in steps 1 and 3
    Double_t k2Err = bkgDecayConstFloating->getError();
    Double_t k4Err = bkgDecayConstFixedToRegD->getError();
    Double_t k5 = 
      ((k2/(k2Err*k2Err)) + (k4/(k4Err*k4Err)))/((1.0/(k2Err*k2Err)) + (1.0/(k4Err*k4Err)));
    label5 << "frame_muHadMass_" << replacePeskyCharacters(k5) << "_0_0_weightedAvgkFixed";
    fitJPsiPeakAndBackground(outFile, k5, 0.0, 0.0, regCHistRaw, label5.str());

    /*create new histograms with only the resonant components
      nominal: background fixed to weighted average of steps 1 and 3
      up/down systematic variations: background fixed to either of steps 1 or 3*/
    cout << "-------------------------------\n";
    cout << "Step 2\n";
    TH1F* resBkg2 = createResonanceBackgroundHistogram(regCHistFinal, outFile, label2.str(), 
						       ABCDNorm, ABCDNormErr, nNewBins, 
						       newBinEdges);
    cout << "Step 4\n";
    TH1F* resBkg4 = createResonanceBackgroundHistogram(regCHistFinal, outFile, label4.str(), 
						       ABCDNorm, ABCDNormErr, nNewBins, 
						       newBinEdges);
    cout << "Step 5\n";
    TH1F* resBkg5 = createResonanceBackgroundHistogram(regCHistFinal, outFile, label5.str(), 
						       ABCDNorm, ABCDNormErr, nNewBins, 
						       newBinEdges);

    /*create new histogram with only the non-resonant component
      nominal: subtracted J/psi background fixed to weighted average of steps 1 and 3
      up/down systematic variations: subtracted J/psi background fixed to either of steps 1 or 3*/
    createNonResonantBackgroundHistogram(regCHistFinal, outFile, label2.str(), nNewBins, 
					 newBinEdges);
    createNonResonantBackgroundHistogram(regCHistFinal, outFile, label4.str(), nNewBins, 
					 newBinEdges);
    createNonResonantBackgroundHistogram(regCHistFinal, outFile, label5.str(), nNewBins, 
					 newBinEdges);

    //plot all resonance background histograms on the same axes for comparison
    TCanvas resBkgSystCanvas("resBkgSystCanvas", "", 600, 600);
    setCanvasOptions(resBkgSystCanvas, 0, 0, 0);
    setHistogramOptions(resBkg2, kRed, 0.7, 20, 1.0, "m_{#mu+had} (GeV)", "");
    setHistogramOptions(resBkg4, kBlue, 0.7, 20, 1.0, "m_{#mu+had} (GeV)", "");
    setHistogramOptions(resBkg5, kBlack, 0.7, 20, 1.0, "m_{#mu+had} (GeV)", "");
    resBkg5->Draw("E");
    resBkg2->Draw("HISTSAME");
    resBkg4->Draw("HISTSAME");
    Double_t legendX1 = 0.4;
    Double_t legendX2 = 0.8;
    if (nNewBins == -1) {
      resBkg5->GetXaxis()->SetRange(7, 16);
      resBkg2->GetXaxis()->SetRange(7, 16);
      resBkg4->GetXaxis()->SetRange(7, 16);
      legendX1 = 0.2;
      legendX2 = 0.6;
    }
    TLegend legend(legendX1, 0.6, legendX2, 0.9);
    setLegendOptions(legend, ("#splitline{#splitline{CMS 19.7 fb^{-1}}{J/#psi bkg. estimate}}{" + 
			      HLTPath + "}").c_str());
    legend.AddEntry(resBkg2, "Bkg. from region C sideband fit", "l");
    legend.AddEntry(resBkg4, "Bkg. from region D fit", "l");
    legend.AddEntry(resBkg5, "Nominal", "lp");
    legend.Draw();
    outFile.cd();
    resBkgSystCanvas.Write();
  }
  else if (MTBin == "_lowMT") {
    /*in the low-MT bin, do not fit and subtract any peak in region C before preparing resonant 
      and non-resonant background histograms*/

    //create new histogram with only the upsilon component
    createUpsilonBackgroundHistogram(regCHistFinal, outFile, "noFit_weightedAvgkFixed", 
				     nNewBins, newBinEdges);

    //create new histogram with only the non-resonant component
    createNonResonantBackgroundHistogram(regCHistFinal, outFile, "noFit_weightedAvgkFixed", 
					 nNewBins, newBinEdges);
  }

  //close files
  if (outFile.IsOpen()) outFile.Close();
  else cerr << "Error: file " << outFile.GetName() << " is not open.\n";
  regCDataFile.Close();
}
