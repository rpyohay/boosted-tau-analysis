#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TProfile.h"
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
void plotNice(const string& inputFileName, const map<string, pair<string, string> >& effHistMap, 
	      const map<string, vector<string> >& binLabelMap, 
	      const map<string, string>& hist1DMap, const string& outputFileName, 
	      const string& savePath)
{
  string fnName("const string& inputFileName, const map<string, pair<string, string> >& ");
  fnName+="effHistMap, const map<string, vector<string> >& binLabelMap, ";
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

  //make efficiency plots
  plotEfficiency(in, out, effHistMap, binLabelMap, savePath);

  //make 1D histogram plots
  plot1DHistograms(in, out, hist1DMap, binLabelMap, savePath);

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
	   const unsigned int size, const bool dataMC)
{
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    if (dataMC) {
      outputCanvases.push_back(new TCanvas(iCanvasName->c_str(), "", 600, 900));
      outputCanvases[outputCanvases.size() - 1]->Divide(1, 2);
      outputCanvases[outputCanvases.size() - 1]->cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
      setCanvasOptions(*outputCanvases[outputCanvases.size() - 1]->cd(1), 1, setLogY, 0);
      setCanvasMargins(*outputCanvases[outputCanvases.size() - 1]->cd(1), 0.2, 0.2, 0.2, 0.2);
      outputCanvases[outputCanvases.size() - 1]->cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
      setCanvasOptions(*outputCanvases[outputCanvases.size() - 1]->cd(2), 1, 0, 0);
      setCanvasMargins(*outputCanvases[outputCanvases.size() - 1]->cd(2), 0.2, 0.2, 0.2, 0.2);
    }
    else {
      outputCanvases.push_back(new TCanvas(iCanvasName->c_str(), "", 600, 600));
      setCanvasOptions(*outputCanvases[outputCanvases.size() - 1], 1, setLogY, 0);
      setCanvasMargins(*outputCanvases[outputCanvases.size() - 1], 0.2, 0.2, 0.2, 0.2);
    }
    legends.push_back(new TLegend(0.4, 0.55, 0.8, 0.75));
    string stackName(*iCanvasName);
    stackName.replace(stackName.find("Canvas"), 6, "Stack");
    stacks.push_back(new THStack(stackName.c_str(), ""));
    if (legendHeaders.size() > canvasIndex) setLegendOptions(*legends[legends.size() - 1], 
							     legendHeaders[canvasIndex].c_str());
    hists.push_back(vector<T*>(size, NULL));
  }
}

template<typename T>
void setup(const vector<string>& canvasNames, vector<TCanvas*>& outputCanvases, 
	   vector<T*>& hists)
{
  vector<TLegend*> legends;
  vector<THStack*> stacks;
  vector<string> legendHeaders;
  vector<vector<TH1F*> > dummyHists;
  setup(canvasNames, outputCanvases, false, legends, stacks, legendHeaders, dummyHists, 1, false);
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) { hists.push_back(NULL); }
}

// void setup(const vector<string>& canvasNames, vector<TCanvas*>& outputCanvases, 
// 	   vector<TH2F*>& hists)
// {
//   setup(canvasNames, outputCanvases, false, vector<TLegend*>(), vector<THStack*>(), 
// 	vector<string>(), vector<vector<TH1F*> >(), 1);
//   for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
//        iCanvasName != canvasNames.end(); ++iCanvasName) { hists.push_back(NULL); }
// }

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
    }
  }
}

void write(vector<TCanvas*>& outputCanvases)
{
  for (vector<TCanvas*>::iterator iOutputCanvas = outputCanvases.begin(); 
       iOutputCanvas != outputCanvases.end(); ++iOutputCanvas) {
    (*iOutputCanvas)->Write();
//     (*iOutputCanvas)->SaveAs((string((*iOutputCanvas)->GetName()) + "_nonIso.pdf").c_str());
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

// void scaleAndAdd(const vector<string>& canvasNames, TFile* in, const vector<string>& graphNames, 
// 		 const float weight, vector<TH2F*>& hists, const unsigned int fileIndex)
// {
//   for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
//        iCanvasName != canvasNames.end(); ++iCanvasName) {
//     const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
//     TCanvas* pCanvas;
//     in->GetObject(iCanvasName->c_str(), pCanvas);
//     TH2F* pHist = (TH2F*)pCanvas->GetPrimitive(graphNames[canvasIndex].c_str());
//     pHist->Scale(weight);
//     if (fileIndex == 0) hists[canvasIndex] = pHist;
//     else hists[canvasIndex]->Add(pHist);
//   }
// }

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
    return;
  }
  TFile outStream(outputFileName.c_str(), "RECREATE");
  vector<TFile*> inputStreams;
  vector<TCanvas*> outputCanvases;
  vector<TLegend*> legends;
  vector<THStack*> stacks;
  vector<vector<TH1F*> > hists;
  setup(canvasNames, outputCanvases, setLogY, legends, stacks, legendHeaders, hists, 
	inputFiles.size(), dataMC);
  bool data = true;
  for (vector<string>::const_iterator iInputFile = inputFiles.begin(); 
       iInputFile != inputFiles.end(); ++iInputFile) {
    const unsigned int fileIndex = iInputFile - inputFiles.begin();
    if ((fileIndex == 0) && (iInputFile->find("Wh1") != string::npos)) data = false;
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
	if (string(pHist->GetName()) == "muHadMass"/*"tauHadIso"*/) {
	  cout << *iInputFile << endl;
	  cout << "Integral: " << pHist->Integral(5, -1)/*pHist->Integral(0, -1)*/ << endl;
	}
	string histName(pHist->GetName());
	if (histName == "jet_pt_etacut") pHist->GetXaxis()->SetTitle("p_{T} (GeV)");
	if (histName == "jet_eta") pHist->GetXaxis()->SetTitle("#eta");
	if (histName == "jet_phi") pHist->GetXaxis()->SetTitle("#phi");
	if (histName == "jet_mass_etacut") pHist->GetXaxis()->SetTitle("m (GeV)");
	if (histName == "jet_ptmj_etacut") {
	  pHist->GetXaxis()->SetTitle("#frac{p_{T}}{m}");
	}
// 	if (string(pHist->GetName()) == "jet_ptmj_etacut") {
// 	  cout << "Integral(second jet pT/m > 13): " << pHist->Integral(14, -1) << endl;
// 	}
	if (setLogY) pHist->GetYaxis()->SetRangeUser(0.1, 10000.0);
	string legendStyle("l");
	if (drawStack) {
// 	  pHist->SetFillStyle(0);
// 	  pHist->SetFillColor(0);
	  pHist->SetFillStyle(1001);
	  if (((fileIndex == 0) || (fileIndex == 1)) && !data) pHist->SetFillColor(0);
	  else pHist->SetFillColor(colors[fileIndex]);
	  if (fileIndex == 0) {
	    if (data) legendStyle = "lp";
	    else legendStyle = "l";
	  }
	  else if ((fileIndex == 1) && !data) legendStyle = "l";
	  else legendStyle = "f";
	}
	hists[canvasIndex][fileIndex] = pHist;
	legends[canvasIndex]->
	  AddEntry(pHist, legendEntries[fileIndex].c_str(), legendStyle.c_str());
	if ((data && (fileIndex != 0)) || (!data && (fileIndex != 0) && (fileIndex != 1))) {
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
      for (vector<string>::const_iterator iInputFile = inputFiles.begin(); 
	   iInputFile != inputFiles.end(); ++iInputFile) {
	const unsigned int fileIndex = iInputFile - inputFiles.begin();
	if (hists[canvasIndex][fileIndex] != NULL) hists[canvasIndex][fileIndex]->Draw("HISTSAME");
      }
    }
    else {
      if (setLogY) {
	Double_t histMin = stacks[canvasIndex]->GetMinimum();
	Double_t axisMin = 1.0;
	int exponent = 0;
	if (histMin != 0.0) {
	  while (histMin < 1.0) {
	    histMin*=10;
	    --exponent;
	  }
	  while (histMin >= 10.0) {
	    histMin/=10;
	    ++exponent;
	  }
	}
	if (exponent < 0) {
	  exponent+=(2*exponent);
	  for (int i = 0; i < exponent; ++i) { axisMin/=10; }
	}
	else for (int i = 0; i < exponent; ++i) { axisMin*=10; }
// 	stacks[canvasIndex]->SetMinimum(/*axisMin == 0.0 ? */0.1/* : axisMin*/);
// 	stacks[canvasIndex]->SetMaximum(10000.0);
	stacks[canvasIndex]->SetMinimum(1.0);
	stacks[canvasIndex]->SetMaximum(10000000.0);
      }
      stacks[canvasIndex]->Draw();
      TList* stackedHists = stacks[canvasIndex]->GetHists();
      TH1F* hist = (TH1F*)stackedHists->First();
      stacks[canvasIndex]->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
      stacks[canvasIndex]->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
      string drawOpt("SAME");
      if (!data) drawOpt = "HISTSAME";
      hists[canvasIndex][0]->Draw(drawOpt.c_str());
      if (!data) hists[canvasIndex][1]->Draw(drawOpt.c_str());
      TH1F* stackSumHist = NULL;
      for (Int_t i = 0; i < stackedHists->GetEntries(); ++i) {
	TH1F* stackHist = (TH1F*)stackedHists->At(i)->Clone();
	stackHist->SetTitle("");
	if (i == 0) stackSumHist = stackHist;
	else stackSumHist->Add(stackHist);
      }
      TH1F* dataHist = (TH1F*)hists[canvasIndex][0]->Clone();
      stackSumHist->Add(dataHist, -1.0);
      stackSumHist->Divide(dataHist);
      setHistogramOptions(stackSumHist, kBlack, 0.7, 20, 1.0, 
			  stackSumHist->GetXaxis()->GetTitle(), 
			  "#frac{N_{MC} - N_{data}}{N_{data}}");
      stackSumHist->GetYaxis()->SetRangeUser(-1.0, 1.0);
      if (dataMC) {
	outputCanvases[canvasIndex]->cd(dataMC ? 2 : 0);
	stackSumHist->Draw();
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
  if ((inputFiles.size() != weights.size()) || (canvasNames1D.size() != graphNames1D.size()) || 
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
  setup(canvasNames2D, outputCanvases2D, hists2D);
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

//test data-driven background estimation method with MC on a given variable
void addClosurePlot(TFile& sigVsBkgIso20InvFbStream, const string& var, const string& unit, 
                    TFile& bkgNonIsoStream, const float nonIsoScale, 
                    const int normRegionLowerBin, const int normRegionUpperBin, TFile& outStream)
{
  //top level declarations
  string canvasName(var + "Canvas");
  string stackName(var + "Stack");

  //get plots of signal and background MC, isolated tau sample, 20 fb^-1 normalization
  TCanvas* canvasIso = NULL;
  sigVsBkgIso20InvFbStream.GetObject(canvasName.c_str(), canvasIso);

  //get the signal histogram
  TH1F* histSig = NULL;
  THStack* stackBkgIso = NULL;
  if (canvasIso != NULL) {
    histSig = (TH1F*)canvasIso->GetPrimitive(var.c_str());
    setHistogramOptions(histSig, kBlack, 0.7, 20, 1.0, unit.c_str(), "");
    histSig->GetYaxis()->SetRangeUser(0.1, 10000.0);

    /*get the background stack histogram in the signal region (isolated taus) and convert it to a 
      TH1*/
    stackBkgIso = (THStack*)canvasIso->GetPrimitive(stackName.c_str());
  }
  else {
    cerr << "Error opening canvas " << canvasName << " from file ";
    cerr << sigVsBkgIso20InvFbStream.GetName() << ".\n";
    return;
  }
  TList* stackedHistsIso = NULL;
  if (stackBkgIso != NULL) stackedHistsIso = stackBkgIso->GetHists();
  else {
    cerr << "Error opening stack " << stackName << " from canvas " << canvasName << " from file ";
    cerr << sigVsBkgIso20InvFbStream.GetName() << ".\n";
    return;
  }
  TH1F* histBkgIso = NULL;
  if (stackedHistsIso != NULL) {
    for (Int_t i = 0; i < stackedHistsIso->GetEntries(); ++i) {
      TH1F* stackHist = (TH1F*)stackedHistsIso->At(i)->Clone();
      if (i == 0) histBkgIso = stackHist;
      else histBkgIso->Add(stackHist);
    }
    setHistogramOptions(histBkgIso, kBlue, 0.7, 21, 1.0, unit.c_str(), "");
    histBkgIso->GetYaxis()->SetRangeUser(0.1, 10000.0);
  }
  else {
    cerr << "Error opening histogram list from stack " << stackName << " from canvas ";
    cerr << canvasName << " from file " << sigVsBkgIso20InvFbStream.GetName() << ".\n";
    return;
  }

  //get plots of background MC, non-isolated tau sample, 2.5 fb^-1 normalization
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
  if (stackedHistsNonIso != NULL) {
    for (Int_t i = 0; i < stackedHistsNonIso->GetEntries(); ++i) {
      TH1F* stackHist = (TH1F*)stackedHistsNonIso->At(i)->Clone();
      if (i == 0) histBkgNonIso = stackHist;
      else histBkgNonIso->Add(stackHist);
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
  TCanvas outCanvas(canvasName.c_str(), "", 600, 600);
  setCanvasOptions(outCanvas, 1, 1, 0);
  setCanvasMargins(outCanvas, 0.2, 0.2, 0.2, 0.2);
  outCanvas.cd();
  Double_t histBkgIsoErr = -1.0;
  Double_t histBkgNonIsoErr = -1.0;
  histSig->Draw("HIST");
  cout << histSig->Integral(5, -1) << endl;
  histBkgIso->Draw("HISTE");
  cout << histBkgIso->IntegralAndError(5, -1, histBkgIsoErr) << "+/-";
  cout << histBkgIsoErr << endl;
  histBkgNonIso->Draw("HISTESAME");
  cout << histBkgNonIso->IntegralAndError(5, -1, histBkgNonIsoErr) << "+/-";
  cout << histBkgNonIsoErr << endl;
  outCanvas.Write();
}

//test data-driven background estimation method with MC
void makeMCClosurePlots(const string& sigVsBkgIso20InvFbFileName, const vector<string>& vars, 
			const vector<string>& units, const string& bkgNonIsoFileName, 
			const float nonIsoScale, const vector<int>& normRegionLowerBins, 
			const vector<int>& normRegionUpperBins, const string& outputFileName)
{
  //open files
  TFile sigVsBkgIso20InvFbStream(sigVsBkgIso20InvFbFileName.c_str());
  TFile bkgNonIsoStream(bkgNonIsoFileName.c_str());
  TFile outStream(outputFileName.c_str(), "RECREATE");
  if (sigVsBkgIso20InvFbStream.IsOpen() && bkgNonIsoStream.IsOpen() && outStream.IsOpen()) {
    for (vector<string>::const_iterator iVar = vars.begin(); iVar != vars.end(); ++iVar) {
      const unsigned int varIndex = iVar - vars.begin();
      addClosurePlot(sigVsBkgIso20InvFbStream, *iVar, units[varIndex], bkgNonIsoStream, 
		     nonIsoScale, normRegionLowerBins[varIndex], normRegionUpperBins[varIndex], 
		     outStream);
    }
  }
  else {
    cerr << "Error opening files " << sigVsBkgIso20InvFbFileName << " or " << bkgNonIsoFileName;
    cerr << " or " << outputFileName << ".\n";
    return;
  }

  //write to file
  outStream.Write();
  outStream.Close();
  sigVsBkgIso20InvFbStream.Close();
  bkgNonIsoStream.Close();
}

//structs defining the draw options for each sample
struct drawOptions {
  string sampleName;
  unsigned int sampleNumber;
  Color_t markerColor;
  Style_t markerStyle;
  Color_t lineColor;
  Color_t fillColor;
  Style_t fillStyle;
} Wh1, gg, WW, ZZ, WZ, WNJets, singleTop, tt, DrellYan;

//set the draw options for each sample and store it all in global scope
void setDrawOptions(drawOptions& sampleDrawOptions, const string& sampleName, 
		    unsigned int sampleNumber, Color_t markerColor, Style_t markerStyle, 
		    Color_t lineColor, Color_t fillColor, Style_t fillStyle)
{
  sampleDrawOptions.sampleName = sampleName;
  sampleDrawOptions.sampleNumber = sampleNumber;
  sampleDrawOptions.markerColor = markerColor;
  sampleDrawOptions.markerStyle = markerStyle;
  sampleDrawOptions.lineColor = lineColor;
  sampleDrawOptions.fillColor = fillColor;
  sampleDrawOptions.fillStyle = fillStyle;
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
		  const string& option, TFile& outStream)
{
  //top level declarations
  string canvasName(var + "Canvas");
  string stackName(var + "Stack");
  string canvasReweightErrSqName(var + "ReweightErrSqCanvas");

  //get plots of signal and background MC and data-driven QCD estimate, isolated tau sample
  TCanvas* canvasIsoSigBkg = NULL;
  isoSigBkgFile.first->GetObject(canvasName.c_str(), canvasIsoSigBkg);

  //get the signal histograms
  vector<TH1F*> isoSig(2, NULL);
  THStack* isoBkg = NULL;
  TLegend legendBkgSep(0.35, 0.55, 0.75, 0.75);
  TLegend legendBkgMain5(0.35, 0.55, 0.75, 0.75);
  TLegend legendBkgAll(0.35, 0.55, 0.75, 0.75);
  setLegendOptions(legendBkgSep, "CMS 2.5 fb^{-1}");
  setLegendOptions(legendBkgMain5, "CMS 2.5 fb^{-1}");
  setLegendOptions(legendBkgAll, "CMS 2.5 fb^{-1}");
  if (canvasIsoSigBkg != NULL) {
    TList* sigs = canvasIsoSigBkg->GetListOfPrimitives();
    for (vector<TH1F*>::iterator iIsoSig = isoSig.begin(); iIsoSig != isoSig.end(); 
	 ++iIsoSig) {
      const unsigned int i = iIsoSig - isoSig.begin();
      *iIsoSig = (TH1F*)sigs->At(i + 2)->Clone();
      (*iIsoSig)->GetYaxis()->SetRangeUser(0.01, 10000.0);
    }
    legendBkgSep.AddEntry(isoSig[0], "Wh_{1}", "l");
    legendBkgSep.AddEntry(isoSig[1], "gg fusion", "l");
    legendBkgMain5.AddEntry(isoSig[0], "Wh_{1}", "l");
    legendBkgMain5.AddEntry(isoSig[1], "gg fusion", "l");
    legendBkgAll.AddEntry(isoSig[0], "Wh_{1}", "l");
    legendBkgAll.AddEntry(isoSig[1], "gg fusion", "l");
    setHistogramOptions(isoSig[0], kSpring - 1, 0.7, 20, isoSigBkgFile.second, unit.c_str(), "");
    setHistogramOptions(isoSig[1], kAzure + 1, 0.7, 20, isoSigBkgFile.second, unit.c_str(), "");

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
      if (i == 0) {
	isoBkgAllHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgAllHist->Scale(isoSigBkgFile.second);
	isoBkgDibosonHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgDibosonHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "WW (MC)", "f");
      }
      else isoBkgAllHist->Add(stackHist);
      if (i == 1) legendBkgSep.AddEntry(stackHist, "ZZ (MC)", "f");
      if (i == 2) legendBkgSep.AddEntry(stackHist, "WZ (MC)", "f");
      if ((i == 1) || (i == 2)) isoBkgDibosonHist->Add(stackHist);
      if (i == 3) {
	isoBkgWNJetsHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgWNJetsHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "W + #geq1 jet (MC)", "f");
      }
      if (i == 4) {
	isoBkgTopHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgTopHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "t/#bar{t} (MC)", "f");
      }
      if (i == 5) {
	isoBkgTopHist->Add(stackHist);
	legendBkgSep.AddEntry(stackHist, "t#bar{t} + jets (MC)", "f");
      }
      if (i == 6) {
	isoBkgDrellYanHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgDrellYanHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "Drell-Yan + jets (MC)", "f");
      }
      if (i == 7) {
	isoBkgQCDHist = (TH1F*)isoBkgHists->At(i)->Clone();
	isoBkgQCDHist->Scale(isoSigBkgFile.second);
	legendBkgSep.AddEntry(stackHist, "QCD (data)", "f");
      }
    }
//     setHistogramOptions(isoBkgAllHist, kBlue, 0.7, 21, 1.0, unit.c_str(), "");
    isoBkgAllHist->GetYaxis()->SetRangeUser(0.01, 10000.0);
    isoBkgAll.Add(isoBkgAllHist, "HIST");
    isoBkgMain5.Add(isoBkgDibosonHist, "HIST");
    isoBkgMain5.Add(isoBkgWNJetsHist, "HIST");
    isoBkgMain5.Add(isoBkgTopHist, "HIST");
    isoBkgMain5.Add(isoBkgDrellYanHist, "HIST");
//     isoBkgMain5.Add(isoBkgQCDHist, "HIST");
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
      (TH1F*)canvasNonIsoDataReweightErrSq->cd(1)->GetPrimitive(var.c_str())->Clone();
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
    legendBkgSep.AddEntry(isoData, "Data", "p");
    legendBkgMain5.AddEntry(isoData, "Data", "p");
    legendBkgAll.AddEntry(isoData, "Data", "p");
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

  //debug
  cerr << "Norm bin 1\n";
  cerr << "----------\n";
  cerr << "NIso = " << isoData->GetBinContent(1) << ", errIso = " << isoData->GetBinError(1);
  cerr << "\n";
  cerr << "NNonIso = " << nonIsoData->GetBinContent(1) << ", errNonIso = ";
  cerr << nonIsoData->GetBinError(1) << "\n";
  cerr << "reweight error squared = ";
  cerr << nonIsoDataReweightErrSq->GetBinError(1)*nonIsoDataReweightErrSq->GetBinError(1) << endl;
  cerr << endl;
  cerr << "Norm bin 2\n";
  cerr << "----------\n";
  cerr << "NIso = " << isoData->GetBinContent(2) << ", errIso = " << isoData->GetBinError(2);
  cerr << "\n";
  cerr << "NNonIso = " << nonIsoData->GetBinContent(2) << ", errNonIso = ";
  cerr << nonIsoData->GetBinError(2) << "\n";
  cerr << "reweight error squared = ";
  cerr << nonIsoDataReweightErrSq->GetBinError(2)*nonIsoDataReweightErrSq->GetBinError(2) << endl;
  cerr << endl;
  cerr << "Search bin 1\n";
  cerr << "----------\n";
  cerr << "NIso = " << isoData->GetBinContent(4) << ", errIso = " << isoData->GetBinError(4);
  cerr << "\n";
  cerr << "NNonIso = " << nonIsoData->GetBinContent(4) << ", errNonIso = ";
  cerr << nonIsoData->GetBinError(4) << "\n";
  cerr << "reweight error squared = ";
  cerr << nonIsoDataReweightErrSq->GetBinError(4)*nonIsoDataReweightErrSq->GetBinError(4) << endl;
  cerr << endl;
  cerr << "norm = " << norm << endl << endl;


  /*calculate the statistical error on the background prediction from the non-isolated data, 
    including the term from the error on the normalization factor*/
  const TH1* nonIsoDataPtrCast = dynamic_cast<const TH1*>(nonIsoData);
  const TH1* nonIsoDataReweightErrSqPtrCast = dynamic_cast<const TH1*>(nonIsoDataReweightErrSq);
  vector<Double_t> nonIsoDataStatErrSq;
  for (Int_t iBin = 1; iBin <= nonIsoData->GetNbinsX(); ++iBin) {
//     nonIsoDataStatErrSq.
//       push_back(bkgErrSqFromNorm(nonIsoDataPtrCast->GetBinContent(iBin), normErrSq) + 
// 		bkgErrSqFromStats(norm, nonIsoDataPtrCast->GetBinError(iBin)) + 
// 		bkgErrSqFromReweight(norm, reweightErrSq));
    nonIsoDataStatErrSq.
      push_back(bkgErrSq(nonIsoDataPtrCast, nonIsoDataReweightErrSqPtrCast, iBin, norm, 
			 normErrSq(dynamic_cast<const TH1*>(isoData), nonIsoDataPtrCast, 
				   nonIsoDataReweightErrSqPtrCast, normRegionLowerBin, 
				   normRegionUpperBin, norm)));
//434401.4766 + 1.200894736 + 0.021129675
  }

  //debug
  cerr << endl << "nonIsoDataStatErrSq[3] = " << nonIsoDataStatErrSq[3] << endl;

  //normalize non-isolated data histogram to isolated data in signal-depleted region
  nonIsoData->Scale(norm);

  //debug
  cerr << "Search bin 1 NNonIso = " << nonIsoData->GetBinContent(4) << endl;

  //set statistical error in each bin of the non-isolated data histogram
  for (Int_t iBin = 1; iBin <= nonIsoData->GetNbinsX(); ++iBin) {
    nonIsoData->SetBinError(iBin, sqrt(nonIsoDataStatErrSq[iBin - 1]));
  }

  //debug
  cerr << "Search bin 1 errNonIso = " << nonIsoData->GetBinError(4) << endl;

  //write to file
  outStream.cd();
  TCanvas outCanvas(canvasName.c_str(), "", 600, 600);
  setCanvasOptions(outCanvas, 1, 1, 0);
  setCanvasMargins(outCanvas, 0.2, 0.2, 0.2, 0.2);
  outCanvas.cd();
  if (option == "separate") {
    isoBkgSep.Draw();
    isoBkgSep.SetMinimum(0.01);
    isoBkgSep.SetMaximum(10000.0);
  }
  else if (option == "main 5") {
    isoBkgMain5.Draw();
    isoBkgMain5.SetMinimum(0.01);
    isoBkgMain5.SetMaximum(10000.0);
  }
  else if (option == "combined") {
    isoBkgAll.Draw();
    isoBkgAll.GetHistogram()->GetYaxis()->SetRangeUser(0.01, 10000.0);
  }
  nonIsoData->Draw("HISTESAME");
  for (vector<TH1F*>::iterator iIsoSig = isoSig.begin(); iIsoSig != isoSig.end(); 
       ++iIsoSig) {
    (*iIsoSig)->Draw("HISTSAME");
  }
  isoData->Draw("ESAME");
  if (option == "separate") legendBkgSep.Draw();
  else if (option == "main 5") legendBkgMain5.Draw();
  else if (option == "combined") legendBkgAll.Draw();
  outCanvas.Write();
}

//create a file of properly formatted final plots
void makeFinalPlot(const pair<string, float>& isoMC, const string& isoDataFileName, 
		   const pair<string, float>& nonIsoData, const vector<string>& vars, 
		   const vector<string>& units, const vector<int>& normRegionLowerBins, 
		   const vector<int>& normRegionUpperBins, const string& outputFileName, 
		   const string& option// , const vector<Color_t>& colors, 
// 		   const vector<Style_t>& styles
		   )
{
//   //set the draw options for each sample and store it all in global scope
//   setDrawOptions(Wh1, "Wh_{1}", 2, kBlack, 20, kBlack, 0, 0);
//   setDrawOptions(gg, "Wh_{1}", 2, kBlack, 20, kBlack, 0, 0);
//   setDrawOptions(WW, "Wh_{1}", 2, kBlack, 20, kBlack, 0, 0);
//   setDrawOptions(ZZ, "Wh_{1}", 2, kBlack, 20, kBlack, 0, 0);
//   setDrawOptions(WZ, "Wh_{1}", 2, kBlack, 20, kBlack, 0, 0);
//   setDrawOptions(WNJets, "Wh_{1}", 2, kBlack, 20, kBlack, 0, 0);
//   setDrawOptions(singleTop, "Wh_{1}", 2, kBlack, 20, kBlack, 0, 0);
//   setDrawOptions(tt, "Wh_{1}", 2, kBlack, 20, kBlack, 0, 0);
//   setDrawOptions(DrellYan, "Wh_{1}", 2, kBlack, 20, kBlack, 0, 0);

  //open files
  pair<TFile*, float> isoSigBkgFile(new TFile(isoMC.first.c_str()), isoMC.second);
  TFile isoDataFile(isoDataFileName.c_str());
  pair<TFile*, float> nonIsoDataFile(new TFile(nonIsoData.first.c_str()), nonIsoData.second);
  TFile outStream(outputFileName.c_str(), "RECREATE");
  if (isoSigBkgFile.first->IsOpen() && isoDataFile.IsOpen() && nonIsoDataFile.first->IsOpen() && 
      outStream.IsOpen()) {

    //add plots
    for (vector<string>::const_iterator iVar = vars.begin(); iVar != vars.end(); ++iVar) {
      const unsigned int varIndex = iVar - vars.begin();
      addFinalPlot(isoSigBkgFile, isoDataFile, nonIsoDataFile, *iVar, units[varIndex], 
		   normRegionLowerBins[varIndex], normRegionUpperBins[varIndex], option, 
		   outStream);
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
