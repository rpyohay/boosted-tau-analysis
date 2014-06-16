#include <string>
#include <map>
#include <iostream>
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
	cout << "    mistagEffVsPTAndEta_->SetBinContent(" << iPTBin << ", " << iAbsEtaBin << ", " << effHist->GetBinContent(iPTBin, iAbsEtaBin) << ");\n";
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
	      const string& savePath)
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
  cout << inputFileName << endl;

  //open output file
  TFile out(outputFileName.c_str(), "RECREATE");
  if (!out.IsOpen()) {
    cerr << errorCannotOpenFile(fnName, outputFileName);
    in.Close();
    return;
  }

  //make efficiency plots
  plotEfficiency(in, out, effHistMap1D, binLabelMap, savePath);
  plot2DEfficiency(in, out, effHistMap2D, savePath);

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
    stackName.replace(stackName.find("Canvas"), 6, "Stack");
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
	  cout << "m > 4 GeV: " << pHist->Integral(5, -1) << endl;
	}
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
      stacks[canvasIndex]->Draw("e");
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
	stackSumHist->Draw("e");
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
    //histDiff.push_back(<TH1F*>);
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
		  histDiff[canvasIndex]->Add(hists[canvasIndex][fileIndex], -1.); // subtract from data
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

//get Region A QCD histograms
void drawQCDRegionAHistograms(const string& outputFileA, 
			      const string& inputFileNameB, 
			      const string& inputFileNameC, 
			      const string& inputFileNameD, 
			      const vector<string>& canvasNames, 
			      const vector<string>& graphNames, 
			      const vector<string>& legendHeaders, 
			      const vector<Color_t>& colors, 
			      const vector<Style_t>& styles, 
			      const vector<string>& legendEntries, 
			      const vector<float>& weights, const bool setLogY, 
			      const bool dataMC)
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
  TFile inputFileC(inputFileNameC.c_str(), "read");
  TFile inputFileD(inputFileNameD.c_str(), "read");
  vector<TCanvas*> outputCanvases;
  vector<TLegend*> legends;
  vector<THStack*> stacks;
  vector<vector<TH1F*> > hists;
  setup(canvasNames, outputCanvases, setLogY, legends, stacks, legendHeaders, hists, 3, dataMC, false);
  
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) { // loop over canvases
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    TCanvas* pCanvasB;
    TCanvas* pCanvasC;
    TCanvas* pCanvasD;
    inputFileB.GetObject(iCanvasName->c_str(), pCanvasB);
    inputFileC.GetObject(iCanvasName->c_str(), pCanvasC);
    inputFileD.GetObject(iCanvasName->c_str(), pCanvasD);
    TH1F* pHistB = NULL;
    TH1F* pHistC = NULL;
    TH1F* pHistD = NULL;
    //string canvasNameB(iCanvasName->c_str());
    //canvasNameB+="_1";
    //TCanvas* pCanvasB1 = (TCanvas*)pCanvasB->GetPrimitive(canvasNameB.c_str());
    pHistB = (TH1F*)pCanvasB->GetPrimitive(graphNames[canvasIndex].c_str())->Clone();
    pHistC = (TH1F*)pCanvasC->GetPrimitive(graphNames[canvasIndex].c_str())->Clone();
    pHistD = (TH1F*)pCanvasD->GetPrimitive(graphNames[canvasIndex].c_str())->Clone();
    pHistB->Scale(pHistC->Integral(0, -1)/pHistD->Integral(0, -1));
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
    outputCanvases[canvasIndex]->cd(dataMC ? 1 : 0);
    //outputCanvases[canvasIndex]->cd(0);
    pHistB->Draw();
  } // loop over canvases

  outStream.cd();
  write(outputCanvases);
  outStream.Write();
  outStream.Close();
  inputFileB.Close();
  inputFileC.Close();
  inputFileD.Close();
  deleteObjects(legends);
  deleteObjects(stacks);
  deleteObjects(outputCanvases);

} // end routine

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
  TH1F* histSig2 = NULL;
  THStack* stackBkgIso = NULL;
  if (canvasIso != NULL) {
    TList* sigs = canvasIso->GetListOfPrimitives();
    histSig1 = (TH1F*)sigs->At(2)->Clone();
    histSig2 = (TH1F*)sigs->At(3)->Clone();
    histSig1->GetYaxis()->SetRangeUser(0.1, 10000.0);
    histSig2->GetYaxis()->SetRangeUser(0.1, 10000.0);

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
  TCanvas outCanvas(canvasName.c_str(), "", 600, 900);
  outCanvas.Divide(1, 2);
  outCanvas.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  setCanvasOptions(*outCanvas.cd(1), 1, 1, 0);
  outCanvas.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  setCanvasOptions(*outCanvas.cd(2), 1, 0, 0);
  outCanvas.cd(1);
//   Double_t histBkgIsoErr = -1.0;
//   Double_t histBkgNonIsoErr = -1.0;
//   histSig1->Draw("HISTE");
//   cout << histSig1->Integral(5, -1) << endl;
//   histSig2->Draw("HISTESAME");
  histBkgIso->Draw("HISTE");
//   cout << histBkgIso->IntegralAndError(5, -1, histBkgIsoErr) << "+/-";
//   cout << histBkgIsoErr << endl;
  histBkgNonIso->Draw("HISTESAME");
//   cout << histBkgNonIso->IntegralAndError(5, -1, histBkgNonIsoErr) << "+/-";
//   cout << histBkgNonIsoErr << endl;
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
			 const string& units, 
			 const int normRegionLowerBin, 
			 const int normRegionUpperBin, const string& outputFileName)
{

  TFile outStream(outputFileName.c_str(), "RECREATE");
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
    cout << "QCDVSMC pCanvas: " << pCanvas << endl;
    TH1F* pHist = (TH1F*)pCanvas->GetPrimitive(var.c_str());
    cout << "QCDVsMC pHist: " << pHist << endl;
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
  //  hists[0]->Scale(hists[1]->Integral(normRegionLowerBin,normRegionUpperBin)/hists[0]->Integral(normRegionLowerBin,normRegionUpperBin));
  hists[0]->Scale(1./hists[0]->Integral());
  hists[1]->Scale(1./hists[1]->Integral());
  hists[1]->GetXaxis()->SetTitle(units.c_str());

  TLegend *leg = new TLegend(0.35, 0.55, 0.75, 0.75, "");
  leg->AddEntry(hists[0], "QCD (data-driven)", "lp");
  leg->AddEntry(hists[1], "EWK+TOP+DY (MC)", "lp");

  //write to file
  outStream.cd();
  TCanvas outCanvas(canvasName.c_str(), "", 600, 600);
  hists[0]->Draw("HISTE");
  hists[1]->Draw("HISTESAME");
  leg->Draw();
  outCanvas.Write();
  outStream.Write();
  outStream.Close();
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
  string varReweightErrSq(var + "ReweightErrSq");

  //get plots of signal and background MC and data-driven QCD estimate, isolated tau sample
  TCanvas* canvasIsoSigBkg = NULL;
  isoSigBkgFile.first->GetObject(canvasName.c_str(), canvasIsoSigBkg);

  //get the signal histograms
  vector<TH1F*> isoSig(2, NULL);
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
      for (Int_t iBin = 5; iBin <= (hist->GetNbinsX() + 1); ++iBin) {
	sum+=hist->GetBinContent(iBin);
	err+=(hist->GetBinError(iBin)*hist->GetBinError(iBin));
      }
    }
    cout << "Region A MC + data-driven QCD, m > 4: " << sum << " +/- " << sqrt(err) <<  endl;
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
  Double_t statErr;
  nonIsoData->IntegralAndError(5, -1, statErr);
  cout << "Region B data, m > 4: " << nonIsoData->Integral(5, -1) << " +/- " << statErr << endl;
  for (vector<TH1F*>::iterator iIsoSig = isoSig.begin(); iIsoSig != isoSig.end(); 
       ++iIsoSig) {
    (*iIsoSig)->Draw("HISTSAME");
    cout << "Region A signal " << iIsoSig - isoSig.begin() << ", m > 4: ";
    cout << (*iIsoSig)->Integral(5, -1) << endl;
  }
  isoData->Draw("ESAME");
  cout << "Region A data, m < 2: " << isoData->Integral(1, 2) << endl;
  if (option == "separate") legendBkgSep.Draw();
  else if (option == "main 5") legendBkgMain5.Draw();
  else if (option == "combined") legendBkgAll.Draw();
  outCanvas.cd(2);
  TH1F* nonIsoDataMinusIsoBkgAll = (TH1F*)nonIsoData->Clone();
  nonIsoDataMinusIsoBkgAll->Add(isoBkgAllHist, -1.0);
  nonIsoDataMinusIsoBkgAll->Divide(nonIsoData);
  nonIsoDataMinusIsoBkgAll->GetYaxis()->SetTitle("#frac{Data (B) - MC (A)}{Data (B)}");
  nonIsoDataMinusIsoBkgAll->GetYaxis()->SetRangeUser(-1.0, 1.0);
  nonIsoDataMinusIsoBkgAll->Draw();
  TH1F* nonIsoDataMinusIsoData = (TH1F*)nonIsoData->Clone();
  nonIsoDataMinusIsoData->Add(isoData, -1.0);
  nonIsoDataMinusIsoData->Divide(nonIsoData);
  nonIsoDataMinusIsoData->GetYaxis()->SetTitle("#frac{Data (B) - Data (A)}{Data (B)}");
  nonIsoDataMinusIsoData->GetYaxis()->SetRangeUser(-1.0, 1.0);
  nonIsoDataMinusIsoData->SetMarkerColor(kBlack);
  nonIsoDataMinusIsoData->SetLineColor(kBlack);
//   nonIsoDataMinusIsoData->Draw("SAME");
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
    obj = (T*)canvas->cd(pad)->GetPrimitive(objName.c_str())->Clone();
  }
  return obj;
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
//   string selNames[8] = {"Baseline", "+ jet", "+#slash{E}_{T}", "+ jet + #slash{E}_{T}", "+M_{T}", "+M_{T} + jet", "+ jet (40 GeV), "+ jet (20 GeV)"};
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
