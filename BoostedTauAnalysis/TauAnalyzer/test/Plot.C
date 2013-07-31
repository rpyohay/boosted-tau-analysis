#include <string>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TProfile.h"

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
void setCanvasOptions(TCanvas& canvas, const Int_t grid, const Int_t logY, const Int_t logZ)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
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

    //rebin
    cerr << "Rebinning\n";
    numeratorHist->Rebin(2);
    denominatorHist->Rebin(2);

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
	   const unsigned int size)
{
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    outputCanvases.push_back(new TCanvas(iCanvasName->c_str(), "", 600, 600));
    setCanvasOptions(*outputCanvases[outputCanvases.size() - 1], 1, setLogY, 0);
    setCanvasMargins(*outputCanvases[outputCanvases.size() - 1], 0.2, 0.2, 0.2, 0.2);
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
  setup(canvasNames, outputCanvases, false, legends, stacks, legendHeaders, dummyHists, 1);
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
		 const float weight, vector<T*>& hists, const unsigned int fileIndex)
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
      hists[canvasIndex]->GetZaxis()->SetRangeUser(0.00001, 1.0);
      hists[canvasIndex]->Draw("COLZ");
      TProfile* profileX = 
	hists[canvasIndex]->ProfileX((string(hists[canvasIndex]->GetName()) + "_pfx").c_str());
      profileX->SetLineWidth(3);
      profileX->Draw("HISTSAME");
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
					   const bool drawStack)
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
	inputFiles.size());
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
			    weight, pHist->GetXaxis()->GetTitle(), 
			    pHist->GetYaxis()->GetTitle());
	if (setLogY) pHist->GetYaxis()->SetRangeUser(0.1, 10000.0);
	if (drawStack) {
// 	  pHist->SetFillStyle(0);
// 	  pHist->SetFillColor(0);
	  pHist->SetFillStyle(1001);
	  pHist->SetFillColor(colors[fileIndex]);
	}
	hists[canvasIndex][fileIndex] = pHist;
	legends[canvasIndex]->AddEntry(pHist, legendEntries[fileIndex].c_str(), "l");
// 	if (fileIndex == (inputFiles.size() - 1)) stacks[canvasIndex]->Add(pHist, "HIST");
// 	else stacks[canvasIndex]->Add(pHist, "HISTE");
	if (fileIndex != 0) stacks[canvasIndex]->Add(pHist, "HIST");
      }
      outStream.cd();
      outputCanvases[canvasIndex]->cd();
      if (fileIndex == 0) { if (pGraph != NULL) pGraph->Draw("AP"); }
      else if (pGraph != NULL) pGraph->Draw("PSAME");
    }
  }
  if (!drawStack) {
    for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
	 iCanvasName != canvasNames.end(); ++iCanvasName) {
      const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
      outStream.cd();
      outputCanvases[canvasIndex]->cd();
      TH1F* pHistWithMaxMaxBin = histWithLargestBinContent(hists[canvasIndex]);
      if (pHistWithMaxMaxBin != NULL) pHistWithMaxMaxBin->Draw();
    }
  }
  for (vector<string>::const_iterator iCanvasName = canvasNames.begin(); 
       iCanvasName != canvasNames.end(); ++iCanvasName) {
    const unsigned int canvasIndex = iCanvasName - canvasNames.begin();
    outStream.cd();
    outputCanvases[canvasIndex]->cd();
    if (!drawStack) {
      for (vector<string>::const_iterator iInputFile = inputFiles.begin(); 
	   iInputFile != inputFiles.end(); ++iInputFile) {
	const unsigned int fileIndex = iInputFile - inputFiles.begin();
	if (hists[canvasIndex][fileIndex] != NULL) hists[canvasIndex][fileIndex]->Draw("SAME");
      }
    }
    else {
      if (setLogY) {
	stacks[canvasIndex]->SetMinimum(0.1);
	stacks[canvasIndex]->SetMaximum(10000.0);
      }
      stacks[canvasIndex]->Draw();
      TH1F* hist = (TH1F*)stacks[canvasIndex]->GetHists()->First();
      stacks[canvasIndex]->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
      stacks[canvasIndex]->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
      hists[canvasIndex][0]->Draw("SAME");
    }
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

//hadd histograms drawn on canvases
void haddCanvases(const string& outputFileName, const vector<string>& inputFiles, 
		  const vector<float>& weights, const vector<string>& canvasNames1D, 
		  const vector<string>& graphNames1D, const vector<string>& canvasNames2D, 
		  const vector<string>& graphNames2D)
{
  if ((inputFiles.size() != weights.size()) || (canvasNames1D.size() != graphNames1D.size()) || 
      (canvasNames2D.size() != graphNames2D.size())) {
    cerr << "Error: vector size mismatch.\n";
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
  for (vector<string>::const_iterator iInputFile = inputFiles.begin(); 
       iInputFile != inputFiles.end(); ++iInputFile) {
    const unsigned int fileIndex = iInputFile - inputFiles.begin();
    inputStreams.push_back(new TFile(iInputFile->c_str()));
    scaleAndAdd(canvasNames1D, inputStreams[inputStreams.size() - 1], graphNames1D, 
		weights[fileIndex], hists1D, fileIndex);
    scaleAndAdd(canvasNames2D, inputStreams[inputStreams.size() - 1], graphNames2D, 
		weights[fileIndex], hists2D, fileIndex);
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
