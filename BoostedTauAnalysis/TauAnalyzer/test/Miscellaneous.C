#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"

void compareVersions(const vector<vector<string> >& versions)
{
  const string histogramName("muHadMass");
  const string canvasName(histogramName + "Canvas");
  vector<vector<vector<Float_t> > > muHadMass;
  for (vector<vector<string> >::const_iterator iVersion = versions.begin(); 
       iVersion != versions.end(); ++iVersion) {
    vector<vector<Float_t> > muHadMassThisVersion;
    for (vector<string>::const_iterator iFile = iVersion->begin(); iFile != iVersion->end(); 
	 ++iFile) {
      vector<Float_t> muHadMassThisFile;
      TFile file(iFile->c_str());
      if (file.IsOpen()) {
	TCanvas* canvas = NULL;
	file.GetObject(canvasName.c_str(), canvas);
	if (canvas != NULL) {
	  TH1F* histogram = NULL;
	  histogram = (TH1F*)canvas->GetPrimitive(histogramName.c_str());
	  if (histogram != NULL) {
	    for (Int_t iBin = 0; iBin <= (histogram->GetNbinsX() + 1); ++iBin) {
	      muHadMassThisFile.push_back(histogram->GetBinContent(iBin));
	    }
	  }
	  else cerr << "Null histogram pointer\n";
	}
	else cerr << "Null canvas pointer\n";
	file.Close();
      }
      else cerr << "File not open\n";
      muHadMassThisVersion.push_back(muHadMassThisFile);
    }
    muHadMass.push_back(muHadMassThisVersion);
  }
  cout << "v86 v87 % difference\n";
  for (unsigned int iFile = 0; iFile < versions[0].size(); ++iFile) {
    cout << versions[0][iFile] << endl;
    for (unsigned int iBin = 0; iBin < muHadMass[0][0].size(); ++iBin) {
      if (muHadMass[0][iFile][iBin] != muHadMass[1][iFile][iBin]) {
	cout << "Bin " << iBin << ": " << muHadMass[0][iFile][iBin] << " ";
	cout << muHadMass[1][iFile][iBin] << " ";
	cout << ((muHadMass[0][iFile][iBin] - muHadMass[1][iFile][iBin])/
		 muHadMass[0][iFile][iBin])*100 << endl;
      }
    }
  }
}
