void analyze3MuonEvents(const string& version)
{
  gROOT->Reset();

  vector<string> MTBins;
  MTBins.push_back("_lowMT");
  MTBins.push_back("_highMT");

  vector<string> tauIso;
  tauIso.push_back("Iso");
  tauIso.push_back("NonIso");

  vector<string> muIso;
  muIso.push_back("data");
  muIso.push_back("nonIsoWData");

  vector<string> methods;
  methods.push_back("ShareTrack");
  methods.push_back("SoftMu");
  methods.push_back("SoftMu5GeV");
  methods.push_back("SoftMu15GeV");
  methods.push_back("SoftMu20GeV");

  for (vector<string>::const_iterator iMTBin = MTBins.begin(); iMTBin != MTBins.end(); ++iMTBin) {

    for (vector<string>::const_iterator iTauIso = tauIso.begin(); iTauIso != tauIso.end(); 
	 ++iTauIso) {

      for (vector<string>::const_iterator iMuIso = muIso.begin(); iMuIso != muIso.end(); 
	   ++iMuIso) {

	string prefix("");
	if (*iMuIso == "nonIsoWData") prefix = "nonIsoW_";
	TFile file(("/data1/yohay/" + *iMuIso + "/analysis/" + prefix + "muHad" + *iTauIso + 
		    "Analysis" + *iMTBin + "_SingleMu_" + version + ".root").c_str());

	for (vector<string>::const_iterator iMethod = methods.begin(); iMethod != methods.end(); 
	     ++iMethod) {

	  TCanvas* canvas = NULL;
	  file.GetObject(("muHadMass3Mu" + *iMethod + "Canvas").c_str(), canvas);
	  if (canvas != NULL) {
	    TH1F* hist = (TH1F*)canvas->GetPrimitive(("muHadMass3Mu" + *iMethod).c_str());
	    if (hist != NULL) {
	      cout << *iMTBin << " " << *iTauIso << " " << *iMuIso << " " << *iMethod << " ";
	      cout << hist->Integral(0, 4);
	      if ((*iTauIso == "Iso") && (*iMuIso == "data")) {
		cout << endl;
	      }
	      else cout << " (" << hist->Integral(5, -1) << ")\n";
	    }
	  }
	}

	file.Close();
      }
    }
  }
}
