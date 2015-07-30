void savePlots(const string& outputVersion, const string& a1Mass, const string& MTBin, 
	       const bool narrowBins, const string& HLTPath)
{
  //initial
  gROOT->Reset();

  //get CMSSW path
  const char* CMSSWPathCString = gSystem->Getenv("CMSSW_BASE");
  if (!CMSSWPathCString) {
    CMSSWPathCString = "";
    cout << "Error: environment variable CMSSW_BASE is not set.  ";
    cout << "Please run cmsenv from within your CMSSW project area.\n";
  }
  string CMSSWPathCPPString(CMSSWPathCString);

  //load
  string macroPath(CMSSWPathCPPString + "/src/BoostedTauAnalysis/TauAnalyzer/test/");
  gROOT->ProcessLine("#include <utility>");
  gSystem->Load((macroPath + "STLDictionary.so").c_str());
  gSystem->Load((macroPath + "Miscellaneous_C.so").c_str());
  gSystem->Load((macroPath + "Error_C.so").c_str());
  gSystem->Load((macroPath + "Plot_C.so").c_str());

  //ignore warnings
  gErrorIgnoreLevel = kError;

  //needed so vector<Color_t> and vector<Style_t> work
  vector<short> dummy;

  //space-saving constant definitions
  string user(gSystem->GetFromPipe("whoami").Data());
  string userDir(1, user[0]);
  const string analysisFilePath("/data1/" + user + "/");
  const string fileExt(".root");
  const string tag19p7InvFb("_19p7fb-1");
  const string tag1("_normalizedTo1");

  //version tag
  const string outputVTag("_" + outputVersion);

  //for error bands
  gStyle->SetErrorX(0.5);

//   const string FranTestDir("/myNMSSMAnalysis/CMSSW_5_3_11/src/BoostedTauAnalysis/TauAnalyzer/test/plotsToSave/");
  //save directory
  const string saveDir("/afs/cern.ch/user/" + userDir + "/" + user + 
		       "/AN-13-254/notes/AN-13-254/trunk/figures/");

  // //HLT path
  // string HLTPath("_IsoMu24_eta2p1");
  // if (outputVersion == "v199") HLTPath = "_Mu40_eta2p1";

  if (/*(outputVersion != "v199") && (outputVersion != "v194")*/!narrowBins) {

    const string sigVsBkgFileName(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + 
				  MTBin + a1Mass + tag19p7InvFb + outputVTag + fileExt);
    TFile sigVsBkgFile(sigVsBkgFileName.c_str());

    //AN-13-254 Sec. 3.5 Fig. 3
    //AN-13-254 Sec. 3.5 Fig. 4
    if (outputVersion == "v202") {
      saveCanvasAsPDF("muHadMassCanvas", "muHadMass" + MTBin + "_beforeOSSF", saveDir, 
		      sigVsBkgFile);
    }

    //AN-13-254 Sec. 3.6 Fig. 5
    if (outputVersion == "v203") {
      saveCanvasAsPDF("muHadChargeCanvas", "muHadCharge" + MTBin + "_beforeSSSF", saveDir, 
		      sigVsBkgFile);
    }

    //AN-13-254 Sec. 3.7 Fig. 6
    if (outputVersion == "v204") {
      saveCanvasAsPDF("bTagDiscrimCanvas", "sigVsBkg_csv_regA" + MTBin + "_v61", saveDir, 
		      sigVsBkgFile);
    }

    //AN-13-254 Sec. 3.10 Fig. 9
    if (outputVersion == "v205") {
      saveCanvasAsPDF("WMuMTCanvas", "sigVsBkg_MT_regA_v62", saveDir, sigVsBkgFile);
    }

    //AN-13-254 Sec. 3.5 Fig. 3
    //AN-13-254 Sec. 3.5 Fig. 4
    //AN-13-254 Sec. 3.9 Fig. 7
    //AN-13-254 Sec. 3.9 Fig. 8
    if (outputVersion == "v206") {
      saveCanvasAsPDF("muHadMassCanvas", "muHadMass" + MTBin + "_afterOSSF", saveDir, 
		      sigVsBkgFile);
      saveCanvasAsPDF("muPVdzCanvas", "sigVsBkg_dztaumu_regA" + MTBin + "_v60", saveDir, 
		      sigVsBkgFile);
      saveCanvasAsPDF("hadPVdzCanvas", "sigVsBkg_dztauhad_regA" + MTBin + "_v60", saveDir, 
		      sigVsBkgFile);
    }

    sigVsBkgFile.Close();

    // //AN-13-254 Sec. 5.5 Fig. 46
    // //AN-13-254 Sec. 5.8 Fig. 53
    // const string finalPlotsFileName(analysisFilePath + "results/final" + MTBin + a1Mass + 
    // 				    outputVTag + fileExt);
    // TFile finalPlotsFile(finalPlotsFileName.c_str());
    // saveCanvasAsPDF("jetFakeBkgSystCanvas", "jetFakeBkgSystCanvas" + MTBin, saveDir, 
    // 		    finalPlotsFile);
    // saveCanvasAsPDF("muHadMassCanvas;1", "muHadMassCanvas_final_a9" + MTBin, saveDir, 
    // 		    finalPlotsFile);
    // finalPlotsFile.Close();

    // //AN-13-254 Sec. 5.6 Fig. 47
    // const string regBDataMinusMCVsRegDDataFileName(analysisFilePath + 
    // 						   "results/regBDataMinusMCVsRegDData" + MTBin + 
    // 						   tag19p7InvFb + outputVTag + fileExt);
    // TFile regBDataMinusMCVsRegDDataFile(regBDataMinusMCVsRegDDataFileName.c_str());
    // saveCanvasAsPDF("muHadMassCanvas", "muHadMassCanvas_regBDataMinusMCVsRegDData" + MTBin, 
    // 		    saveDir, regBDataMinusMCVsRegDDataFile);
    // regBDataMinusMCVsRegDDataFile.Close();

    // //AN-13-254 Sec. 5.7 Fig. 48
    // TH1F upsilonBkgVsNormRegion("upsilonBkgVsNormRegion", 
    // 				";Normalization window (GeV);Pred. #Upsilon bkg.", 3, -0.5, 2.5);
    // TAxis* upsilonBkgVsNormRegionXAxis = upsilonBkgVsNormRegion.GetXaxis();
    // upsilonBkgVsNormRegionXAxis->SetBinLabel(1, "<2");
    // upsilonBkgVsNormRegionXAxis->SetBinLabel(2, "2-4");
    // upsilonBkgVsNormRegionXAxis->SetBinLabel(3, "#geq4");
    // if (MTBin == "_lowMT") {
    //   upsilonBkgVsNormRegion.SetBinContent(1, 3.5);
    //   upsilonBkgVsNormRegion.SetBinError(1, 1.64);
    //   upsilonBkgVsNormRegion.SetBinContent(2, 1.47);
    //   upsilonBkgVsNormRegion.SetBinError(2, 1.38);
    //   upsilonBkgVsNormRegion.SetBinContent(3, 3.98);
    //   upsilonBkgVsNormRegion.SetBinError(3, 1.83);
    // }
    // else if (MTBin == "_highMT") {
    //   upsilonBkgVsNormRegion.SetBinContent(1, 0.913);
    //   upsilonBkgVsNormRegion.SetBinError(1, 0.971);
    //   upsilonBkgVsNormRegion.SetBinContent(2, 1.16);
    //   upsilonBkgVsNormRegion.SetBinError(2, 1.33);
    //   upsilonBkgVsNormRegion.SetBinContent(3, 1.23);
    //   upsilonBkgVsNormRegion.SetBinError(3, 1.27);
    // }
    // TCanvas upsilonBkgVsNormRegionCanvas("upsilonBkgVsNormRegionCanvas", "", 600, 600);
    // upsilonBkgVsNormRegionCanvas.SetFillStyle(0);
    // upsilonBkgVsNormRegionCanvas.SetFillColor(0);
    // upsilonBkgVsNormRegion.Draw();
    // upsilonBkgVsNormRegionCanvas.SaveAs((saveDir + 
    // 					 string(upsilonBkgVsNormRegionCanvas.GetName()) + MTBin + 
    // 					 ".pdf").c_str());
  }
  // else {

  //   //AN-13-254 Sec. 5.7 Fig. 49
  //   const string regCDataVsRegDDataFileName(analysisFilePath + "results/regCDataVsRegDData" + 
  // 					    MTBin + tag19p7InvFb + outputVTag + fileExt);
  //   TFile regCDataVsRegDDataFile(regCDataVsRegDDataFileName.c_str());
  //   saveCanvasAsPDF("muHadMassCanvas", "muHadMassCanvas_regCDataVsRegDData_" + HLTPath + MTBin, 
  // 		    saveDir, regCDataVsRegDDataFile);
  //   regCDataVsRegDDataFile.Close();

  //   //AN-13-254 Sec. 5.7 Figs. 50-52
  //   if (MTBin == "_highMT") {
  //     const string resBkgFileName(analysisFilePath + "results/resBkg_noRebin" + MTBin + 
  // 				  tag19p7InvFb + outputVTag + fileExt);
  //     TFile resBkgFile(resBkgFileName.c_str());
  //     TIter iKey(resBkgFile.GetListOfKeys());
  //     TKey* key;
  //     while (key = (TKey*)iKey()) {
  // 	string objName(key->GetName());
  // 	if ((objName.find("_0_0_weightedAvgkFixed") != string::npos) || 
  // 	    (objName.find("_0_0_regCkFixed") != string::npos) || 
  // 	    (objName.find("_0_0_regDkFixed") != string::npos)) {
  // 	  TCanvas* obj = NULL;
  // 	  resBkgFile.GetObject(objName.c_str(), obj);
  // 	  if (obj != NULL) {
  // 	    obj->SaveAs((saveDir + objName + "_" + HLTPath + ".pdf").c_str());
  // 	  }
  // 	}
  //     }
  //     saveCanvasAsPDF("frame_muHadMass_regD", "frame_muHadMass_regD_" + HLTPath, saveDir, 
  // 		      resBkgFile);
  //     saveCanvasAsPDF("resBkgSystCanvas", "resBkgSystCanvas_" + HLTPath, saveDir, resBkgFile);
  //     resBkgFile.Close();
  //   }
  // }
}
