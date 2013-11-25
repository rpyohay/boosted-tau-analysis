void formatPlots(const string& inputVersion, const string& outputVersion)
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
  gROOT->LoadMacro((macroPath + "Plot.C++").c_str());
//   gSystem->Load((macroPath + "Plot_C.so").c_str());

  //needed so vector<Color_t> and vector<Style_t> work
  vector<short> dummy;

  //set up canvas and graph names
  vector<string> canvasNames1D;
  canvasNames1D.push_back("mu1PTCanvas");
  canvasNames1D.push_back("mu2PTCanvas");
  canvasNames1D.push_back("dimuonPTCanvas");
  canvasNames1D.push_back("mu1EtaCanvas");
  canvasNames1D.push_back("mu2EtaCanvas");
  canvasNames1D.push_back("dimuonEtaCanvas");
  canvasNames1D.push_back("mu1PhiCanvas");
  canvasNames1D.push_back("mu2PhiCanvas");
  canvasNames1D.push_back("dimuonPhiCanvas");
  canvasNames1D.push_back("mu1METMTCanvas");
  canvasNames1D.push_back("mu2METMTCanvas");
  canvasNames1D.push_back("dimuonMETMTCanvas");
  canvasNames1D.push_back("dPhiMu1METCanvas");
  canvasNames1D.push_back("dPhiMu2METCanvas");
  canvasNames1D.push_back("dPhiDimuonMETCanvas");
  canvasNames1D.push_back("mu1IsoCanvas");
  canvasNames1D.push_back("mu2IsoCanvas");
  canvasNames1D.push_back("dEtaMu1Mu2Canvas");
  canvasNames1D.push_back("dPhiMu1Mu2Canvas");
  canvasNames1D.push_back("dRMu1Mu2Canvas");
  canvasNames1D.push_back("METCanvas");
  canvasNames1D.push_back("nGoodVtxCanvas");
  canvasNames1D.push_back("mDimuonCanvas");
  vector<string> canvasNames2D;
  vector<string> graphNames1D;
  graphNames1D.push_back("mu1PT");
  graphNames1D.push_back("mu2PT");
  graphNames1D.push_back("dimuonPT");
  graphNames1D.push_back("mu1Eta");
  graphNames1D.push_back("mu2Eta");
  graphNames1D.push_back("dimuonEta");
  graphNames1D.push_back("mu1Phi");
  graphNames1D.push_back("mu2Phi");
  graphNames1D.push_back("dimuonPhi");
  graphNames1D.push_back("mu1METMT");
  graphNames1D.push_back("mu2METMT");
  graphNames1D.push_back("dimuonMETMT");
  graphNames1D.push_back("dPhiMu1MET");
  graphNames1D.push_back("dPhiMu2MET");
  graphNames1D.push_back("dPhiDimuonMET");
  graphNames1D.push_back("mu1Iso");
  graphNames1D.push_back("mu2Iso");
  graphNames1D.push_back("dEtaMu1Mu2");
  graphNames1D.push_back("dPhiMu1Mu2");
  graphNames1D.push_back("dRMu1Mu2");
  graphNames1D.push_back("MET");
  graphNames1D.push_back("nGoodVtx");
  graphNames1D.push_back("mDimuon");
  vector<string> graphNames2D;

  //set up plot style options
  vector<string> legendHeaders19p7InvFb(canvasNames1D.size(), "Normalized to 19.7 fb^{-1}");
  vector<string> legendHeaders19p7InvFbDYJetsToLL(canvasNames1D.size(), 
						  "Drell-Yan + jets normalized to 19.7 fb^{-1}");
  vector<Color_t> colors;
  colors.push_back(kBlack);
  colors.push_back(/*kRed*/kAzure + 1);
  colors.push_back(/*kBlue*/kOrange + 1);
  colors.push_back(/*kMagenta*/kGreen - 2);
  colors.push_back(/*kGreen*/kMagenta + 2);
  colors.push_back(kCyan + 2);
  colors.push_back(kRed + 2);
  colors.push_back(kSpring + 4);
  colors.push_back(kViolet - 7);
  colors.push_back(kYellow);
  colors.push_back(kBlue + 1);
  colors.push_back(kGray + 2);
  colors.push_back(kMagenta - 2);
  colors.push_back(kGreen + 3);
  colors.push_back(kRed);
  vector<Style_t> styles;
  styles.push_back(20);
  styles.push_back(21);
  styles.push_back(22);
  styles.push_back(23);
  styles.push_back(24);
  styles.push_back(25);
  styles.push_back(26);
  styles.push_back(27);
  styles.push_back(28);
  styles.push_back(29);
  styles.push_back(30);
  styles.push_back(31);
  styles.push_back(32);
  styles.push_back(33);
  styles.push_back(34);
  vector<string> legendEntriesMCData;
  legendEntriesMCData.push_back("Data 19.7 fb^{-1}");
  legendEntriesMCData.push_back("W + jets");
  legendEntriesMCData.push_back("t#bar{t} + jets");
  legendEntriesMCData.push_back("Drell-Yan + jets");
  const bool setLogY = true;
  const bool setLinY = false;
  const bool drawStack = true;
  const bool drawSame = false;
  const bool dataMC = true;
  const bool sigBkg = false;

  //best available weights according to Dropbox spreadsheet
  vector<float> weightsMCData;
  weightsMCData.push_back(1.0); //data (int. lumi. = 19.7 fb^-1)
  weightsMCData.push_back(40.1864153331496); //W+jets
  weightsMCData.push_back(3.54908787697385); //tt+jets
  weightsMCData.push_back(1.0); //DYJetsToLLRelXSecWeights already weighted to 19.7 fb^-1
  vector<float> DYJetsToLLRelXSecWeights; //weighted to 19.7 fb^-1
  DYJetsToLLRelXSecWeights.push_back(9.10904899920841); //(10 < m < 50) GeV using ttH cross section
  DYJetsToLLRelXSecWeights.push_back(2.26675101231954); //m > 50 GeV using SM@8TeV Twiki

  //space-saving constant definitions
  string user(gSystem->GetFromPipe("whoami").Data());
  const string analysisFilePath("/data1/" + user + "/");
  const string fileExt(".root");
  const string tag19p7InvFb("_19p7fb-1");

  //version tags
  const string outputVTag("_" + outputVersion);
  const string dataVTag("_" + inputVersion);
  const string WJetsToLNuVTag("_" + inputVersion);
  const string TTJetsVTag("_" + inputVersion);
  const string DYJetsToLLVTag("_" + inputVersion);

  //hadd data samples from different eras
  string dataPrefix(analysisFilePath + "data/analysis/DrellYanAnalysis_SingleMu");
  string dataSuffix(dataVTag + fileExt);
  string dataHaddOutputFile(dataPrefix + dataSuffix);
  vector<string> dataHaddInputFiles;
  vector<string> runEras;
  runEras.push_back("_Run2012A");
  runEras.push_back("_Run2012B");
  runEras.push_back("_Run2012C");
  runEras.push_back("_Run2012D");
  for (vector<string>::const_iterator iRunEra = runEras.begin(); iRunEra != runEras.end(); 
       ++iRunEra) {
    stringstream dataName;
    dataName << dataPrefix << *iRunEra << dataSuffix;
    dataHaddInputFiles.push_back(dataName.str());
  }
  haddCanvases(dataHaddOutputFile, dataHaddInputFiles, vector<float>(4, 1.0), 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //hadd Drell-Yan+jets ml+l- binned samples
  string DYJetsToLLPrefix(analysisFilePath + "DYJetsToLL/analysis/DrellYanAnalysis_DYJetsToLL");
  string DYJetsToLLSuffix(DYJetsToLLVTag + fileExt);
  string DYJetsToLLHaddOutputFile(DYJetsToLLPrefix + DYJetsToLLSuffix);
  vector<string> DYJetsToLLHaddInputFiles;
  vector<string> massBins;
  massBins.push_back("_M-10To50");
  massBins.push_back("_M-50");
  for (vector<string>::const_iterator iMassBin = massBins.begin(); iMassBin != massBins.end(); 
       ++iMassBin) {
    stringstream DYJetsToLLName;
    DYJetsToLLName << DYJetsToLLPrefix << *iMassBin << DYJetsToLLSuffix;
    DYJetsToLLHaddInputFiles.push_back(DYJetsToLLName.str());
  }
  haddCanvases(DYJetsToLLHaddOutputFile, DYJetsToLLHaddInputFiles, 
	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //compare data to MC in control region
  string dataVsMCOutputFile19p7InvFb(analysisFilePath + "results/dataVsMC_DrellYanAnalysis" + 
				     tag19p7InvFb + outputVTag + fileExt);
  vector<string> dataVsMCInputFiles;
  dataVsMCInputFiles.push_back(dataHaddOutputFile);
  dataVsMCInputFiles.push_back(analysisFilePath + 
			       "WJetsToLNu/analysis/DrellYanAnalysis_WJetsToLNu" + 
			       WJetsToLNuVTag + fileExt);
  dataVsMCInputFiles.push_back(analysisFilePath + "TTJets/analysis/DrellYanAnalysis_TTJets" + 
			       TTJetsVTag + fileExt);
  dataVsMCInputFiles.push_back(DYJetsToLLHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(dataVsMCOutputFile19p7InvFb, dataVsMCInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders19p7InvFb, 
					colors, styles, legendEntriesMCData, 
					weightsMCData, setLogY, drawStack, dataMC);
}
