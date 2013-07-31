{
  //load
  gROOT->Reset();
  gROOT->ProcessLine("#include <utility>");
  string macroPath("/afs/cern.ch/user/y/yohay/CMSSW_5_3_3/src/BoostedTauAnalysis/");
  macroPath+="TauAnalyzer/test/";
  gSystem->Load((macroPath + "STLDictionary.so").c_str());
  gROOT->LoadMacro((macroPath + "Plot.C++").c_str());
//   gSystem->Load((macroPath + "Plot_C.so").c_str());

  //needed so vector<Color_t> and vector<Style_t> work
  vector<short> dummy;

  //unit strings
  string unitPTTau("Reco #tau p_{T} (GeV)");
  string unitPTMu("Reco #mu p_{T} (GeV)");
  string unitEtaTau("Reco #tau #eta");
  string unitEtaMu("Reco #mu #eta");
  string unitDR("#DeltaR(visible gen #tau, gen #mu)");
  string noUnit("");

  //map of bin labels for certain efficiency plots
  vector<string> binLabels;
  binLabels.push_back("kNull");
  binLabels.push_back("kOneProng0PiZero");
  binLabels.push_back("kOneProng1PiZero");
  binLabels.push_back("kOneProng2PiZero");
  binLabels.push_back("kOneProng3PiZero");
  binLabels.push_back("kOneProngNPiZero");
  binLabels.push_back("kTwoProng0PiZero");
  binLabels.push_back("kTwoProng1PiZero");
  binLabels.push_back("kTwoProng2PiZero");
  binLabels.push_back("kTwoProng3PiZero");
  binLabels.push_back("kTwoProngNPiZero");
  binLabels.push_back("kThreeProng0PiZero");
  binLabels.push_back("kThreeProng1PiZero");
  binLabels.push_back("kThreeProng2PiZero");
  binLabels.push_back("kThreeProng3PiZero");
  binLabels.push_back("kThreeProngNPiZero");
  binLabels.push_back("kRareDecayMode");
  map<string, vector<string> > binLabelMap;
  binLabelMap["muHadGenDecayMode"] = binLabels;
  binLabelMap["muHadCorrectRecoDecayModeGenDecayMode"] = binLabels;
  binLabelMap["muHadRecoDecayMode"] = binLabels;
  binLabelMap["muHadGen1ProngRecoDecayMode"] = binLabels;
  binLabelMap["muHadGen1Prong1Pi0RecoDecayMode"] = binLabels;
  binLabelMap["muHadGen3ProngRecoDecayMode"] = binLabels;

  //map of inputs to efficiency histograms
  map<string, pair<string, string> > effHistMap;
  effHistMap["numeratorPT"] = make_pair(string("denominatorPT"), unitPTMu);
  map<string, pair<string, string> > effHistMapMu;
  effHistMapMu["numeratorPT"] = make_pair(string("denominatorPT"), unitPTMu);
  effHistMapMu["numeratorEta"] = make_pair(string("denominatorEta"), unitEtaMu);
  map<string, pair<string, string> > effHistMapTau;
  effHistMapTau["numeratorPT"] = make_pair(string("denominatorPT"), unitPTTau);
  effHistMapTau["numeratorEta"] = make_pair(string("denominatorEta"), unitEtaTau);

  //map of inputs to 1D histograms
  map<string, string> hist1DMap;
  hist1DMap["numeratorPT"] = unitPTMu;
  hist1DMap["denominatorPT"] = unitPTMu;
  map<string, string> hist1DMapMu;
  hist1DMapMu["numeratorPT"] = unitPTMu;
  hist1DMapMu["denominatorPT"] = unitPTMu;
  hist1DMapMu["numeratorEta"] = unitEtaMu;
  hist1DMapMu["denominatorEta"] = unitEtaMu;
  map<string, string> hist1DMapTau;
  hist1DMapTau["numeratorPT"] = unitPTTau;
  hist1DMapTau["denominatorPT"] = unitPTTau;
  hist1DMapTau["numeratorEta"] = unitEtaTau;
  hist1DMapTau["denominatorEta"] = unitEtaTau;

//   vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > > histMap;
//   histMap.push_back(make_pair(make_pair(), make_pair()));

//   map<pair<string, string>, vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > > > 
//     canvasMap;
//   canvasMap[make_pair(string("muHadGen1ProngRecoDecayMode"), noUnit)] = ;
//   canvasMap[make_pair(string("muHadGen1Prong1Pi0RecoDecayMode"), noUnit)] = ;
//   canvasMap[make_pair(string("muHadGen3ProngRecoDecayMode"), noUnit)] = ;

  //set up canvas and graph names
  vector<string> canvasNames1D;
  canvasNames1D.push_back("hadTauAssociatedMuMultiplicityCanvas");
  canvasNames1D.push_back("muHadMassCanvas");
  canvasNames1D.push_back("muHadChargeCanvas");
  canvasNames1D.push_back("METCanvas");
  canvasNames1D.push_back("WMuMTCanvas");
  canvasNames1D.push_back("tauMuMTCanvas");
  canvasNames1D.push_back("dPhiWMuMETCanvas");
  canvasNames1D.push_back("dPhiTauMuMETCanvas");
  canvasNames1D.push_back("tauMuTauHadJetHTCanvas");
  canvasNames1D.push_back("diJetHTCanvas");
  canvasNames1D.push_back("jetTauJetHTCanvas");
  canvasNames1D.push_back("tauMuTauHadJetWMuHTCanvas");
  canvasNames1D.push_back("diJetWMuHTCanvas");
  canvasNames1D.push_back("jetTauJetWMuHTCanvas");
  canvasNames1D.push_back("tauMuPTCanvas");
  canvasNames1D.push_back("tauHadPTCanvas");
  canvasNames1D.push_back("tauHadIsoCanvas");
  canvasNames1D.push_back("tauHadEtaCanvas");
  canvasNames1D.push_back("softMuPTOverMuHadMassCanvas");
  canvasNames1D.push_back("muHadPTOverMuHadMassCanvas");
  canvasNames1D.push_back("dRSoftMuNearestGenMuHistCanvas");
  canvasNames1D.push_back("muHadPTCanvas");
  canvasNames1D.push_back("muHadMultiplicityCanvas");
//   canvasNames1D.push_back("nGoodVtxCanvas");
  vector<string> canvasNames2D;
  canvasNames2D.push_back("muHadMassVsDRSoftMuTauCanvas");
  canvasNames2D.push_back("tauHadIsoVsSoftMuPTCanvas");
  canvasNames2D.push_back("cleanedJetPTVsCleanedTauPTCanvas");
  canvasNames2D.push_back("uncleanedJetPTVsCleanedTauPTCanvas");
  canvasNames2D.push_back("muHadMassVsSoftMuPTCanvas");
  canvasNames2D.push_back("genMuExistsVsSoftMuNearestMuPropertiesCanvas");
  canvasNames2D.push_back("muHadMassVsTauHadEtaCanvas");
  canvasNames2D.push_back("muHadMassVsSoftMuEtaCanvas");
  canvasNames2D.push_back("muHadMassVsTauHadIsoCanvas");
  canvasNames2D.push_back("muHadMassVsTauHadPTCanvas");
  canvasNames2D.push_back("tauHadIsoVsEtaCanvas");
  canvasNames2D.push_back("tauHadEtaVsSoftMuEtaCanvas");
  canvasNames2D.push_back("dEtaTauHadSoftMuVsDPhiTauHadSoftMuCanvas");
  canvasNames2D.push_back("tauHadPTOverMuHadMassVsTauHadIsoCanvas");
  canvasNames2D.push_back("softMuPTOverMuHadMassVsTauHadIsoCanvas");
  canvasNames2D.push_back("avgTauHadSoftMuPTOverMuHadMassVsTauHadIsoCanvas");
  canvasNames2D.push_back("muHadPTOverMuHadMassVsTauHadIsoCanvas");
  canvasNames2D.push_back("softMuPTVsTauHadPTCanvas");
  vector<string> graphNames1D;
  graphNames1D.push_back("hadTauAssociatedMuMultiplicity");
  graphNames1D.push_back("muHadMass");
  graphNames1D.push_back("muHadCharge");
  graphNames1D.push_back("MET");
  graphNames1D.push_back("WMuMT");
  graphNames1D.push_back("tauMuMT");
  graphNames1D.push_back("dPhiWMuMET");
  graphNames1D.push_back("dPhiTauMuMET");
  graphNames1D.push_back("tauMuTauHadJetHT");
  graphNames1D.push_back("diJetHT");
  graphNames1D.push_back("jetTauJetHT");
  graphNames1D.push_back("tauMuTauHadJetWMuHT");
  graphNames1D.push_back("diJetWMuHT");
  graphNames1D.push_back("jetTauJetWMuHT");
  graphNames1D.push_back("tauMuPT");
  graphNames1D.push_back("tauHadPT");
  graphNames1D.push_back("tauHadIso");
  graphNames1D.push_back("tauHadEta");
  graphNames1D.push_back("softMuPTOverMuHadMass");
  graphNames1D.push_back("muHadPTOverMuHadMass");
  graphNames1D.push_back("dRSoftMuNearestGenMuHist");
  graphNames1D.push_back("muHadPT");
  graphNames1D.push_back("muHadMultiplicity");
//   graphNames1D.push_back("nGoodVtx");
  vector<string> graphNames2D;
  graphNames2D.push_back("muHadMassVsDRSoftMuTau");
  graphNames2D.push_back("tauHadIsoVsSoftMuPT");
  graphNames2D.push_back("cleanedJetPTVsCleanedTauPT");
  graphNames2D.push_back("uncleanedJetPTVsCleanedTauPT");
  graphNames2D.push_back("muHadMassVsSoftMuPT");
  graphNames2D.push_back("genMuExistsVsSoftMuNearestMuProperties");
  graphNames2D.push_back("muHadMassVsTauHadEta");
  graphNames2D.push_back("muHadMassVsSoftMuEta");
  graphNames2D.push_back("muHadMassVsTauHadIso");
  graphNames2D.push_back("muHadMassVsTauHadPT");
  graphNames2D.push_back("tauHadIsoVsEta");
  graphNames2D.push_back("tauHadEtaVsSoftMuEta");
  graphNames2D.push_back("dEtaTauHadSoftMuVsDPhiTauHadSoftMu");
  graphNames2D.push_back("tauHadPTOverMuHadMassVsTauHadIso");
  graphNames2D.push_back("softMuPTOverMuHadMassVsTauHadIso");
  graphNames2D.push_back("avgTauHadSoftMuPTOverMuHadMassVsTauHadIso");
  graphNames2D.push_back("muHadPTOverMuHadMassVsTauHadIso");
  graphNames2D.push_back("softMuPTVsTauHadPT");

  //set up plot style options
  vector<string> legendHeaders20InvFb(canvasNames1D.size(), "Normalized to 20 fb^{-1}");
  vector<string> legendHeaders2p5InvFb(canvasNames1D.size(), "Normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbDYJetsToLL(canvasNames1D.size(), 
						 "Drell-Yan + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbTTJets(canvasNames1D.size(), 
					     "t#bar{t} + jets normalized to 2.5 fb^{-1}");
  vector<string> 
    legendHeaders2p5InvFbWNJetsToLNu(canvasNames1D.size(), 
				     "W#rightarrow#mu#nu + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders1(canvasNames1D.size(), "Normalized to 1");
  vector<string> legendHeaders1DYJetsToLL(canvasNames1D.size(), 
					  "Drell-Yan + jets normalized to 1");
  vector<string> legendHeaders1TTJets(canvasNames1D.size(), "t#bar{t} + jets normalized to 1");
  vector<string> legendHeaders1WNJetsToLNu(canvasNames1D.size(), 
					   "W#rightarrow#mu#nu + jets normalized to 1");
  vector<Color_t> colors;
  colors.push_back(kBlack);
  colors.push_back(/*kRed*/kAzure + 1);
  colors.push_back(/*kBlue*/kOrange + 1);
  colors.push_back(/*kMagenta*/kGreen - 2);
  colors.push_back(/*kGreen*/kMagenta + 2);
  colors.push_back(kCyan + 2);
  colors.push_back(kRed + 2);
  colors.push_back(kSpring + 4);
  vector<Style_t> styles;
  styles.push_back(20);
  styles.push_back(21);
  styles.push_back(22);
  styles.push_back(23);
  styles.push_back(24);
  styles.push_back(25);
  styles.push_back(26);
  styles.push_back(27);
  vector<string> legendEntriesSigBkgInd;
  legendEntriesSigBkgInd.push_back("W#rightarrow#mu#nu + h_{1}");
  legendEntriesSigBkgInd.push_back("Drell-Yan + jets (10 < m_{l^{+}l^{-}} < 50) GeV");
  legendEntriesSigBkgInd.push_back("Drell-Yan + jets m_{l^{+}l^{-}} > 50 GeV");
  legendEntriesSigBkgInd.push_back("t#bar{t} + jets");
  legendEntriesSigBkgInd.push_back("W#rightarrow#mu#nu + 1 jet");
  legendEntriesSigBkgInd.push_back("W#rightarrow#mu#nu + 2 jets");
  legendEntriesSigBkgInd.push_back("W#rightarrow#mu#nu + 3 jets");
  legendEntriesSigBkgInd.push_back("W#rightarrow#mu#nu + 4 jets");
  vector<string> legendEntriesMCDataInd(legendEntriesSigBkgInd);
  legendEntriesMCDataInd[0] = "Data 2.5 fb^{-1}";
  vector<string> legendEntriesSigBkg;
  legendEntriesSigBkg.push_back("W#rightarrow#mu#nu + h_{1}");
  legendEntriesSigBkg.push_back("Drell-Yan + jets");
  legendEntriesSigBkg.push_back("t#bar{t} + jets");
  legendEntriesSigBkg.push_back("W#rightarrow#mu#nu + jets");
  vector<string> legendEntriesMCData(legendEntriesSigBkg);
  legendEntriesMCData[0] = "Data 2.5 fb^{-1}";
  vector<string> legendEntriesSearchVsControl;
  legendEntriesSearchVsControl.push_back("Isolated #tau leptons");
  legendEntriesSearchVsControl.push_back("Non-isolated #tau leptons");
  const bool setLogY = true;
  const bool setLinY = false;
  const bool drawStack = true;
  const bool drawSame = false;

  //weights (sig. figs are probably wrong)
  const float Wh1Weight20InvFb = 0.0681;
  const float Wh1Weight2p5InvFb = 0.00851; //20 fb^-1 weight * (2.5/20)
  vector<float> weights1(8, 0.0);
  vector<float> weightsSigBkgInd;
  weightsSigBkgInd.push_back(Wh1Weight20InvFb);
  weightsSigBkgInd.push_back(5.841109917); /*Drell-Yan + jets (10 < ml+l- < 50) GeV weighted to 
					     20 fb^-1*/
  weightsSigBkgInd.push_back(1.936998118); /*Drell-Yan + jets ml+l- > 50 GeV weighted to 
					     20 fb^-1*/
  weightsSigBkgInd.push_back(1.99738713); //tt+jets weighted to 20 fb^-1
  weightsSigBkgInd.push_back(4.57); //W-->munu + 1 jet weighted to 20 fb^-1
  weightsSigBkgInd.push_back(1.01); //W-->munu + 2 jets weighted to 20 fb^-1
  weightsSigBkgInd.push_back(0.654); //W-->munu + 3 jets weighted to 20 fb^-1
  weightsSigBkgInd.push_back(0.313); //W-->munu + 4 jets weighted to 20 fb^-1
  vector<float> weightsMCDataInd;
  weightsMCDataInd.push_back(1.0); //data
  weightsMCDataInd.push_back(0.7301387396); /*Drell-Yan + jets (10 < ml+l- < 50) GeV weighted to 
					      2.5 fb^-1
					      11050 pb(PREP xs) * 2500 pb^-1/
					      37835275(no. events processed)*/
  weightsMCDataInd.push_back(0.2421247648); /*Drell-Yan + jets ml+l- > 50 GeV weighted to 
					      2.5 fb^-1
					      2950 pb(PREP xs) * 2500 pb^-1/
					      30459503(no. events processed)*/
  weightsMCDataInd.push_back(0.2496733913); //tt+jets weighted to 2.5 fb^-1
  weightsMCDataInd.push_back(0.57125); //W-->munu + 1 jet weighted to 2.5 fb^-1
  weightsMCDataInd.push_back(0.12625); //W-->munu + 2 jets weighted to 2.5 fb^-1
  weightsMCDataInd.push_back(0.08175); //W-->munu + 3 jets weighted to 2.5 fb^-1
  weightsMCDataInd.push_back(0.039125); //W-->munu + 4 jets weighted to 2.5 fb^-1
  vector<float> weightsSigBkg;
  weightsSigBkg.push_back(Wh1Weight20InvFb);
  weightsSigBkg.push_back(8.0); /*DYJetsToLLRelXSecWeights already weighted to 2.5 fb^-1 ==> 
				  multiply by 20/2.5 to get overall weight for 20 fb^-1*/
  weightsSigBkg.push_back(1.99738713); //tt+jets weighted to 20 fb^-1
  weightsSigBkg.push_back(6524.691358); //W-->munu + jets weighted to 20 fb^-1
  vector<float> weightsMCData;
  weightsMCData.push_back(1.0); //data
  weightsMCData.push_back(1.0); /*DYJetsToLLRelXSecWeights already weighted to 2.5 fb^-1 ==> 
				  multiply by 1.0 to get overall weight for 2.5 fb^-1*/
  weightsMCData.push_back(0.2496733913); /*tt+jets weighted to 2.5 fb^-1
					   (136.3 pb(PREP xs) * 2500 pb^-1)/
					   1364783(no. events processed)*/
  weightsMCData.push_back(815.5864198); //W-->munu + jets 20 fb^-1 weight * (2.5/20)
  vector<float> WNJetsToLNuRelXSecWeights;
  WNJetsToLNuRelXSecWeights.push_back(0.00070041554);
  WNJetsToLNuRelXSecWeights.push_back(0.00015423878);
  WNJetsToLNuRelXSecWeights.push_back(0.00010017824);
  WNJetsToLNuRelXSecWeights.push_back(0.000047972683);
  vector<float> DYJetsToLLRelXSecWeights;
  DYJetsToLLRelXSecWeights.push_back(0.7301387396); //weighted to 2.5 fb^-1
  DYJetsToLLRelXSecWeights.push_back(0.2421247648); //weighted to 2.5 fb^-1

  //space-saving constant definitions
  const string analysisFilePath("/data1/yohay/");
  const string fileExt(".root");
  const string tag20InvFb("_20fb-1");
  const string tag2p5InvFb("_2p5fb-1");
  const string tag1("_normalizedTo1");

  //version tags
  const string outputVTag("_v34");
  const string dataVTag("_v32");
  const string sigVTag("_v25");
  const string WNJetsToLNuVTag("_v33");
  const string TTJetsVTag("_v33");
  const string DYJetsToLLVTag("_v33");

  //hadd data samples from different eras
//   string dataIsoPrefix(analysisFilePath + "data/analysis/muHadIsoAnalysis_SingleMu");
//   string dataIsoSuffix(dataVTag + fileExt);
//   string dataIsoHaddOutputFile(dataIsoPrefix + dataIsoSuffix); //BLINDED!!!
  string dataNonIsoPrefix(analysisFilePath + "data/analysis/muHadNonIsoAnalysis_SingleMu");
  string dataNonIsoSuffix(dataVTag + fileExt);
  string dataNonIsoHaddOutputFile(dataNonIsoPrefix + dataNonIsoSuffix);
//   vector<string> dataIsoHaddInputFiles; //BLINDED!!!
  vector<string> dataNonIsoHaddInputFiles;
  vector<string> runEras;
  runEras.push_back("_Run2012A");
  runEras.push_back("_Run2012B");
  runEras.push_back("_Run2012C");
  runEras.push_back("_Run2012D");
  for (vector<string>::const_iterator iRunEra = runEras.begin(); iRunEra != runEras.end(); 
       ++iRunEra) {
//     stringstream dataIsoName;
//     dataIsoName << dataIsoPrefix << *iRunEra << dataIsoSuffix; //BLINDED!!!
//     dataIsoHaddInputFiles.push_back(dataIsoName.str());
    stringstream dataNonIsoName;
    dataNonIsoName << dataNonIsoPrefix << *iRunEra << dataNonIsoSuffix;
    dataNonIsoHaddInputFiles.push_back(dataNonIsoName.str());
  }
//   haddCanvases(dataIsoHaddOutputFile, dataIsoHaddInputFiles, vector<float>(4, 1.0), 
// 	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D); //BLINDED!!!
  haddCanvases(dataNonIsoHaddOutputFile, dataNonIsoHaddInputFiles, vector<float>(4, 1.0), 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //hadd Drell-Yan+jets ml+l- binned samples
  string DYJetsToLLIsoPrefix(analysisFilePath + "DYJetsToLL/analysis/muHadIsoAnalysis_DYJetsToLL");
  string DYJetsToLLIsoSuffix(DYJetsToLLVTag + fileExt);
  string DYJetsToLLIsoHaddOutputFile(DYJetsToLLIsoPrefix + DYJetsToLLIsoSuffix);
  string DYJetsToLLNonIsoPrefix(analysisFilePath + 
				"DYJetsToLL/analysis/muHadNonIsoAnalysis_DYJetsToLL");
  string DYJetsToLLNonIsoSuffix(DYJetsToLLVTag + fileExt);
  string DYJetsToLLNonIsoHaddOutputFile(DYJetsToLLNonIsoPrefix + DYJetsToLLNonIsoSuffix);
  vector<string> DYJetsToLLIsoHaddInputFiles;
  vector<string> DYJetsToLLNonIsoHaddInputFiles;
  vector<string> massBins;
  massBins.push_back("_M-10To50");
  massBins.push_back("_M-50");
  for (vector<string>::const_iterator iMassBin = massBins.begin(); iMassBin != massBins.end(); 
       ++iMassBin) {
    stringstream DYJetsToLLIsoName;
    DYJetsToLLIsoName << DYJetsToLLIsoPrefix << *iMassBin << DYJetsToLLIsoSuffix;
    DYJetsToLLIsoHaddInputFiles.push_back(DYJetsToLLIsoName.str());
    stringstream DYJetsToLLNonIsoName;
    DYJetsToLLNonIsoName << DYJetsToLLNonIsoPrefix << *iMassBin << DYJetsToLLNonIsoSuffix;
    DYJetsToLLNonIsoHaddInputFiles.push_back(DYJetsToLLNonIsoName.str());
  }
//   haddCanvases(DYJetsToLLIsoHaddOutputFile, DYJetsToLLIsoHaddInputFiles, DYJetsToLLRelXSecWeights, 
// 	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
//   haddCanvases(DYJetsToLLNonIsoHaddOutputFile, DYJetsToLLNonIsoHaddInputFiles, 
// 	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //hadd W+jets Njets binned samples
  string WNJetsToLNuIsoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadIsoAnalysis_W");
  string WNJetsToLNuNonIsoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadNonIsoAnalysis_W");
  string WNJetsToLNuAllTauPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadAnalysis_W");
  string WNJetsToLNuSuffix("JetsToLNu" + WNJetsToLNuVTag + fileExt);
  string WNJetsToLNuIsoHaddOutputFile(WNJetsToLNuIsoPrefix + "N" + WNJetsToLNuSuffix);
  string WNJetsToLNuNonIsoHaddOutputFile(WNJetsToLNuNonIsoPrefix + "N" + WNJetsToLNuSuffix);
  string WNJetsToLNuAllTauHaddOutputFile(WNJetsToLNuAllTauPrefix + "N" + WNJetsToLNuSuffix);
  vector<string> WNJetsToLNuIsoHaddInputFiles;
  vector<string> WNJetsToLNuNonIsoHaddInputFiles;
  vector<string> WNJetsToLNuAllTauHaddInputFiles;
  for (unsigned int iNJets = 1; iNJets <= 4; ++iNJets) {
    stringstream WNJetsToLNuIsoName;
    WNJetsToLNuIsoName << WNJetsToLNuIsoPrefix << iNJets << WNJetsToLNuSuffix;
    WNJetsToLNuIsoHaddInputFiles.push_back(WNJetsToLNuIsoName.str());
    stringstream WNJetsToLNuNonIsoName;
    WNJetsToLNuNonIsoName << WNJetsToLNuNonIsoPrefix << iNJets << WNJetsToLNuSuffix;
    WNJetsToLNuNonIsoHaddInputFiles.push_back(WNJetsToLNuNonIsoName.str());
    stringstream WNJetsToLNuAllTauName;
    WNJetsToLNuAllTauName << WNJetsToLNuAllTauPrefix << iNJets << WNJetsToLNuSuffix;
    WNJetsToLNuAllTauHaddInputFiles.push_back(WNJetsToLNuAllTauName.str());
  }
//   haddCanvases(WNJetsToLNuIsoHaddOutputFile, WNJetsToLNuIsoHaddInputFiles, 
// 	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
// 	       graphNames2D);
//   haddCanvases(WNJetsToLNuNonIsoHaddOutputFile, WNJetsToLNuNonIsoHaddInputFiles, 
// 	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
// 	       graphNames2D);
//   haddCanvases(WNJetsToLNuAllTauHaddOutputFile, WNJetsToLNuAllTauHaddInputFiles, 
// 	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
// 	       graphNames2D);

  //compare MC signal to background
  string sigVsBkgOutputFile20InvFb(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + 
				   tag20InvFb + outputVTag + fileExt);
  string sigVsBkgOutputFile1(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + tag1 + 
			     outputVTag + fileExt);
  vector<string> sigVsBkgInputFiles;
  sigVsBkgInputFiles.push_back(analysisFilePath + "Wh1_Medium/muHadIsoAnalysis_Wh1" + sigVTag + 
			       fileExt);
  sigVsBkgInputFiles.push_back(DYJetsToLLIsoHaddOutputFile);
//   for (vector<string>::const_iterator iFile = DYJetsToLLNonIsoHaddInputFiles.begin(); 
//        iFile != DYJetsToLLNonIsoHaddInputFiles.end(); ++iFile) {
//     sigVsBkgInputFiles.push_back(*iFile);
//   }
  sigVsBkgInputFiles.push_back(analysisFilePath + "TTJets/analysis/muHadIsoAnalysis_TTJets" + 
			       TTJetsVTag + fileExt);
  sigVsBkgInputFiles.push_back(WNJetsToLNuIsoHaddOutputFile);
//   for (vector<string>::const_iterator iFile = WNJetsToLNuNonIsoHaddInputFiles.begin(); 
//        iFile != WNJetsToLNuNonIsoHaddInputFiles.end(); ++iFile) {
//     sigVsBkgInputFiles.push_back(*iFile);
//   }
//   drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile1, sigVsBkgInputFiles, 
// 					canvasNames1D, graphNames1D, legendHeaders1, 
// 					colors, styles, legendEntriesSigBkgInd, weights1, 
// 					setLinY, drawSame);
  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile20InvFb, sigVsBkgInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders20InvFb, 
					colors, styles, legendEntriesSigBkg, weightsSigBkg, 
					setLogY, drawStack);
  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile1, sigVsBkgInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders1, colors, 
					styles, legendEntriesSigBkg, weights1, setLinY, drawSame);

  //compare data to MC in control region
  string dataVsMCOutputFile2p5InvFb(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + 
				    tag2p5InvFb + outputVTag + fileExt);
  vector<string> dataVsMCInputFiles;
  dataVsMCInputFiles.push_back(dataNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(DYJetsToLLNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(analysisFilePath + "TTJets/analysis/muHadNonIsoAnalysis_TTJets" + 
			       TTJetsVTag + fileExt);
  dataVsMCInputFiles.push_back(WNJetsToLNuNonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(dataVsMCOutputFile2p5InvFb, dataVsMCInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders2p5InvFb, 
					colors, styles, legendEntriesMCData, 
					weightsMCData, setLogY, drawStack);

  //compare Drell-Yan+jets search sample to control sample
  string DYJetsToLLSearchVsControlOutputFile(analysisFilePath + 
					     "DYJetsToLL/analysis/isoVsNonIsoTaus" + tag1 + 
					     outputVTag + fileExt);
  vector<string> DYJetsToLLSearchVsControlInputFiles;
  DYJetsToLLSearchVsControlInputFiles.push_back(DYJetsToLLIsoHaddOutputFile);
  DYJetsToLLSearchVsControlInputFiles.push_back(DYJetsToLLNonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(DYJetsToLLSearchVsControlOutputFile, 
					DYJetsToLLSearchVsControlInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders1DYJetsToLL, colors, styles, 
					legendEntriesSearchVsControl, weights1, setLinY, drawSame);

  //compare tt+jets search sample to control sample
  string TTJetsSearchVsControlOutputFile(analysisFilePath + "TTJets/analysis/isoVsNonIsoTaus" + 
					 tag1 + outputVTag + fileExt);
  vector<string> TTJetsSearchVsControlInputFiles;
  TTJetsSearchVsControlInputFiles.push_back(analysisFilePath + 
					    "TTJets/analysis/muHadIsoAnalysis_TTJets" + 
					    TTJetsVTag + fileExt);
  TTJetsSearchVsControlInputFiles.push_back(analysisFilePath + 
					    "TTJets/analysis/muHadNonIsoAnalysis_TTJets" + 
					    TTJetsVTag + fileExt);
  drawMultipleEfficiencyGraphsOn1Canvas(TTJetsSearchVsControlOutputFile, 
					TTJetsSearchVsControlInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders1TTJets, colors, styles, 
					legendEntriesSearchVsControl, weights1, setLinY, drawSame);

  //compare W+jets search sample to control sample
  string WNJetsToLNuSearchVsControlOutputFile(analysisFilePath + 
					      "WNJetsToLNu/analysis/isoVsNonIsoTaus" + tag1 + 
					      outputVTag + fileExt);
  vector<string> WNJetsToLNuSearchVsControlInputFiles;
  WNJetsToLNuSearchVsControlInputFiles.push_back(WNJetsToLNuIsoHaddOutputFile);
  WNJetsToLNuSearchVsControlInputFiles.push_back(WNJetsToLNuNonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(WNJetsToLNuSearchVsControlOutputFile, 
					WNJetsToLNuSearchVsControlInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders1WNJetsToLNu, colors, styles, 
					legendEntriesSearchVsControl, weights1, setLinY, drawSame);

//   //compare dR(W muon, soft muon) for events with mu+had mass > 0 and > 2
//   mergePlotsIn1File(analysisFilePath + "Wh1_Medium/muHadAnalysisV8.root", 
// 		    analysisFilePath + "Wh1_Medium/dRWMuSoftMu_muHadMassGt0VsGt2_WMuonExcluded.root");

// //space-saving constant definitions
//   const string analysisFilePath(analysisFilePath + "WJets/WJets_tau_analysis");
//   const string fileExt(".root");

//   //make individual efficiency plots for signal and Z-->mumu samples
//   vector<string> effInputFiles;
//   vector<string> comparisonInputFiles;
//   effInputFiles.push_back(analysisFilePath + controlSample + objTag + controlGenSel + trigger);
//   effInputFiles.push_back(analysisFilePath + signalSample + objTag + leg + signalGenSel + trigger);
//   for (vector<string>::const_iterator iFile = effInputFiles.begin(); iFile != effInputFiles.end(); 
//        ++iFile) {
// //     comparisonInputFiles.push_back(*iFile + effTag + fileExt);
//     comparisonInputFiles.push_back(*iFile + fileExt);
// //     plotNice(*iFile + fileExt, jetAnalysisEffMap, binLabelMap, jetAnalysisMap, 
// // 	     *iFile + effTag + fileExt, "noPDF");
//   }

//   //compare Z-->mumu and signal muon efficiency
//   vector<string> canvasNames;
// //   canvasNames.push_back("eff_numeratorPT_over_denominatorPT");
// //   canvasNames.push_back("eff_numeratorEta_over_denominatorEta");
//   canvasNames.push_back("muEnergyFractionCanvas");
//   vector<string> graphNames;
// //   graphNames.push_back("divide_numeratorPT_by_denominatorPT");
// //   graphNames.push_back("divide_numeratorEta_by_denominatorEta");
//   graphNames.push_back("muEnergyFraction");
//   Color_t colors[2] = {kBlack, kRed};
//   Style_t styles[2] = {20, 21};
// //   drawMultipleEfficiencyGraphsOn1Canvas(analysisFilePath + "eff" + objTag + outputTag + 
// // 					fileExt, comparisonInputFiles, canvasNames, graphNames, 
// // 					colors, styles);
}
