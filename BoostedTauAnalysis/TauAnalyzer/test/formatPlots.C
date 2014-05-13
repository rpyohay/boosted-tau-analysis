void formatPlots(const string& inputVersion, const string& outputVersion, 
		 const bool compile, const bool doNoHPSIsoCut = false)
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
  if (compile) {
    gROOT->LoadMacro((macroPath + "Miscellaneous.C++").c_str());
    gROOT->LoadMacro((macroPath + "Error.C++").c_str());
    gROOT->LoadMacro((macroPath + "Plot.C++").c_str());
  }
  else {
    gSystem->Load((macroPath + "Miscellaneous_C.so").c_str());
    gSystem->Load((macroPath + "Error_C.so").c_str());
    gSystem->Load((macroPath + "Plot_C.so").c_str());
  }

  //needed so vector<Color_t> and vector<Style_t> work
  vector<short> dummy;

  //set up canvas and graph names and blinded bins for data
  vector<string> canvasNames1D;
  canvasNames1D.push_back("hadTauAssociatedMuMultiplicityCanvas");
  canvasNames1D.push_back("muHadMassCanvas");
  canvasNames1D.push_back("muHadMass1ProngCanvas");
  canvasNames1D.push_back("muHadMass1Prong1Pi0Canvas");
  canvasNames1D.push_back("muHadMass1Prong2Pi0Canvas");
  canvasNames1D.push_back("muHadMass3ProngCanvas");
  canvasNames1D.push_back("muHadMassReweightErrSqCanvas");
  canvasNames1D.push_back("muHadChargeCanvas");
  canvasNames1D.push_back("METCanvas");
  canvasNames1D.push_back("bTagDiscrimCanvas");
  canvasNames1D.push_back("WMuIsoCanvas");
  canvasNames1D.push_back("WMuMTCanvas");
  canvasNames1D.push_back("tauMuMTCanvas");
  canvasNames1D.push_back("tauHadMTCanvas");
  canvasNames1D.push_back("dPhiWMuMETCanvas");
  canvasNames1D.push_back("dPhiTauMuMETCanvas");
  canvasNames1D.push_back("tauMuTauHadJetHTCanvas");
  canvasNames1D.push_back("diJetHTCanvas");
  canvasNames1D.push_back("jetTauJetHTCanvas");
  canvasNames1D.push_back("tauMuTauHadJetWMuHTCanvas");
  canvasNames1D.push_back("tauMuTauHadJetWMuMETHTCanvas");
  canvasNames1D.push_back("diJetWMuHTCanvas");
  canvasNames1D.push_back("jetTauJetWMuHTCanvas");
  canvasNames1D.push_back("tauMuPTCanvas");
  canvasNames1D.push_back("tauHadPTCanvas");
  canvasNames1D.push_back("tauHadPT1ProngCanvas");
  canvasNames1D.push_back("tauHadPT1Prong1Pi0Canvas");
  canvasNames1D.push_back("tauHadPT1Prong2Pi0Canvas");
  canvasNames1D.push_back("tauHadPT3ProngCanvas");
  canvasNames1D.push_back("tauHadIsoCanvas");
  canvasNames1D.push_back("WMuLeptonRelIsoCanvas");
  canvasNames1D.push_back("tauHadEtaCanvas");
  canvasNames1D.push_back("softMuPTOverMuHadMassCanvas");
  canvasNames1D.push_back("muHadPTOverMuHadMassCanvas");
  canvasNames1D.push_back("dRSoftMuNearestGenMuHistCanvas");
  canvasNames1D.push_back("muHadPTCanvas");
  canvasNames1D.push_back("muHadMultiplicityCanvas");
  canvasNames1D.push_back("nGoodVtxCanvas");
  canvasNames1D.push_back("mWMuTauMuCanvas");
  canvasNames1D.push_back("PDGIDNearestStatus1GenParticleToSoftMuCanvas");
  canvasNames1D.push_back("PDGIDMuSrcCanvas");
  canvasNames1D.push_back("mSecondJetMuHadCanvas");
  canvasNames1D.push_back("dPhiMuHadSecondJetCanvas");
  canvasNames1D.push_back("muHadUncleanedJetPTRankCanvas");
  canvasNames1D.push_back("nAddlHardMuonsCanvas");
  canvasNames1D.push_back("nAddlJetsPTGeq0Canvas");
  canvasNames1D.push_back("nAddlJetsPTGeq20Canvas");
  canvasNames1D.push_back("nAddlJetsPTGeq40Canvas");
  canvasNames1D.push_back("tauHadDecayModeCanvas");
  canvasNames1D.push_back("dRSoftMuNearestGenZOrTTMuCanvas");
  canvasNames1D.push_back("dRSoftMuNearestGenZOrTTMuFSRCanvas");
  canvasNames1D.push_back("jet_pt_etacutCanvas");
  canvasNames1D.push_back("jet_etaCanvas");
  canvasNames1D.push_back("jet_phiCanvas");
  canvasNames1D.push_back("jet_mass_etacutCanvas");
  canvasNames1D.push_back("jet_ptmj_etacutCanvas");
  canvasNames1D.push_back("muHad_t3t1Canvas");
  canvasNames1D.push_back("muHad_t2t1Canvas");
  canvasNames1D.push_back("muHad_t3t1_pT1020Canvas");
  canvasNames1D.push_back("muHad_t3t1_pT2030Canvas");
  canvasNames1D.push_back("muHad_t3t1_pT3040Canvas");
  canvasNames1D.push_back("muHad_t3t1_pT4050Canvas");
  canvasNames1D.push_back("muHad_t3t1_pT50UpCanvas");
  canvasNames1D.push_back("muHad_t3t1_0JetsCanvas");
  canvasNames1D.push_back("muHad_t3t1_1JetsCanvas");
  canvasNames1D.push_back("muHad_t3t1_2JetsCanvas");
  canvasNames1D.push_back("muHad_t3t1_3JetsCanvas");
  canvasNames1D.push_back("muHad_t3t1_4JetsCanvas");
  canvasNames1D.push_back("muHad_t3t1_5JetsCanvas");
  canvasNames1D.push_back("muHad_t3t1_MoreJetsCanvas");
  canvasNames1D.push_back("muHadNchtrk_0_Canvas");
  canvasNames1D.push_back("muHadNchtrk_1_Canvas");
  canvasNames1D.push_back("muHadNchtrk_10_Canvas");
  canvasNames1D.push_back("muHadNchtrk_30_Canvas");
  canvasNames1D.push_back("secondNchtrk_0_Canvas");
  canvasNames1D.push_back("secondNchtrk_1_Canvas");
  canvasNames1D.push_back("secondNchtrk_10_Canvas");
  canvasNames1D.push_back("secondNchtrk_30_Canvas");
  canvasNames1D.push_back("dPhiWMuSoftMuCanvas");
  canvasNames1D.push_back("dPhiWMuSoftMuWithCutCanvas");
  canvasNames1D.push_back("dPhiWMuSecJetCanvas");
  canvasNames1D.push_back("WMu3rdTightMuChargeProductCanvas");
  canvasNames1D.push_back("tauHadPhotonEnergyFractionCanvas");
  canvasNames1D.push_back("dThetaPhotonOtherTauConstituentsCanvas");
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
  canvasNames2D.push_back("WMuIsoVsTauHadIsoCanvas");
  canvasNames2D.push_back("softMuPTVsTauHadPTCanvas");
  canvasNames2D.push_back("muHadPTOverMuHadMassVsMWMuSoftMuCanvas");
  canvasNames2D.push_back("softMuPFPTVsRECOPTCanvas");
  canvasNames2D.push_back("softMuPFEtaVsRECOEtaCanvas");
  canvasNames2D.push_back("mSecondJetMuHadVsMuHadMassCanvas");
  canvasNames2D.push_back("muHadMassVsSecondJetMassCanvas");
  canvasNames2D.push_back("muHadPTOverMVsSecondJetPTOverMCanvas");
  canvasNames2D.push_back("muHadMassVsSecondJetPTOverMCanvas");
  canvasNames2D.push_back("tauMuTauHadJetWMuHTVsMETCanvas");
  canvasNames2D.push_back("mWMuTauMuVsMWMuTauMuTauHadCanvas");
  canvasNames2D.push_back("mWMuTauMuVsMWMuTauMuTauHadGenFSRCanvas");
  canvasNames2D.push_back("muHadMassVsMWMuTauMuTauHadCanvas");
  canvasNames2D.push_back("WMuMTVsMETCanvas");
  canvasNames2D.push_back("muHad_t3t1VsptmjCanvas");
  canvasNames2D.push_back("muHad_t3t1VsDecayModeCanvas");
  canvasNames2D.push_back("muHadMassVsNAddlJetsCanvas");
  canvasNames2D.push_back("muHadMassVsCSVScoreCanvas");
  canvasNames2D.push_back("tauHadJetEnergyFractionVsTauHadIsoCanvas");
  canvasNames2D.push_back("tauHadCleanedJetEnergyFractionVsTauHadIsoCanvas");
  vector<string> graphNames1D;
  graphNames1D.push_back("hadTauAssociatedMuMultiplicity");
  graphNames1D.push_back("muHadMass");
  graphNames1D.push_back("muHadMass1Prong");
  graphNames1D.push_back("muHadMass1Prong1Pi0");
  graphNames1D.push_back("muHadMass1Prong2Pi0");
  graphNames1D.push_back("muHadMass3Prong");
  graphNames1D.push_back("muHadMassReweightErrSq");
  graphNames1D.push_back("muHadCharge");
  graphNames1D.push_back("MET");
  graphNames1D.push_back("bTagDiscrim");
  graphNames1D.push_back("WMuIso");
  graphNames1D.push_back("WMuMT");
  graphNames1D.push_back("tauMuMT");
  graphNames1D.push_back("tauHadMT");
  graphNames1D.push_back("dPhiWMuMET");
  graphNames1D.push_back("dPhiTauMuMET");
  graphNames1D.push_back("tauMuTauHadJetHT");
  graphNames1D.push_back("diJetHT");
  graphNames1D.push_back("jetTauJetHT");
  graphNames1D.push_back("tauMuTauHadJetWMuHT");
  graphNames1D.push_back("tauMuTauHadJetWMuMETHT");
  graphNames1D.push_back("diJetWMuHT");
  graphNames1D.push_back("jetTauJetWMuHT");
  graphNames1D.push_back("tauMuPT");
  graphNames1D.push_back("tauHadPT");
  graphNames1D.push_back("tauHadPT1Prong");
  graphNames1D.push_back("tauHadPT1Prong1Pi0");
  graphNames1D.push_back("tauHadPT1Prong2Pi0");
  graphNames1D.push_back("tauHadPT3Prong");
  graphNames1D.push_back("tauHadIso");
  graphNames1D.push_back("WMuLeptonRelIso");
  graphNames1D.push_back("tauHadEta");
  graphNames1D.push_back("softMuPTOverMuHadMass");
  graphNames1D.push_back("muHadPTOverMuHadMass");
  graphNames1D.push_back("dRSoftMuNearestGenMuHist");
  graphNames1D.push_back("muHadPT");
  graphNames1D.push_back("muHadMultiplicity");
  graphNames1D.push_back("nGoodVtx");
  graphNames1D.push_back("mWMuTauMu");
  graphNames1D.push_back("PDGIDNearestStatus1GenParticleToSoftMu");
  graphNames1D.push_back("PDGIDMuSrc");
  graphNames1D.push_back("mSecondJetMuHad");
  graphNames1D.push_back("dPhiMuHadSecondJet");
  graphNames1D.push_back("muHadUncleanedJetPTRank");
  graphNames1D.push_back("nAddlHardMuons");
  graphNames1D.push_back("nAddlJetsPTGeq0");
  graphNames1D.push_back("nAddlJetsPTGeq20");
  graphNames1D.push_back("nAddlJetsPTGeq40");
  graphNames1D.push_back("tauHadDecayMode");
  graphNames1D.push_back("dRSoftMuNearestGenZOrTTMu");
  graphNames1D.push_back("dRSoftMuNearestGenZOrTTMuFSR");
  graphNames1D.push_back("jet_pt_etacut");
  graphNames1D.push_back("jet_eta");
  graphNames1D.push_back("jet_phi");
  graphNames1D.push_back("jet_mass_etacut");
  graphNames1D.push_back("jet_ptmj_etacut");
  graphNames1D.push_back("muHad_t3t1");
  graphNames1D.push_back("muHad_t2t1");
  graphNames1D.push_back("muHad_t3t1_pT1020");
  graphNames1D.push_back("muHad_t3t1_pT2030");
  graphNames1D.push_back("muHad_t3t1_pT3040");
  graphNames1D.push_back("muHad_t3t1_pT4050");
  graphNames1D.push_back("muHad_t3t1_pT50Up");
  graphNames1D.push_back("muHad_t3t1_0Jets");
  graphNames1D.push_back("muHad_t3t1_1Jets");
  graphNames1D.push_back("muHad_t3t1_2Jets");
  graphNames1D.push_back("muHad_t3t1_3Jets");
  graphNames1D.push_back("muHad_t3t1_4Jets");
  graphNames1D.push_back("muHad_t3t1_5Jets");
  graphNames1D.push_back("muHad_t3t1_MoreJets");
  graphNames1D.push_back("muHad_Nchtrk_0");
  graphNames1D.push_back("muHad_Nchtrk_1");
  graphNames1D.push_back("muHad_Nchtrk_10");
  graphNames1D.push_back("muHad_Nchtrk_30");
  graphNames1D.push_back("second_Nchtrk_0");
  graphNames1D.push_back("second_Nchtrk_1");
  graphNames1D.push_back("second_Nchtrk_10");
  graphNames1D.push_back("second_Nchtrk_30");
  graphNames1D.push_back("dPhiWMuSoftMu");
  graphNames1D.push_back("dPhiWMuSoftMu_withCut");
  graphNames1D.push_back("dPhiWMuSecJet");
  graphNames1D.push_back("WMu3rdTightMuChargeProduct");
  graphNames1D.push_back("tauHadPhotonEnergyFraction");
  graphNames1D.push_back("dThetaPhotonOtherTauConstituents");
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
  graphNames2D.push_back("WMuIsoVsTauHadIso");
  graphNames2D.push_back("softMuPTVsTauHadPT");
  graphNames2D.push_back("muHadPTOverMuHadMassVsMWMuSoftMu");
  graphNames2D.push_back("softMuPFPTVsRECOPT");
  graphNames2D.push_back("softMuPFEtaVsRECOEta");
  graphNames2D.push_back("mSecondJetMuHadVsMuHadMass");
  graphNames2D.push_back("muHadMassVsSecondJetMass");
  graphNames2D.push_back("muHadPTOverMVsSecondJetPTOverM");
  graphNames2D.push_back("muHadMassVsSecondJetPTOverM");
  graphNames2D.push_back("tauMuTauHadJetWMuHTVsMET");
  graphNames2D.push_back("mWMuTauMuVsMWMuTauMuTauHad");
  graphNames2D.push_back("mWMuTauMuVsMWMuTauMuTauHadGenFSR");
  graphNames2D.push_back("muHadMassVsMWMuTauMuTauHad");
  graphNames2D.push_back("WMuMTVsMET");
  graphNames2D.push_back("muHad_t3t1Vsptmj");
  graphNames2D.push_back("muHad_t3t1VsDecayMode");
  graphNames2D.push_back("muHadMassVsNAddlJets");
  graphNames2D.push_back("muHadMassVsCSVScore");
  graphNames2D.push_back("tauHadJetEnergyFractionVsTauHadIso");
  graphNames2D.push_back("tauHadCleanedJetEnergyFractionVsTauHadIso");
  vector<Int_t> nullBlindLow(canvasNames1D.size(), 0);
  vector<Int_t> nullBlindHigh(canvasNames1D.size(), -2);
  vector<Int_t> dataBlindLow(canvasNames1D.size(), 0);
  dataBlindLow[1] = 3;
  dataBlindLow[2] = 3;
  dataBlindLow[3] = 3;
  dataBlindLow[4] = 3;
  dataBlindLow[5] = 3;

  vector<Int_t> dataBlindHigh(canvasNames1D.size(), -1);

  //set up plot style options
  vector<string> legendHeaders19p7InvFb(canvasNames1D.size(), "Normalized to 19.7 fb^{-1}");
  vector<string> legendHeaders1(canvasNames1D.size(), "Normalized to 1");
  vector<string> legendHeaders1QCD(canvasNames1D.size(), "QCD #mu-enriched normalized to 1");
  vector<string> legendHeaders1QCDB(canvasNames1D.size(), "QCD b-enriched normalized to 1");
  vector<string> legendHeaders1QCDBMu(canvasNames1D.size(), 
				      "QCD b#rightarrow#mu-enriched normalized to 1");
  vector<string> legendHeaders1DYJetsToLL(canvasNames1D.size(), 
					  "Drell-Yan + jets normalized to 1");
  vector<string> legendHeaders1TTJets(canvasNames1D.size(), "t#bar{t} + jets normalized to 1");
  vector<string> legendHeaders1T(canvasNames1D.size(), "t/#bar{t} normalized to 1");
  vector<string> legendHeaders1WNJetsToLNu(canvasNames1D.size(), "W + #geq1 jet normalized to 1");
  vector<string> legendHeaders1Wbb(canvasNames1D.size(), "W + b#bar{b} normalized to 1");
  vector<string> legendHeaders1WJetsToLNu(canvasNames1D.size(), "W + jets normalized to 1");
  vector<string> legendHeaders1WZ(canvasNames1D.size(), "WZ normalized to 1");
  vector<string> legendHeaders1ZZ(canvasNames1D.size(), "ZZ normalized to 1");
  vector<string> legendHeaders1WW(canvasNames1D.size(), "WW normalized to 1");
  vector<string> legendHeaders1NonIsoWData(canvasNames1D.size(), 
					   "Non-isolated W data normalized to 1");
  vector<string> legendHeaders1SinglePhotonData(canvasNames1D.size(), 
						"Single photon data normalized to 1");
  vector<Color_t> colors;
  colors.push_back(kBlack);
  colors.push_back(kAzure + 1);
  colors.push_back(kOrange + 1);
  colors.push_back(kGreen - 2);
  colors.push_back(kMagenta + 2);
  colors.push_back(kCyan + 2);
  colors.push_back(kRed + 2);
  colors.push_back(kYellow);
  colors.push_back(kViolet - 7);
  colors.push_back(kSpring + 4);
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
  vector<string> legendEntriesSigBkg;
  legendEntriesSigBkg.push_back("Wh_{1}");
  legendEntriesSigBkg.push_back("gg fusion");
//   legendEntriesSigBkg.push_back("Data 19.7 fb^{-1} region A");
//   legendEntriesSigBkg.push_back("QCD");
//   legendEntriesSigBkg.push_back("QCDB");
//   legendEntriesSigBkg.push_back("QCDBMu");
  legendEntriesSigBkg.push_back("Drell-Yan + jets");
  legendEntriesSigBkg.push_back("t#bar{t} + jets");
  legendEntriesSigBkg.push_back("t/#bar{t}");
  legendEntriesSigBkg.push_back("W + #geq1 jet");
//   legendEntriesSigBkg.push_back("W + b#bar{b}");
//   legendEntriesSigBkg.push_back("W + jets");
  legendEntriesSigBkg.push_back("WZ");
  legendEntriesSigBkg.push_back("ZZ");
  legendEntriesSigBkg.push_back("WW");
  std::reverse(legendEntriesSigBkg.begin() + 2, legendEntriesSigBkg.end());
  vector<string> legendEntriesSigBkgQCDFromData(legendEntriesSigBkg);
  legendEntriesSigBkgQCDFromData.push_back("QCD (from data)");
  vector<string> legendEntriesMCData(legendEntriesSigBkg);
  legendEntriesMCData[0] = "Data 19.7 fb^{-1}";
  legendEntriesMCData.erase(legendEntriesMCData.begin() + 1);
  vector<string> legendEntriesMCDataQCDFromData(legendEntriesMCData);
  legendEntriesMCDataQCDFromData.push_back("QCD (from data)");
  vector<string> legendEntriesSearchVsControl;
  legendEntriesSearchVsControl.push_back("Isolated #tau leptons");
  legendEntriesSearchVsControl.push_back("Non-isolated #tau leptons");
  const bool setLogY = true;
  const bool setLinY = false;
  const bool drawStack = true;
  const bool drawSame = false;
  const bool dataMC = true;
  const bool sigBkg = false;

  //best available weights according to Dropbox spreadsheet
  const float Wh1Weight19p7InvFb = 0.0058903;
  const float ggWeight19p7InvFb = 1.95295217378794;
  vector<float> weights1(15, 0.0);
  vector<float> weightsSigBkg;
  weightsSigBkg.push_back(1.0); //Wh1 already weighted to 19.7 fb^-1
  weightsSigBkg.push_back(1.0); //gg already weighted to 19.7 fb^-1
//   weightsSigBkg.push_back(1.0);
//   weightsSigBkg.push_back(1.0); //QCDRelXSecWeights already weighted to 19.7 fb^-1
//   weightsSigBkg.push_back(1.0); //QCDBRelXSecWeights already weighted to 19.7 fb^-1
//   weightsSigBkg.push_back(1.0); //QCDBMuRelXSecWeights already weighted to 19.7 fb^-1
  weightsSigBkg.push_back(1.0); //DYJetsToLLRelXSecWeights already weighted to 19.7 fb^-1
  weightsSigBkg.push_back(1.0); //tt + jets already weighted to 19.7 fb^-1
  weightsSigBkg.push_back(1.0); //single top already weighted to 19.7 fb^-1
  weightsSigBkg.push_back(1.0); //W + >=1 jet already weighted to 19.7 fb^-1
//   weightsSigBkg.push_back(1.0); //W + bb already weighted to 19.7 fb^-1
//   weightsSigBkg.push_back(1.0); //W + jets already weighted to 19.7 fb^-1
  weightsSigBkg.push_back(1.0); //WZ already weighted to 19.7 fb^-1
  weightsSigBkg.push_back(1.0); //ZZ already weighted to 19.7 fb^-1
  weightsSigBkg.push_back(1.0); //WW already weighted to 19.7 fb^-1
  std::reverse(weightsSigBkg.begin() + 2, weightsSigBkg.end());
  vector<float> weightsSigBkgQCDFromData(weightsSigBkg);
  weightsSigBkgQCDFromData.push_back(1.0); //QCD estimate from data already weighted to 19.7 fb^-1
  vector<float> weightsMCData;
  weightsMCData.push_back(1.0); //data (int. lumi. = 19.7 fb^-1)
//   weightsMCData.push_back(1.0);
//   weightsMCData.push_back(1.0); //QCDRelXSecWeights already weighted to 19.7 fb^-1
//   weightsMCData.push_back(1.0); //QCDBRelXSecWeights already weighted to 19.7 fb^-1
//   weightsMCData.push_back(1.0); //QCDBMuRelXSecWeights already weighted to 19.7 fb^-1
  weightsMCData.push_back(1.0); //DYJetsToLLRelXSecWeights already weighted to 19.7 fb^-1
  weightsMCData.push_back(1.0); //tt + jets already weighted to 19.7 fb^-1
  weightsMCData.push_back(1.0); //single top already weighted to 19.7 fb^-1
  weightsMCData.push_back(1.0); //W + >=1 jet already weighted to 19.7 fb^-1
//   weightsMCData.push_back(1.0); //W + bb already weighted to 19.7 fb^-1
//   weightsMCData.push_back(1.0); //W + jets already weighted to 19.7 fb^-1
  weightsMCData.push_back(1.0); //WZ already weighted to 19.7 fb^-1
  weightsMCData.push_back(1.0); //ZZ already weighted to 19.7 fb^-1
  weightsMCData.push_back(1.0); //WW already weighted to 19.7 fb^-1
  std::reverse(weightsMCData.begin() + 1, weightsMCData.end());
  vector<float> weightsMCDataQCDFromData(weightsMCData);
  weightsMCDataQCDFromData.push_back(1.0); //QCD estimate from data already weighted to 19.7 fb^-1
  vector<float> QCDRelXSecWeights;  //weighted to 19.7 fb^-1 using PREP cross sections
  QCDRelXSecWeights.push_back(4330.24221789241); // Pt 20-30
  QCDRelXSecWeights.push_back(1661.46760576197); // Pt 30-50
  QCDRelXSecWeights.push_back(334.859498535006); // Pt 50-80
  QCDRelXSecWeights.push_back(86.2492128172084); // Pt 80-120
  QCDRelXSecWeights.push_back(17.2948414684422); // Pt 120-170
  QCDRelXSecWeights.push_back(5.90683539273479); // Pt 170-300
  QCDRelXSecWeights.push_back(0.381825191984792); // Pt 300-470
  QCDRelXSecWeights.push_back(0.061429134916651); // Pt 470-600
  QCDRelXSecWeights.push_back(0.0136978188679245); // Pt 600-800
  QCDRelXSecWeights.push_back(0.00176856029171443); // Pt 800-1000
  QCDRelXSecWeights.push_back(0.000431775067953546); // Pt 1000
  vector<float> QCDBRelXSecWeights; //weighted to 19.7 fb^-1 using PREP cross sections
  QCDBRelXSecWeights.push_back(265005.9602); // Pt 15-30
  QCDBRelXSecWeights.push_back(21786.76565); // Pt 30-50
  QCDBRelXSecWeights.push_back(7109.004289); // Pt 50-150
  QCDBRelXSecWeights.push_back(352.2478877); // Pt 150
  vector<float> QCDBMuRelXSecWeights; /*weighted to 19.7 fb^-1 using cross sections from J. 
					Antonelli's e-mail*/
  QCDBMuRelXSecWeights.push_back(282.0); // Pt 15-30
  QCDBMuRelXSecWeights.push_back(163.0); // Pt 30-50
  QCDBMuRelXSecWeights.push_back(91.2); // Pt 50-150
  QCDBMuRelXSecWeights.push_back(2.22); // Pt 150
  vector<float> DYJetsToLLRelXSecWeights; //weighted to 19.7 fb^-1
  DYJetsToLLRelXSecWeights.push_back(7.65500977593); /*(10 < m < 50) GeV using ttH cross section*/
  DYJetsToLLRelXSecWeights.push_back(2.26606084150487); //m > 50 GeV using SM@8TeV Twiki
  vector<float> TRelXSecWeights; //weighted to 19.7 fb^-1 using SM@8TeV Twiki
  TRelXSecWeights.push_back(0.287208465885267); //t s-channel
  TRelXSecWeights.push_back(0.34672); //tbar s-channel
  TRelXSecWeights.push_back(0.321286023155796); //t t-channel
  TRelXSecWeights.push_back(0.312541342130939); //tbar t-channel
  vector<float> WNJetsToLNuRelXSecWeights; //weighted to 19.7 fb^-1 using PREP cross sections
  WNJetsToLNuRelXSecWeights.push_back(5.61973075498071); //W + 1 jet
  WNJetsToLNuRelXSecWeights.push_back(1.22833526483929); //W + 2 jets
  WNJetsToLNuRelXSecWeights.push_back(0.819050332546934); //W + 3 jets
  WNJetsToLNuRelXSecWeights.push_back(0.31501621894905); //W + 4 jets

  //space-saving constant definitions
  string user(gSystem->GetFromPipe("whoami").Data());
  const string analysisFilePath("/data1/" + user + "/");
  const string fileExt(".root");
  const string tag19p7InvFb("_19p7fb-1");
  const string tag1("_normalizedTo1");

  //version tags
  const string outputVTag("_" + outputVersion);
  const string dataVTag("_" + inputVersion);
  const string nonIsoWDataVTag("_" + inputVersion);
  const string SinglePhotonDataVTag("_" + inputVersion);
  const string Wh1SigVTag("_" + inputVersion);
  const string ggSigVTag("_" + inputVersion);
  const string QCDVTag("_" + inputVersion);
  const string QCDBVTag("_" + inputVersion);
  const string QCDBMuVTag("_" + inputVersion);
  const string DYJetsToLLVTag("_" + inputVersion);
  const string TTJetsVTag("_" + inputVersion);
  const string TVTag("_" + inputVersion);
  const string WNJetsToLNuVTag("_" + inputVersion);
  const string WbbVTag("_" + inputVersion);
  const string WJetsToLNuVTag("_" + inputVersion);
  const string WZVTag("_" + inputVersion);
  const string ZZVTag("_" + inputVersion);
  const string WWVTag("_" + inputVersion);

  cout << "Begin hadding...\n";

  //hadd data samples from different eras
  cout << "...data\n";
  string dataSuffix(dataVTag + fileExt);
  string dataIsoPrefix(analysisFilePath + "data/analysis/muHadIsoAnalysis_SingleMu");
  string dataIsoHaddOutputFile(dataIsoPrefix + dataSuffix); //BLINDED!!!
  string dataNonIsoPrefix(analysisFilePath + "data/analysis/muHadNonIsoAnalysis_SingleMu");
  string dataNonIsoHaddOutputFile(dataNonIsoPrefix + dataSuffix);
  string dataNonIsoReweightPrefix(analysisFilePath + 
				  "data/analysis/muHadNonIsoReweightAnalysis_SingleMu");
  string dataNonIsoReweightHaddOutputFile(dataNonIsoReweightPrefix + dataSuffix);
  string dataAllPrefix(analysisFilePath + "data/analysis/muHadAnalysis_SingleMu");
  string dataAllHaddOutputFile(dataAllPrefix + dataSuffix); //BLINDED!!!
  vector<string> dataIsoHaddInputFiles; //BLINDED!!!
  vector<string> dataNonIsoHaddInputFiles;
  vector<string> dataNonIsoReweightHaddInputFiles;
  vector<string> dataAllHaddInputFiles; //BLINDED!!!
  vector<string> runEras;
  runEras.push_back("_Run2012A");
  runEras.push_back("_Run2012B");
  runEras.push_back("_Run2012C");
  runEras.push_back("_Run2012D");
  vector<string> subJobs;
  subJobs.push_back("_0");
  subJobs.push_back("_1");
  subJobs.push_back("_2");
  subJobs.push_back("_3");
  subJobs.push_back("_4");
  subJobs.push_back("_5");
  subJobs.push_back("_6");
  subJobs.push_back("_7");
  for (vector<string>::const_iterator iRunEra = runEras.begin(); iRunEra != runEras.end(); 
       ++iRunEra) {
    for (vector<string>::const_iterator iSubJob = subJobs.begin(); iSubJob != subJobs.end(); 
	 ++iSubJob) {
      stringstream dataIsoName;
      dataIsoName << dataIsoPrefix << *iRunEra << *iSubJob << dataSuffix; //BLINDED!!!
      dataIsoHaddInputFiles.push_back(dataIsoName.str());
      stringstream dataNonIsoName;
      dataNonIsoName << dataNonIsoPrefix << *iRunEra << *iSubJob << dataSuffix;
      dataNonIsoHaddInputFiles.push_back(dataNonIsoName.str());
      stringstream dataNonIsoReweightName;
      dataNonIsoReweightName << dataNonIsoReweightPrefix << *iRunEra << *iSubJob << dataSuffix;
      dataNonIsoReweightHaddInputFiles.push_back(dataNonIsoReweightName.str());
      stringstream dataAllName;
      dataAllName << dataAllPrefix << *iRunEra << *iSubJob << dataSuffix;
      dataAllHaddInputFiles.push_back(dataAllName.str());
    }
  }
  haddCanvases(dataIsoHaddOutputFile, dataIsoHaddInputFiles, vector<float>(32, 1.0), 
  	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, dataBlindLow, 
  	       dataBlindHigh); //BLINDED!!!
  haddCanvases(dataNonIsoHaddOutputFile, dataNonIsoHaddInputFiles, vector<float>(32, 1.0), 
  	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
  	       nullBlindHigh);
  haddCanvases(dataNonIsoReweightHaddOutputFile, dataNonIsoReweightHaddInputFiles, 
  	       vector<float>(32, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
  	       nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(dataAllHaddOutputFile, dataAllHaddInputFiles, 
		 vector<float>(32, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
		 dataBlindLow, dataBlindHigh); //BLINDED!!!
  }

  //hadd non-isolated W data samples from different eras
  cout << "...non-isolated W data\n";
  string nonIsoWDataSuffix(nonIsoWDataVTag + fileExt);
  string nonIsoWDataIsoPrefix(analysisFilePath + 
			      "nonIsoWData/analysis/nonIsoW_muHadIsoAnalysis_SingleMu");
  string nonIsoWDataIsoHaddOutputFile(nonIsoWDataIsoPrefix + nonIsoWDataSuffix);
  string nonIsoWDataNonIsoPrefix(analysisFilePath + 
				 "nonIsoWData/analysis/nonIsoW_muHadNonIsoAnalysis_SingleMu");
  string nonIsoWDataNonIsoHaddOutputFile(nonIsoWDataNonIsoPrefix + nonIsoWDataSuffix);
  string nonIsoWDataNonIsoReweightPrefix(analysisFilePath + "nonIsoWData/analysis/nonIsoW_muHadNonIsoReweightAnalysis_SingleMu");
  string nonIsoWDataNonIsoReweightHaddOutputFile(nonIsoWDataNonIsoReweightPrefix + 
						 nonIsoWDataSuffix);
  string nonIsoWDataAllPrefix(analysisFilePath + 
			      "nonIsoWData/analysis/nonIsoW_muHadAnalysis_SingleMu");
  string nonIsoWDataAllHaddOutputFile(nonIsoWDataAllPrefix + nonIsoWDataSuffix);
  vector<string> nonIsoWDataIsoHaddInputFiles;
  vector<string> nonIsoWDataNonIsoHaddInputFiles;
  vector<string> nonIsoWDataNonIsoReweightHaddInputFiles;
  vector<string> nonIsoWDataAllHaddInputFiles;
  for (vector<string>::const_iterator iRunEra = runEras.begin(); iRunEra != runEras.end(); 
       ++iRunEra) {
    stringstream nonIsoWDataIsoName;
    nonIsoWDataIsoName << nonIsoWDataIsoPrefix << *iRunEra << nonIsoWDataSuffix;
    nonIsoWDataIsoHaddInputFiles.push_back(nonIsoWDataIsoName.str());
    stringstream nonIsoWDataNonIsoName;
    nonIsoWDataNonIsoName << nonIsoWDataNonIsoPrefix << *iRunEra << nonIsoWDataSuffix;
    nonIsoWDataNonIsoHaddInputFiles.push_back(nonIsoWDataNonIsoName.str());
    stringstream nonIsoWDataNonIsoReweightName;
    nonIsoWDataNonIsoReweightName << nonIsoWDataNonIsoReweightPrefix << *iRunEra;
    nonIsoWDataNonIsoReweightName << nonIsoWDataSuffix;
    nonIsoWDataNonIsoReweightHaddInputFiles.push_back(nonIsoWDataNonIsoReweightName.str());
    stringstream nonIsoWDataAllName;
    nonIsoWDataAllName << nonIsoWDataAllPrefix << *iRunEra << nonIsoWDataSuffix;
    nonIsoWDataAllHaddInputFiles.push_back(nonIsoWDataAllName.str());
  }
  haddCanvases(nonIsoWDataIsoHaddOutputFile, nonIsoWDataIsoHaddInputFiles, vector<float>(4, 1.0), 
  	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
  	       nullBlindHigh);
  haddCanvases(nonIsoWDataNonIsoHaddOutputFile, nonIsoWDataNonIsoHaddInputFiles, 
  	       vector<float>(4, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
  	       nullBlindLow, nullBlindHigh);
  haddCanvases(nonIsoWDataNonIsoReweightHaddOutputFile, nonIsoWDataNonIsoReweightHaddInputFiles, 
  	       vector<float>(4, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
  	       nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(nonIsoWDataAllHaddOutputFile, nonIsoWDataAllHaddInputFiles, 
		 vector<float>(4, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
		 nullBlindLow, nullBlindHigh);
  }
  
  //hadd single photon data samples
  cout << "...single photon data\n";
  string SinglePhotonDataSuffix(SinglePhotonDataVTag + fileExt);
  string SinglePhotonDataIsoPrefix(analysisFilePath + "SinglePhotonParkedData/analysis/muHadIsoAnalysis_SinglePhotonParked_Run2012D");
  string SinglePhotonDataIsoHaddOutputFile(SinglePhotonDataIsoPrefix + SinglePhotonDataSuffix);
  string SinglePhotonDataNonIsoPrefix(analysisFilePath + "SinglePhotonParkedData/analysis/muHadNonIsoAnalysis_SinglePhotonParked_Run2012D");
  string 
    SinglePhotonDataNonIsoHaddOutputFile(SinglePhotonDataNonIsoPrefix + SinglePhotonDataSuffix);
  string SinglePhotonDataNonIsoReweightPrefix(analysisFilePath + "SinglePhotonParkedData/analysis/muHadNonIsoReweightAnalysis_SinglePhotonParked_Run2012D");
  string SinglePhotonDataNonIsoReweightHaddOutputFile(SinglePhotonDataNonIsoReweightPrefix + 
						      SinglePhotonDataSuffix);
  string SinglePhotonDataAllPrefix(analysisFilePath + "SinglePhotonParkedData/analysis/muHadAnalysis_SinglePhotonParked_Run2012D");
  string SinglePhotonDataAllHaddOutputFile(SinglePhotonDataAllPrefix + SinglePhotonDataSuffix);
  vector<string> SinglePhotonDataIsoHaddInputFiles;
  vector<string> SinglePhotonDataNonIsoHaddInputFiles;
  vector<string> SinglePhotonDataNonIsoReweightHaddInputFiles;
  vector<string> SinglePhotonDataAllHaddInputFiles;
  for (unsigned int iJob = 0; iJob < 12; ++iJob) {
    stringstream SinglePhotonDataIsoName;
    SinglePhotonDataIsoName << SinglePhotonDataIsoPrefix << "_" << iJob << SinglePhotonDataSuffix;
    SinglePhotonDataIsoHaddInputFiles.push_back(SinglePhotonDataIsoName.str());
    stringstream SinglePhotonDataNonIsoName;
    SinglePhotonDataNonIsoName << SinglePhotonDataNonIsoPrefix << "_" << iJob;
    SinglePhotonDataNonIsoName << SinglePhotonDataSuffix;
    SinglePhotonDataNonIsoHaddInputFiles.push_back(SinglePhotonDataNonIsoName.str());
    stringstream SinglePhotonDataNonIsoReweightName;
    SinglePhotonDataNonIsoReweightName << SinglePhotonDataNonIsoReweightPrefix << "_" << iJob;
    SinglePhotonDataNonIsoReweightName << SinglePhotonDataSuffix;
    SinglePhotonDataNonIsoReweightHaddInputFiles.
      push_back(SinglePhotonDataNonIsoReweightName.str());
    stringstream SinglePhotonDataAllName;
    SinglePhotonDataAllName << SinglePhotonDataAllPrefix << "_" << iJob << SinglePhotonDataSuffix;
    SinglePhotonDataAllHaddInputFiles.push_back(SinglePhotonDataAllName.str());
  }
  haddCanvases(SinglePhotonDataIsoHaddOutputFile, SinglePhotonDataIsoHaddInputFiles, 
  	       vector<float>(12, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
  	       nullBlindLow, nullBlindHigh);
  haddCanvases(SinglePhotonDataNonIsoHaddOutputFile, SinglePhotonDataNonIsoHaddInputFiles, 
  	       vector<float>(12, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
  	       nullBlindLow, nullBlindHigh);
  haddCanvases(SinglePhotonDataNonIsoReweightHaddOutputFile, 
  	       SinglePhotonDataNonIsoReweightHaddInputFiles, vector<float>(12, 1.0), 
  	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
  	       nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(SinglePhotonDataAllHaddOutputFile, SinglePhotonDataAllHaddInputFiles, 
		 vector<float>(12, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
		 nullBlindLow, nullBlindHigh);
  }

  //"hadd" Wh1 sample just to get the formatting of the 2D plots the same
  cout << "...Wh1\n";
  string Wh1Suffix(Wh1SigVTag + fileExt);
  string Wh1IsoPrefix(analysisFilePath + "Wh1_Medium/muHadIsoAnalysis_Wh1");
  string Wh1IsoHaddOutputFile(Wh1IsoPrefix + "_hadd" + Wh1Suffix);
  string Wh1AllPrefix(analysisFilePath + "Wh1_Medium/muHadAnalysis_Wh1");
  string Wh1AllHaddOutputFile(Wh1AllPrefix + "_hadd" + Wh1Suffix);
  vector<string> Wh1IsoHaddInputFiles;
  vector<string> Wh1AllHaddInputFiles;
  stringstream Wh1IsoName;
  Wh1IsoName << Wh1IsoPrefix << Wh1Suffix;
  Wh1IsoHaddInputFiles.push_back(Wh1IsoName.str());
  stringstream Wh1AllName;
  Wh1AllName << Wh1AllPrefix << Wh1Suffix;
  Wh1AllHaddInputFiles.push_back(Wh1AllName.str());
  haddCanvases(Wh1IsoHaddOutputFile, Wh1IsoHaddInputFiles, vector<float>(1, Wh1Weight19p7InvFb), 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
	       nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(Wh1AllHaddOutputFile, Wh1AllHaddInputFiles, vector<float>(1, Wh1Weight19p7InvFb), 
		 canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
		 nullBlindHigh);
  }

  //"hadd" gg sample just to get the formatting of the 2D plots the same
  cout << "...gg fusion\n";
  string ggSuffix(ggSigVTag + fileExt);
  string ggIsoPrefix(analysisFilePath + "gg/muHadIsoAnalysis_gg");
  string ggIsoHaddOutputFile(ggIsoPrefix + "_hadd" + ggSuffix);
  string ggAllPrefix(analysisFilePath + "gg/muHadAnalysis_gg");
  string ggAllHaddOutputFile(ggAllPrefix + "_hadd" + ggSuffix);
  vector<string> ggIsoHaddInputFiles;
  vector<string> ggAllHaddInputFiles;
  stringstream ggIsoName;
  ggIsoName << ggIsoPrefix << ggSuffix;
  ggIsoHaddInputFiles.push_back(ggIsoName.str());
  stringstream ggAllName;
  ggAllName << ggAllPrefix << ggSuffix;
  ggAllHaddInputFiles.push_back(ggAllName.str());
  haddCanvases(ggIsoHaddOutputFile, ggIsoHaddInputFiles, vector<float>(1, ggWeight19p7InvFb), 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
	       nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(ggAllHaddOutputFile, ggAllHaddInputFiles, vector<float>(1, ggWeight19p7InvFb), 
		 canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
		 nullBlindHigh);
  }

  //hadd QCD Mu-enriched Pt-binned samples
  cout << "...muon-enriched QCD\n";
  string QCDSuffix(QCDVTag + fileExt);
  string QCDIsoPrefix(analysisFilePath + "QCD/analysis/muHadIsoAnalysis_QCD");
  string QCDIsoHaddOutputFile(QCDIsoPrefix + QCDSuffix);
  string QCDNonIsoPrefix(analysisFilePath + "QCD/analysis/muHadNonIsoAnalysis_QCD");
  string QCDNonIsoHaddOutputFile(QCDNonIsoPrefix + QCDSuffix);
  string QCDNonIsoReweightPrefix(analysisFilePath + 
				 "QCD/analysis/muHadNonIsoReweightAnalysis_QCD");
  string QCDNonIsoReweightHaddOutputFile(QCDNonIsoReweightPrefix + QCDSuffix);
  string QCDAllPrefix(analysisFilePath + "QCD/analysis/muHadAnalysis_QCD");
  string QCDAllHaddOutputFile(QCDAllPrefix + QCDSuffix);
  vector<string> QCDIsoHaddInputFiles;
  vector<string> QCDNonIsoHaddInputFiles;
  vector<string> QCDNonIsoReweightHaddInputFiles;
  vector<string> QCDAllHaddInputFiles;
  vector<string> ptBins;
  ptBins.push_back("_Pt-20to30");
  ptBins.push_back("_Pt-30to50");
  ptBins.push_back("_Pt-50to80");
  ptBins.push_back("_Pt-80to120");
  ptBins.push_back("_Pt-120to170");
  ptBins.push_back("_Pt-170to300");
  ptBins.push_back("_Pt-300to470");
  ptBins.push_back("_Pt-470to600");
  ptBins.push_back("_Pt-600to800");
  ptBins.push_back("_Pt-800to1000");
  ptBins.push_back("_Pt-1000");
   for (vector<string>::const_iterator iPtBin = ptBins.begin(); iPtBin != ptBins.end(); 
       ++iPtBin) {
    stringstream QCDIsoName;
    QCDIsoName << QCDIsoPrefix << *iPtBin << QCDSuffix;
    QCDIsoHaddInputFiles.push_back(QCDIsoName.str());
    stringstream QCDNonIsoName;
    QCDNonIsoName << QCDNonIsoPrefix << *iPtBin << QCDSuffix;
    QCDNonIsoHaddInputFiles.push_back(QCDNonIsoName.str());
    stringstream QCDNonIsoReweightName;
    QCDNonIsoReweightName << QCDNonIsoReweightPrefix << *iPtBin << QCDSuffix;
    QCDNonIsoReweightHaddInputFiles.push_back(QCDNonIsoReweightName.str());
    stringstream QCDAllName;
    QCDAllName << QCDAllPrefix << *iPtBin << QCDSuffix;
    QCDAllHaddInputFiles.push_back(QCDAllName.str());
  }
  haddCanvases(QCDIsoHaddOutputFile, QCDIsoHaddInputFiles, QCDRelXSecWeights, canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(QCDNonIsoHaddOutputFile, QCDNonIsoHaddInputFiles, QCDRelXSecWeights, canvasNames1D, 
  	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(QCDNonIsoReweightHaddOutputFile, QCDNonIsoReweightHaddInputFiles, 
  	       QCDRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
  	       nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(QCDAllHaddOutputFile, QCDAllHaddInputFiles, QCDRelXSecWeights, canvasNames1D, 
		 graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //hadd QCD b-enriched Pt-binned samples
  cout << "...b-enriched QCD\n";
  string QCDBSuffix(QCDBVTag + fileExt);
  string QCDBIsoPrefix(analysisFilePath + "QCDB/analysis/muHadIsoAnalysis_QCDB");
  string QCDBIsoHaddOutputFile(QCDBIsoPrefix + QCDBSuffix);
  string QCDBNonIsoPrefix(analysisFilePath + "QCDB/analysis/muHadNonIsoAnalysis_QCDB");
  string QCDBNonIsoHaddOutputFile(QCDBNonIsoPrefix + QCDBSuffix);
  string QCDBNonIsoReweightPrefix(analysisFilePath + 
				  "QCDB/analysis/muHadNonIsoReweightAnalysis_QCDB");
  string QCDBNonIsoReweightHaddOutputFile(QCDBNonIsoReweightPrefix + QCDBSuffix);
  string QCDBAllPrefix(analysisFilePath + "QCDB/analysis/muHadAnalysis_QCDB");
  string QCDBAllHaddOutputFile(QCDBAllPrefix + QCDBSuffix);
  vector<string> QCDBIsoHaddInputFiles;
  vector<string> QCDBNonIsoHaddInputFiles;
  vector<string> QCDBNonIsoReweightHaddInputFiles;
  vector<string> QCDBAllHaddInputFiles;
  vector<string> ptBinsB;
  ptBinsB.push_back("_Pt-15To30");
  ptBinsB.push_back("_Pt-30To50");
  ptBinsB.push_back("_Pt-50To150");
  ptBinsB.push_back("_Pt-150");
  for (vector<string>::const_iterator iPtBin = ptBinsB.begin(); iPtBin != ptBinsB.end(); 
       ++iPtBin) {
    stringstream QCDBIsoName;
    QCDBIsoName << QCDBIsoPrefix << *iPtBin << QCDBSuffix;
    QCDBIsoHaddInputFiles.push_back(QCDBIsoName.str());
    stringstream QCDBNonIsoName;
    QCDBNonIsoName << QCDBNonIsoPrefix << *iPtBin << QCDBSuffix;
    QCDBNonIsoHaddInputFiles.push_back(QCDBNonIsoName.str());
    stringstream QCDBNonIsoReweightName;
    QCDBNonIsoReweightName << QCDBNonIsoReweightPrefix << *iPtBin << QCDBSuffix;
    QCDBNonIsoReweightHaddInputFiles.push_back(QCDBNonIsoReweightName.str());
    stringstream QCDBAllName;
    QCDBAllName << QCDBAllPrefix << *iPtBin << QCDBSuffix;
    QCDBAllHaddInputFiles.push_back(QCDBAllName.str());
  }
  haddCanvases(QCDBIsoHaddOutputFile, QCDBIsoHaddInputFiles, QCDBRelXSecWeights, canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(QCDBNonIsoHaddOutputFile, QCDBNonIsoHaddInputFiles, QCDBRelXSecWeights, 
  	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
  	       nullBlindHigh);
  haddCanvases(QCDBNonIsoReweightHaddOutputFile, QCDBNonIsoReweightHaddInputFiles, 
  	       QCDBRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
  	       nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(QCDBAllHaddOutputFile, QCDBAllHaddInputFiles, QCDBRelXSecWeights, canvasNames1D, 
		 graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //hadd QCD bToMu-enriched Pt-binned samples
  cout << "...b-to-muon-enriched QCD\n";
  string QCDBMuSuffix(QCDBMuVTag + fileExt);
  string QCDBMuIsoPrefix(analysisFilePath + "QCDBMu/analysis/muHadIsoAnalysis_QCDBMu");
  string QCDBMuIsoHaddOutputFile(QCDBMuIsoPrefix + QCDBMuSuffix);
  string QCDBMuNonIsoPrefix(analysisFilePath + "QCDBMu/analysis/muHadNonIsoAnalysis_QCDBMu");
  string QCDBMuNonIsoHaddOutputFile(QCDBMuNonIsoPrefix + QCDBMuSuffix);
  string QCDBMuNonIsoReweightPrefix(analysisFilePath + 
				    "QCDBMu/analysis/muHadNonIsoReweightAnalysis_QCDBMu");
  string QCDBMuNonIsoReweightHaddOutputFile(QCDBMuNonIsoReweightPrefix + QCDBMuSuffix);
  string QCDBMuAllPrefix(analysisFilePath + "QCDBMu/analysis/muHadAnalysis_QCDBMu");
  string QCDBMuAllHaddOutputFile(QCDBMuAllPrefix + QCDBMuSuffix);
  vector<string> QCDBMuIsoHaddInputFiles;
  vector<string> QCDBMuNonIsoHaddInputFiles;
  vector<string> QCDBMuNonIsoReweightHaddInputFiles;
  vector<string> QCDBMuAllHaddInputFiles;
  vector<string> ptBinsBMu;
  ptBinsBMu.push_back("_pt15to30");
  ptBinsBMu.push_back("_pt30to50");
  ptBinsBMu.push_back("_pt50to150");
  ptBinsBMu.push_back("_pt150");
  for (vector<string>::const_iterator iPtBin = ptBinsBMu.begin(); iPtBin != ptBinsBMu.end(); 
       ++iPtBin) {
    stringstream QCDBMuIsoName;
    QCDBMuIsoName << QCDBMuIsoPrefix << *iPtBin << QCDBMuSuffix;
    QCDBMuIsoHaddInputFiles.push_back(QCDBMuIsoName.str());
    stringstream QCDBMuNonIsoName;
    QCDBMuNonIsoName << QCDBMuNonIsoPrefix << *iPtBin << QCDBMuSuffix;
    QCDBMuNonIsoHaddInputFiles.push_back(QCDBMuNonIsoName.str());
    stringstream QCDBMuNonIsoReweightName;
    QCDBMuNonIsoReweightName << QCDBMuNonIsoReweightPrefix << *iPtBin << QCDBMuSuffix;
    QCDBMuNonIsoReweightHaddInputFiles.push_back(QCDBMuNonIsoReweightName.str());
    stringstream QCDBMuAllName;
    QCDBMuAllName << QCDBMuAllPrefix << *iPtBin << QCDBMuSuffix;
    QCDBMuAllHaddInputFiles.push_back(QCDBMuAllName.str());
  }
  haddCanvases(QCDBMuIsoHaddOutputFile, QCDBMuIsoHaddInputFiles, QCDBMuRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
	       nullBlindHigh);
  haddCanvases(QCDBMuNonIsoHaddOutputFile, QCDBMuNonIsoHaddInputFiles, QCDBMuRelXSecWeights, 
  	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
  	       nullBlindHigh);
  haddCanvases(QCDBMuNonIsoReweightHaddOutputFile, QCDBMuNonIsoReweightHaddInputFiles, 
  	       QCDBMuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
  	       nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(QCDBMuAllHaddOutputFile, QCDBMuAllHaddInputFiles, QCDBMuRelXSecWeights, 
		 canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
		 nullBlindHigh);
  }

  //hadd Drell-Yan+jets ml+l- binned samples
  cout << "...Drell-Yan\n";
  string DYJetsToLLSuffix(DYJetsToLLVTag + fileExt);
  string DYJetsToLLIsoPrefix(analysisFilePath + "DYJetsToLL/analysis/muHadIsoAnalysis_DYJetsToLL");
  string DYJetsToLLIsoHaddOutputFile(DYJetsToLLIsoPrefix + DYJetsToLLSuffix);
  string DYJetsToLLNonIsoPrefix(analysisFilePath + 
				"DYJetsToLL/analysis/muHadNonIsoAnalysis_DYJetsToLL");
  string DYJetsToLLNonIsoHaddOutputFile(DYJetsToLLNonIsoPrefix + DYJetsToLLSuffix);
  string 
    DYJetsToLLNonIsoReweightPrefix(analysisFilePath + 
				   "DYJetsToLL/analysis/muHadNonIsoReweightAnalysis_DYJetsToLL");
  string DYJetsToLLNonIsoReweightHaddOutputFile(DYJetsToLLNonIsoReweightPrefix + DYJetsToLLSuffix);
  string DYJetsToLLAllPrefix(analysisFilePath + "DYJetsToLL/analysis/muHadAnalysis_DYJetsToLL");
  string DYJetsToLLAllHaddOutputFile(DYJetsToLLAllPrefix + DYJetsToLLSuffix);
  vector<string> DYJetsToLLIsoHaddInputFiles;
  vector<string> DYJetsToLLNonIsoHaddInputFiles;
  vector<string> DYJetsToLLNonIsoReweightHaddInputFiles;
  vector<string> DYJetsToLLAllHaddInputFiles;
  vector<string> massBins;
  massBins.push_back("_M-10To50");
  massBins.push_back("_M-50");
  for (vector<string>::const_iterator iMassBin = massBins.begin(); iMassBin != massBins.end(); 
       ++iMassBin) {
    stringstream DYJetsToLLIsoName;
    DYJetsToLLIsoName << DYJetsToLLIsoPrefix << *iMassBin << DYJetsToLLSuffix;
    DYJetsToLLIsoHaddInputFiles.push_back(DYJetsToLLIsoName.str());
    stringstream DYJetsToLLNonIsoName;
    DYJetsToLLNonIsoName << DYJetsToLLNonIsoPrefix << *iMassBin << DYJetsToLLSuffix;
    DYJetsToLLNonIsoHaddInputFiles.push_back(DYJetsToLLNonIsoName.str());
    stringstream DYJetsToLLNonIsoReweightName;
    DYJetsToLLNonIsoReweightName << DYJetsToLLNonIsoReweightPrefix << *iMassBin;
    DYJetsToLLNonIsoReweightName << DYJetsToLLSuffix;
    DYJetsToLLNonIsoReweightHaddInputFiles.push_back(DYJetsToLLNonIsoReweightName.str());
    stringstream DYJetsToLLAllName;
    DYJetsToLLAllName << DYJetsToLLAllPrefix << *iMassBin << DYJetsToLLSuffix;
    DYJetsToLLAllHaddInputFiles.push_back(DYJetsToLLAllName.str());
  }
  haddCanvases(DYJetsToLLIsoHaddOutputFile, DYJetsToLLIsoHaddInputFiles, DYJetsToLLRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
	       nullBlindHigh);
  haddCanvases(DYJetsToLLNonIsoHaddOutputFile, DYJetsToLLNonIsoHaddInputFiles, 
  	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
  	       nullBlindLow, nullBlindHigh);
  haddCanvases(DYJetsToLLNonIsoReweightHaddOutputFile, DYJetsToLLNonIsoReweightHaddInputFiles, 
  	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
  	       nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(DYJetsToLLAllHaddOutputFile, DYJetsToLLAllHaddInputFiles, 
		 DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
		 graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //"hadd" ttbar sample just to get the formatting of the 2D plots the same
  cout << "...ttbar\n";
  string TTJetsSuffix(TTJetsVTag + fileExt);
  string TTJetsIsoPrefix(analysisFilePath + "TTJets/analysis/muHadIsoAnalysis_TTJets");
  string TTJetsIsoHaddOutputFile(TTJetsIsoPrefix + "_hadd" + TTJetsSuffix);
  string TTJetsNonIsoPrefix(analysisFilePath + "TTJets/analysis/muHadNonIsoAnalysis_TTJets");
  string TTJetsNonIsoHaddOutputFile(TTJetsNonIsoPrefix + "_hadd" + TTJetsSuffix);
  string TTJetsNonIsoReweightPrefix(analysisFilePath + 
				    "TTJets/analysis/muHadNonIsoReweightAnalysis_TTJets");
  string TTJetsNonIsoReweightHaddOutputFile(TTJetsNonIsoReweightPrefix + "_hadd" + TTJetsSuffix);
  string TTJetsAllPrefix(analysisFilePath + "TTJets/analysis/muHadAnalysis_TTJets");
  string TTJetsAllHaddOutputFile(TTJetsAllPrefix + "_hadd" + TTJetsSuffix);
  vector<string> TTJetsIsoHaddInputFiles;
  vector<string> TTJetsNonIsoHaddInputFiles;
  vector<string> TTJetsNonIsoReweightHaddInputFiles;
  vector<string> TTJetsAllHaddInputFiles;
  stringstream TTJetsIsoName;
  TTJetsIsoName << TTJetsIsoPrefix << TTJetsSuffix;
  TTJetsIsoHaddInputFiles.push_back(TTJetsIsoName.str());
  stringstream TTJetsNonIsoName;
  TTJetsNonIsoName << TTJetsNonIsoPrefix << TTJetsSuffix;
  TTJetsNonIsoHaddInputFiles.push_back(TTJetsNonIsoName.str());
  stringstream TTJetsNonIsoReweightName;
  TTJetsNonIsoReweightName << TTJetsNonIsoReweightPrefix << TTJetsSuffix;
  TTJetsNonIsoReweightHaddInputFiles.push_back(TTJetsNonIsoReweightName.str());
  stringstream TTJetsAllName;
  TTJetsAllName << TTJetsAllPrefix << TTJetsSuffix;
  TTJetsAllHaddInputFiles.push_back(TTJetsAllName.str());
  haddCanvases(TTJetsIsoHaddOutputFile, TTJetsIsoHaddInputFiles, 
	       vector<float>(1, 3.54800726562391), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(TTJetsNonIsoHaddOutputFile, TTJetsNonIsoHaddInputFiles, 
  	       vector<float>(1, 3.54800726562391), canvasNames1D, graphNames1D, 
  	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(TTJetsNonIsoReweightHaddOutputFile, TTJetsNonIsoReweightHaddInputFiles, 
  	       vector<float>(1, 3.54800726562391), canvasNames1D, graphNames1D, 
  	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(TTJetsAllHaddOutputFile, TTJetsAllHaddInputFiles, 
		 vector<float>(1, 3.54800726562391), canvasNames1D, graphNames1D, 
		 canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //hadd single top samples
  cout << "...single top\n";
  string TSuffix(TVTag + fileExt);
  string TIsoPrefix(analysisFilePath + "SingleTop/analysis/muHadIsoAnalysis_T");
  string TIsoHaddOutputFile(TIsoPrefix + TSuffix);
  string TNonIsoPrefix(analysisFilePath + "SingleTop/analysis/muHadNonIsoAnalysis_T");
  string TNonIsoHaddOutputFile(TNonIsoPrefix + TSuffix);
  string TNonIsoReweightPrefix(analysisFilePath + 
			       "SingleTop/analysis/muHadNonIsoReweightAnalysis_T");
  string TNonIsoReweightHaddOutputFile(TNonIsoReweightPrefix + TSuffix);
  string TAllPrefix(analysisFilePath + "SingleTop/analysis/muHadAnalysis_T");
  string TAllHaddOutputFile(TAllPrefix + TSuffix);
  vector<string> TIsoHaddInputFiles;
  vector<string> TNonIsoHaddInputFiles;
  vector<string> TNonIsoReweightHaddInputFiles;
  vector<string> TAllHaddInputFiles;
  vector<string> singleTopSamples;
  singleTopSamples.push_back("_s-channel");
  singleTopSamples.push_back("bar_s-channel");
  singleTopSamples.push_back("_t-channel");
  singleTopSamples.push_back("bar_t-channel");
  for (vector<string>::const_iterator iSample = singleTopSamples.begin(); 
       iSample != singleTopSamples.end(); ++iSample) {
    stringstream TIsoName;
    TIsoName << TIsoPrefix << *iSample << TSuffix;
    TIsoHaddInputFiles.push_back(TIsoName.str());
    stringstream TNonIsoName;
    TNonIsoName << TNonIsoPrefix << *iSample << TSuffix;
    TNonIsoHaddInputFiles.push_back(TNonIsoName.str());
    stringstream TNonIsoReweightName;
    TNonIsoReweightName << TNonIsoReweightPrefix << *iSample << TSuffix;
    TNonIsoReweightHaddInputFiles.push_back(TNonIsoReweightName.str());
    stringstream TAllName;
    TAllName << TAllPrefix << *iSample << TSuffix;
    TAllHaddInputFiles.push_back(TAllName.str());
  }
  haddCanvases(TIsoHaddOutputFile, TIsoHaddInputFiles, TRelXSecWeights, canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(TNonIsoHaddOutputFile, TNonIsoHaddInputFiles, TRelXSecWeights, canvasNames1D, 
  	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(TNonIsoReweightHaddOutputFile, TNonIsoReweightHaddInputFiles, TRelXSecWeights, 
  	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
  	       nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(TAllHaddOutputFile, TAllHaddInputFiles, TRelXSecWeights, canvasNames1D, 
		 graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //hadd W+>=1 jet samples
  cout << "...W+>=1 jet\n";
  string WNJetsToLNuSuffix("JetsToLNu" + WNJetsToLNuVTag + fileExt);
  string WNJetsToLNuIsoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadIsoAnalysis_W");
  string WNJetsToLNuIsoHaddOutputFile(WNJetsToLNuIsoPrefix + "N" + WNJetsToLNuSuffix);
  string WNJetsToLNuNonIsoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadNonIsoAnalysis_W");
  string WNJetsToLNuNonIsoHaddOutputFile(WNJetsToLNuNonIsoPrefix + "N" + WNJetsToLNuSuffix);
  string WNJetsToLNuNonIsoReweightPrefix(analysisFilePath + 
					 "WNJetsToLNu/analysis/muHadNonIsoReweightAnalysis_W");
  string WNJetsToLNuNonIsoReweightHaddOutputFile(WNJetsToLNuNonIsoReweightPrefix + "N" + 
						 WNJetsToLNuSuffix);
  string WNJetsToLNuAllTauPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadAnalysis_W");
  string WNJetsToLNuAllTauHaddOutputFile(WNJetsToLNuAllTauPrefix + "N" + WNJetsToLNuSuffix);
  vector<string> WNJetsToLNuIsoHaddInputFiles;
  vector<string> WNJetsToLNuNonIsoHaddInputFiles;
  vector<string> WNJetsToLNuNonIsoReweightHaddInputFiles;
  vector<string> WNJetsToLNuAllTauHaddInputFiles;
  for (unsigned int iNJets = 1; iNJets <= 4; ++iNJets) {
    stringstream WNJetsToLNuIsoName;
    WNJetsToLNuIsoName << WNJetsToLNuIsoPrefix << iNJets << WNJetsToLNuSuffix;
    WNJetsToLNuIsoHaddInputFiles.push_back(WNJetsToLNuIsoName.str());
    stringstream WNJetsToLNuNonIsoName;
    WNJetsToLNuNonIsoName << WNJetsToLNuNonIsoPrefix << iNJets << WNJetsToLNuSuffix;
    WNJetsToLNuNonIsoHaddInputFiles.push_back(WNJetsToLNuNonIsoName.str());
    stringstream WNJetsToLNuNonIsoReweightName;
    WNJetsToLNuNonIsoReweightName << WNJetsToLNuNonIsoReweightPrefix << iNJets;
    WNJetsToLNuNonIsoReweightName << WNJetsToLNuSuffix;
    WNJetsToLNuNonIsoReweightHaddInputFiles.push_back(WNJetsToLNuNonIsoReweightName.str());
    stringstream WNJetsToLNuAllTauName;
    WNJetsToLNuAllTauName << WNJetsToLNuAllTauPrefix << iNJets << WNJetsToLNuSuffix;
    WNJetsToLNuAllTauHaddInputFiles.push_back(WNJetsToLNuAllTauName.str());
  }
  haddCanvases(WNJetsToLNuIsoHaddOutputFile, WNJetsToLNuIsoHaddInputFiles, 
	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
	       graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WNJetsToLNuNonIsoHaddOutputFile, WNJetsToLNuNonIsoHaddInputFiles, 
  	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
  	       graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WNJetsToLNuNonIsoReweightHaddOutputFile, WNJetsToLNuNonIsoReweightHaddInputFiles, 
  	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
  	       graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(WNJetsToLNuAllTauHaddOutputFile, WNJetsToLNuAllTauHaddInputFiles, 
		 WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
		 graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //"hadd" W + bbbar sample just to get the formatting of the 2D plots the same
  cout << "...Wbb\n";
  string WbbSuffix(WbbVTag + fileExt);
  string WbbIsoPrefix(analysisFilePath + "Wbb/analysis/muHadIsoAnalysis_Wbb");
  string WbbIsoHaddOutputFile(WbbIsoPrefix + "_hadd" + WbbSuffix);
  string WbbNonIsoPrefix(analysisFilePath + "Wbb/analysis/muHadNonIsoAnalysis_Wbb");
  string WbbNonIsoHaddOutputFile(WbbNonIsoPrefix + "_hadd" + WbbSuffix);
  string 
    WbbNonIsoReweightPrefix(analysisFilePath + "Wbb/analysis/muHadNonIsoReweightAnalysis_Wbb");
  string WbbNonIsoReweightHaddOutputFile(WbbNonIsoReweightPrefix + "_hadd" + WbbSuffix);
  string WbbAllPrefix(analysisFilePath + "Wbb/analysis/muHadAnalysis_Wbb");
  string WbbAllHaddOutputFile(WbbAllPrefix + "_hadd" + WbbSuffix);
  vector<string> WbbIsoHaddInputFiles;
  vector<string> WbbNonIsoHaddInputFiles;
  vector<string> WbbNonIsoReweightHaddInputFiles;
  vector<string> WbbAllHaddInputFiles;
  stringstream WbbIsoName;
  WbbIsoName << WbbIsoPrefix << WbbSuffix;
  WbbIsoHaddInputFiles.push_back(WbbIsoName.str());
  stringstream WbbNonIsoName;
  WbbNonIsoName << WbbNonIsoPrefix << WbbSuffix;
  WbbNonIsoHaddInputFiles.push_back(WbbNonIsoName.str());
  stringstream WbbNonIsoReweightName;
  WbbNonIsoReweightName << WbbNonIsoReweightPrefix << WbbSuffix;
  WbbNonIsoReweightHaddInputFiles.push_back(WbbNonIsoReweightName.str());
  stringstream WbbAllName;
  WbbAllName << WbbAllPrefix << WbbSuffix;
  WbbAllHaddInputFiles.push_back(WbbAllName.str());
  haddCanvases(WbbIsoHaddOutputFile, WbbIsoHaddInputFiles, vector<float>(1, 0.11667483305847), 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
	       nullBlindHigh);
  haddCanvases(WbbNonIsoHaddOutputFile, WbbNonIsoHaddInputFiles, 
  	       vector<float>(1, 0.11667483305847), canvasNames1D, graphNames1D, canvasNames2D, 
  	       graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WbbNonIsoReweightHaddOutputFile, WbbNonIsoReweightHaddInputFiles, 
  	       vector<float>(1, 0.11667483305847), canvasNames1D, graphNames1D, canvasNames2D, 
  	       graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(WbbAllHaddOutputFile, WbbAllHaddInputFiles, vector<float>(1, 0.11667483305847), 
		 canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
		 nullBlindHigh);
  }

  //"hadd" W+jets sample just to get the formatting of the 2D plots the same
  cout << "...W+jets\n";
  string WJetsToLNuSuffix(WJetsToLNuVTag + fileExt);
  string WJetsToLNuIsoPrefix(analysisFilePath + "WJetsToLNu/analysis/muHadIsoAnalysis_WJetsToLNu");
  string WJetsToLNuIsoHaddOutputFile(WJetsToLNuIsoPrefix + "_hadd" + WJetsToLNuSuffix);
  string WJetsToLNuNonIsoPrefix(analysisFilePath + 
				"WJetsToLNu/analysis/muHadNonIsoAnalysis_WJetsToLNu");
  string WJetsToLNuNonIsoHaddOutputFile(WJetsToLNuNonIsoPrefix + "_hadd" + WJetsToLNuSuffix);
  string 
    WJetsToLNuNonIsoReweightPrefix(analysisFilePath + 
				   "WJetsToLNu/analysis/muHadNonIsoReweightAnalysis_WJetsToLNu");
  string WJetsToLNuNonIsoReweightHaddOutputFile(WJetsToLNuNonIsoReweightPrefix + "_hadd" + 
						WJetsToLNuSuffix);
  string WJetsToLNuAllPrefix(analysisFilePath + "WJetsToLNu/analysis/muHadAnalysis_WJetsToLNu");
  string WJetsToLNuAllHaddOutputFile(WJetsToLNuAllPrefix + "_hadd" + WJetsToLNuSuffix);
  vector<string> WJetsToLNuIsoHaddInputFiles;
  vector<string> WJetsToLNuNonIsoHaddInputFiles;
  vector<string> WJetsToLNuNonIsoReweightHaddInputFiles;
  vector<string> WJetsToLNuAllHaddInputFiles;
  stringstream WJetsToLNuIsoName;
  WJetsToLNuIsoName << WJetsToLNuIsoPrefix << WJetsToLNuSuffix;
  WJetsToLNuIsoHaddInputFiles.push_back(WJetsToLNuIsoName.str());
  stringstream WJetsToLNuNonIsoName;
  WJetsToLNuNonIsoName << WJetsToLNuNonIsoPrefix << WJetsToLNuSuffix;
  WJetsToLNuNonIsoHaddInputFiles.push_back(WJetsToLNuNonIsoName.str());
  stringstream WJetsToLNuNonIsoReweightName;
  WJetsToLNuNonIsoReweightName << WJetsToLNuNonIsoReweightPrefix << WJetsToLNuSuffix;
  WJetsToLNuNonIsoReweightHaddInputFiles.push_back(WJetsToLNuNonIsoReweightName.str());
  stringstream WJetsToLNuAllName;
  WJetsToLNuAllName << WJetsToLNuAllPrefix << WJetsToLNuSuffix;
  WJetsToLNuAllHaddInputFiles.push_back(WJetsToLNuAllName.str());
  haddCanvases(WJetsToLNuIsoHaddOutputFile, WJetsToLNuIsoHaddInputFiles, 
	       vector<float>(1, 40.174179542426), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WJetsToLNuNonIsoHaddOutputFile, WJetsToLNuNonIsoHaddInputFiles, 
  	       vector<float>(1, 40.174179542426), canvasNames1D, graphNames1D, 
  	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WJetsToLNuNonIsoReweightHaddOutputFile, WJetsToLNuNonIsoReweightHaddInputFiles, 
  	       vector<float>(1, 40.174179542426), canvasNames1D, graphNames1D, 
  	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(WJetsToLNuAllHaddOutputFile, WJetsToLNuAllHaddInputFiles, 
		 vector<float>(1, 40.174179542426), canvasNames1D, graphNames1D, 
		 canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //"hadd" WZ sample just to get the formatting of the 2D plots the same
  cout << "...WZ\n";
  string WZSuffix(WZVTag + fileExt);
  string WZIsoPrefix(analysisFilePath + "WZ/analysis/muHadIsoAnalysis_WZ");
  string WZIsoHaddOutputFile(WZIsoPrefix + "_hadd" + WZSuffix);
  string WZNonIsoPrefix(analysisFilePath + "WZ/analysis/muHadNonIsoAnalysis_WZ");
  string WZNonIsoHaddOutputFile(WZNonIsoPrefix + "_hadd" + WZSuffix);
  string WZNonIsoReweightPrefix(analysisFilePath + "WZ/analysis/muHadNonIsoReweightAnalysis_WZ");
  string WZNonIsoReweightHaddOutputFile(WZNonIsoReweightPrefix + "_hadd" + WZSuffix);
  string WZAllPrefix(analysisFilePath + "WZ/analysis/muHadAnalysis_WZ");
  string WZAllHaddOutputFile(WZAllPrefix + "_hadd" + WZSuffix);
  vector<string> WZIsoHaddInputFiles;
  vector<string> WZNonIsoHaddInputFiles;
  vector<string> WZNonIsoReweightHaddInputFiles;
  vector<string> WZAllHaddInputFiles;
  stringstream WZIsoName;
  WZIsoName << WZIsoPrefix << WZSuffix;
  WZIsoHaddInputFiles.push_back(WZIsoName.str());
  stringstream WZNonIsoName;
  WZNonIsoName << WZNonIsoPrefix << WZSuffix;
  WZNonIsoHaddInputFiles.push_back(WZNonIsoName.str());
  stringstream WZNonIsoReweightName;
  WZNonIsoReweightName << WZNonIsoReweightPrefix << WZSuffix;
  WZNonIsoReweightHaddInputFiles.push_back(WZNonIsoReweightName.str());
  stringstream WZAllName;
  WZAllName << WZAllPrefix << WZSuffix;
  WZAllHaddInputFiles.push_back(WZAllName.str());
  haddCanvases(WZIsoHaddOutputFile, WZIsoHaddInputFiles, 
	       vector<float>(1, 0.0667569497737973), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WZNonIsoHaddOutputFile, WZNonIsoHaddInputFiles, 
  	       vector<float>(1, 0.0667569497737973), canvasNames1D, graphNames1D, 
  	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WZNonIsoReweightHaddOutputFile, WZNonIsoReweightHaddInputFiles, 
  	       vector<float>(1, 0.0667569497737973), canvasNames1D, graphNames1D, 
  	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(WZAllHaddOutputFile, WZAllHaddInputFiles, 
		 vector<float>(1, 0.0667569497737973), canvasNames1D, graphNames1D, 
		 canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //"hadd" ZZ sample just to get the formatting of the 2D plots the same
  cout << "...ZZ\n";
  string ZZSuffix(ZZVTag + fileExt);
  string ZZIsoPrefix(analysisFilePath + "ZZ/analysis/muHadIsoAnalysis_ZZ");
  string ZZIsoHaddOutputFile(ZZIsoPrefix + "_hadd" + ZZSuffix);
  string ZZNonIsoPrefix(analysisFilePath + "ZZ/analysis/muHadNonIsoAnalysis_ZZ");
  string ZZNonIsoHaddOutputFile(ZZNonIsoPrefix + "_hadd" + ZZSuffix);
  string ZZNonIsoReweightPrefix(analysisFilePath + "ZZ/analysis/muHadNonIsoReweightAnalysis_ZZ");
  string ZZNonIsoReweightHaddOutputFile(ZZNonIsoReweightPrefix + "_hadd" + ZZSuffix);
  string ZZAllPrefix(analysisFilePath + "ZZ/analysis/muHadAnalysis_ZZ");
  string ZZAllHaddOutputFile(ZZAllPrefix + "_hadd" + ZZSuffix);
  vector<string> ZZIsoHaddInputFiles;
  vector<string> ZZNonIsoHaddInputFiles;
  vector<string> ZZNonIsoReweightHaddInputFiles;
  vector<string> ZZAllHaddInputFiles;
  stringstream ZZIsoName;
  ZZIsoName << ZZIsoPrefix << ZZSuffix;
  ZZIsoHaddInputFiles.push_back(ZZIsoName.str());
  stringstream ZZNonIsoName;
  ZZNonIsoName << ZZNonIsoPrefix << ZZSuffix;
  ZZNonIsoHaddInputFiles.push_back(ZZNonIsoName.str());
  stringstream ZZNonIsoReweightName;
  ZZNonIsoReweightName << ZZNonIsoReweightPrefix << ZZSuffix;
  ZZNonIsoReweightHaddInputFiles.push_back(ZZNonIsoReweightName.str());
  stringstream ZZAllName;
  ZZAllName << ZZAllPrefix << ZZSuffix;
  ZZAllHaddInputFiles.push_back(ZZAllName.str());
  haddCanvases(ZZIsoHaddOutputFile, ZZIsoHaddInputFiles, 
	       vector<float>(1, 0.0377509625152821), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(ZZNonIsoHaddOutputFile, ZZNonIsoHaddInputFiles, 
  	       vector<float>(1, 0.0377509625152821), canvasNames1D, graphNames1D, 
  	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(ZZNonIsoReweightHaddOutputFile, ZZNonIsoReweightHaddInputFiles, 
  	       vector<float>(1, 0.0377509625152821), canvasNames1D, graphNames1D, 
  	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(ZZAllHaddOutputFile, ZZAllHaddInputFiles, 
		 vector<float>(1, 0.0377509625152821), canvasNames1D, graphNames1D, 
		 canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //"hadd" WW sample just to get the formatting of the 2D plots the same
  cout << "...WW\n";
  string WWSuffix(WWVTag + fileExt);
  string WWIsoPrefix(analysisFilePath + "WW/analysis/muHadIsoAnalysis_WW");
  string WWIsoHaddOutputFile(WWIsoPrefix + "_hadd" + WWSuffix);
  string WWNonIsoPrefix(analysisFilePath + "WW/analysis/muHadNonIsoAnalysis_WW");
  string WWNonIsoHaddOutputFile(WWNonIsoPrefix + "_hadd" + WWSuffix);
  string WWNonIsoReweightPrefix(analysisFilePath + "WW/analysis/muHadNonIsoReweightAnalysis_WW");
  string WWNonIsoReweightHaddOutputFile(WWNonIsoReweightPrefix + "_hadd" + WWSuffix);
  string WWAllPrefix(analysisFilePath + "WW/analysis/muHadAnalysis_WW");
  string WWAllHaddOutputFile(WWAllPrefix + "_hadd" + WWSuffix);
  vector<string> WWIsoHaddInputFiles;
  vector<string> WWNonIsoHaddInputFiles;
  vector<string> WWNonIsoReweightHaddInputFiles;
  vector<string> WWAllHaddInputFiles;
  stringstream WWIsoName;
  WWIsoName << WWIsoPrefix << WWSuffix;
  WWIsoHaddInputFiles.push_back(WWIsoName.str());
  stringstream WWNonIsoName;
  WWNonIsoName << WWNonIsoPrefix << WWSuffix;
  WWNonIsoHaddInputFiles.push_back(WWNonIsoName.str());
  stringstream WWNonIsoReweightName;
  WWNonIsoReweightName << WWNonIsoReweightPrefix << WWSuffix;
  WWNonIsoReweightHaddInputFiles.push_back(WWNonIsoReweightName.str());
  stringstream WWAllName;
  WWAllName << WWAllPrefix << WWSuffix;
  WWAllHaddInputFiles.push_back(WWAllName.str());
  haddCanvases(WWIsoHaddOutputFile, WWIsoHaddInputFiles, 
	       vector<float>(1, 0.124167251024691), canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WWNonIsoHaddOutputFile, WWNonIsoHaddInputFiles, 
  	       vector<float>(1, 0.124167251024691), canvasNames1D, 
  	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WWNonIsoReweightHaddOutputFile, WWNonIsoReweightHaddInputFiles, 
  	       vector<float>(1, 0.124167251024691), canvasNames1D, 
  	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(WWAllHaddOutputFile, WWAllHaddInputFiles, 
		 vector<float>(1, 0.124167251024691), canvasNames1D, 
		 graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }

  cout << endl;

  //compare MC signal to background
  string sigVsBkgOutputFile(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + 
			    tag19p7InvFb + outputVTag + fileExt);
  string sigVsBkgOutputFile1(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + tag1 + 
			     outputVTag + fileExt);
  string sigVsBkgOutputFileNoHPSIsoCut(analysisFilePath + "results/sigVsBkg_muHadAnalysis" + 
				       tag19p7InvFb + outputVTag + fileExt);
  string sigVsBkgOutputFileNoHPSIsoCutNorm1(analysisFilePath + "results/sigVsBkg_muHadAnalysis" + 
					    tag1 + outputVTag + fileExt);
  vector<string> sigVsBkgInputFiles;
  sigVsBkgInputFiles.push_back(Wh1IsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(ggIsoHaddOutputFile);
//   sigVsBkgInputFiles.push_back(dataIsoHaddOutputFile);
//   sigVsBkgInputFiles.push_back(QCDIsoHaddOutputFile);
//   sigVsBkgInputFiles.push_back(QCDBIsoHaddOutputFile);
//   sigVsBkgInputFiles.push_back(QCDBMuIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(DYJetsToLLIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(TTJetsIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(TIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(WNJetsToLNuIsoHaddOutputFile);
//   sigVsBkgInputFiles.push_back(WbbIsoHaddOutputFile);
//   sigVsBkgInputFiles.push_back(WJetsToLNuIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(WZIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(ZZIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(WWIsoHaddOutputFile);
  std::reverse(sigVsBkgInputFiles.begin() + 2, sigVsBkgInputFiles.end());
  vector<string> sigVsBkgNoHPSIsoCutInputFiles;
  sigVsBkgNoHPSIsoCutInputFiles.push_back(Wh1AllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(ggAllHaddOutputFile);
//   sigVsBkgNoHPSIsoCutInputFiles.push_back(QCDAllHaddOutputFile);
//   sigVsBkgNoHPSIsoCutInputFiles.push_back(QCDBAllHaddOutputFile);
//   sigVsBkgNoHPSIsoCutInputFiles.push_back(QCDBMuAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(DYJetsToLLAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(TTJetsAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(TAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(WNJetsToLNuAllTauHaddOutputFile);
//   sigVsBkgNoHPSIsoCutInputFiles.push_back(WbbAllHaddOutputFile);
//   sigVsBkgNoHPSIsoCutInputFiles.push_back(WJetsToLNuAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(WZAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(ZZAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(WWAllHaddOutputFile);
  std::reverse(sigVsBkgNoHPSIsoCutInputFiles.begin() + 2, sigVsBkgNoHPSIsoCutInputFiles.end());
  cout << "Plot signal vs. background normalized to data luminosity\n---\n";
  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile, sigVsBkgInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders19p7InvFb, 
					colors, styles, legendEntriesSigBkg, weightsSigBkg, 
					setLogY, drawStack, sigBkg);
  cout << "\nPlot signal vs. background normalized to 1\n---\n";
  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile1, sigVsBkgInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders1, colors, 
					styles, legendEntriesSigBkg, weights1, setLinY, drawSame, 
					sigBkg);
  if (doNoHPSIsoCut) {
    cout << "\nPlot signal vs. background normalized to data luminosity, ";
    cout << "no cut on tau isolation\n---\n";
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFileNoHPSIsoCut, 
					  sigVsBkgNoHPSIsoCutInputFiles, canvasNames1D, 
					  graphNames1D, legendHeaders19p7InvFb, colors, styles, 
					  legendEntriesSigBkg, weightsSigBkg, setLogY, drawStack, 
					  sigBkg);
    cout << "\nPlot signal vs. background normalized to 1, no cut on tau isolation\n---\n";
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFileNoHPSIsoCutNorm1, 
					  sigVsBkgNoHPSIsoCutInputFiles, canvasNames1D, 
					  graphNames1D, legendHeaders1, colors, styles, 
					  legendEntriesSigBkg, weights1, setLinY, drawSame, 
					  sigBkg);
  }

  //compare data to MC in control region and compute data - MC for data-driven QCD shape
  string dataVsMCOutputFile(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + 
			    tag19p7InvFb + outputVTag + fileExt);
  string dataVsMCOutputDiff(analysisFilePath + "results/dataVsMC_muHadNonIsoDifference" + 
			    tag19p7InvFb + outputVTag + fileExt);
  string dataVsMCReweightOutputFile(analysisFilePath + 
				    "results/dataVsMC_muHadNonIsoReweightAnalysis" + 
				    tag19p7InvFb + outputVTag + fileExt);
  vector<string> dataVsMCInputFiles;
  dataVsMCInputFiles.push_back(dataNonIsoHaddOutputFile);
//   dataVsMCInputFiles.push_back(QCDNonIsoHaddOutputFile);
//   dataVsMCInputFiles.push_back(QCDBNonIsoHaddOutputFile);
//   dataVsMCInputFiles.push_back(QCDBMuNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(DYJetsToLLNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(TTJetsNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(TNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(WNJetsToLNuNonIsoHaddOutputFile);
//   dataVsMCInputFiles.push_back(WbbNonIsoHaddOutputFile);
//   dataVsMCInputFiles.push_back(WJetsToLNuNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(WZNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(ZZNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(WWNonIsoHaddOutputFile);
  std::reverse(dataVsMCInputFiles.begin() + 1, dataVsMCInputFiles.end());
  vector<string> dataVsMCReweightInputFiles;
  dataVsMCReweightInputFiles.push_back(dataNonIsoReweightHaddOutputFile);
//   dataVsMCReweightInputFiles.push_back(QCDNonIsoReweightHaddOutputFile);
//   dataVsMCReweightInputFiles.push_back(QCDBNonIsoReweightHaddOutputFile);
//   dataVsMCReweightInputFiles.push_back(QCDBMuNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(DYJetsToLLNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(TTJetsNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(TNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(WNJetsToLNuNonIsoReweightHaddOutputFile);
//   dataVsMCReweightInputFiles.push_back(WbbNonIsoReweightHaddOutputFile);
//   dataVsMCReweightInputFiles.push_back(WJetsToLNuNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(WZNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(ZZNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(WWNonIsoReweightHaddOutputFile);
  std::reverse(dataVsMCReweightInputFiles.begin() + 1, dataVsMCReweightInputFiles.end());
  cout << "\nPlot data vs. MC normalized to data luminosity\n---\n";
  drawMultipleEfficiencyGraphsOn1Canvas(dataVsMCOutputFile, dataVsMCInputFiles, 
  					canvasNames1D, graphNames1D, legendHeaders19p7InvFb, 
  					colors, styles, legendEntriesMCData, 
  					weightsMCData, setLogY, drawStack, dataMC);
  cout << "\nPlot data vs. MC normalized to 1\n---\n";
  drawMultipleEfficiencyGraphsOn1Canvas(dataVsMCReweightOutputFile, 
  					dataVsMCReweightInputFiles, canvasNames1D, graphNames1D, 
  					legendHeaders19p7InvFb, colors, styles, 
  					legendEntriesMCData, weightsMCData, setLogY, drawStack, 
  					dataMC);
  cout << "\nPlot data minus MC normalized to data luminosity\n---\n";
  drawDifferenceGraphsOn1Canvas(dataVsMCOutputDiff, dataVsMCInputFiles, 
  				canvasNames1D, graphNames1D, legendHeaders19p7InvFb, colors, 
  				styles, legendEntriesMCData, weightsMCData, setLogY, sigBkg);

  //compute data-driven QCD estimate in signal (i.e. isolated W muon + isolated tau) region
  string outputFileNameA(analysisFilePath + "results/dataVsMC_RegionAQCDEstimate" + dataVTag + 
			 fileExt);
  string inputFileNameB(nonIsoWDataIsoHaddOutputFile);
  string inputFileNameC(dataVsMCOutputDiff);
  string inputFileNameD(nonIsoWDataNonIsoHaddOutputFile);
  cout << "\nPlot data-driven QCD estimate for region A\n---\n";
  drawQCDRegionAHistograms(outputFileNameA,inputFileNameB,inputFileNameC,
  			   inputFileNameD,canvasNames1D, graphNames1D,
  			   legendHeaders19p7InvFb,colors, styles, legendEntriesMCData,
  			   weightsMCData, setLogY, sigBkg);

  //compute data-driven QCD estimate in control (i.e. isolated W muon + non-isolated tau) region
  string outputFileNameB(analysisFilePath + "results/dataVsMC_RegionBQCDEstimate" + dataVTag + 
			 fileExt);
  cout << "\nPlot data-driven QCD estimate for region B\n---\n";
  drawQCDRegionAHistograms(outputFileNameB,inputFileNameD,inputFileNameC,
  			   inputFileNameD,canvasNames1D, graphNames1D,
  			   legendHeaders19p7InvFb,colors, styles, legendEntriesMCData,
  			   weightsMCData, setLogY, sigBkg);

  //compare MC signal + data-driven QCD to background
  string sigVsBkgQCDFromDataOutputFile(analysisFilePath + 
				       "results/sigVsBkgQCDFromData_muHadIsoAnalysis" + 
				       tag19p7InvFb + outputVTag + fileExt);
  string sigVsBkgQCDFromDataOutputFile1(analysisFilePath + 
					"results/sigVsBkgQCDFromData_muHadIsoAnalysis" + tag1 + 
					outputVTag + fileExt);
  vector<string> sigVsBkgQCDFromDataInputFiles(sigVsBkgInputFiles);
  sigVsBkgQCDFromDataInputFiles.push_back(outputFileNameA);
  cout << "\nPlot signal vs. background with data-driven QCD estimate ";
  cout << "normalized to data luminosity\n---\n";
  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgQCDFromDataOutputFile, 
					sigVsBkgQCDFromDataInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders19p7InvFb, colors, styles, 
					legendEntriesSigBkgQCDFromData, weightsSigBkgQCDFromData, 
					setLogY, drawStack, sigBkg);
  cout << "\nPlot signal vs. background with data-driven QCD estimate normalized to 1\n---\n";
  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgQCDFromDataOutputFile1, 
					sigVsBkgQCDFromDataInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders1, colors, styles, 
					legendEntriesSigBkgQCDFromData, weights1, setLinY, 
					drawSame, sigBkg);

  //compare data to MC + data-driven QCD in control region
  string dataVsMCQCDFromDataOutputFile(analysisFilePath + 
				       "results/dataVsMCQCDFromData_muHadNonIsoAnalysis" + 
				       tag19p7InvFb + outputVTag + fileExt);
  string dataVsMCQCDFromDataReweightOutputFile
    (analysisFilePath + "results/dataVsMCQCDFromData_muHadNonIsoReweightAnalysis" + tag1 + 
     outputVTag + fileExt);
  vector<string> dataVsMCQCDFromDataInputFiles(dataVsMCInputFiles);
  dataVsMCQCDFromDataInputFiles.push_back(outputFileNameB);
  cout << "\nPlot data vs. MC with data-driven QCD estimate ";
  cout << "normalized to data luminosity\n---\n";
  drawMultipleEfficiencyGraphsOn1Canvas(dataVsMCQCDFromDataOutputFile, 
  					dataVsMCQCDFromDataInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders19p7InvFb, colors, styles, 
  					legendEntriesMCDataQCDFromData, weightsMCDataQCDFromData, 
  					setLogY, drawStack, dataMC);

  cout << "\nBegin region A vs. region B plots, sample by sample...\n\n";
  
  //compare QCD search sample to control sample
  cout << "...muon-enriched QCD\n";
  string QCDSearchVsControlOutputFile(analysisFilePath + "QCD/analysis/isoVsNonIsoTaus" + tag1 + 
  				      outputVTag + fileExt);
  string QCDSearchVsControlReweightOutputFile = 
    smartReplace(QCDSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> QCDSearchVsControlInputFiles;
  QCDSearchVsControlInputFiles.push_back(QCDIsoHaddOutputFile);
  QCDSearchVsControlInputFiles.push_back(QCDNonIsoHaddOutputFile);
  vector<string> QCDSearchVsControlReweightInputFiles(QCDSearchVsControlInputFiles);
  QCDSearchVsControlReweightInputFiles[1] = QCDNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(QCDSearchVsControlOutputFile, 
  					QCDSearchVsControlInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1QCD, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(QCDSearchVsControlReweightOutputFile, 
  					QCDSearchVsControlReweightInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1QCD, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);

  //compare QCDB search sample to control sample
  cout << "...b-enriched QCD\n";
  string QCDBSearchVsControlOutputFile(analysisFilePath + "QCDB/analysis/isoVsNonIsoTaus" + tag1 + 
  				       outputVTag + fileExt);
  string QCDBSearchVsControlReweightOutputFile = 
    smartReplace(QCDBSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> QCDBSearchVsControlInputFiles;
  QCDBSearchVsControlInputFiles.push_back(QCDBIsoHaddOutputFile);
  QCDBSearchVsControlInputFiles.push_back(QCDBNonIsoHaddOutputFile);
  vector<string> QCDBSearchVsControlReweightInputFiles(QCDBSearchVsControlInputFiles);
  QCDBSearchVsControlReweightInputFiles[1] = QCDBNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(QCDBSearchVsControlOutputFile, 
  					QCDBSearchVsControlInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1QCDB, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(QCDBSearchVsControlReweightOutputFile, 
  					QCDBSearchVsControlReweightInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1QCDB, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);

  //compare QCDBMu search sample to control sample
  cout << "...b-to-muon-enriched QCD\n";
  string QCDBMuSearchVsControlOutputFile(analysisFilePath + "QCDBMu/analysis/isoVsNonIsoTaus" + 
  					 tag1 + outputVTag + fileExt);
  string QCDBMuSearchVsControlReweightOutputFile = 
    smartReplace(QCDBMuSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> QCDBMuSearchVsControlInputFiles;
  QCDBMuSearchVsControlInputFiles.push_back(QCDBMuIsoHaddOutputFile);
  QCDBMuSearchVsControlInputFiles.push_back(QCDBMuNonIsoHaddOutputFile);
  vector<string> QCDBMuSearchVsControlReweightInputFiles(QCDBMuSearchVsControlInputFiles);
  QCDBMuSearchVsControlReweightInputFiles[1] = QCDBMuNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(QCDBMuSearchVsControlOutputFile, 
  					QCDBMuSearchVsControlInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1QCDBMu, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(QCDBMuSearchVsControlReweightOutputFile, 
  					QCDBMuSearchVsControlReweightInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1QCDBMu, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);

  //compare Drell-Yan+jets search sample to control sample
  cout << "...Drell-Yan\n";
  string DYJetsToLLSearchVsControlOutputFile(analysisFilePath + 
  					     "DYJetsToLL/analysis/isoVsNonIsoTaus" + tag1 + 
  					     outputVTag + fileExt);
  string DYJetsToLLSearchVsControlReweightOutputFile = 
    smartReplace(DYJetsToLLSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> DYJetsToLLSearchVsControlInputFiles;
  DYJetsToLLSearchVsControlInputFiles.push_back(DYJetsToLLIsoHaddOutputFile);
  DYJetsToLLSearchVsControlInputFiles.push_back(DYJetsToLLNonIsoHaddOutputFile);
  vector<string> DYJetsToLLSearchVsControlReweightInputFiles(DYJetsToLLSearchVsControlInputFiles);
  DYJetsToLLSearchVsControlReweightInputFiles[1] = DYJetsToLLNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(DYJetsToLLSearchVsControlOutputFile, 
  					DYJetsToLLSearchVsControlInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1DYJetsToLL, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(DYJetsToLLSearchVsControlReweightOutputFile, 
  					DYJetsToLLSearchVsControlReweightInputFiles, 
  					canvasNames1D, graphNames1D, legendHeaders1DYJetsToLL, 
  					colors, styles, legendEntriesSearchVsControl, weights1, 
  					setLinY, drawSame, dataMC);

  //compare tt+jets search sample to control sample
  cout << "...ttbar\n";
  string TTJetsSearchVsControlOutputFile(analysisFilePath + "TTJets/analysis/isoVsNonIsoTaus" + 
  					 tag1 + outputVTag + fileExt);
  string TTJetsSearchVsControlReweightOutputFile = 
    smartReplace(TTJetsSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> TTJetsSearchVsControlInputFiles;
  TTJetsSearchVsControlInputFiles.push_back(TTJetsIsoHaddOutputFile);
  TTJetsSearchVsControlInputFiles.push_back(TTJetsNonIsoHaddOutputFile);
  vector<string> TTJetsSearchVsControlReweightInputFiles(TTJetsSearchVsControlInputFiles);
  TTJetsSearchVsControlReweightInputFiles[1] = TTJetsNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(TTJetsSearchVsControlOutputFile, 
  					TTJetsSearchVsControlInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1TTJets, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(TTJetsSearchVsControlReweightOutputFile, 
  					TTJetsSearchVsControlReweightInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1TTJets, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);

  //compare single top search sample to control sample
  cout << "...single top\n";
  string TSearchVsControlOutputFile(analysisFilePath + "SingleTop/analysis/isoVsNonIsoTaus" + 
  				    tag1 + outputVTag + fileExt);
  string TSearchVsControlReweightOutputFile = 
    smartReplace(TSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> TSearchVsControlInputFiles;
  TSearchVsControlInputFiles.push_back(TIsoHaddOutputFile);
  TSearchVsControlInputFiles.push_back(TNonIsoHaddOutputFile);
  vector<string> TSearchVsControlReweightInputFiles(TSearchVsControlInputFiles);
  TSearchVsControlReweightInputFiles[1] = TNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(TSearchVsControlOutputFile, TSearchVsControlInputFiles, 
  					canvasNames1D, graphNames1D, legendHeaders1T, colors, 
  					styles, legendEntriesSearchVsControl, weights1, setLinY, 
  					drawSame, dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(TSearchVsControlReweightOutputFile, 
  					TSearchVsControlReweightInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1T, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);

  //compare W+>=1 jet search sample to control sample
  cout << "...W+>=1 jet\n";
  string WNJetsToLNuSearchVsControlOutputFile(analysisFilePath + 
  					      "WNJetsToLNu/analysis/isoVsNonIsoTaus" + tag1 + 
  					      outputVTag + fileExt);
  string WNJetsToLNuSearchVsControlReweightOutputFile = 
    smartReplace(WNJetsToLNuSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> WNJetsToLNuSearchVsControlInputFiles;
  WNJetsToLNuSearchVsControlInputFiles.push_back(WNJetsToLNuIsoHaddOutputFile);
  WNJetsToLNuSearchVsControlInputFiles.push_back(WNJetsToLNuNonIsoHaddOutputFile);
  vector<string> 
    WNJetsToLNuSearchVsControlReweightInputFiles(WNJetsToLNuSearchVsControlInputFiles);
  WNJetsToLNuSearchVsControlReweightInputFiles[1] = WNJetsToLNuNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(WNJetsToLNuSearchVsControlOutputFile, 
  					WNJetsToLNuSearchVsControlInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1WNJetsToLNu, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(WNJetsToLNuSearchVsControlReweightOutputFile, 
  					WNJetsToLNuSearchVsControlReweightInputFiles, 
  					canvasNames1D, graphNames1D, legendHeaders1WNJetsToLNu, 
  					colors, styles, legendEntriesSearchVsControl, weights1, 
  					setLinY, drawSame, dataMC);

  //compare W+bbbar search sample to control sample
  cout << "...Wbb\n";
  string WbbSearchVsControlOutputFile(analysisFilePath + "Wbb/analysis/isoVsNonIsoTaus" + tag1 + 
  				      outputVTag + fileExt);
  string WbbSearchVsControlReweightOutputFile = 
    smartReplace(WbbSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> WbbSearchVsControlInputFiles;
  WbbSearchVsControlInputFiles.push_back(WbbIsoHaddOutputFile);
  WbbSearchVsControlInputFiles.push_back(WbbNonIsoHaddOutputFile);
  vector<string> WbbSearchVsControlReweightInputFiles(WbbSearchVsControlInputFiles);
  WbbSearchVsControlReweightInputFiles[1] = WbbNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(WbbSearchVsControlOutputFile, 
  					WbbSearchVsControlInputFiles, canvasNames1D, graphNames1D, 
  					legendHeaders1Wbb, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(WbbSearchVsControlReweightOutputFile, 
  					WbbSearchVsControlReweightInputFiles, 
  					canvasNames1D, graphNames1D, legendHeaders1Wbb, colors, 
  					styles, legendEntriesSearchVsControl, weights1, setLinY, 
  					drawSame, dataMC);

  //compare W+jets jet search sample to control sample
  cout << "...W+jets\n";
  string WJetsToLNuSearchVsControlOutputFile(analysisFilePath + 
  					     "WJetsToLNu/analysis/isoVsNonIsoTaus" + tag1 + 
  					     outputVTag + fileExt);
  string WJetsToLNuSearchVsControlReweightOutputFile = 
    smartReplace(WJetsToLNuSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> WJetsToLNuSearchVsControlInputFiles;
  WJetsToLNuSearchVsControlInputFiles.push_back(WJetsToLNuIsoHaddOutputFile);
  WJetsToLNuSearchVsControlInputFiles.push_back(WJetsToLNuNonIsoHaddOutputFile);
  vector<string> WJetsToLNuSearchVsControlReweightInputFiles(WJetsToLNuSearchVsControlInputFiles);
  WJetsToLNuSearchVsControlReweightInputFiles[1] = WJetsToLNuNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(WJetsToLNuSearchVsControlOutputFile, 
  					WJetsToLNuSearchVsControlInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1WJetsToLNu, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(WJetsToLNuSearchVsControlReweightOutputFile, 
  					WJetsToLNuSearchVsControlReweightInputFiles, 
  					canvasNames1D, graphNames1D, legendHeaders1WJetsToLNu, 
  					colors, styles, legendEntriesSearchVsControl, weights1, 
  					setLinY, drawSame, dataMC);

  //compare WZ search sample to control sample
  cout << "...WZ\n";
  string WZSearchVsControlOutputFile(analysisFilePath + "WZ/analysis/isoVsNonIsoTaus" + tag1 + 
  				     outputVTag + fileExt);
  string WZSearchVsControlReweightOutputFile = 
    smartReplace(WZSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> WZSearchVsControlInputFiles;
  WZSearchVsControlInputFiles.push_back(WZIsoHaddOutputFile);
  WZSearchVsControlInputFiles.push_back(WZNonIsoHaddOutputFile);
  vector<string> WZSearchVsControlReweightInputFiles(WZSearchVsControlInputFiles);
  WZSearchVsControlReweightInputFiles[1] = WZNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(WZSearchVsControlOutputFile, WZSearchVsControlInputFiles, 
  					canvasNames1D, graphNames1D, legendHeaders1WZ, colors, 
  					styles, legendEntriesSearchVsControl, weights1, setLinY, 
  					drawSame, dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(WZSearchVsControlReweightOutputFile, 
  					WZSearchVsControlReweightInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1WZ, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);

  //compare ZZ search sample to control sample
  cout << "...ZZ\n";
  string ZZSearchVsControlOutputFile(analysisFilePath + "ZZ/analysis/isoVsNonIsoTaus" + tag1 + 
  				     outputVTag + fileExt);
  string ZZSearchVsControlReweightOutputFile = 
    smartReplace(ZZSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> ZZSearchVsControlInputFiles;
  ZZSearchVsControlInputFiles.push_back(ZZIsoHaddOutputFile);
  ZZSearchVsControlInputFiles.push_back(ZZNonIsoHaddOutputFile);
  vector<string> ZZSearchVsControlReweightInputFiles(ZZSearchVsControlInputFiles);
  ZZSearchVsControlReweightInputFiles[1] = ZZNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(ZZSearchVsControlOutputFile, ZZSearchVsControlInputFiles, 
  					canvasNames1D, graphNames1D, legendHeaders1ZZ, colors, 
  					styles, legendEntriesSearchVsControl, weights1, setLinY, 
  					drawSame, dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(ZZSearchVsControlReweightOutputFile, 
  					ZZSearchVsControlReweightInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1ZZ, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);

  //compare WW search sample to control sample
  cout << "...WW\n";
  string WWSearchVsControlOutputFile(analysisFilePath + "WW/analysis/isoVsNonIsoTaus" + tag1 + 
  				     outputVTag + fileExt);
  string WWSearchVsControlReweightOutputFile = 
    smartReplace(WWSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> WWSearchVsControlInputFiles;
  WWSearchVsControlInputFiles.push_back(WWIsoHaddOutputFile);
  WWSearchVsControlInputFiles.push_back(WWNonIsoHaddOutputFile);
  vector<string> WWSearchVsControlReweightInputFiles(WWSearchVsControlInputFiles);
  WWSearchVsControlReweightInputFiles[1] = WWNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(WWSearchVsControlOutputFile, WWSearchVsControlInputFiles, 
  					canvasNames1D, graphNames1D, legendHeaders1WW, colors, 
  					styles, legendEntriesSearchVsControl, weights1, setLinY, 
  					drawSame, dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(WWSearchVsControlReweightOutputFile, 
  					WWSearchVsControlReweightInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1WW, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);

  //compare region C to region D
  cout << "...non-isolated W data\n";
  string nonIsoWDataSearchVsControlOutputFile(analysisFilePath + 
  					      "nonIsoWData/analysis/isoVsNonIsoTaus" + tag1 + 
  					      outputVTag + fileExt);
  string nonIsoWDataSearchVsControlReweightOutputFile = 
    smartReplace(nonIsoWDataSearchVsControlOutputFile, "NonIso", "NonIsoReweight");
  vector<string> nonIsoWDataSearchVsControlInputFiles;
  nonIsoWDataSearchVsControlInputFiles.push_back(nonIsoWDataIsoHaddOutputFile);
  nonIsoWDataSearchVsControlInputFiles.push_back(nonIsoWDataNonIsoHaddOutputFile);
  vector<string> 
    nonIsoWDataSearchVsControlReweightInputFiles(nonIsoWDataSearchVsControlInputFiles);
  nonIsoWDataSearchVsControlReweightInputFiles[1] = nonIsoWDataNonIsoReweightHaddOutputFile;
  cout << "...without reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(nonIsoWDataSearchVsControlOutputFile, 
  					nonIsoWDataSearchVsControlInputFiles, canvasNames1D, 
  					graphNames1D, legendHeaders1NonIsoWData, colors, styles, 
  					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
  					dataMC);
  cout << "...with reweighting\n";
  drawMultipleEfficiencyGraphsOn1Canvas(nonIsoWDataSearchVsControlReweightOutputFile, 
  					nonIsoWDataSearchVsControlReweightInputFiles, 
  					canvasNames1D, graphNames1D, legendHeaders1NonIsoWData, 
  					colors, styles, legendEntriesSearchVsControl, weights1, 
  					setLinY, drawSame, dataMC);

  cout << "---\nMaking final plots\n";

  //make the final plot showing all background methods, signals, data, and errors
  makeFinalPlot(pair<string, float>(sigVsBkgQCDFromDataOutputFile, 1.0), 
		dataIsoHaddOutputFile, 
		pair<string, float>(dataVsMCOutputFile, 1.0), 
		vector<string>(1, "muHadMass"), vector<string>(1, "m_{#mu+had} (GeV)"), 
		vector<int>(1, 1), vector<int>(1, 2), 
		analysisFilePath + "results/final" + outputVTag + fileExt, "main 5");

//   //print the hadronic tau pT weights and their statistical errors
//   printWeightsAndErrors(nonIsoWDataIsoHaddOutputFile, dataNonIsoHaddOutputFile);

  //MC closure plots
  vector<string> vars;
  vars.push_back("muHadMass");
  vars.push_back("tauHadPT");
  vars.push_back("tauHadEta");
  vars.push_back("bTagDiscrim");
  vars.push_back("MET");
  vector<string> units;
  units.push_back("m_{#mu+had} (GeV)");
  units.push_back("p_{T} (GeV)");
  units.push_back("#eta");
  units.push_back("CSV discriminant");
  units.push_back("#slash{E}_{T} (GeV)");
  vector<int> normRegionLowerBins;
  normRegionLowerBins.push_back(1);
  normRegionLowerBins.push_back(1);
  normRegionLowerBins.push_back(1);
  normRegionLowerBins.push_back(50);
  normRegionLowerBins.push_back(7);
  vector<int> normRegionUpperBins;
  normRegionUpperBins.push_back(2);
  normRegionUpperBins.push_back(3);
  normRegionUpperBins.push_back(46);
  normRegionUpperBins.push_back(50);
  normRegionUpperBins.push_back(9);
  makeMCClosurePlots(sigVsBkgOutputFile, vars, units, dataVsMCOutputFile, 1.0, 
  		     normRegionLowerBins, normRegionUpperBins, 
  		     analysisFilePath + "results/MC_closure_" + outputVersion + fileExt);

  //make plots of hadronic tau pT to support reweighting
  vector<string> fileNames;
  fileNames.push_back(sigVsBkgOutputFile); //region A
  fileNames.push_back(dataVsMCOutputFile); //region B
  fileNames.push_back(nonIsoWDataIsoHaddOutputFile); //region C
  vector<pair<Color_t, Color_t> > colorPairs;
  colorPairs.push_back(pair<Color_t, Color_t>(kBlack, kBlack));
  colorPairs.push_back(pair<Color_t, Color_t>(kRed, kBlue));
  colorPairs.push_back(pair<Color_t, Color_t>(kMagenta, kBlack));
  vector<pair<Style_t, Style_t> > stylePairs;
  stylePairs.push_back(pair<Style_t, Style_t>(20, 20));
  stylePairs.push_back(pair<Style_t, Style_t>(21, 22));
  stylePairs.push_back(pair<Style_t, Style_t>(23, 20));
  plotTauHadPT(fileNames, colorPairs, stylePairs, 
  	       analysisFilePath + "results/tauHadPT" + outputVTag + fileExt);

  //plot ratio of region C and B data tau pT spectra and fit
  cout << "---\nFitting for weights\n";
  calculateTauPTWeightsFromFit(nonIsoWDataIsoHaddOutputFile, dataNonIsoHaddOutputFile, 
  			       analysisFilePath + "results/tauHadPTWeights" + outputVTag + 
  			       fileExt);

  cout << "---\nCalculating fake rates\n";

  //calculate jet-->tau fake rate in data
  const string 
    dataFakeRateFileName(analysisFilePath + "results/fake_rate_data" + outputVTag + fileExt);
  plotFakeRate(SinglePhotonDataIsoHaddOutputFile, SinglePhotonDataNonIsoHaddOutputFile, 
  	       dataFakeRateFileName);

  //calculate jet-->tau fake rate in Drell-Yan MC
  const string DYJetsToLLFakeRateFileName(analysisFilePath + "results/fake_rate_DYJetsToLL" + 
  					  outputVTag + fileExt);
  plotFakeRate(DYJetsToLLIsoHaddOutputFile, DYJetsToLLNonIsoHaddOutputFile, 
  	       DYJetsToLLFakeRateFileName);

  //calculate jet-->tau fake rate in W+jets MC
  const string WNJetsToLNuFakeRateFileName(analysisFilePath + "results/fake_rate_WNJetsToLNu" + 
  					   outputVTag + fileExt);
  plotFakeRate(WNJetsToLNuIsoHaddOutputFile, WNJetsToLNuNonIsoHaddOutputFile, 
  	       WNJetsToLNuFakeRateFileName);

  //calculate jet-->tau fake rate in ttbar MC
  const string 
    TTJetsFakeRateFileName(analysisFilePath + "results/fake_rate_TTJets" + outputVTag + fileExt);
  plotFakeRate(TTJetsIsoHaddOutputFile, TTJetsNonIsoHaddOutputFile, TTJetsFakeRateFileName);

  //calculate jet-->tau fake rate in all MC
  const string 
    MCFakeRateFileName(analysisFilePath + "results/fake_rate_MC" + outputVTag + fileExt);
  plotFakeRate(sigVsBkgOutputFile, dataVsMCOutputFile, MCFakeRateFileName, true);

  //calculate ratio of jet-->tau fake rate in data and Drell-Yan MC
  plotFakeRateRatio(dataFakeRateFileName, DYJetsToLLFakeRateFileName, analysisFilePath + 
  		    "results/fake_rate_ratio_DYJetsToLL" + outputVTag + fileExt);

  //calculate ratio of jet-->tau fake rate in data and W+jets MC
  plotFakeRateRatio(dataFakeRateFileName, WNJetsToLNuFakeRateFileName, analysisFilePath + 
  		    "results/fake_rate_ratio_WNJetsToLNu" + outputVTag + fileExt);

  //calculate ratio of jet-->tau fake rate in data and ttbar MC
  plotFakeRateRatio(dataFakeRateFileName, TTJetsFakeRateFileName, analysisFilePath + 
  		    "results/fake_rate_ratio_TTJets" + outputVTag + fileExt);

  //calculate ratio of jet-->tau fake rate in data and MC
  plotFakeRateRatio(dataFakeRateFileName, MCFakeRateFileName, analysisFilePath + 
  		    "results/fake_rate_ratio_MC" + outputVTag + fileExt);

  // //compare the same plot from 2 versions of the analysis
  // vector<string> fileNamesForComparison1;
  // fileNamesForComparison1.push_back(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + 
  // 				    tag19p7InvFb + "_v149" + fileExt);
  // fileNamesForComparison1.push_back(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + 
  // 				    tag19p7InvFb + "_v149" + fileExt);
  // fileNamesForComparison1.push_back(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + 
  // 				    tag19p7InvFb + "_v149" + fileExt);
  // fileNamesForComparison1.push_back(nonIsoWDataIsoPrefix + "_v149" + fileExt);
  // fileNamesForComparison1.push_back(nonIsoWDataNonIsoPrefix + "_v149" + fileExt);
  // vector<string> fileNamesForComparison2;
  // fileNamesForComparison2.push_back(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + 
  // 				    tag19p7InvFb + "_v152" + fileExt);
  // fileNamesForComparison2.push_back(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + 
  // 				    tag19p7InvFb + "_v152" + fileExt);
  // fileNamesForComparison2.push_back(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + 
  // 				    tag19p7InvFb + "_v152" + fileExt);
  // fileNamesForComparison2.push_back(nonIsoWDataIsoPrefix + "_v152" + fileExt);
  // fileNamesForComparison2.push_back(nonIsoWDataNonIsoPrefix + "_v152" + fileExt);
  // vector<string> outputCanvasTagsForComparison;
  // outputCanvasTagsForComparison.push_back("_regionBData");
  // outputCanvasTagsForComparison.push_back("_regionBMC");
  // outputCanvasTagsForComparison.push_back("_regionAMC");
  // outputCanvasTagsForComparison.push_back("_regionCData");
  // outputCanvasTagsForComparison.push_back("_regionDData");
  // vector<bool> stack;
  // stack.push_back(false);
  // stack.push_back(true);
  // stack.push_back(true);
  // stack.push_back(false);
  // stack.push_back(false);
  // vector<unsigned int> pad;
  // pad.push_back(1);
  // pad.push_back(1);
  // pad.push_back(1);
  // pad.push_back(0);
  // pad.push_back(0);
  // compare2Versions(fileNamesForComparison1, fileNamesForComparison2, analysisFilePath + 
  // 		   "results/comparison_v149_v152" + fileExt, outputCanvasTagsForComparison, stack, 
  // 		   pad);
}
