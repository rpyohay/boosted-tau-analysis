void formatPlots(const string& inputVersion, const string& outputVersion, 
		 const string& rawVersion, const string& reweightVersion)
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
  canvasNames1D.push_back("tauMuTauHadJetWMuMETHTCanvas");
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
  canvasNames2D.push_back("muHad_t3t1VsptmjCanvas");
  canvasNames2D.push_back("muHad_t3t1VsDecayModeCanvas");
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
  graphNames1D.push_back("tauMuTauHadJetWMuMETHT");
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
  graphNames2D.push_back("muHad_t3t1Vsptmj");
  graphNames2D.push_back("muHad_t3t1VsDecayMode");

  //set up plot style options
  vector<string> legendHeaders20InvFb(canvasNames1D.size(), "Normalized to 20 fb^{-1}");
  vector<string> legendHeaders2p5InvFb(canvasNames1D.size(), "Normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbQCD(canvasNames1D.size(), 
					  "QCD #mu-enriched normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbQCDB(canvasNames1D.size(), 
					   "QCD b-enriched normalized to 2.5 fb^{-1}");
  vector<string> 
    legendHeaders2p5InvFbQCDBMu(canvasNames1D.size(), 
				"QCD b#rightarrow#mu-enriched normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbDYJetsToLL(canvasNames1D.size(), 
						 "Drell-Yan + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbTTJets(canvasNames1D.size(), 
					     "t#bar{t} + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbT(canvasNames1D.size(), 
					"t/#bar{t} normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbWNJetsToLNu(canvasNames1D.size(), 
						  "W + #geq1 jet normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbWJetsToLNu(canvasNames1D.size(), 
						 "W + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbWZ(canvasNames1D.size(), "WZ normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbZZ(canvasNames1D.size(), "ZZ normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbWW(canvasNames1D.size(), "WW normalized to 2.5 fb^{-1}");
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
  vector<string> legendHeaders1WJetsToLNu(canvasNames1D.size(), "W + jets normalized to 1");
  vector<string> legendHeaders1WZ(canvasNames1D.size(), "WZ normalized to 1");
  vector<string> legendHeaders1ZZ(canvasNames1D.size(), "ZZ normalized to 1");
  vector<string> legendHeaders1WW(canvasNames1D.size(), "WW normalized to 1");
  vector<Color_t> colors;
  colors.push_back(kBlack);
  colors.push_back(kAzure + 1);
  colors.push_back(kOrange + 1);
  colors.push_back(kGreen - 2);
  colors.push_back(kMagenta + 2);
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
//   colors.push_back(2);
//   colors.push_back(3);
//   colors.push_back(4);
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
  legendEntriesSigBkg.push_back("QCD");
  legendEntriesSigBkg.push_back("QCDB");
  legendEntriesSigBkg.push_back("QCDBMu");
  legendEntriesSigBkg.push_back("Drell-Yan + jets");
  legendEntriesSigBkg.push_back("t#bar{t} + jets");
  legendEntriesSigBkg.push_back("t/#bar{t}");
  legendEntriesSigBkg.push_back("W + #geq1 jet");
//   legendEntriesSigBkg.push_back("W + jets");
  legendEntriesSigBkg.push_back("WZ");
  legendEntriesSigBkg.push_back("ZZ");
  legendEntriesSigBkg.push_back("WW");
  std::reverse(legendEntriesSigBkg.begin() + 1, legendEntriesSigBkg.end());
  vector<string> legendEntriesMCData(legendEntriesSigBkg);
  legendEntriesMCData[0] = "Data 2.5 fb^{-1}";
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
  const float Wh1Weight20InvFb = 0.07208;
  const float Wh1Weight2p5InvFb = 0.00901;
  vector<float> weights1(15, 0.0);
  vector<float> weightsSigBkg;
  weightsSigBkg.push_back(Wh1Weight20InvFb);
  weightsSigBkg.push_back(8.0); /*QCDXSecWeights already weighted to 2.5 fb^-1 ==>
				  multiply by 20/2.5 to get overall weight for 20 fb^-1*/
  weightsSigBkg.push_back(8.0); /*QCDBXSecWeights already weighted to 2.5 fb^-1 ==>
				  multiply by 20/2.5 to get overall weight for 20 fb^-1*/
  weightsSigBkg.push_back(8.0); /*QCDBMuXSecWeights already weighted to 2.5 fb^-1 ==>
				  multiply by 20/2.5 to get overall weight for 20 fb^-1*/
  weightsSigBkg.push_back(1.0); //DYJetsToLLRelXSecWeights already weighted to 20 fb^-1
  weightsSigBkg.push_back(3.60203783312072); //tt + jets
  weightsSigBkg.push_back(1.0); //single top already weighted to 20 fb^-1
  weightsSigBkg.push_back(1.0); //W + >=1 jet already weighted to 20 fb^-1
//   weightsSigBkg.push_back(40.7859690786051); //W + jets
  weightsSigBkg.push_back(0.067773553069845); //WZ
  weightsSigBkg.push_back(0.0383258502693219); //ZZ
  weightsSigBkg.push_back(0.126058122867706); //WW
  std::reverse(weightsSigBkg.begin() + 1, weightsSigBkg.end());
  vector<float> weightsMCData;
  weightsMCData.push_back(1.0); //data (int. lumi. = 2.5 fb^-1)
  weightsMCData.push_back(1.0); /*QCDRelXSecWeights already weighted to 2.5 fb^-1 ==>
				  multiply by 1.0 to get overall weight for 2.5 fb^-1*/
  weightsMCData.push_back(1.0); /*QCDBRelXSecWeights already weighted to 2.5 fb^-1 ==>
				  multiply by 1.0 to get overall weight for 2.5 fb^-1*/
  weightsMCData.push_back(1.0); /*QCDBMuRelXSecWeights already weighted to 2.5 fb^-1 ==>
				  multiply by 1.0 to get overall weight for 2.5 fb^-1*/
  weightsMCData.push_back(0.125); /*DYJetsToLLRelXSecWeights already weighted to 20 fb^-1 ==> 
				    multiply by 2.5/20 to get overall weight for 2.5 fb^-1*/
  weightsMCData.push_back(0.45025472914009); //tt + jets
  weightsMCData.push_back(0.125); /*single top already weighted to 20 fb^-1 ==> 
				    multiply by 2.5/20 to get overall weight for 2.5 fb^-1*/
  weightsMCData.push_back(0.125); /*W + >=1 jet already weighted to 20 fb^-1 ==> 
				    multiply by 2.5/20 to get overall weight for 2.5 fb^-1*/
//   weightsMCData.push_back(5.09824613482563); //W + jets
  weightsMCData.push_back(0.00847169413373063); //WZ
  weightsMCData.push_back(0.00479073128366524); //ZZ
  weightsMCData.push_back(0.0157572653584633); //WW
  std::reverse(weightsMCData.begin() + 1, weightsMCData.end());
  vector<float> QCDRelXSecWeights;  //weighted to 20 fb^-1 using PREP cross sections
  QCDRelXSecWeights.push_back(4396.184992); // Pt 20-30
  QCDRelXSecWeights.push_back(1686.769144); // Pt 30-50
  QCDRelXSecWeights.push_back(339.9588816); // Pt 50-80
  QCDRelXSecWeights.push_back(87.5626528); // Pt 80-120
  QCDRelXSecWeights.push_back(17.55821472); // Pt 120-170
  QCDRelXSecWeights.push_back(5.9967872); // Pt 170-300
  QCDRelXSecWeights.push_back(0.38763976); // Pt 300-470
  QCDRelXSecWeights.push_back(0.06236464); // Pt 470-600
  QCDRelXSecWeights.push_back(0.0130624); // Pt 600-800
  QCDRelXSecWeights.push_back(0.00179552); // Pt 800-1000
  QCDRelXSecWeights.push_back(2.602e-10); // Pt 1000
  vector<float> QCDBRelXSecWeights; //weighted to 20 fb^-1 using PREP cross sections
  QCDBRelXSecWeights.push_back(269041.584); // Pt 15-30
  QCDBRelXSecWeights.push_back(22118.54381); // Pt 30-50
  QCDBRelXSecWeights.push_back(7217.263238); // Pt 50-150
  QCDBRelXSecWeights.push_back(357.6120687); // Pt 150
  vector<float> QCDBMuRelXSecWeights; //weighted to 20 fb^-1 using PREP cross sections
  QCDBMuRelXSecWeights.push_back(0.); // Pt 15-30
  QCDBMuRelXSecWeights.push_back(0.); // Pt 30-50
  QCDBMuRelXSecWeights.push_back(0.); // Pt 50-150
  QCDBMuRelXSecWeights.push_back(0.); // Pt 150
  vector<float> DYJetsToLLRelXSecWeights; //weighted to 20 fb^-1
  DYJetsToLLRelXSecWeights.push_back(7.77158352886295); /*(10 < m < 50) GeV using ttH cross 
							  section*/
  DYJetsToLLRelXSecWeights.push_back(2.30056938223844); //m > 50 GeV using SM@8TeV Twiki
  vector<float> TRelXSecWeights; //weighted to 20 fb^-1 using SM@8TeV Twiki
  TRelXSecWeights.push_back(0.291582198868292); //t s-channel
  TRelXSecWeights.push_back(0.352); //tbar s-channel
  TRelXSecWeights.push_back(0.326178703711468); //t t-channel
  TRelXSecWeights.push_back(0.317300854955268); //tbar t-channel
  vector<float> WNJetsToLNuRelXSecWeights; //weighted to 20 fb^-1 using PREP cross sections
  WNJetsToLNuRelXSecWeights.push_back(5.70531041114793); //W + 1 jet
  WNJetsToLNuRelXSecWeights.push_back(1.24704087800944); //W + 2 jets
  WNJetsToLNuRelXSecWeights.push_back(0.831523180250694); //W + 3 jets
  WNJetsToLNuRelXSecWeights.push_back(0.319813420252842); //W + 4 jets

  //space-saving constant definitions
  string user(gSystem->GetFromPipe("whoami").Data());
  const string analysisFilePath("/data1/" + user + "/");
  const string fileExt(".root");
  const string tag20InvFb("_20fb-1");
  const string tag2p5InvFb("_2p5fb-1");
  const string tag1("_normalizedTo1");

  //version tags
  const string outputVTag("_" + outputVersion);
  const string dataVTag("_" + inputVersion);
  const string sigVTag("_" + inputVersion);
  const string QCDVTag("_" + inputVersion);
  const string QCDBVTag("_" + inputVersion);
  const string QCDBMuVTag("_" + inputVersion);
  const string DYJetsToLLVTag("_" + inputVersion);
  const string TTJetsVTag("_" + inputVersion);
  const string TVTag("_" + inputVersion);
  const string WNJetsToLNuVTag("_" + inputVersion);
  const string WJetsToLNuVTag("_" + inputVersion);
  const string WZVTag("_" + inputVersion);
  const string ZZVTag("_" + inputVersion);
  const string WWVTag("_" + inputVersion);

  //hadd data samples from different eras
  string dataSuffix(dataVTag + fileExt);
//   string dataIsoPrefix(analysisFilePath + "data/analysis/muHadIsoAnalysis_SingleMu");
//   string dataIsoHaddOutputFile(dataIsoPrefix + dataSuffix); //BLINDED!!!
  string dataNonIsoPrefix(analysisFilePath + "data/analysis/muHadNonIsoAnalysis_SingleMu");
  string dataNonIsoHaddOutputFile(dataNonIsoPrefix + dataSuffix);
  string dataNonIsoReweightPrefix(analysisFilePath + "data/analysis/muHadNonIsoReweightAnalysis_SingleMu");
  string dataNonIsoReweightHaddOutputFile(dataNonIsoReweightPrefix + dataSuffix);
//   vector<string> dataIsoHaddInputFiles; //BLINDED!!!
  vector<string> dataNonIsoHaddInputFiles;
  vector<string> dataNonIsoReweightHaddInputFiles;
  vector<string> runEras;
  runEras.push_back("_Run2012A");
  runEras.push_back("_Run2012B");
  runEras.push_back("_Run2012C");
  runEras.push_back("_Run2012D");
  for (vector<string>::const_iterator iRunEra = runEras.begin(); iRunEra != runEras.end(); 
       ++iRunEra) {
//     stringstream dataIsoName;
//     dataIsoName << dataIsoPrefix << *iRunEra << dataSuffix; //BLINDED!!!
//     dataIsoHaddInputFiles.push_back(dataIsoName.str());
    stringstream dataNonIsoName;
    dataNonIsoName << dataNonIsoPrefix << *iRunEra << dataSuffix;
    dataNonIsoHaddInputFiles.push_back(dataNonIsoName.str());
    stringstream dataNonIsoReweightName;
    dataNonIsoReweightName << dataNonIsoReweightPrefix << *iRunEra << dataSuffix;
    dataNonIsoReweightHaddInputFiles.push_back(dataNonIsoReweightName.str());
  }
//   haddCanvases(dataIsoHaddOutputFile, dataIsoHaddInputFiles, vector<float>(4, 1.0), 
// 	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D); //BLINDED!!!
  haddCanvases(dataNonIsoHaddOutputFile, dataNonIsoHaddInputFiles, vector<float>(4, 1.0), 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(dataNonIsoReweightHaddOutputFile, dataNonIsoReweightHaddInputFiles, 
	       vector<float>(4, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //"hadd" Wh1 sample just to get the formatting of the 2D plots the same
  string Wh1IsoPrefix(analysisFilePath + "Wh1_Medium/muHadIsoAnalysis_Wh1");
  string Wh1Suffix(sigVTag + fileExt);
  string Wh1IsoHaddOutputFile(Wh1IsoPrefix + "_hadd" + Wh1Suffix);
  vector<string> Wh1IsoHaddInputFiles;
  stringstream Wh1IsoName;
  Wh1IsoName << Wh1IsoPrefix << Wh1Suffix;
  Wh1IsoHaddInputFiles.push_back(Wh1IsoName.str());
  haddCanvases(Wh1IsoHaddOutputFile, Wh1IsoHaddInputFiles, vector<float>(1, Wh1Weight20InvFb), 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //hadd QCD Mu-enriched Pt-binned samples
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
  haddCanvases(QCDIsoHaddOutputFile, QCDIsoHaddInputFiles, QCDRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(QCDNonIsoHaddOutputFile, QCDNonIsoHaddInputFiles, QCDRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(QCDNonIsoReweightHaddOutputFile, QCDNonIsoReweightHaddInputFiles, 
	       QCDRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
//   haddCanvases(QCDAllHaddOutputFile, QCDAllHaddInputFiles, QCDRelXSecWeights, 
// 	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //hadd QCD b-enriched Pt-binned samples
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
  haddCanvases(QCDBIsoHaddOutputFile, QCDBIsoHaddInputFiles, QCDBRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(QCDBNonIsoHaddOutputFile, QCDBNonIsoHaddInputFiles, QCDBRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(QCDBNonIsoReweightHaddOutputFile, QCDBNonIsoReweightHaddInputFiles, 
	       QCDBRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
//   haddCanvases(QCDBAllHaddOutputFile, QCDBAllHaddInputFiles, QCDBRelXSecWeights, 
// 	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //hadd QCD bToMu-enriched Pt-binned samples
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
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(QCDBMuNonIsoHaddOutputFile, QCDBMuNonIsoHaddInputFiles, QCDBMuRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(QCDBMuNonIsoReweightHaddOutputFile, QCDBMuNonIsoReweightHaddInputFiles, 
	       QCDBMuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
//   haddCanvases(QCDBMuAllHaddOutputFile, QCDBMuAllHaddInputFiles, QCDBMuRelXSecWeights, 
// 	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //hadd Drell-Yan+jets ml+l- binned samples
  string DYJetsToLLSuffix(DYJetsToLLVTag + fileExt);
  string DYJetsToLLIsoPrefix(analysisFilePath + "DYJetsToLL/analysis/muHadIsoAnalysis_DYJetsToLL");
  string DYJetsToLLIsoHaddOutputFile(DYJetsToLLIsoPrefix + DYJetsToLLSuffix);
  string DYJetsToLLNonIsoPrefix(analysisFilePath + 
				"DYJetsToLL/analysis/muHadNonIsoAnalysis_DYJetsToLL");
  string DYJetsToLLNonIsoHaddOutputFile(DYJetsToLLNonIsoPrefix + DYJetsToLLSuffix);
  string DYJetsToLLNonIsoReweightPrefix(analysisFilePath + 
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
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(DYJetsToLLNonIsoHaddOutputFile, DYJetsToLLNonIsoHaddInputFiles, 
	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(DYJetsToLLNonIsoReweightHaddOutputFile, DYJetsToLLNonIsoReweightHaddInputFiles, 
	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
//   haddCanvases(DYJetsToLLAllHaddOutputFile, DYJetsToLLAllHaddInputFiles, 
// 	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //"hadd" ttbar sample just to get the formatting of the 2D plots the same
  string TTJetsSuffix(sigVTag + fileExt);
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
	       vector<float>(1, 3.60203783312072), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(TTJetsNonIsoHaddOutputFile, TTJetsNonIsoHaddInputFiles, 
	       vector<float>(1, 3.60203783312072), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(TTJetsNonIsoReweightHaddOutputFile, TTJetsNonIsoReweightHaddInputFiles, 
	       vector<float>(1, 3.60203783312072), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
//   haddCanvases(TTJetsAllHaddOutputFile, TTJetsAllHaddInputFiles, 
// 	       vector<float>(1, 3.60203783312072), canvasNames1D, graphNames1D, 
// 	       canvasNames2D, graphNames2D);

  //hadd single top samples
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
	       graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(TNonIsoHaddOutputFile, TNonIsoHaddInputFiles, TRelXSecWeights, canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(TNonIsoReweightHaddOutputFile, TNonIsoReweightHaddInputFiles, TRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
//   haddCanvases(TAllHaddOutputFile, TAllHaddInputFiles, TRelXSecWeights, canvasNames1D, 
// 	       graphNames1D, canvasNames2D, graphNames2D);

  //hadd W+>=1 jet samples
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
	       graphNames2D);
  haddCanvases(WNJetsToLNuNonIsoHaddOutputFile, WNJetsToLNuNonIsoHaddInputFiles, 
	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
	       graphNames2D);
  haddCanvases(WNJetsToLNuNonIsoReweightHaddOutputFile, WNJetsToLNuNonIsoReweightHaddInputFiles, 
	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
	       graphNames2D);
//   haddCanvases(WNJetsToLNuAllTauHaddOutputFile, WNJetsToLNuAllTauHaddInputFiles, 
// 	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
// 	       graphNames2D);

  //"hadd" W+jets sample just to get the formatting of the 2D plots the same
  string WJetsToLNuSuffix(sigVTag + fileExt);
  string WJetsToLNuIsoPrefix(analysisFilePath + "WJetsToLNu/analysis/muHadIsoAnalysis_WJetsToLNu");
  string WJetsToLNuIsoHaddOutputFile(WJetsToLNuIsoPrefix + "_hadd" + WJetsToLNuSuffix);
  string WJetsToLNuNonIsoPrefix(analysisFilePath + 
				"WJetsToLNu/analysis/muHadNonIsoAnalysis_WJetsToLNu");
  string WJetsToLNuNonIsoHaddOutputFile(WJetsToLNuNonIsoPrefix + "_hadd" + WJetsToLNuSuffix);
  string WJetsToLNuNonIsoReweightPrefix(analysisFilePath + 
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
	       vector<float>(1, 5.09824613482563), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(WJetsToLNuNonIsoHaddOutputFile, WJetsToLNuNonIsoHaddInputFiles, 
	       vector<float>(1, 5.09824613482563), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(WJetsToLNuNonIsoReweightHaddOutputFile, WJetsToLNuNonIsoReweightHaddInputFiles, 
	       vector<float>(1, 5.09824613482563), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
//   haddCanvases(WJetsToLNuAllHaddOutputFile, WJetsToLNuAllHaddInputFiles, 
// 	       vector<float>(1, 5.09824613482563), canvasNames1D, graphNames1D, 
// 	       canvasNames2D, graphNames2D);

  //"hadd" WZ sample just to get the formatting of the 2D plots the same
  string WZSuffix(sigVTag + fileExt);
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
	       vector<float>(1, 0.067773553069845), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(WZNonIsoHaddOutputFile, WZNonIsoHaddInputFiles, 
	       vector<float>(1, 0.067773553069845), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(WZNonIsoReweightHaddOutputFile, WZNonIsoReweightHaddInputFiles, 
	       vector<float>(1, 0.067773553069845), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
//   haddCanvases(WZAllHaddOutputFile, WZAllHaddInputFiles, 
// 	       vector<float>(1, 0.067773553069845), canvasNames1D, graphNames1D, 
// 	       canvasNames2D, graphNames2D);

  //"hadd" ZZ sample just to get the formatting of the 2D plots the same
  string ZZSuffix(sigVTag + fileExt);
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
	       vector<float>(1, 0.0383258502693219), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(ZZNonIsoHaddOutputFile, ZZNonIsoHaddInputFiles, 
	       vector<float>(1, 0.0383258502693219), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(ZZNonIsoReweightHaddOutputFile, ZZNonIsoReweightHaddInputFiles, 
	       vector<float>(1, 0.0383258502693219), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
//   haddCanvases(ZZAllHaddOutputFile, ZZAllHaddInputFiles, 
// 	       vector<float>(1, 0.0383258502693219), canvasNames1D, graphNames1D, 
// 	       canvasNames2D, graphNames2D);

  //"hadd" WW sample just to get the formatting of the 2D plots the same
  string WWSuffix(sigVTag + fileExt);
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
	       vector<float>(1, 0.126058122867706), canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(WWNonIsoHaddOutputFile, WWNonIsoHaddInputFiles, 
	       vector<float>(1, 0.126058122867706), canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(WWNonIsoReweightHaddOutputFile, WWNonIsoReweightHaddInputFiles, 
	       vector<float>(1, 0.126058122867706), canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D);
//   haddCanvases(WWAllHaddOutputFile, WWAllHaddInputFiles, 
// 	       vector<float>(1, 0.126058122867706), canvasNames1D, 
// 	       graphNames1D, canvasNames2D, graphNames2D);

  //compare MC signal to background
  string sigVsBkgOutputFile20InvFb(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + 
				   tag20InvFb + outputVTag + fileExt);
  string sigVsBkgOutputFile1(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + tag1 + 
			     outputVTag + fileExt);
  string sigVsBkgOutputFileNoHPSIsoCut(analysisFilePath + "results/sigVsBkg_muHadAnalysis" + 
				       tag20InvFb + outputVTag + fileExt);
  string sigVsBkgOutputFileNoHPSIsoCutNorm1(analysisFilePath + "results/sigVsBkg_muHadAnalysis" + 
					    tag1 + outputVTag + fileExt);
  vector<string> sigVsBkgInputFiles;
  sigVsBkgInputFiles.push_back(analysisFilePath + "Wh1_Medium/muHadIsoAnalysis_Wh1" + sigVTag + 
			       fileExt);
  sigVsBkgInputFiles.push_back(QCDIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(QCDBIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(QCDBMuIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(DYJetsToLLIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(analysisFilePath + "TTJets/analysis/muHadIsoAnalysis_TTJets" + 
			       TTJetsVTag + fileExt);
  sigVsBkgInputFiles.push_back(TIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(WNJetsToLNuIsoHaddOutputFile);
//   sigVsBkgInputFiles.push_back(analysisFilePath + 
// 			       "WJetsToLNu/analysis/muHadIsoAnalysis_WJetsToLNu" + 
// 			       WJetsToLNuVTag + fileExt);
  sigVsBkgInputFiles.push_back(analysisFilePath + "WZ/analysis/muHadIsoAnalysis_WZ" + WZVTag + 
			       fileExt);
  sigVsBkgInputFiles.push_back(analysisFilePath + "ZZ/analysis/muHadIsoAnalysis_ZZ" + ZZVTag + 
			       fileExt);
  sigVsBkgInputFiles.push_back(analysisFilePath + "WW/analysis/muHadIsoAnalysis_WW" + WWVTag + 
			       fileExt);
  std::reverse(sigVsBkgInputFiles.begin() + 1, sigVsBkgInputFiles.end());
  vector<string> sigVsBkgNoHPSIsoCutInputFiles;
  sigVsBkgNoHPSIsoCutInputFiles.push_back(analysisFilePath + "Wh1_Medium/muHadAnalysis_Wh1" + 
					  sigVTag + fileExt);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(QCDAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(QCDBAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(QCDBMuAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(DYJetsToLLAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(analysisFilePath + 
					  "TTJets/analysis/muHadAnalysis_TTJets" + TTJetsVTag + 
					  fileExt);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(TAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(WNJetsToLNuAllTauHaddOutputFile);
//   sigVsBkgNoHPSIsoCutInputFiles.push_back(analysisFilePath + 
// 			       "WJetsToLNu/analysis/muHadAnalysis_WJetsToLNu" + 
// 			       WJetsToLNuVTag + fileExt);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(analysisFilePath + "WZ/analysis/muHadAnalysis_WZ" + 
					  WZVTag + fileExt);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(analysisFilePath + "ZZ/analysis/muHadAnalysis_ZZ" + 
					  ZZVTag + fileExt);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(analysisFilePath + "WW/analysis/muHadAnalysis_WW" + 
					  WWVTag + fileExt);
  std::reverse(sigVsBkgNoHPSIsoCutInputFiles.begin() + 1, sigVsBkgNoHPSIsoCutInputFiles.end());
  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile20InvFb, sigVsBkgInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders20InvFb, 
					colors, styles, legendEntriesSigBkg, weightsSigBkg, 
					setLogY, drawStack, sigBkg);
  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile1, sigVsBkgInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders1, colors, 
					styles, legendEntriesSigBkg, weights1, setLinY, drawSame, 
					sigBkg);
//   drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFileNoHPSIsoCut, 
// 					sigVsBkgNoHPSIsoCutInputFiles, canvasNames1D, 
// 					graphNames1D, legendHeaders20InvFb, colors, styles, 
// 					legendEntriesSigBkg, weightsSigBkg, setLogY, drawStack, 
// 					sigBkg);
//   drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFileNoHPSIsoCutNorm1, 
// 					sigVsBkgNoHPSIsoCutInputFiles, canvasNames1D, 
// 					graphNames1D, legendHeaders1, colors, styles, 
// 					legendEntriesSigBkg, weights1, setLinY, drawSame, 
// 					sigBkg);

  //compare data to MC in control region
  string dataVsMCOutputFile2p5InvFb(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + 
				    tag2p5InvFb + outputVTag + fileExt);
  string dataVsMCReweightOutputFile2p5InvFb(analysisFilePath + 
					    "results/dataVsMC_muHadNonIsoReweightAnalysis" + 
					    tag2p5InvFb + outputVTag + fileExt);
  vector<string> dataVsMCInputFiles;
  dataVsMCInputFiles.push_back(dataNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(QCDNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(QCDBNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(QCDBMuNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(DYJetsToLLNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(analysisFilePath + "TTJets/analysis/muHadNonIsoAnalysis_TTJets" + 
			       TTJetsVTag + fileExt);
  dataVsMCInputFiles.push_back(TNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(WNJetsToLNuNonIsoHaddOutputFile);
//   dataVsMCInputFiles.push_back(analysisFilePath + 
// 			       "WJetsToLNu/analysis/muHadNonIsoAnalysis_WJetsToLNu" + 
// 			       WJetsToLNuVTag + fileExt);
  dataVsMCInputFiles.push_back(analysisFilePath + "WZ/analysis/muHadNonIsoAnalysis_WZ" + WZVTag + 
			       fileExt);
  dataVsMCInputFiles.push_back(analysisFilePath + "ZZ/analysis/muHadNonIsoAnalysis_ZZ" + ZZVTag + 
			       fileExt);
  dataVsMCInputFiles.push_back(analysisFilePath + "WW/analysis/muHadNonIsoAnalysis_WW" + WWVTag + 
			       fileExt);
  std::reverse(dataVsMCInputFiles.begin() + 1, dataVsMCInputFiles.end());
  vector<string> dataVsMCReweightInputFiles;
  dataVsMCReweightInputFiles.push_back(dataNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(QCDNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(QCDBNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(QCDBMuNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(DYJetsToLLNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(analysisFilePath + 
				       "TTJets/analysis/muHadNonIsoReweightAnalysis_TTJets" + 
				       TTJetsVTag + fileExt);
  dataVsMCReweightInputFiles.push_back(TNonIsoReweightHaddOutputFile);
  dataVsMCReweightInputFiles.push_back(WNJetsToLNuNonIsoReweightHaddOutputFile);
//   dataVsMCReweightInputFiles.
//     push_back(analysisFilePath + "WJetsToLNu/analysis/muHadNonIsoReweightAnalysis_WJetsToLNu" + 
// 	      WJetsToLNuVTag + fileExt);
  dataVsMCReweightInputFiles.push_back(analysisFilePath + 
				       "WZ/analysis/muHadNonIsoReweightAnalysis_WZ" + WZVTag + 
				       fileExt);
  dataVsMCReweightInputFiles.push_back(analysisFilePath + 
				       "ZZ/analysis/muHadNonIsoReweightAnalysis_ZZ" + ZZVTag + 
				       fileExt);
  dataVsMCReweightInputFiles.push_back(analysisFilePath + 
				       "WW/analysis/muHadNonIsoReweightAnalysis_WW" + WWVTag + 
				       fileExt);
  std::reverse(dataVsMCReweightInputFiles.begin() + 1, dataVsMCReweightInputFiles.end());
  drawMultipleEfficiencyGraphsOn1Canvas(dataVsMCOutputFile2p5InvFb, dataVsMCInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders2p5InvFb, 
					colors, styles, legendEntriesMCData, 
					weightsMCData, setLogY, drawStack, dataMC);
  drawMultipleEfficiencyGraphsOn1Canvas(dataVsMCReweightOutputFile2p5InvFb, 
					dataVsMCReweightInputFiles, canvasNames1D, graphNames1D, 
					legendHeaders2p5InvFb, colors, styles, 
					legendEntriesMCData, weightsMCData, setLogY, drawStack, 
					dataMC);

  //compare QCD search sample to control sample
  string QCDSearchVsControlOutputFile(analysisFilePath + "QCD/analysis/isoVsNonIsoTaus" + tag1 + 
				      outputVTag + fileExt);
  vector<string> QCDSearchVsControlInputFiles;
  QCDSearchVsControlInputFiles.push_back(QCDIsoHaddOutputFile);
  QCDSearchVsControlInputFiles.push_back(QCDNonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(QCDSearchVsControlOutputFile, 
					QCDSearchVsControlInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders1QCD, colors, styles, 
					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
					sigBkg);

  //compare QCDB search sample to control sample
  string QCDBSearchVsControlOutputFile(analysisFilePath + "QCDB/analysis/isoVsNonIsoTaus" + tag1 + 
				       outputVTag + fileExt);
  vector<string> QCDBSearchVsControlInputFiles;
  QCDBSearchVsControlInputFiles.push_back(QCDBIsoHaddOutputFile);
  QCDBSearchVsControlInputFiles.push_back(QCDBNonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(QCDBSearchVsControlOutputFile, 
					QCDBSearchVsControlInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders1QCDB, colors, styles, 
					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
					sigBkg);

  //compare QCDBMu search sample to control sample
  string QCDBMuSearchVsControlOutputFile(analysisFilePath + "QCDBMu/analysis/isoVsNonIsoTaus" + 
					 tag1 + outputVTag + fileExt);
  vector<string> QCDBMuSearchVsControlInputFiles;
  QCDBMuSearchVsControlInputFiles.push_back(QCDBMuIsoHaddOutputFile);
  QCDBMuSearchVsControlInputFiles.push_back(QCDBMuNonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(QCDBMuSearchVsControlOutputFile, 
					QCDBMuSearchVsControlInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders1QCDBMu, colors, styles, 
					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
					sigBkg);

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
					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
					sigBkg);

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
					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
					sigBkg);

  //compare single top search sample to control sample
  string TSearchVsControlOutputFile(analysisFilePath + "SingleTop/analysis/isoVsNonIsoTaus" + 
				    tag1 + outputVTag + fileExt);
  vector<string> TSearchVsControlInputFiles;
  TSearchVsControlInputFiles.push_back(TIsoHaddOutputFile);
  TSearchVsControlInputFiles.push_back(TNonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(TSearchVsControlOutputFile, TSearchVsControlInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders1T, colors, 
					styles, legendEntriesSearchVsControl, weights1, setLinY, 
					drawSame, sigBkg);

  //compare W+>=1 jet search sample to control sample
  string WNJetsToLNuSearchVsControlOutputFile(analysisFilePath + 
					      "WNJetsToLNu/analysis/isoVsNonIsoTaus" + tag1 + 
					      outputVTag + fileExt);
  vector<string> WNJetsToLNuSearchVsControlInputFiles;
  WNJetsToLNuSearchVsControlInputFiles.push_back(WNJetsToLNuIsoHaddOutputFile);
  WNJetsToLNuSearchVsControlInputFiles.push_back(WNJetsToLNuNonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(WNJetsToLNuSearchVsControlOutputFile, 
					WNJetsToLNuSearchVsControlInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders1WNJetsToLNu, colors, styles, 
					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
					sigBkg);

  //compare W+jets jet search sample to control sample
  string WJetsToLNuSearchVsControlOutputFile(analysisFilePath + 
					     "WJetsToLNu/analysis/isoVsNonIsoTaus" + tag1 + 
					     outputVTag + fileExt);
  vector<string> WJetsToLNuSearchVsControlInputFiles;
  WJetsToLNuSearchVsControlInputFiles.push_back(WJetsToLNuIsoHaddOutputFile);
  WJetsToLNuSearchVsControlInputFiles.push_back(WJetsToLNuNonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(WJetsToLNuSearchVsControlOutputFile, 
					WJetsToLNuSearchVsControlInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders1WJetsToLNu, colors, styles, 
					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
					sigBkg);

  //compare WZ search sample to control sample
  string WZSearchVsControlOutputFile(analysisFilePath + "WZ/analysis/isoVsNonIsoTaus" + tag1 + 
				     outputVTag + fileExt);
  vector<string> WZSearchVsControlInputFiles;
  WZSearchVsControlInputFiles.push_back(analysisFilePath + "WZ/analysis/muHadIsoAnalysis_WZ" + 
					WZVTag + fileExt);
  WZSearchVsControlInputFiles.push_back(analysisFilePath + "WZ/analysis/muHadNonIsoAnalysis_WZ" + 
					WZVTag + fileExt);
  drawMultipleEfficiencyGraphsOn1Canvas(WZSearchVsControlOutputFile, WZSearchVsControlInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders1WZ, colors, 
					styles, legendEntriesSearchVsControl, weights1, setLinY, 
					drawSame, sigBkg);

  //compare ZZ search sample to control sample
  string ZZSearchVsControlOutputFile(analysisFilePath + "ZZ/analysis/isoVsNonIsoTaus" + tag1 + 
				     outputVTag + fileExt);
  vector<string> ZZSearchVsControlInputFiles;
  ZZSearchVsControlInputFiles.push_back(analysisFilePath + "ZZ/analysis/muHadIsoAnalysis_ZZ" + 
					ZZVTag + fileExt);
  ZZSearchVsControlInputFiles.push_back(analysisFilePath + "ZZ/analysis/muHadNonIsoAnalysis_ZZ" + 
					ZZVTag + fileExt);
  drawMultipleEfficiencyGraphsOn1Canvas(ZZSearchVsControlOutputFile, ZZSearchVsControlInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders1ZZ, colors, 
					styles, legendEntriesSearchVsControl, weights1, setLinY, 
					drawSame, sigBkg);

  //compare WW search sample to control sample
  string WWSearchVsControlOutputFile(analysisFilePath + "WW/analysis/isoVsNonIsoTaus" + tag1 + 
				     outputVTag + fileExt);
  vector<string> WWSearchVsControlInputFiles;
  WWSearchVsControlInputFiles.push_back(analysisFilePath + "WW/analysis/muHadIsoAnalysis_WW" + 
					WWVTag + fileExt);
  WWSearchVsControlInputFiles.push_back(analysisFilePath + "WW/analysis/muHadNonIsoAnalysis_WW" + 
					WWVTag + fileExt);
  drawMultipleEfficiencyGraphsOn1Canvas(WWSearchVsControlOutputFile, WWSearchVsControlInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders1WW, colors, 
					styles, legendEntriesSearchVsControl, weights1, setLinY, 
					drawSame, sigBkg);

  //MC closure plots
  vector<string> vars;
  vars.push_back("muHadMass");
  vars.push_back("muHadPTOverMuHadMass");
  vars.push_back("tauHadPT");
  vars.push_back("tauMuPT");
  vector<string> units;
  units.push_back("m_{#mu+had} (GeV)");
  units.push_back("p_{T}^{#mu+had}/m^{#mu+had}");
  units.push_back("p_{T} (GeV)");
  units.push_back("p_{T} (GeV)");
  vector<int> normRegionLowerBins;
  normRegionLowerBins.push_back(1);
  normRegionLowerBins.push_back(1);
  normRegionLowerBins.push_back(1);
  normRegionLowerBins.push_back(1);
  vector<int> normRegionUpperBins;
  normRegionUpperBins.push_back(2);
  normRegionUpperBins.push_back(7);
  normRegionUpperBins.push_back(3);
  normRegionUpperBins.push_back(3);
  makeMCClosurePlots(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + tag20InvFb + "_" + 
		     rawVersion + fileExt, vars, units, analysisFilePath + 
		     "results/dataVsMC_muHadNonIsoReweightAnalysis" + tag2p5InvFb + "_" + 
		     reweightVersion + fileExt, 20.0/2.5, normRegionLowerBins, 
		     normRegionUpperBins, analysisFilePath + "results/MC_closure_" + 
		     reweightVersion + fileExt);
}
