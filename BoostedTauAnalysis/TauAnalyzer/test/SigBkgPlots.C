{
  //load
  gROOT->Reset();
  gROOT->ProcessLine("#include <utility>");
  string macroPath("/afs/cern.ch/user/f/friccita/myOtherOther533Area/CMSSW_5_3_3/src/BoostedTauAnalysis/");
  macroPath+="TauAnalyzer/test/";
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
  //canvasNames1D.push_back("tauHadMTCanvas");
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
  canvasNames2D.push_back("muHad_t3t1VsptmjCanvas");
  canvasNames2D.push_back("muHad_t3t1VsDecayModeCanvas");
  vector<string> graphNames1D;
  graphNames1D.push_back("hadTauAssociatedMuMultiplicity");
  graphNames1D.push_back("muHadMass");
  graphNames1D.push_back("muHadCharge");
  graphNames1D.push_back("MET");
  graphNames1D.push_back("WMuMT");
  graphNames1D.push_back("tauMuMT");
  //graphNames1D.push_back("tauHadMT");
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
  graphNames2D.push_back("muHad_t3t1Vsptmj");
  graphNames2D.push_back("muHad_t3t1VsDecayMode");
  //set up plot style options
  vector<string> legendHeaders20InvFb(canvasNames1D.size(), "Normalized to 20 fb^{-1}");
  vector<string> legendHeaders2p5InvFb(canvasNames1D.size(), "Normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbQCD(canvasNames1D.size(), 
					  "QCD data-driven normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbDYJetsToLL(canvasNames1D.size(), 
						 "Drell-Yan + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbTTJets(canvasNames1D.size(), 
					     "t#bar{t} + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbT(canvasNames1D.size(), 
					"t/#bar{t} normalized to 2.5 fb^{-1}");
  vector<string> 
    legendHeaders2p5InvFbWNJetsToLNu(canvasNames1D.size(), "W + jets normalized to 2.5 fb^{-1}");
  //vector<string> 
  //legendHeaders2p5InvFbWJetsToLNu(canvasNames1D.size(), "W + jets inclusive normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbWZ(canvasNames1D.size(), "WZ normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbZZ(canvasNames1D.size(), "ZZ normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders1(canvasNames1D.size(), "Normalized to 1");
  vector<string> legendHeaders1QCD(canvasNames1D.size(), "QCD data-driven normalized to 1");
  vector<string> legendHeaders1DYJetsToLL(canvasNames1D.size(), 
					  "Drell-Yan + jets normalized to 1");
  vector<string> legendHeaders1TTJets(canvasNames1D.size(), "t#bar{t} + jets normalized to 1");
  vector<string> legendHeaders1T(canvasNames1D.size(), "t/#bar{t} normalized to 1");
  vector<string> legendHeaders1WNJetsToLNu(canvasNames1D.size(), "W + jets normalized to 1");
  //vector<string> legendHeaders1WJetsToLNu(canvasNames1D.size(), "W + jets inclusive normalized to 1");
  vector<string> legendHeaders1WZ(canvasNames1D.size(), "WZ normalized to 1");
  vector<string> legendHeaders1ZZ(canvasNames1D.size(), "ZZ normalized to 1");
  vector<Color_t> colors;
  colors.push_back(kBlack);
  colors.push_back(/*kRed*/kAzure + 1);
  colors.push_back(/*kBlue*/kOrange + 1);
  colors.push_back(/*kMagenta*/kGreen - 2);
  colors.push_back(/*kGreen*/kMagenta + 2);
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
  legendEntriesSigBkg.push_back("QCDData"); //data-driven
  legendEntriesSigBkg.push_back("Drell-Yan + jets");
  legendEntriesSigBkg.push_back("t#bar{t} + jets");
  legendEntriesSigBkg.push_back("t/#bar{t}");
  legendEntriesSigBkg.push_back("W + jets");
   //  legendEntriesSigBkg.push_back("W + jets inclusive");
  legendEntriesSigBkg.push_back("WZ");
  legendEntriesSigBkg.push_back("ZZ");
//   legendEntriesSigBkg.push_back("WW");
  std::reverse(legendEntriesSigBkg.begin() + 1, legendEntriesSigBkg.end());
  const bool setLogY = true;
  const bool setLinY = false;
  const bool drawStack = true;
  const bool drawSame = false;
  const bool dataMC = false;
  const bool sigBkg = true;

  //weights (sig. figs are probably wrong)
  //first number in parentheses is the PREP weight
  //second number in parentheses is the best available weight
  const float Wh1Weight20InvFb = 0.07208; /*(0.3604 pb(Pythia LO xs) * 20000 pb^-1)/
  100000(no. events processed)*/
  const float Wh1Weight2p5InvFb = 0.00901; //20 fb^-1 weight * (2.5/20)
  vector<float> weights1(15, 0.0);
  vector<float> weightsSigBkg;
  weightsSigBkg.push_back(Wh1Weight20InvFb);
  weightsSigBkg.push_back(1.0); /*QCDXSecWeights already weighted to 20 fb^-1*/
  weightsSigBkg.push_back(1.0); /*DYJetsToLLRelXSecWeights already weighted to 20 fb^-1 ==> 
 multiply by 1.0 to get overall weight for 20 fb^-1*/
  weightsSigBkg.push_back(/*1.99738713*/3.60203783312072); //tt+jets weighted to 20 fb^-1
  weightsSigBkg.push_back(1.0); //t already weighted to 20 fb^-1
  weightsSigBkg.push_back(1.0); //W + jets already weighted to 20 fb^-1
  //  weightsSigBkg.push_back(/*33.055892185598*/40.785969078605063); /*W + jets inclusive weighted to 20 fb^-1
//      (37509 pb(8 TeV Twiki xs) * 20000 pb^-1)/
//      18393090(no. events processed)*/
  weightsSigBkg.push_back(/*0.0257747659*/0.067773553069845); //WZ weighted to 20 fb^-1
  weightsSigBkg.push_back(/*0.0112802265*/0.0383258502693219); //ZZ weighted to 20 fb^-1
//   weightsSigBkg.push_back(/*0.0772605403111639*/0.126058122867706); //WW weighted to 20 fb^-1
  std::reverse(weightsSigBkg.begin() + 1, weightsSigBkg.end());

  vector<float> WNJetsToLNuRelXSecWeights;
  WNJetsToLNuRelXSecWeights.push_back(/*3.38697439995285*/5.7053104111); /*W + 1 jet weighted to 20 fb^-1*/
  WNJetsToLNuRelXSecWeights.push_back(/*0.461044066257386*/1.2984609364); /*W + 2 jets weighted to 20 fb^-1*/
  WNJetsToLNuRelXSecWeights.push_back(/*0.164178582464271*/0.81547009579); /*W + 3 jets weighted to 20 fb^-1*/
  WNJetsToLNuRelXSecWeights.push_back(/*0.0314134490360502*/0.3198134203); /*W + 4 jets weighted to 20 fb^-1*/
  vector<float> DYJetsToLLRelXSecWeights;
  //  DYJetsToLLRelXSecWeights.push_back(/*0.7301387396*/0.148111781928372); /*(10 < m < 50) GeV 
  //  weighted to 2.5 fb^-1*/
  DYJetsToLLRelXSecWeights.push_back(/*0.7301387396*//*0.97144794111*/7.77158352888); /*(10 < m < 50) GeV 
  weighted to (2.5)20 fb^-1*/
  DYJetsToLLRelXSecWeights.push_back(/*0.2421247648*//*0.287571172779805*/2.30056938223844); /*m > 50 GeV 
  weighted to (2.5)20 fb^-1*/
  vector<float> TRelXSecWeights;
  TRelXSecWeights.push_back(0.2169556203); //t s-channel weighted to 20 fb^-1
  TRelXSecWeights.push_back(0.314); //tbar s-channel weighted to 20 fb^-1
  TRelXSecWeights.push_back(0.2718155864); //t t-channel weighted to 20 fb^-1
  TRelXSecWeights.push_back(0.2583883184); //tbar t-channel weighted to 20 fb^-1

  //space-saving constant definitions
  const string analysisFilePath("/data1/friccita/");
  const string fileExt(".root");
  const string tag20InvFb("_20fb-1");
  const string tag2p5InvFb("_2p5fb-1");
  const string tag1("_normalizedTo1");

  //version tags
  const string outputVTag("_v12");
  const string dataVTag("_v12");
  const string sigVTag("_v12");
  const string WNJetsToLNuVTag("_v12");
  const string WJetsVTag("_v12");
  const string TTJetsVTag("_v12");
  const string TVTag("_v12");
  const string DYJetsToLLVTag("_v12");
  const string WZVTag("_v12");
  const string ZZVTag("_v12");
  const string WWVTag("_v12");
  const string QCDVTag("_v12");


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
  haddCanvases(DYJetsToLLIsoHaddOutputFile, DYJetsToLLIsoHaddInputFiles, DYJetsToLLRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(DYJetsToLLNonIsoHaddOutputFile, DYJetsToLLNonIsoHaddInputFiles, 
	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //data-driven QCD estimate
  string QCDIsoHaddOutputFile = "/data1/friccita/results/dataVsMC_RegionAQCDEstimate";
  QCDIsoHaddOutputFile+=dataVTag;
  QCDIsoHaddOutputFile+=".root";

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
  haddCanvases(WNJetsToLNuIsoHaddOutputFile, WNJetsToLNuIsoHaddInputFiles, 
	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
	       graphNames2D);
  haddCanvases(WNJetsToLNuNonIsoHaddOutputFile, WNJetsToLNuNonIsoHaddInputFiles, 
	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
	       graphNames2D);
//   haddCanvases(WNJetsToLNuAllTauHaddOutputFile, WNJetsToLNuAllTauHaddInputFiles, 
// 	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
// 	       graphNames2D);

  //hadd single top samples
  string TIsoPrefix(analysisFilePath + "SingleTop/analysis/muHadIsoAnalysis_T");
  string TNonIsoPrefix(analysisFilePath + "SingleTop/analysis/muHadNonIsoAnalysis_T");
  string TSuffix(TVTag + fileExt);
  string TIsoHaddOutputFile(TIsoPrefix + TSuffix);
  string TNonIsoHaddOutputFile(TNonIsoPrefix + TSuffix);
  vector<string> TIsoHaddInputFiles;
  vector<string> TNonIsoHaddInputFiles;
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
  }
  haddCanvases(TIsoHaddOutputFile, TIsoHaddInputFiles, TRelXSecWeights, canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(TNonIsoHaddOutputFile, TNonIsoHaddInputFiles, TRelXSecWeights, canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D);



  //compare MC signal to background
  string sigVsBkgOutputFile20InvFb(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis_QCDDataDriven" + 
				   tag20InvFb + outputVTag + fileExt);
  string sigVsBkgOutputFile1(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis_QCDDataDriven" + tag1 + 
			     outputVTag + fileExt);
  vector<string> sigVsBkgInputFiles;
  sigVsBkgInputFiles.push_back(analysisFilePath + "Wh1_Medium/muHadIsoAnalysis_Wh1" + sigVTag + 
			       fileExt);
  sigVsBkgInputFiles.push_back(QCDIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(DYJetsToLLIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(analysisFilePath + "TTJets/analysis/muHadIsoAnalysis_TTJets" + 
			       TTJetsVTag + fileExt);
  sigVsBkgInputFiles.push_back(TIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(WNJetsToLNuIsoHaddOutputFile);
  //  sigVsBkgInputFiles.push_back(analysisFilePath + "WJetsToLNu/analysis/muHadIsoAnalysis_WJetsToLNu" + WJetsVTag + 
  //			       fileExt);
  sigVsBkgInputFiles.push_back(analysisFilePath + "WZ/analysis/muHadIsoAnalysis_WZ" + WZVTag + 
			       fileExt);
  sigVsBkgInputFiles.push_back(analysisFilePath + "ZZ/analysis/muHadIsoAnalysis_ZZ" + ZZVTag + 
			       fileExt);
//   sigVsBkgInputFiles.push_back(analysisFilePath + "WW/analysis/muHadIsoAnalysis_WW" + WWVTag + 
// 			       fileExt);
  std::reverse(sigVsBkgInputFiles.begin() + 1, sigVsBkgInputFiles.end());

  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile20InvFb, sigVsBkgInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders20InvFb, 
					colors, styles, legendEntriesSigBkg, weightsSigBkg, 
					setLogY, drawStack, dataMC);
  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile1, sigVsBkgInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders1, colors, 
					styles, legendEntriesSigBkg, weights1, setLinY, drawSame, 
					dataMC);
  
}
