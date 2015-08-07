//REGION A DATA 2D HISTOGRAMS ARE NOT BLINDED!!!  BEWARE!!!

void formatSigPlots(const string& inputVersion, const string& outputVersion, 
		    const string& versionNarrow, const string& uncTag, const string& a1Mass, 
		    const string& MTBin, const unsigned int firstBinToBlind, 
		    const bool doNoHPSIsoCut = false)
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

  //set up canvas and graph names and blinded bins for data
  vector<string> canvasNames1D;
  canvasNames1D.push_back("hadTauAssociatedMuMultiplicityCanvas");
  canvasNames1D.push_back("muHadMassCanvas");
  canvasNames1D.push_back("muHadMassZMuCanvas");
  canvasNames1D.push_back("muHadMassZTauMuCanvas");
  canvasNames1D.push_back("muHadMassOtherTauMuCanvas");
  canvasNames1D.push_back("muHadMassMuElseCanvas");
  canvasNames1D.push_back("muHadMassNotMuElseCanvas");
  canvasNames1D.push_back("muHadMass1ProngCanvas");
  canvasNames1D.push_back("muHadMass1Prong1Pi0Canvas");
  canvasNames1D.push_back("muHadMass1Prong2Pi0Canvas");
  canvasNames1D.push_back("muHadMass3ProngCanvas");
  canvasNames1D.push_back("muHadMass3MuShareTrackCanvas");
  canvasNames1D.push_back("muHadMass3MuSoftMuCanvas");
  canvasNames1D.push_back("muHadMass3MuSoftMu5GeVCanvas");
  canvasNames1D.push_back("muHadMass3MuSoftMu15GeVCanvas");
  canvasNames1D.push_back("muHadMass3MuSoftMu20GeVCanvas");
  canvasNames1D.push_back("muHadMassReweightErrSqCanvas");
  canvasNames1D.push_back("muHadChargeCanvas");
  canvasNames1D.push_back("muHadDdxyCanvas");
  canvasNames1D.push_back("muHadDdzCanvas");
  canvasNames1D.push_back("WMuPVdzCanvas");
  canvasNames1D.push_back("muPVdzCanvas");
  canvasNames1D.push_back("hadPVdzCanvas");
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
  canvasNames1D.push_back("dRSoftMuTauHadCanvas");
  canvasNames1D.push_back("HPTCanvas");
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
  canvasNames1D.push_back("nAddlJetsPTGeq20Canvas");
  canvasNames1D.push_back("nAddlJetsPTGeq30Canvas");
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
  canvasNames1D.push_back("muHad_t1t2Canvas");
  canvasNames1D.push_back("muHad_t2t3Canvas");
  canvasNames1D.push_back("muHad_t3t4Canvas");
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
  canvasNames1D.push_back("dRWMuSoftMuCanvas");
  canvasNames1D.push_back("dRWMuTauHadCanvas");
  canvasNames1D.push_back("dRTauMuTauHadCanvas");
  canvasNames1D.push_back("dRWMuTauMuTauHadCanvas");
  canvasNames1D.push_back("dPhiWMuSoftMuCanvas");
  canvasNames1D.push_back("dPhiWMuSoftMuWithCutCanvas");
  canvasNames1D.push_back("dPhiWMuSecJetCanvas");
  canvasNames1D.push_back("WMu3rdTightMuChargeProductCanvas");
  canvasNames1D.push_back("tauHadPhotonEnergyFractionCanvas");
  canvasNames1D.push_back("dThetaPhotonOtherTauConstituentsCanvas");
  canvasNames1D.push_back("hardestCorrJetEtaCanvas");
  canvasNames1D.push_back("WMuPTCanvas");
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
  canvasNames2D.push_back("muHadMassVsWMuMTCanvas");
  canvasNames2D.push_back("tauHadJetEnergyFractionVsTauHadIsoCanvas");
  canvasNames2D.push_back("tauHadCleanedJetEnergyFractionVsTauHadIsoCanvas");
  canvasNames2D.push_back("mu1ProngMassVsMu1Prong1Pi0MassCanvas");
  canvasNames2D.push_back("dPhiTauMuMETVsMuHadMassCanvas");
  vector<string> graphNames1D;
  graphNames1D.push_back("hadTauAssociatedMuMultiplicity");
  graphNames1D.push_back("muHadMass");
  graphNames1D.push_back("muHadMassZMu");
  graphNames1D.push_back("muHadMassZTauMu");
  graphNames1D.push_back("muHadMassOtherTauMu");
  graphNames1D.push_back("muHadMassMuElse");
  graphNames1D.push_back("muHadMassNotMuElse");
  graphNames1D.push_back("muHadMass1Prong");
  graphNames1D.push_back("muHadMass1Prong1Pi0");
  graphNames1D.push_back("muHadMass1Prong2Pi0");
  graphNames1D.push_back("muHadMass3Prong");
  graphNames1D.push_back("muHadMass3MuShareTrack");
  graphNames1D.push_back("muHadMass3MuSoftMu");
  graphNames1D.push_back("muHadMass3MuSoftMu5GeV");
  graphNames1D.push_back("muHadMass3MuSoftMu15GeV");
  graphNames1D.push_back("muHadMass3MuSoftMu20GeV");
  graphNames1D.push_back("muHadMassReweightErrSq");
  graphNames1D.push_back("muHadCharge");
  graphNames1D.push_back("muHadDdxy");
  graphNames1D.push_back("muHadDdz");
  graphNames1D.push_back("WMuPVdz");
  graphNames1D.push_back("muPVdz");
  graphNames1D.push_back("hadPVdz");
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
  graphNames1D.push_back("dRSoftMuTauHad");
  graphNames1D.push_back("HPT");
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
  graphNames1D.push_back("nAddlJetsPTGeq20");
  graphNames1D.push_back("nAddlJetsPTGeq30");
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
  graphNames1D.push_back("muHad_t1t2");
  graphNames1D.push_back("muHad_t2t3");
  graphNames1D.push_back("muHad_t3t4");
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
  graphNames1D.push_back("dRWMuSoftMu");
  graphNames1D.push_back("dRWMuTauHad");
  graphNames1D.push_back("dRTauMuTauHad");
  graphNames1D.push_back("dRWMuTauMuTauHad");
  graphNames1D.push_back("dPhiWMuSoftMu");
  graphNames1D.push_back("dPhiWMuSoftMu_withCut");
  graphNames1D.push_back("dPhiWMuSecJet");
  graphNames1D.push_back("WMu3rdTightMuChargeProduct");
  graphNames1D.push_back("tauHadPhotonEnergyFraction");
  graphNames1D.push_back("dThetaPhotonOtherTauConstituents");
  graphNames1D.push_back("hardestCorrJetEta");
  graphNames1D.push_back("WMuPT");
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
  graphNames2D.push_back("muHadMassVsWMuMT");
  graphNames2D.push_back("tauHadJetEnergyFractionVsTauHadIso");
  graphNames2D.push_back("tauHadCleanedJetEnergyFractionVsTauHadIso");
  graphNames2D.push_back("mu1ProngMassVsMu1Prong1Pi0Mass");
  graphNames2D.push_back("dPhiTauMuMETVsMuHadMass");
  vector<Int_t> nullBlindLow(canvasNames1D.size(), 0);
  vector<Int_t> nullBlindHigh(canvasNames1D.size(), -2);

  //flag for pseudoscalar mass
  const bool ma9GeV = (a1Mass == "_a9");

  //set up plot style options
  vector<string> legendHeaders19p7InvFb(canvasNames1D.size(), "Normalized to 19.7 fb^{-1}");
  vector<string> legendHeaders1(canvasNames1D.size(), "Normalized to 1");
  vector<Color_t> colorsSigBkg;
  colorsSigBkg.push_back(kSpring - 1); //WH, data
  colorsSigBkg.push_back(kBlue); //ggH
  if (ma9GeV) {
    colorsSigBkg.push_back(kSpring - 7); //ZH
    colorsSigBkg.push_back(kMagenta); //VBF
  }
  colorsSigBkg.push_back(kMagenta + 2); //WW
  colorsSigBkg.push_back(kCyan + 2); //ZZ
  colorsSigBkg.push_back(kRed + 2); //WZ
  colorsSigBkg.push_back(kYellow); //W + jets
  colorsSigBkg.push_back(kViolet - 7); //single top
  colorsSigBkg.push_back(kSpring + 4); //tt
  colorsSigBkg.push_back(kBlue + 1); //Drell-Yan
  colorsSigBkg.push_back(kGray + 2); //QCD
  colorsSigBkg.push_back(kMagenta - 2);
  colorsSigBkg.push_back(kGreen + 3);
  colorsSigBkg.push_back(kRed);
  vector<Style_t> styles;
  styles.push_back(20); //data
  styles.push_back(21); //WW
  styles.push_back(22); //ZZ
  styles.push_back(23); //WZ
  styles.push_back(24); //W + jets
  styles.push_back(25); //single top
  styles.push_back(26); //tt
  styles.push_back(27); //Drell-Yan
  styles.push_back(28); //QCD
  styles.push_back(29);
  styles.push_back(30);
  styles.push_back(31);
  styles.push_back(32);
  styles.push_back(33);
  styles.push_back(34);
  vector<string> legendEntriesSigBkg;
  legendEntriesSigBkg.push_back("WH");
  legendEntriesSigBkg.push_back("ggH");
  if (ma9GeV) {
    legendEntriesSigBkg.push_back("ZH");
    legendEntriesSigBkg.push_back("VBF");
  }
  legendEntriesSigBkg.push_back("WW");
  legendEntriesSigBkg.push_back("ZZ");
  legendEntriesSigBkg.push_back("WZ");
  legendEntriesSigBkg.push_back("W + #geq1 jet");
  legendEntriesSigBkg.push_back("t/#bar{t}");
  legendEntriesSigBkg.push_back("t#bar{t} + jets");
  legendEntriesSigBkg.push_back("Drell-Yan + jets");
  vector<string> legendEntriesSigBkgQCDFromData(legendEntriesSigBkg);
  legendEntriesSigBkgQCDFromData.push_back("QCD (from data)");
  const bool setLogY = true;
  const bool setLinY = false;
  const bool drawStack = true;
  const bool drawSame = false;
  const bool dataMC = true;
  const bool sigBkg = false;

  //best available weights according to Dropbox spreadsheet
  //const float Wh1a5Weight19p7InvFb = 0.00558166666666667; // Pythia xsec
  const float Wh1a5Weight19p7InvFb = 0.01507435331999; // SM xsec
  //const float Wh1a7Weight19p7InvFb = 0.00545988813757409; // Pythia xsec
  const float Wh1a7Weight19p7InvFb = 0.01510077968443; // SM xsec
  //const float Wh1a9Weight19p7InvFb = 0.00591; // Pythia xsec
  const float Wh1a9Weight19p7InvFb = 0.01507435332; // SM xsec
  //const float Wh1a11Weight19p7InvFb = 0.0000913508064516129; // Pythia xsec
  const float Wh1a11Weight19p7InvFb = 0.01519592068411; // SM xsec
  //const float Wh1a13Weight19p7InvFb = 0.0000570729270729271; // Pythia xsec
  const float Wh1a13Weight19p7InvFb = 0.01505930553923; // SM xsec
  //const float Wh1a15Weight19p7InvFb = 0.0000493486973947896; // Pythia xsec
  const float Wh1a15Weight19p7InvFb = 0.01510456244299; // SM xsec
  //const float gga5Weight19p7InvFb = 0.0552668358395599; // Pythia xsec
  const float gga5Weight19p7InvFb = 0.14925861038; // SM xsec
  //const float gga7Weight19p7InvFb = 0.0542857677477334; // Pythia xsec
  const float gga7Weight19p7InvFb = 0.15014179743; // SM xsec
  //const float gga9Weight19p7InvFb = 0.0584748299845453; // Pythia xsec
  const float gga9Weight19p7InvFb = 0.14914894205; // SM xsec
  //const float gga11Weight19p7InvFb = 0.000907280100622663; // Pythia xsec
  const float gga11Weight19p7InvFb = 0.15092320455; // SM xsec
  //const float gga13Weight19p7InvFb = 0.000567340969680438; // Pythia xsec
  const float gga13Weight19p7InvFb = 0.14969890127; // SM xsec
  //const float gga15Weight19p7InvFb = 0.00048980947568; // Pythia xsec
  const float gga15Weight19p7InvFb = 0.14992002225; // SM xsec
  const float ZHa9Weight19p7InvFb = 0.00277590622762097; /*SM ggH cross section + 100% 
							   BR(H-->aa-->4tau)*/
  const float VBFa9Weight19p7InvFb = 0.0132638419622011; /*SM ggH cross section + 100% 
							   BR(H-->aa-->4tau)*/
  float Wh1Weight19p7InvFb = 1.0;
  float ggWeight19p7InvFb = 1.0;
  float ZHWeight19p7InvFb = 1.0;
  float VBFWeight19p7InvFb = 1.0;
  if (a1Mass == "_a5") {
    Wh1Weight19p7InvFb = Wh1a5Weight19p7InvFb;
    ggWeight19p7InvFb = gga5Weight19p7InvFb;
  }
  else if (a1Mass == "_a7") {
    Wh1Weight19p7InvFb = Wh1a7Weight19p7InvFb;
    ggWeight19p7InvFb = gga7Weight19p7InvFb;
  }
  else if (ma9GeV) {
    Wh1Weight19p7InvFb = Wh1a9Weight19p7InvFb;
    ggWeight19p7InvFb = gga9Weight19p7InvFb;
    ZHWeight19p7InvFb = ZHa9Weight19p7InvFb;
    VBFWeight19p7InvFb = VBFa9Weight19p7InvFb;
  }
  else if (a1Mass == "_a11") {
    Wh1Weight19p7InvFb = Wh1a11Weight19p7InvFb;
    ggWeight19p7InvFb = gga11Weight19p7InvFb;
  }
  else if (a1Mass == "_a13") {
    Wh1Weight19p7InvFb = Wh1a13Weight19p7InvFb;
    ggWeight19p7InvFb = gga13Weight19p7InvFb;
  }
  else if (a1Mass == "_a15") {
    Wh1Weight19p7InvFb = Wh1a15Weight19p7InvFb;
    ggWeight19p7InvFb = gga15Weight19p7InvFb;
  }
  else cout << "Unrecognized a1 mass, weight defaulting to 1\n";

  //weights for getting VBF estimates from ggH samples
  const float VBFa5WeightFromggH = 0.01248931;
  const float VBFa7WeightFromggH = 0.0125632114;
  const float VBFa9WeightFromggH = 0.012480134;
  const float VBFa11WeightFromggH = 0.012628596;
  const float VBFa13WeightFromggH = 0.01256152;
  const float VBFa15WeightFromggH = 0.012544654;
  //weights for getting ZH estimates from WH samples
  const float ZHa5WeightFromWH = 0.002781679*1.1;
  const float ZHa7WeightFromWH = 0.002786556*1.1;
  const float ZHa9WeightFromWH = 0.002781679*1.1;
  const float ZHa11WeightFromWH = 0.0028041123*1.1;
  const float ZHa13WeightFromWH = 0.0027789*1.1;
  const float ZHa15WeightFromWH = 0.002787254*1.1;

  float VBFWeightFromggH = 1.0;
  float ZHWeightFromWH = 1.0;
  if (a1Mass == "_a5") {
    VBFWeightFromggH = VBFa5WeightFromggH;
    ZHWeightFromWH = ZHa5WeightFromWH;
  }
  else if (a1Mass == "_a7") {
    VBFWeightFromggH = VBFa7WeightFromggH;
    ZHWeightFromWH = ZHa7WeightFromWH;
  }
  else if (a1Mass == "_a9") {
    VBFWeightFromggH = VBFa9WeightFromggH;
    ZHWeightFromWH = ZHa9WeightFromWH;
  }
  else if (a1Mass == "_a11") {
    VBFWeightFromggH = VBFa11WeightFromggH;
    ZHWeightFromWH = ZHa11WeightFromWH;
  }
  else if (a1Mass == "_a13") {
    VBFWeightFromggH = VBFa13WeightFromggH;
    ZHWeightFromWH = ZHa13WeightFromWH;
  }
  else if (a1Mass == "_a15") {
    VBFWeightFromggH = VBFa15WeightFromggH;
    ZHWeightFromWH = ZHa15WeightFromWH;
  }

  vector<float> weights1(15, 0.0);
  vector<float> weightsSigBkg(11, 1.0); //MC samples already properly weighted during hadding
  vector<float> weightsSigBkgQCDFromData(12, 1.0);

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
  const string Wh1SigVTag("_" + inputVersion);
  const string ggSigVTag("_" + inputVersion);
  const string ZHSigVTag("_" + inputVersion);
  const string VBFSigVTag("_" + inputVersion);
  const string DYJetsToLLVTag("_" + inputVersion);
  const string TTJetsVTag("_" + inputVersion);
  const string TVTag("_" + inputVersion);
  const string WNJetsToLNuVTag("_" + inputVersion);
  const string WZVTag("_" + inputVersion);
  const string ZZVTag("_" + inputVersion);
  const string WWVTag("_" + inputVersion);
  const string narrowBinsVTag("_" + versionNarrow);

  cout << "Begin hadding...\n";

  //"hadd" Wh1 sample just to get the formatting of the 2D plots the same
  cout << "...Wh1\n";
  string Wh1Suffix(Wh1SigVTag + fileExt);
  string Wh1IsoPrefix(analysisFilePath + "Wh1_Medium/muHadIsoAnalysis" + MTBin + uncTag + "_Wh1" + 
		      a1Mass);
  string Wh1IsoHaddOutputFile(Wh1IsoPrefix + "_hadd" + Wh1Suffix);
  string Wh1AllPrefix(analysisFilePath + "Wh1_Medium/muHadAnalysis" + MTBin + uncTag + "_Wh1" + 
		      a1Mass);
  string Wh1AllHaddOutputFile(Wh1AllPrefix + "_hadd" + Wh1Suffix);
  vector<string> Wh1IsoHaddInputFiles;
  vector<string> Wh1AllHaddInputFiles;
  stringstream Wh1IsoName;
  Wh1IsoName << Wh1IsoPrefix << Wh1Suffix;
  Wh1IsoHaddInputFiles.push_back(Wh1IsoName.str());
  stringstream Wh1AllName;
  Wh1AllName << Wh1AllPrefix << Wh1Suffix;
  Wh1AllHaddInputFiles.push_back(Wh1AllName.str());
  haddCanvases(Wh1IsoHaddOutputFile, Wh1IsoHaddInputFiles, 
	       vector<float>(1, Wh1Weight19p7InvFb), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(Wh1AllHaddOutputFile, Wh1AllHaddInputFiles, 
		 vector<float>(1, Wh1Weight19p7InvFb), canvasNames1D, graphNames1D, 
		 canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //"hadd" gg sample just to get the formatting of the 2D plots the same
  cout << "...gg fusion\n";
  string ggSuffix(ggSigVTag + fileExt);
  string ggIsoPrefix(analysisFilePath + "gg/muHadIsoAnalysis" + MTBin + uncTag + "_gg" + a1Mass);
  string ggIsoHaddOutputFile(ggIsoPrefix + "_hadd" + ggSuffix);
  string ggAllPrefix(analysisFilePath + "gg/muHadAnalysis" + MTBin + uncTag + "_gg" + a1Mass);
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

  //get VBF estimates from hadded ggH sample
   cout << "...VBF estimates from ggH\n";
  string VBFEstSuffix(VBFSigVTag + fileExt);
  string 
    VBFIsoEstPrefix(analysisFilePath + "VBF/muHadIsoAnalysis" + MTBin + uncTag + "_VBFEst" + a1Mass);
  string VBFIsoEstHaddOutputFile(VBFIsoEstPrefix + "_hadd" + VBFEstSuffix);
  GetVBFZHEstimate(VBFIsoEstHaddOutputFile, ggIsoHaddOutputFile, ggWeight19p7InvFb, VBFWeightFromggH);

  //get ZH estimates from hadded WH sample
   cout << "...ZH estimates from WH\n";
  string ZHEstSuffix(ZHSigVTag + fileExt);
  string 
    ZHIsoEstPrefix(analysisFilePath + "ZH/muHadIsoAnalysis" + MTBin + uncTag + "_ZHEst" + a1Mass);
  string ZHIsoEstHaddOutputFile(ZHIsoEstPrefix + "_hadd" + ZHEstSuffix);
  GetVBFZHEstimate(ZHIsoEstHaddOutputFile, Wh1IsoHaddOutputFile, Wh1Weight19p7InvFb, ZHWeightFromWH);

  //"hadd" ZH sample just to get the formatting of the 2D plots the same
  string ZHSuffix(ZHSigVTag + fileExt);
  string ZHIsoPrefix(analysisFilePath + "ZH/muHadIsoAnalysis" + MTBin + uncTag + "_ZH" + a1Mass);
  string ZHIsoHaddOutputFile(ZHIsoPrefix + "_hadd" + ZHSuffix);
  string ZHAllPrefix(analysisFilePath + "ZH/muHadAnalysis" + MTBin + uncTag + "_ZH" + a1Mass);
  string ZHAllHaddOutputFile(ZHAllPrefix + "_hadd" + ZHSuffix);
  vector<string> ZHIsoHaddInputFiles;
  vector<string> ZHAllHaddInputFiles;
  stringstream ZHIsoName;
  ZHIsoName << ZHIsoPrefix << ZHSuffix;
  ZHIsoHaddInputFiles.push_back(ZHIsoName.str());
  stringstream ZHAllName;
  ZHAllName << ZHAllPrefix << ZHSuffix;
  ZHAllHaddInputFiles.push_back(ZHAllName.str());
  if (ma9GeV)
    {
      cout << "...ZH\n";
      haddCanvases(ZHIsoHaddOutputFile, ZHIsoHaddInputFiles, vector<float>(1, ZHWeight19p7InvFb), 
		   canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
		   nullBlindHigh);
      if (doNoHPSIsoCut) {
	haddCanvases(ZHAllHaddOutputFile, ZHAllHaddInputFiles, vector<float>(1, ZHWeight19p7InvFb), 
		     canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
		     nullBlindHigh);
      }
    }

  //"hadd" VBF sample just to get the formatting of the 2D plots the same
  string VBFSuffix(VBFSigVTag + fileExt);
  string 
    VBFIsoPrefix(analysisFilePath + "VBF/muHadIsoAnalysis" + MTBin + uncTag + "_VBF" + a1Mass);
  string VBFIsoHaddOutputFile(VBFIsoPrefix + "_hadd" + VBFSuffix);
  string VBFAllPrefix(analysisFilePath + "VBF/muHadAnalysis" + MTBin + uncTag + "_VBF" + a1Mass);
  string VBFAllHaddOutputFile(VBFAllPrefix + "_hadd" + VBFSuffix);
  vector<string> VBFIsoHaddInputFiles;
  vector<string> VBFAllHaddInputFiles;
  stringstream VBFIsoName;
  VBFIsoName << VBFIsoPrefix << VBFSuffix;
  VBFIsoHaddInputFiles.push_back(VBFIsoName.str());
  stringstream VBFAllName;
  VBFAllName << VBFAllPrefix << VBFSuffix;
  VBFAllHaddInputFiles.push_back(VBFAllName.str());
  if (ma9GeV)
    {
      cout << "...VBF\n";
      haddCanvases(VBFIsoHaddOutputFile, VBFIsoHaddInputFiles, vector<float>(1, VBFWeight19p7InvFb), 
		   canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
		   nullBlindHigh);
      if (doNoHPSIsoCut) {
	haddCanvases(VBFAllHaddOutputFile, VBFAllHaddInputFiles, vector<float>(1, VBFWeight19p7InvFb), 
		     canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, 
		     nullBlindHigh);
      }
    }
  
  cout << endl;

  //data region A files
  string dataSuffix(dataVTag + fileExt);
  string dataIsoPrefix(analysisFilePath + "data/analysis/muHadIsoAnalysis" + MTBin + "_SingleMu");
  string dataIsoHaddOutputFile(dataIsoPrefix + dataSuffix); //BLINDED!!!

  //DY region A files
  string DYJetsToLLSuffix(DYJetsToLLVTag + fileExt);
  string DYJetsToLLIsoPrefix(analysisFilePath + "DYJetsToLL/analysis/muHadIsoAnalysis" + MTBin + 
			     "_DYJetsToLL");
  string DYJetsToLLIsoHaddOutputFile(DYJetsToLLIsoPrefix + DYJetsToLLSuffix);
  string DYJetsToLLAllPrefix(analysisFilePath + "DYJetsToLL/analysis/muHadAnalysis" + MTBin + 
			     "_DYJetsToLL");
  string DYJetsToLLAllHaddOutputFile(DYJetsToLLAllPrefix + DYJetsToLLSuffix);

  //ttbar region A files
  string TTJetsSuffix(TTJetsVTag + fileExt);
  string TTJetsIsoPrefix(analysisFilePath + "TTJets/analysis/muHadIsoAnalysis" + MTBin + 
			 "_TTJets");
  string TTJetsIsoHaddOutputFile(TTJetsIsoPrefix + "_hadd" + TTJetsSuffix);
  string TTJetsAllPrefix(analysisFilePath + "TTJets/analysis/muHadAnalysis" + MTBin + "_TTJets");
  string TTJetsAllHaddOutputFile(TTJetsAllPrefix + "_hadd" + TTJetsSuffix);

  //single top region A files
  string TSuffix(TVTag + fileExt);
  string TIsoPrefix(analysisFilePath + "SingleTop/analysis/muHadIsoAnalysis" + MTBin + "_T");
  string TIsoHaddOutputFile(TIsoPrefix + TSuffix);
  string TAllPrefix(analysisFilePath + "SingleTop/analysis/muHadAnalysis" + MTBin + "_T");
  string TAllHaddOutputFile(TAllPrefix + TSuffix);

  //W + jets region A files
  string WNJetsToLNuSuffix("JetsToLNu" + WNJetsToLNuVTag + fileExt);
  string WNJetsToLNuIsoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadIsoAnalysis" + MTBin + 
			      "_W");
  string WNJetsToLNuIsoHaddOutputFile(WNJetsToLNuIsoPrefix + "N" + WNJetsToLNuSuffix);
  string WNJetsToLNuAllTauPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadAnalysis" + MTBin + 
				 "_W");
  string WNJetsToLNuAllTauHaddOutputFile(WNJetsToLNuAllTauPrefix + "N" + WNJetsToLNuSuffix);

  //WZ region A files
  string WZSuffix(WZVTag + fileExt);
  string WZIsoPrefix(analysisFilePath + "WZ/analysis/muHadIsoAnalysis" + MTBin + "_WZ");
  string WZIsoHaddOutputFile(WZIsoPrefix + "_hadd" + WZSuffix);
  string WZAllPrefix(analysisFilePath + "WZ/analysis/muHadAnalysis" + MTBin + "_WZ");
  string WZAllHaddOutputFile(WZAllPrefix + "_hadd" + WZSuffix);

  //ZZ region A files
  string ZZSuffix(ZZVTag + fileExt);
  string ZZIsoPrefix(analysisFilePath + "ZZ/analysis/muHadIsoAnalysis" + MTBin + "_ZZ");
  string ZZIsoHaddOutputFile(ZZIsoPrefix + "_hadd" + ZZSuffix);
  string ZZAllPrefix(analysisFilePath + "ZZ/analysis/muHadAnalysis" + MTBin + "_ZZ");
  string ZZAllHaddOutputFile(ZZAllPrefix + "_hadd" + ZZSuffix);

  //WW region A files
  string WWSuffix(WWVTag + fileExt);
  string WWIsoPrefix(analysisFilePath + "WW/analysis/muHadIsoAnalysis" + MTBin + "_WW");
  string WWIsoHaddOutputFile(WWIsoPrefix + "_hadd" + WWSuffix);
  string WWAllPrefix(analysisFilePath + "WW/analysis/muHadAnalysis" + MTBin + "_WW");
  string WWAllHaddOutputFile(WWAllPrefix + "_hadd" + WWSuffix);

  //compare MC signal to background
  string sigVsBkgOutputFile(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + MTBin + 
			    uncTag + a1Mass + tag19p7InvFb + outputVTag + fileExt);
  string sigVsBkgOutputFile1(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + MTBin + 
			     tag1 + uncTag + a1Mass + outputVTag + fileExt);
  string sigVsBkgOutputFileNoHPSIsoCut(analysisFilePath + "results/sigVsBkg_muHadAnalysis" + 
				       MTBin + uncTag + a1Mass + tag19p7InvFb + outputVTag + 
				       fileExt);
  string sigVsBkgOutputFileNoHPSIsoCutNorm1(analysisFilePath + "results/sigVsBkg_muHadAnalysis" + 
					    MTBin + uncTag + a1Mass + tag1 + outputVTag + fileExt);
  vector<string> sigVsBkgInputFiles;
  sigVsBkgInputFiles.push_back(Wh1IsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(ggIsoHaddOutputFile);
  if (ma9GeV) {
    sigVsBkgInputFiles.push_back(ZHIsoHaddOutputFile);
    sigVsBkgInputFiles.push_back(VBFIsoHaddOutputFile);
  }
  sigVsBkgInputFiles.push_back(WWIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(ZZIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(WZIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(WNJetsToLNuIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(TIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(TTJetsIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(DYJetsToLLIsoHaddOutputFile);
  vector<string> sigVsBkgNoHPSIsoCutInputFiles;
  sigVsBkgNoHPSIsoCutInputFiles.push_back(Wh1AllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(ggAllHaddOutputFile);
  if (ma9GeV) {
    sigVsBkgNoHPSIsoCutInputFiles.push_back(ZHAllHaddOutputFile);
    sigVsBkgNoHPSIsoCutInputFiles.push_back(VBFAllHaddOutputFile);
  }
  sigVsBkgNoHPSIsoCutInputFiles.push_back(WWAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(ZZAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(WZAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(WNJetsToLNuAllTauHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(TAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(TTJetsAllHaddOutputFile);
  sigVsBkgNoHPSIsoCutInputFiles.push_back(DYJetsToLLAllHaddOutputFile);
  cout << "Plot signal vs. background normalized to data luminosity\n---\n";
  if (uncTag == "") {
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile, sigVsBkgInputFiles, 
					  canvasNames1D, graphNames1D, legendHeaders19p7InvFb, 
					  colorsSigBkg, styles, legendEntriesSigBkg, 
					  weightsSigBkg, setLogY, drawStack, sigBkg);
    cout << "\nPlot signal vs. background normalized to 1\n---\n";
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile1, sigVsBkgInputFiles, 
					  canvasNames1D, graphNames1D, legendHeaders1, 
					  colorsSigBkg, styles, legendEntriesSigBkg, weights1, 
					  setLinY, drawSame, sigBkg);
  }
  else {
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile, sigVsBkgInputFiles, 
					  vector<string>(1, "muHadMassCanvas"), 
					  vector<string>(1, "muHadMass"), legendHeaders19p7InvFb, 
					  colorsSigBkg, styles, legendEntriesSigBkg, 
					  weightsSigBkg, setLogY, drawStack, sigBkg);
    cout << "\nPlot signal vs. background normalized to 1\n---\n";
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile1, sigVsBkgInputFiles, 
					  vector<string>(1, "muHadMassCanvas"), 
					  vector<string>(1, "muHadMass"), legendHeaders1, 
					  colorsSigBkg, styles, legendEntriesSigBkg, weights1, 
					  setLinY, drawSame, sigBkg);
  }
  if (doNoHPSIsoCut) {
    cout << "\nPlot signal vs. background normalized to data luminosity, ";
    cout << "no cut on tau isolation\n---\n";
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFileNoHPSIsoCut, 
					  sigVsBkgNoHPSIsoCutInputFiles, canvasNames1D, 
					  graphNames1D, legendHeaders19p7InvFb, colorsSigBkg, 
					  styles, legendEntriesSigBkg, weightsSigBkg, setLogY, 
					  drawStack, sigBkg);
    cout << "\nPlot signal vs. background normalized to 1, no cut on tau isolation\n---\n";
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFileNoHPSIsoCutNorm1, 
					  sigVsBkgNoHPSIsoCutInputFiles, canvasNames1D, 
					  graphNames1D, legendHeaders1, colorsSigBkg, styles, 
					  legendEntriesSigBkg, weights1, setLinY, drawSame, 
					  sigBkg);
  }

  //compare MC signal + data-driven QCD to background
  string outputFileNameA(analysisFilePath + "results/dataVsMC_RegionAQCDEstimate" + MTBin + 
			 dataVTag + fileExt);
  string sigVsBkgQCDFromDataOutputFile(analysisFilePath + 
				       "results/sigVsBkgQCDFromData_muHadIsoAnalysis" + MTBin + 
				       uncTag + a1Mass + tag19p7InvFb + outputVTag + fileExt);
  string sigVsBkgQCDFromDataOutputFile1(analysisFilePath + 
					"results/sigVsBkgQCDFromData_muHadIsoAnalysis" + MTBin + 
					uncTag + a1Mass + tag1 + outputVTag + fileExt);
  vector<string> sigVsBkgQCDFromDataInputFiles(sigVsBkgInputFiles);
  sigVsBkgQCDFromDataInputFiles.push_back(outputFileNameA);
  cout << "\nPlot signal vs. background with data-driven QCD estimate ";
  cout << "normalized to data luminosity\n---\n";
  if (uncTag == "") {
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgQCDFromDataOutputFile, 
					  sigVsBkgQCDFromDataInputFiles, canvasNames1D, 
					  graphNames1D, legendHeaders19p7InvFb, colorsSigBkg, 
					  styles, legendEntriesSigBkgQCDFromData, 
					  weightsSigBkgQCDFromData, setLogY, drawStack, sigBkg);
    cout << "\nPlot signal vs. background with data-driven QCD estimate normalized to 1\n---\n";
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgQCDFromDataOutputFile1, 
					  sigVsBkgQCDFromDataInputFiles, canvasNames1D, 
					  graphNames1D, legendHeaders1, colorsSigBkg, styles, 
					  legendEntriesSigBkgQCDFromData, weights1, setLinY, 
					  drawSame, sigBkg);
  }
  else {
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgQCDFromDataOutputFile, 
					  sigVsBkgQCDFromDataInputFiles, 
					  vector<string>(1, "muHadMassCanvas"), 
					  vector<string>(1, "muHadMass"), legendHeaders19p7InvFb, 
					  colorsSigBkg, styles, legendEntriesSigBkgQCDFromData, 
					  weightsSigBkgQCDFromData, setLogY, drawStack, sigBkg);
    cout << "\nPlot signal vs. background with data-driven QCD estimate normalized to 1\n---\n";
    drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgQCDFromDataOutputFile1, 
					  sigVsBkgQCDFromDataInputFiles, 
					  vector<string>(1, "muHadMassCanvas"), 
					  vector<string>(1, "muHadMass"), legendHeaders1, 
					  colorsSigBkg, styles, legendEntriesSigBkgQCDFromData, 
					  weights1, setLinY, drawSame, sigBkg);
  }

  cout << "---\nMaking final plots\n";

  //make the final plot showing all background methods, signals, data, and errors
  string resBkgOutputFile(analysisFilePath + "results/resBkg" + MTBin + tag19p7InvFb + 
			  narrowBinsVTag + fileExt);
  string dataVsMCOutputFile(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + MTBin + 
			    tag19p7InvFb + outputVTag + fileExt);
  const string 
    HLTMu40eta2p1RegDInputFile(analysisFilePath + "nonIsoWDataHLTMu40eta2p1/analysis/" + 
			       "nonIsoW_HLT_Mu40_eta2p1_muHadNonIsoAnalysis" + MTBin + 
			       "_SingleMu" + nonIsoWDataVTag + fileExt);
  string nonIsoWDataSuffix(nonIsoWDataVTag + fileExt);
  string nonIsoWDataNonIsoPrefix(analysisFilePath + 
				 "nonIsoWData/analysis/nonIsoW_muHadNonIsoAnalysis" + MTBin + 
				 "_SingleMu");
  string nonIsoWDataNonIsoHaddOutputFile(nonIsoWDataNonIsoPrefix + nonIsoWDataSuffix);
  gStyle->SetErrorX(1);
  makeFinalPlot(pair<string, float>(sigVsBkgQCDFromDataOutputFile, 1.0), 
  		dataIsoHaddOutputFile, 
  		pair<string, float>(dataVsMCOutputFile, 1.0), 
  		pair<string, float>(resBkgOutputFile, 1.0), 
		pair<string, float>(nonIsoWDataNonIsoHaddOutputFile, 1.0), 
		vector<string>(1, "muHadMass"), vector<string>(1, "m_{#mu+had} (GeV)"), 
  		vector<int>(1, 1), vector<int>(1, firstBinToBlind - 1), 
  		analysisFilePath + "results/final" + MTBin + uncTag + a1Mass + outputVTag + 
		fileExt, "main 5", ma9GeV);

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
  normRegionUpperBins.push_back(firstBinToBlind - 1);
  normRegionUpperBins.push_back(3);
  normRegionUpperBins.push_back(46);
  normRegionUpperBins.push_back(50);
  normRegionUpperBins.push_back(9);
  if (uncTag == "") {
    makeMCClosurePlots(sigVsBkgOutputFile, vars, units, dataVsMCOutputFile, 1.0, 
		       normRegionLowerBins, normRegionUpperBins, 
		       analysisFilePath + "results/MC_closure" + MTBin + uncTag + outputVTag + 
		       fileExt);
  }

//   //compare the same plot from 2 versions of the analysis
//   string nonIsoWDataIsoPrefix(analysisFilePath + 
// 			      "nonIsoWData/analysis/nonIsoW_muHadIsoAnalysis" + MTBin + 
// 			      "_SingleMu");
//   string nonIsoWDataNonIsoPrefix(analysisFilePath + 
// 				 "nonIsoWData/analysis/nonIsoW_muHadNonIsoAnalysis" + MTBin + 
// 				 "_SingleMu");
//   vector<string> fileNamesForComparison1;
//   fileNamesForComparison1.push_back(analysisFilePath + 
// 				    "results/dataVsMCQCDFromData_muHadNonIsoAnalysis" + MTBin +
//   				    tag19p7InvFb + "_v35" + fileExt);
//   fileNamesForComparison1.push_back(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + 
//   				    tag19p7InvFb + "_v149" + fileExt);
//   fileNamesForComparison1.push_back(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + 
// 				    MTBin + uncTag + a1Mass + tag19p7InvFb + "_v149" + fileExt);
//   fileNamesForComparison1.push_back(nonIsoWDataIsoPrefix + "_v149" + fileExt);
//   fileNamesForComparison1.push_back(nonIsoWDataNonIsoPrefix + "_v149" + fileExt);
//   vector<string> fileNamesForComparison2;
//   fileNamesForComparison2.push_back(analysisFilePath + 
// 				    "results/dataVsMCQCDFromData_muHadNonIsoAnalysis" + MTBin +
// 				    tag19p7InvFb + "_v41" + fileExt);
//   fileNamesForComparison2.push_back(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + 
//   				    tag19p7InvFb + "_v152" + fileExt);
//   fileNamesForComparison2.push_back(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + 
//   				    MTBin + uncTag + a1Mass + tag19p7InvFb + "_v152" + fileExt);
//   fileNamesForComparison2.push_back(nonIsoWDataIsoPrefix + "_v152" + fileExt);
//   fileNamesForComparison2.push_back(nonIsoWDataNonIsoPrefix + "_v152" + fileExt);
//   vector<string> outputCanvasTagsForComparison;
//   outputCanvasTagsForComparison.push_back("_regionBData");
//   outputCanvasTagsForComparison.push_back("_regionBMC");
//   outputCanvasTagsForComparison.push_back("_regionAMC");
//   outputCanvasTagsForComparison.push_back("_regionCData");
//   outputCanvasTagsForComparison.push_back("_regionDData");
//   vector<bool> stack;
//   stack.push_back(false);
//   stack.push_back(true);
//   stack.push_back(true);
//   stack.push_back(false);
//   stack.push_back(false);
//   vector<unsigned int> pad;
//   pad.push_back(1);
//   pad.push_back(1);
//   pad.push_back(1);
//   pad.push_back(0);
//   pad.push_back(0);
//   compare2Versions(fileNamesForComparison1, fileNamesForComparison2, analysisFilePath + 
// 		   "results/comparison_v35_v41" + fileExt, outputCanvasTagsForComparison, stack, 
//     		   pad);
}
