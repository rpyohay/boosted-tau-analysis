//REGION A DATA 2D HISTOGRAMS ARE NOT BLINDED!!!  BEWARE!!!

void formatBoostedResonanceDataBkgPlots(const string& inputVersion, const string& outputVersion, 
					const string& MTBin, const unsigned int firstBinToBlind, 
					const string& HLTPath, const bool doNoHPSIsoCut = false)
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

  //blind the data mu+X mass histograms above 2 GeV
  vector<Int_t> nullBlindLow(canvasNames1D.size(), 0);
  vector<Int_t> nullBlindHigh(canvasNames1D.size(), -2);
  vector<Int_t> dataBlindLow(canvasNames1D.size(), 0);
  for (unsigned int iPlot = 1; iPlot <= 11; ++iPlot) { dataBlindLow[iPlot] = firstBinToBlind; }
  vector<Int_t> dataBlindHigh(canvasNames1D.size(), -1);

  //set up plot style options
  vector<string> legendHeaders19p7InvFb(canvasNames1D.size(), "Normalized to 19.7 fb^{-1}");
  vector<string> legendHeaders1(canvasNames1D.size(), "Normalized to 1");
  vector<string> legendHeaders1DYJetsToLL(canvasNames1D.size(), 
					  "Drell-Yan + jets normalized to 1");
  vector<string> legendHeaders1TTJets(canvasNames1D.size(), "t#bar{t} + jets normalized to 1");
  vector<string> legendHeaders1T(canvasNames1D.size(), "t/#bar{t} normalized to 1");
  vector<string> legendHeaders1WNJetsToLNu(canvasNames1D.size(), "W + #geq1 jet normalized to 1");
  vector<string> legendHeaders1WJetsToLNu(canvasNames1D.size(), "W + jets normalized to 1");
  vector<string> legendHeaders1WZ(canvasNames1D.size(), "WZ normalized to 1");
  vector<string> legendHeaders1ZZ(canvasNames1D.size(), "ZZ normalized to 1");
  vector<string> legendHeaders1WW(canvasNames1D.size(), "WW normalized to 1");
  vector<string> legendHeaders1NonIsoWData(canvasNames1D.size(), 
					   "Non-isolated W data normalized to 1");
  vector<Color_t> colorsMCData;
  colorsMCData.push_back(kBlack); //data
  colorsMCData.push_back(kMagenta + 2); //WW
  colorsMCData.push_back(kCyan + 2); //ZZ
  colorsMCData.push_back(kRed + 2); //WZ
  colorsMCData.push_back(kYellow); //W + jets
  colorsMCData.push_back(kViolet - 7); //single top
  colorsMCData.push_back(kSpring + 4); //tt
  colorsMCData.push_back(kBlue + 1); //Drell-Yan
  colorsMCData.push_back(kGray + 2); //QCD
  colorsMCData.push_back(kMagenta - 2);
  colorsMCData.push_back(kGreen + 3);
  colorsMCData.push_back(kRed);
  vector<Color_t> colorsNonIsoWMCData;
  colorsNonIsoWMCData.push_back(kBlack); //data
  colorsNonIsoWMCData.push_back(kBlue + 1); //Drell-Yan
  colorsNonIsoWMCData.push_back(kSpring + 4); //tt
  colorsNonIsoWMCData.push_back(kYellow); //W + jets
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
  vector<string> legendEntriesMCData;
  legendEntriesMCData.push_back("Data 19.7 fb^{-1}");
  legendEntriesMCData.push_back("WW");
  legendEntriesMCData.push_back("ZZ");
  legendEntriesMCData.push_back("WZ");
  legendEntriesMCData.push_back("W + #geq1 jet");
  legendEntriesMCData.push_back("t/#bar{t}");
  legendEntriesMCData.push_back("t#bar{t} + jets");
  legendEntriesMCData.push_back("Drell-Yan + jets");
  vector<string> legendEntriesMCDataQCDFromData(legendEntriesMCData);
  legendEntriesMCDataQCDFromData.push_back("QCD (from data)");
  vector<string> legendEntriesSearchVsControl;
  legendEntriesSearchVsControl.push_back("Isolated #tau leptons");
  legendEntriesSearchVsControl.push_back("Non-isolated #tau leptons");
  vector<string> legendEntriesNonIsoWMCData;
  legendEntriesNonIsoWMCData.push_back("Data 19.7 fb^{-1}");
  legendEntriesNonIsoWMCData.push_back("Drell-Yan + jets");
  legendEntriesNonIsoWMCData.push_back("t#bar{t} + jets");
  legendEntriesNonIsoWMCData.push_back("W + #geq1 jet");
  const bool setLogY = true;
  const bool setLinY = false;
  const bool drawStack = true;
  const bool drawSame = false;
  const bool dataMC = true;
  const bool sigBkg = false;

  //best available weights according to Dropbox spreadsheet
  vector<float> weights1(15, 0.0);
  vector<float> weightsMCData(8, 1.0); //MC samples already properly weighted during hadding
  vector<float> weightsMCDataQCDFromData(9, 1.0);
  vector<float> weightsNonIsoWMCData(4, 1.0);
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
  //TRelXSecWeights.push_back(0.34672); //tbar s-channel
  //TRelXSecWeights.push_back(0.321286023155796); //t t-channel
  TRelXSecWeights.push_back(0.2477031449); //tbar s-channel
  TRelXSecWeights.push_back(0.2956394066); //t t-channel
  TRelXSecWeights.push_back(0.312541342130939); //tbar t-channel
  vector<float> WNJetsToLNuRelXSecWeights; //weighted to 19.7 fb^-1 using PREP cross sections
  WNJetsToLNuRelXSecWeights.push_back(5.61973075498071); //W + 1 jet
  //  WNJetsToLNuRelXSecWeights.push_back(1.22833526483929); //W + 2 jets
  WNJetsToLNuRelXSecWeights.push_back(1.2211192971); //W + 2 jets
  //  WNJetsToLNuRelXSecWeights.push_back(0.819050332546934); //W + 3 jets
  WNJetsToLNuRelXSecWeights.push_back(0.8032380444); //W + 3 jets
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
  const string DYJetsToLLVTag("_" + inputVersion);
  const string TTJetsVTag("_" + inputVersion);
  const string TVTag("_" + inputVersion);
  const string WNJetsToLNuVTag("_" + inputVersion);
  const string WJetsToLNuVTag("_" + inputVersion);
  const string WZVTag("_" + inputVersion);
  const string ZZVTag("_" + inputVersion);
  const string WWVTag("_" + inputVersion);
  const string NonIsoWDYJetsToLLVTag("_" + inputVersion);
  const string NonIsoWTTJetsVTag("_" + inputVersion);
  const string NonIsoWWNJetsToLNuVTag("_" + inputVersion);

  cout << "Begin hadding...\n";

  //hadd data samples from different eras
  cout << "...data\n";
  string dataSuffix(dataVTag + fileExt);
  string dataIsoPrefix(analysisFilePath + "data/analysis/muHadIsoAnalysis" + MTBin + "_SingleMu");
  string dataIsoHaddOutputFile(dataIsoPrefix + dataSuffix); //BLINDED!!!
  string dataNonIsoPrefix(analysisFilePath + "data/analysis/muHadNonIsoAnalysis" + MTBin + 
			  "_SingleMu");
  string dataNonIsoHaddOutputFile(dataNonIsoPrefix + dataSuffix);
  string dataAllPrefix(analysisFilePath + "data/analysis/muHadAnalysis" + MTBin + "_SingleMu");
  string dataAllHaddOutputFile(dataAllPrefix + dataSuffix); //BLINDED!!!
  vector<string> dataIsoHaddInputFiles; //BLINDED!!!
  vector<string> dataNonIsoHaddInputFiles;
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
      dataIsoName << dataIsoPrefix << *iRunEra << *iSubJob << dataSuffix;
      dataIsoHaddInputFiles.push_back(dataIsoName.str());
      stringstream dataNonIsoName;
      dataNonIsoName << dataNonIsoPrefix << *iRunEra << *iSubJob << dataSuffix;
      dataNonIsoHaddInputFiles.push_back(dataNonIsoName.str());
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
  if (doNoHPSIsoCut) {
    haddCanvases(dataAllHaddOutputFile, dataAllHaddInputFiles, 
		 vector<float>(32, 1.0), canvasNames1D, graphNames1D, canvasNames2D, 
		 graphNames2D, dataBlindLow, dataBlindHigh); //BLINDED!!!
  }

  //hadd non-isolated W data samples from different eras
  cout << "...non-isolated W data\n";
  string nonIsoWDataSuffix(nonIsoWDataVTag + fileExt);
  string nonIsoWDataIsoPrefix(analysisFilePath + 
			      "nonIsoWData/analysis/nonIsoW_muHadIsoAnalysis" + MTBin + 
			      "_SingleMu");
  string nonIsoWDataIsoHaddOutputFile(nonIsoWDataIsoPrefix + nonIsoWDataSuffix);
  string nonIsoWDataNonIsoPrefix(analysisFilePath + 
				 "nonIsoWData/analysis/nonIsoW_muHadNonIsoAnalysis" + MTBin + 
				 "_SingleMu");
  string nonIsoWDataNonIsoHaddOutputFile(nonIsoWDataNonIsoPrefix + nonIsoWDataSuffix);
  string nonIsoWDataAllPrefix(analysisFilePath + 
			      "nonIsoWData/analysis/nonIsoW_muHadAnalysis" + MTBin + "_SingleMu");
  string nonIsoWDataAllHaddOutputFile(nonIsoWDataAllPrefix + nonIsoWDataSuffix);
  vector<string> nonIsoWDataIsoHaddInputFiles;
  vector<string> nonIsoWDataNonIsoHaddInputFiles;
  vector<string> nonIsoWDataAllHaddInputFiles;
  for (vector<string>::const_iterator iRunEra = runEras.begin(); iRunEra != runEras.end(); 
       ++iRunEra) {
    stringstream nonIsoWDataIsoName;
    nonIsoWDataIsoName << nonIsoWDataIsoPrefix << *iRunEra << nonIsoWDataSuffix;
    nonIsoWDataIsoHaddInputFiles.push_back(nonIsoWDataIsoName.str());
    stringstream nonIsoWDataNonIsoName;
    nonIsoWDataNonIsoName << nonIsoWDataNonIsoPrefix << *iRunEra << nonIsoWDataSuffix;
    nonIsoWDataNonIsoHaddInputFiles.push_back(nonIsoWDataNonIsoName.str());
    stringstream nonIsoWDataAllName;
    nonIsoWDataAllName << nonIsoWDataAllPrefix << *iRunEra << nonIsoWDataSuffix;
    nonIsoWDataAllHaddInputFiles.push_back(nonIsoWDataAllName.str());
  }
  haddCanvases(nonIsoWDataIsoHaddOutputFile, nonIsoWDataIsoHaddInputFiles, 
	       vector<float>(4, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
	       nullBlindLow, nullBlindHigh);
  haddCanvases(nonIsoWDataNonIsoHaddOutputFile, nonIsoWDataNonIsoHaddInputFiles, 
	       vector<float>(4, 1.0), canvasNames1D, graphNames1D, canvasNames2D, graphNames2D, 
	       nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(nonIsoWDataAllHaddOutputFile, nonIsoWDataAllHaddInputFiles, 
		 vector<float>(4, 1.0), canvasNames1D, graphNames1D, canvasNames2D, 
		 graphNames2D, nullBlindLow, nullBlindHigh);
  }
    
  //hadd Drell-Yan+jets ml+l- binned samples
  cout << "...Drell-Yan\n";
  string DYJetsToLLSuffix(DYJetsToLLVTag + fileExt);
  string DYJetsToLLIsoPrefix(analysisFilePath + "DYJetsToLL/analysis/muHadIsoAnalysis" + MTBin + 
			     "_DYJetsToLL");
  string DYJetsToLLIsoHaddOutputFile(DYJetsToLLIsoPrefix + DYJetsToLLSuffix);
  string DYJetsToLLNonIsoPrefix(analysisFilePath + 
				"DYJetsToLL/analysis/muHadNonIsoAnalysis" + MTBin + "_DYJetsToLL");
  string DYJetsToLLNonIsoHaddOutputFile(DYJetsToLLNonIsoPrefix + DYJetsToLLSuffix);
  string DYJetsToLLAllPrefix(analysisFilePath + "DYJetsToLL/analysis/muHadAnalysis" + MTBin + 
			     "_DYJetsToLL");
  string DYJetsToLLAllHaddOutputFile(DYJetsToLLAllPrefix + DYJetsToLLSuffix);
  vector<string> DYJetsToLLIsoHaddInputFiles;
  vector<string> DYJetsToLLNonIsoHaddInputFiles;
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
    stringstream DYJetsToLLAllName;
    DYJetsToLLAllName << DYJetsToLLAllPrefix << *iMassBin << DYJetsToLLSuffix;
    DYJetsToLLAllHaddInputFiles.push_back(DYJetsToLLAllName.str());
  }
  haddCanvases(DYJetsToLLIsoHaddOutputFile, DYJetsToLLIsoHaddInputFiles, 
	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
	       graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(DYJetsToLLNonIsoHaddOutputFile, DYJetsToLLNonIsoHaddInputFiles, 
	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
	       graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(DYJetsToLLAllHaddOutputFile, DYJetsToLLAllHaddInputFiles, 
		 DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
		 graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //"hadd" ttbar sample just to get the formatting of the 2D plots the same
  cout << "...ttbar\n";
  string TTJetsSuffix(TTJetsVTag + fileExt);
  string TTJetsIsoPrefix(analysisFilePath + "TTJets/analysis/muHadIsoAnalysis" + MTBin + 
			 "_TTJets");
  string TTJetsIsoHaddOutputFile(TTJetsIsoPrefix + "_hadd" + TTJetsSuffix);
  string TTJetsNonIsoPrefix(analysisFilePath + "TTJets/analysis/muHadNonIsoAnalysis" + MTBin + 
			    "_TTJets");
  string TTJetsNonIsoHaddOutputFile(TTJetsNonIsoPrefix + "_hadd" + TTJetsSuffix);
  string TTJetsAllPrefix(analysisFilePath + "TTJets/analysis/muHadAnalysis" + MTBin + "_TTJets");
  string TTJetsAllHaddOutputFile(TTJetsAllPrefix + "_hadd" + TTJetsSuffix);
  vector<string> TTJetsIsoHaddInputFiles;
  vector<string> TTJetsNonIsoHaddInputFiles;
  vector<string> TTJetsAllHaddInputFiles;
  stringstream TTJetsIsoName;
  TTJetsIsoName << TTJetsIsoPrefix << TTJetsSuffix;
  TTJetsIsoHaddInputFiles.push_back(TTJetsIsoName.str());
  stringstream TTJetsNonIsoName;
  TTJetsNonIsoName << TTJetsNonIsoPrefix << TTJetsSuffix;
  TTJetsNonIsoHaddInputFiles.push_back(TTJetsNonIsoName.str());
  stringstream TTJetsAllName;
  TTJetsAllName << TTJetsAllPrefix << TTJetsSuffix;
  TTJetsAllHaddInputFiles.push_back(TTJetsAllName.str());
  haddCanvases(TTJetsIsoHaddOutputFile, TTJetsIsoHaddInputFiles, 
	       vector<float>(1, 3.54800726562391), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(TTJetsNonIsoHaddOutputFile, TTJetsNonIsoHaddInputFiles, 
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
  string TIsoPrefix(analysisFilePath + "SingleTop/analysis/muHadIsoAnalysis" + MTBin + "_T");
  string TIsoHaddOutputFile(TIsoPrefix + TSuffix);
  string TNonIsoPrefix(analysisFilePath + "SingleTop/analysis/muHadNonIsoAnalysis" + MTBin + "_T");
  string TNonIsoHaddOutputFile(TNonIsoPrefix + TSuffix);
  string TAllPrefix(analysisFilePath + "SingleTop/analysis/muHadAnalysis" + MTBin + "_T");
  string TAllHaddOutputFile(TAllPrefix + TSuffix);
  vector<string> TIsoHaddInputFiles;
  vector<string> TNonIsoHaddInputFiles;
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
    stringstream TAllName;
    TAllName << TAllPrefix << *iSample << TSuffix;
    TAllHaddInputFiles.push_back(TAllName.str());
  }
  haddCanvases(TIsoHaddOutputFile, TIsoHaddInputFiles, TRelXSecWeights, canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(TNonIsoHaddOutputFile, TNonIsoHaddInputFiles, TRelXSecWeights, canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(TAllHaddOutputFile, TAllHaddInputFiles, TRelXSecWeights, canvasNames1D, 
		 graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //hadd W+>=1 jet samples
  cout << "...W+>=1 jet\n";
  string WNJetsToLNuSuffix("JetsToLNu" + WNJetsToLNuVTag + fileExt);
  string WNJetsToLNuIsoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadIsoAnalysis" + MTBin + 
			      "_W");
  string WNJetsToLNuIsoHaddOutputFile(WNJetsToLNuIsoPrefix + "N" + WNJetsToLNuSuffix);
  string WNJetsToLNuNonIsoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadNonIsoAnalysis" + 
				 MTBin + "_W");
  string WNJetsToLNuNonIsoHaddOutputFile(WNJetsToLNuNonIsoPrefix + "N" + WNJetsToLNuSuffix);
  string WNJetsToLNuAllTauPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadAnalysis" + MTBin + 
				 "_W");
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
	       graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WNJetsToLNuNonIsoHaddOutputFile, WNJetsToLNuNonIsoHaddInputFiles, 
	       WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
	       graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(WNJetsToLNuAllTauHaddOutputFile, WNJetsToLNuAllTauHaddInputFiles, 
		 WNJetsToLNuRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, 
		 graphNames2D, nullBlindLow, nullBlindHigh);
  }

  //"hadd" WZ sample just to get the formatting of the 2D plots the same
  cout << "...WZ\n";
  string WZSuffix(WZVTag + fileExt);
  string WZIsoPrefix(analysisFilePath + "WZ/analysis/muHadIsoAnalysis" + MTBin + "_WZ");
  string WZIsoHaddOutputFile(WZIsoPrefix + "_hadd" + WZSuffix);
  string WZNonIsoPrefix(analysisFilePath + "WZ/analysis/muHadNonIsoAnalysis" + MTBin + "_WZ");
  string WZNonIsoHaddOutputFile(WZNonIsoPrefix + "_hadd" + WZSuffix);
  string WZAllPrefix(analysisFilePath + "WZ/analysis/muHadAnalysis" + MTBin + "_WZ");
  string WZAllHaddOutputFile(WZAllPrefix + "_hadd" + WZSuffix);
  vector<string> WZIsoHaddInputFiles;
  vector<string> WZNonIsoHaddInputFiles;
  vector<string> WZAllHaddInputFiles;
  stringstream WZIsoName;
  WZIsoName << WZIsoPrefix << WZSuffix;
  WZIsoHaddInputFiles.push_back(WZIsoName.str());
  stringstream WZNonIsoName;
  WZNonIsoName << WZNonIsoPrefix << WZSuffix;
  WZNonIsoHaddInputFiles.push_back(WZNonIsoName.str());
  stringstream WZAllName;
  WZAllName << WZAllPrefix << WZSuffix;
  WZAllHaddInputFiles.push_back(WZAllName.str());
  haddCanvases(WZIsoHaddOutputFile, WZIsoHaddInputFiles, 
	       vector<float>(1, 0.0667569497737973), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WZNonIsoHaddOutputFile, WZNonIsoHaddInputFiles, 
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
  string ZZIsoPrefix(analysisFilePath + "ZZ/analysis/muHadIsoAnalysis" + MTBin + "_ZZ");
  string ZZIsoHaddOutputFile(ZZIsoPrefix + "_hadd" + ZZSuffix);
  string ZZNonIsoPrefix(analysisFilePath + "ZZ/analysis/muHadNonIsoAnalysis" + MTBin + "_ZZ");
  string ZZNonIsoHaddOutputFile(ZZNonIsoPrefix + "_hadd" + ZZSuffix);
  string ZZAllPrefix(analysisFilePath + "ZZ/analysis/muHadAnalysis" + MTBin + "_ZZ");
  string ZZAllHaddOutputFile(ZZAllPrefix + "_hadd" + ZZSuffix);
  vector<string> ZZIsoHaddInputFiles;
  vector<string> ZZNonIsoHaddInputFiles;
  vector<string> ZZAllHaddInputFiles;
  stringstream ZZIsoName;
  ZZIsoName << ZZIsoPrefix << ZZSuffix;
  ZZIsoHaddInputFiles.push_back(ZZIsoName.str());
  stringstream ZZNonIsoName;
  ZZNonIsoName << ZZNonIsoPrefix << ZZSuffix;
  ZZNonIsoHaddInputFiles.push_back(ZZNonIsoName.str());
  stringstream ZZAllName;
  ZZAllName << ZZAllPrefix << ZZSuffix;
  ZZAllHaddInputFiles.push_back(ZZAllName.str());
  haddCanvases(ZZIsoHaddOutputFile, ZZIsoHaddInputFiles, 
	       vector<float>(1, 0.0377509625152821), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(ZZNonIsoHaddOutputFile, ZZNonIsoHaddInputFiles, 
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
  string WWIsoPrefix(analysisFilePath + "WW/analysis/muHadIsoAnalysis" + MTBin + "_WW");
  string WWIsoHaddOutputFile(WWIsoPrefix + "_hadd" + WWSuffix);
  string WWNonIsoPrefix(analysisFilePath + "WW/analysis/muHadNonIsoAnalysis" + MTBin + "_WW");
  string WWNonIsoHaddOutputFile(WWNonIsoPrefix + "_hadd" + WWSuffix);
  string WWAllPrefix(analysisFilePath + "WW/analysis/muHadAnalysis" + MTBin + "_WW");
  string WWAllHaddOutputFile(WWAllPrefix + "_hadd" + WWSuffix);
  vector<string> WWIsoHaddInputFiles;
  vector<string> WWNonIsoHaddInputFiles;
  vector<string> WWAllHaddInputFiles;
  stringstream WWIsoName;
  WWIsoName << WWIsoPrefix << WWSuffix;
  WWIsoHaddInputFiles.push_back(WWIsoName.str());
  stringstream WWNonIsoName;
  WWNonIsoName << WWNonIsoPrefix << WWSuffix;
  WWNonIsoHaddInputFiles.push_back(WWNonIsoName.str());
  stringstream WWAllName;
  WWAllName << WWAllPrefix << WWSuffix;
  WWAllHaddInputFiles.push_back(WWAllName.str());
  haddCanvases(WWIsoHaddOutputFile, WWIsoHaddInputFiles, 
	       vector<float>(1, 0.124167251024691), canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  haddCanvases(WWNonIsoHaddOutputFile, WWNonIsoHaddInputFiles, 
	       vector<float>(1, 0.124167251024691), canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  if (doNoHPSIsoCut) {
    haddCanvases(WWAllHaddOutputFile, WWAllHaddInputFiles, 
		 vector<float>(1, 0.124167251024691), canvasNames1D, 
		 graphNames1D, canvasNames2D, graphNames2D, nullBlindLow, nullBlindHigh);
  }
  cout << endl;

  //compare data to MC in control region and compute data - MC for data-driven QCD shape
  string dataVsMCOutputFile(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + MTBin + 
			    tag19p7InvFb + outputVTag + fileExt);
  string dataVsMCOutputDiff(analysisFilePath + "results/dataVsMC_muHadNonIsoDifference" + MTBin + 
			    tag19p7InvFb + outputVTag + fileExt);
  vector<string> dataVsMCInputFiles;
  dataVsMCInputFiles.push_back(dataNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(WWNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(ZZNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(WZNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(WNJetsToLNuNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(TNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(TTJetsNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(DYJetsToLLNonIsoHaddOutputFile);
  cout << "\nPlot data vs. MC normalized to data luminosity\n---\n";
  drawMultipleEfficiencyGraphsOn1Canvas(dataVsMCOutputFile, dataVsMCInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders19p7InvFb, 
					colorsMCData, styles, legendEntriesMCData, 
					weightsMCData, setLogY, drawStack, dataMC);
  cout << "\nPlot data minus MC normalized to data luminosity\n---\n";
  drawDifferenceGraphsOn1Canvas(dataVsMCOutputDiff, dataVsMCInputFiles, 
				canvasNames1D, graphNames1D, legendHeaders19p7InvFb, 
				colorsMCData, styles, legendEntriesMCData, weightsMCData, 
				setLogY, sigBkg);

  //conservative estimate of the J/psi-->mumu or upsilon-->mumu/tautau background from region C
  //compute data-driven QCD estimate in signal (i.e. isolated W muon + isolated tau) region: 
  //region C - resonance estimate
  const pair<Int_t, Int_t> normReg(1, 12);
  const string resBkgOutputFileNoRebin(analysisFilePath + "results/resBkg_noRebin" + MTBin + 
				       tag19p7InvFb + outputVTag + fileExt);
  vector<Double_t> nominalBinEdges;
  cout << "\nEstimate boosted resonance background with narrow bins\n---\n";
  estimateBoostedResonanceBackground(resBkgOutputFileNoRebin, dataVsMCOutputDiff, /*first, no 
										    rebinning*/
				     nonIsoWDataIsoHaddOutputFile, 
				     nonIsoWDataNonIsoHaddOutputFile, normReg, HLTPath, MTBin, 
				     nominalBinEdges.size() - 1, nominalBinEdges);
  nominalBinEdges.push_back(0.0);
  nominalBinEdges.push_back(1.0);
  nominalBinEdges.push_back(2.0);
  nominalBinEdges.push_back(3.0);
  nominalBinEdges.push_back(4.0);
  nominalBinEdges.push_back(11.0);
  const string resBkgOutputFile(analysisFilePath + "results/resBkg" + MTBin + tag19p7InvFb + 
				outputVTag + fileExt);
  cout << "\nEstimate boosted resonance background and rebin to match nominal binning\n---\n";
  estimateBoostedResonanceBackground(resBkgOutputFile, dataVsMCOutputDiff, //rebin
				     nonIsoWDataIsoHaddOutputFile, 
				     nonIsoWDataNonIsoHaddOutputFile, normReg, HLTPath, MTBin, 
				     nominalBinEdges.size() - 1, nominalBinEdges);

  //compute data-driven QCD estimate in control (i.e. isolated W muon + non-isolated tau) region
  string outputFileNameB(analysisFilePath + "results/dataVsMC_RegionBQCDEstimate" + MTBin + 
			 dataVTag + fileExt);
  string inputFileNameC(dataVsMCOutputDiff); // Region B
  string inputFileNameD(nonIsoWDataNonIsoHaddOutputFile); // Region D
  const pair<Double_t, Double_t> SFAndErrTotInt = 
    normFactorAndError(pair<string, string>(inputFileNameC, "muHadMass"), 
		       pair<string, string>(inputFileNameD, "muHadMass"), 
		       pair<Int_t, Int_t>(0, -1));
  cout << "\nPlot data-driven QCD estimate for region B\n---\n";
  drawQCDRegionAHistograms(outputFileNameB,inputFileNameD,canvasNames1D, graphNames1D,
			   legendHeaders19p7InvFb,colorsMCData, styles, legendEntriesMCData,
			   weightsMCData, setLogY, sigBkg, SFAndErrTotInt);

  //compare region C data to region D data
  const string variable("muHadMass");
  const string theunit("m_{#mu+X} (GeV)");
  string inputFileNameB(nonIsoWDataIsoHaddOutputFile); // Region C
  vector<string> compRegCDataToRegDDataInputFiles;
  compRegCDataToRegDDataInputFiles.push_back(inputFileNameD);
  compRegCDataToRegDDataInputFiles.push_back(inputFileNameB);
  string compRegCDataToRegDDataOutputFile(analysisFilePath + "results/regCDataVsRegDData" + 
					  MTBin + tag19p7InvFb + outputVTag + fileExt);
  QCDVsMCClosurePlots(compRegCDataToRegDDataInputFiles, variable, theunit, 
		      pair<string, string>("Region D data", "Region C data"), 1, 12, 0.0, 3.75, 
		      compRegCDataToRegDDataOutputFile, "RECREATE");
  QCDVsMCClosurePlots(compRegCDataToRegDDataInputFiles, "muHadMass1Prong", theunit, 
		      pair<string, string>("Region D data", "Region C data"), 1, 12, 0.0, 3.75, 
		      compRegCDataToRegDDataOutputFile, "UPDATE");
  QCDVsMCClosurePlots(compRegCDataToRegDDataInputFiles, "muHadMass1Prong1Pi0", theunit, 
		      pair<string, string>("Region D data", "Region C data"), 1, 12, 0.0, 3.75, 
		      compRegCDataToRegDDataOutputFile, "UPDATE");
  QCDVsMCClosurePlots(compRegCDataToRegDDataInputFiles, "muHadMass1Prong2Pi0", theunit, 
		      pair<string, string>("Region D data", "Region C data"), 1, 12, 0.0, 3.75, 
		      compRegCDataToRegDDataOutputFile, "UPDATE");
  QCDVsMCClosurePlots(compRegCDataToRegDDataInputFiles, "muHadMass3Prong", theunit, 
		      pair<string, string>("Region D data", "Region C data"), 1, 12, 0.0, 3.75, 
		      compRegCDataToRegDDataOutputFile, "UPDATE");
}
