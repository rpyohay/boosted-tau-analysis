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
  				  "QCD Mu-enriched normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbQCDB(canvasNames1D.size(), 
					   "QCD b-enriched normalized to 2.5 fb^{-1}");
  //vector<string> legendHeaders2p5InvFbQCDBMu(canvasNames1D.size(), 
  //				     "QCD bToMu-enriched normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbDYJetsToLL(canvasNames1D.size(), 
						 "Drell-Yan + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbTTJets(canvasNames1D.size(), 
					     "t#bar{t} + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbT(canvasNames1D.size(), 
					"t/#bar{t} normalized to 2.5 fb^{-1}");
  vector<string> 
    legendHeaders2p5InvFbWNJetsToLNu(canvasNames1D.size(), "W + jets normalized to 2.5 fb^{-1}");
  vector<string> 
    legendHeaders2p5InvFbWJetsToLNu(canvasNames1D.size(), "W + jets inclusive normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbWZ(canvasNames1D.size(), "WZ normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbZZ(canvasNames1D.size(), "ZZ normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders1(canvasNames1D.size(), "Normalized to 1");
  vector<string> legendHeaders1QCD(canvasNames1D.size(), 
  			   "QCD Mu-enriched normalized to 1");
  vector<string> legendHeaders1QCDB(canvasNames1D.size(), 
  			    "QCD b-enriched normalized to 1");
  //vector<string> legendHeaders1QCDBMu(canvasNames1D.size(), 
  //				  "QCD bToMu-enriched normalized to 1");
  vector<string> legendHeaders1DYJetsToLL(canvasNames1D.size(), 
					  "Drell-Yan + jets normalized to 1");
  vector<string> legendHeaders1TTJets(canvasNames1D.size(), "t#bar{t} + jets normalized to 1");
  vector<string> legendHeaders1T(canvasNames1D.size(), "t/#bar{t} normalized to 1");
  vector<string> legendHeaders1WNJetsToLNu(canvasNames1D.size(), "W + jets normalized to 1");
  vector<string> legendHeaders1WJetsToLNu(canvasNames1D.size(), "W + jets inclusive normalized to 1");
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
  legendEntriesSigBkg.push_back("QCDMu");
  legendEntriesSigBkg.push_back("QCDB");
  //legendEntriesSigBkg.push_back("QCDBMu");
  legendEntriesSigBkg.push_back("Drell-Yan + jets");
  legendEntriesSigBkg.push_back("t#bar{t} + jets");
  legendEntriesSigBkg.push_back("t/#bar{t}");
  legendEntriesSigBkg.push_back("W + jets");
   //  legendEntriesSigBkg.push_back("W + jets inclusive");
  legendEntriesSigBkg.push_back("WZ");
  legendEntriesSigBkg.push_back("ZZ");
//   legendEntriesSigBkg.push_back("WW");
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

  //weights (sig. figs are probably wrong)
  //first number in parentheses is the PREP weight
  //second number in parentheses is the best available weight
  const float Wh1Weight20InvFb = 0.07208; /*(0.3604 pb(Pythia LO xs) * 20000 pb^-1)/
  100000(no. events processed)*/
  const float Wh1Weight2p5InvFb = 0.00901; //20 fb^-1 weight * (2.5/20)
  vector<float> weights1(15, 0.0);
  vector<float> weightsSigBkg;
  weightsSigBkg.push_back(Wh1Weight20InvFb);
  weightsSigBkg.push_back(1.0); /*QCDXSecWeights already weighted to 20 fb^-1 ==>
  multiply by 1.0 to get overall weight for 20 fb^-1*/
  weightsSigBkg.push_back(1.0); /*QCDBXSecWeights already weighted to 20 fb^-1 ==>
 multiply by 1.0 to get overall weight for 20 fb^-1*/
  weightsSigBkg.push_back(1.0); /*QCDBMuXSecWeights already weighted to 20fb^-1 ==>
 multiply by 1.0 to get overall weight for 20 fb^-1*/
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
  vector<float> weightsMCData;
  weightsMCData.push_back(1.0); //data
  weightsMCData.push_back(1.0); /*QCDRelXSecWeights already weighted to 20 fb^-1 ==>
 multiply by 1.0 to get overall weight for 20 fb^-1*/
  weightsMCData.push_back(1.0); /*QCDBRelXSecWeights already weighted to 20 fb^-1 ==>
 multiply by 1.0 to get overall weight for 20 fb^-1*/
  weightsMCData.push_back(1.0); /*QCDBMuRelXSecWeights already weighted to 20 fb^-1 ==>
 multiply by 1.0 to get overall weight for 20 fb^-1*/
  weightsMCData.push_back(1.0); /*DYJetsToLLRelXSecWeights already weighted to 20 fb^-1 ==> 
 multiply by 1.0 to get overall weight for 20 fb^-1*/
  weightsMCData.push_back(/*0.45025472914009*/3.60203783312072); /*tt+jets weighted to (2.5)20 fb^-1
      (136.3 pb(PREP xs) * 2500 pb^-1)/
      1364783(no. events processed)*/
  weightsMCData.push_back(1.0); /*TRelXSecWeights already weighted to 20 fb^-1 ==> 
   multiply by 1.0 to get overall weight for 20 fb^-1*/
  weightsMCData.push_back(1.0);
  /*WNJetsToLNuRelXSecWeights already weighted to 20 fb^-1 ==> 
   multiply by 1.0 to get overall weight for 20 fb^-1*/
  //   weightsMCData.push_back(1.0); /*W+jets inclusive already weighted to 20 fb^-1
//       20 fb^-1 weight * (2.5/20)*/
  weightsMCData.push_back(/*0.00847169413373063*/0.067773553069845); //WZ weighted to (2.5)20 fb^-1
  weightsMCData.push_back(/*0.00479073128366524*/0.038325850269322); //ZZ weighted to (2.5)20 fb^-1
//   weightsMCData.push_back(/*0.00965756753889549*/0.0157572653584633); //WW weighted to 2.5 fb^-1
  std::reverse(weightsMCData.begin() + 1, weightsMCData.end());
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
  vector<float> QCDRelXSecWeights;
  QCDRelXSecWeights.push_back(/*549.523124*/4396.184992); // Pt 20-30 bin weighted to (2.5)20 fb^-1
  QCDRelXSecWeights.push_back(/*210.846143*/1686.769144); // Pt 30-50 bin weighted to (2.5)20 fb^-1
  QCDRelXSecWeights.push_back(/*42.4948602*/339.9588816); // Pt 50-80 bin weighted to (2.5)20 fb^-1
  QCDRelXSecWeights.push_back(/*10.9453316*/87.5626528); // Pt 80-120 bin weighted to (2.5)20 fb^-1
  QCDRelXSecWeights.push_back(/*2.19477684*/17.55821472); // Pt 120-170 bin weighted to (2.5)20 fb^-1
  QCDRelXSecWeights.push_back(/*0.7495984*/5.9967872); // Pt 170-300 bin weighted to (2.5)20 fb^-1
  QCDRelXSecWeights.push_back(/*0.04845497*/0.38763976); // Pt 300-470 bin weighted to (2.5)20 fb^-1
  QCDRelXSecWeights.push_back(/*0.00779558*/0.06236464); // Pt 470-600 bin weighted to (2.5)20 fb^-1
  QCDRelXSecWeights.push_back(/*0.0016328*/0.0130624); // Pt 600-800 bin weighted to (2.5)20 fb^-1
  QCDRelXSecWeights.push_back(/*0.00022444*/0.00179552); // Pt 800-1000 bin weighted to (2.5)20 fb^-1
  QCDRelXSecWeights.push_back(/*0.000000000032525*/0.0000000002602); // Pt 1000 bin weighted to (2.5)20 fb^-1
  vector<float> QCDBRelXSecWeights;
  QCDBRelXSecWeights.push_back(/*33630.1979938*/269041.5839504); // Pt 15-30 bin weighted to (2.5)20 fb^-1
  QCDBRelXSecWeights.push_back(/*2764.81797576*/22118.54380608); // Pt 30-50 bin weighted to (2.5)20 fb^-1
  QCDBRelXSecWeights.push_back(/*902.157904718*/7217.263237744); // Pt 50-150 bin weighted to (2.5)20 fb^-1
  QCDBRelXSecWeights.push_back(/*44.7015085911*/357.6120687288); // Pt 150 bin weighted to (2.5)20 fb^-1
  vector<float> QCDBMuRelXSecWeights;
  QCDBMuRelXSecWeights.push_back(0.); // Pt 15-30 bin weighted to 2.5 fb^-1
  QCDBMuRelXSecWeights.push_back(0.); // Pt 30-50 bin weighted to 2.5 fb^-1
  QCDBMuRelXSecWeights.push_back(0.); // Pt 50-150 bin weighted to 2.5 fb^-1
  QCDBMuRelXSecWeights.push_back(0.); // Pt 150 bin weighted to 2.5 fb^-1

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
  const string QCDBVTag("_v12");
  const string QCDBMuVTag("_v12");

  //hadd data samples from different eras
  //string dataIsoPrefix(analysisFilePath + "data/analysis/muHadIsoAnalysis_SingleMu");
  //string dataIsoSuffix(dataVTag + fileExt);
  //string dataIsoHaddOutputFile(dataIsoPrefix + dataIsoSuffix); //BLINDED!!!
  string dataNonIsoPrefix(analysisFilePath + "data/analysis/muHadNonIsoAnalysis_SingleMu");
  string dataNonIsoSuffix(dataVTag + fileExt);
  string dataNonIsoHaddOutputFile(dataNonIsoPrefix + dataNonIsoSuffix);
  //vector<string> dataIsoHaddInputFiles; //BLINDED!!!
  vector<string> dataNonIsoHaddInputFiles;
  vector<string> runEras;
  runEras.push_back("_Run2012A");
  runEras.push_back("_Run2012B");
  runEras.push_back("_Run2012C");
  runEras.push_back("_Run2012D");
  for (vector<string>::const_iterator iRunEra = runEras.begin(); iRunEra != runEras.end(); 
       ++iRunEra) {
    //stringstream dataIsoName;
    //dataIsoName << dataIsoPrefix << *iRunEra << dataIsoSuffix; //BLINDED!!!
    //dataIsoHaddInputFiles.push_back(dataIsoName.str());
    stringstream dataNonIsoName;
    dataNonIsoName << dataNonIsoPrefix << *iRunEra << dataNonIsoSuffix;
    dataNonIsoHaddInputFiles.push_back(dataNonIsoName.str());
  }
  //haddCanvases(dataIsoHaddOutputFile, dataIsoHaddInputFiles, vector<float>(4, 1.0), 
  //     canvasNames1D, graphNames1D, canvasNames2D, graphNames2D); //BLINDED!!!
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
  haddCanvases(DYJetsToLLIsoHaddOutputFile, DYJetsToLLIsoHaddInputFiles, DYJetsToLLRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(DYJetsToLLNonIsoHaddOutputFile, DYJetsToLLNonIsoHaddInputFiles, 
	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);


  //hadd QCD Mu-enriched Pt-binned samples
  string QCDIsoPrefix(analysisFilePath + "QCD/analysis/muHadIsoAnalysis_QCD");
  string QCDIsoSuffix(QCDVTag + fileExt);
  string QCDIsoHaddOutputFile(QCDIsoPrefix + QCDIsoSuffix);
  string QCDNonIsoPrefix(analysisFilePath + "QCD/analysis/muHadNonIsoAnalysis_QCD");
  string QCDNonIsoSuffix(QCDVTag + fileExt);
  string QCDNonIsoHaddOutputFile(QCDNonIsoPrefix + QCDNonIsoSuffix);
  vector<string> QCDIsoHaddInputFiles;
  vector<string> QCDNonIsoHaddInputFiles;
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
     QCDIsoName << QCDIsoPrefix << *iPtBin << QCDIsoSuffix;
    QCDIsoHaddInputFiles.push_back(QCDIsoName.str());
    stringstream QCDNonIsoName;
    QCDNonIsoName << QCDNonIsoPrefix << *iPtBin << QCDNonIsoSuffix;
    QCDNonIsoHaddInputFiles.push_back(QCDNonIsoName.str());
    }
  
  haddCanvases(QCDIsoHaddOutputFile, QCDIsoHaddInputFiles, QCDRelXSecWeights, 
         canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(QCDNonIsoHaddOutputFile, QCDNonIsoHaddInputFiles, QCDRelXSecWeights, 
         canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
 
  //hadd QCD b-enriched Pt-binned samples
  string QCDBIsoPrefix(analysisFilePath + "QCDB/analysis/muHadIsoAnalysis_QCDB");
  string QCDBIsoSuffix(QCDBVTag + fileExt);
  string QCDBIsoHaddOutputFile(QCDBIsoPrefix + QCDBIsoSuffix);
  string QCDBNonIsoPrefix(analysisFilePath + "QCDB/analysis/muHadNonIsoAnalysis_QCDB");
  string QCDBNonIsoSuffix(QCDBVTag + fileExt);
  string QCDBNonIsoHaddOutputFile(QCDBNonIsoPrefix + QCDBNonIsoSuffix);
  vector<string> QCDBIsoHaddInputFiles;
  vector<string> QCDBNonIsoHaddInputFiles;
  vector<string> ptBinsB;
  ptBinsB.push_back("_Pt-15To30");
  ptBinsB.push_back("_Pt-30To50");
  ptBinsB.push_back("_Pt-50To150");
  ptBinsB.push_back("_Pt-150");
  for (vector<string>::const_iterator iPtBin = ptBinsB.begin(); iPtBin != ptBinsB.end(); 
       ++iPtBin) {
    stringstream QCDBIsoName;
    QCDBIsoName << QCDBIsoPrefix << *iPtBin << QCDBIsoSuffix;
    QCDBIsoHaddInputFiles.push_back(QCDBIsoName.str());
    stringstream QCDBNonIsoName;
    QCDBNonIsoName << QCDBNonIsoPrefix << *iPtBin << QCDBNonIsoSuffix;
    QCDBNonIsoHaddInputFiles.push_back(QCDBNonIsoName.str());
  }
  haddCanvases(QCDBIsoHaddOutputFile, QCDBIsoHaddInputFiles, QCDBRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(QCDBNonIsoHaddOutputFile, QCDBNonIsoHaddInputFiles, QCDBRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

  //hadd QCD bToMu-enriched Pt-binned samples
  string QCDBMuIsoPrefix(analysisFilePath + "QCDBMu/analysis/muHadIsoAnalysis_QCDBMu");
  string QCDBMuIsoSuffix(QCDBMuVTag + fileExt);
  string QCDBMuIsoHaddOutputFile(QCDBMuIsoPrefix + QCDBMuIsoSuffix);
  string QCDBMuNonIsoPrefix(analysisFilePath + "QCDBMu/analysis/muHadNonIsoAnalysis_QCDBMu");
  string QCDBMuNonIsoSuffix(QCDBMuVTag + fileExt);
  string QCDBMuNonIsoHaddOutputFile(QCDBMuNonIsoPrefix + QCDBMuNonIsoSuffix);
  vector<string> QCDBMuIsoHaddInputFiles;
  vector<string> QCDBMuNonIsoHaddInputFiles;
  vector<string> ptBinsBMu;
  ptBinsBMu.push_back("_pt15to30");
  ptBinsBMu.push_back("_pt30to50");
  ptBinsBMu.push_back("_pt50to150");
  ptBinsBMu.push_back("_pt150");
  for (vector<string>::const_iterator iPtBin = ptBinsBMu.begin(); iPtBin != ptBinsBMu.end(); 
       ++iPtBin) {
    stringstream QCDBMuIsoName;
    QCDBMuIsoName << QCDBMuIsoPrefix << *iPtBin << QCDBMuIsoSuffix;
    QCDBMuIsoHaddInputFiles.push_back(QCDBMuIsoName.str());
    stringstream QCDBMuNonIsoName;
    QCDBMuNonIsoName << QCDBMuNonIsoPrefix << *iPtBin << QCDBMuNonIsoSuffix;
    QCDBMuNonIsoHaddInputFiles.push_back(QCDBMuNonIsoName.str());
  }
  haddCanvases(QCDBMuIsoHaddOutputFile, QCDBMuIsoHaddInputFiles, QCDBMuRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(QCDBMuNonIsoHaddOutputFile, QCDBMuNonIsoHaddInputFiles, QCDBMuRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

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
  string sigVsBkgOutputFile20InvFb(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + 
				   tag20InvFb + outputVTag + fileExt);
  string sigVsBkgOutputFile1(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + tag1 + 
			     outputVTag + fileExt);
  vector<string> sigVsBkgInputFiles;
  sigVsBkgInputFiles.push_back(analysisFilePath + "Wh1_Medium/muHadIsoAnalysis_Wh1" + sigVTag + 
			       fileExt);
  sigVsBkgInputFiles.push_back(QCDIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(QCDBIsoHaddOutputFile);
  //sigVsBkgInputFiles.push_back(QCDBMuIsoHaddOutputFile);
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
					setLogY, drawStack, sigBkg);
  drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile1, sigVsBkgInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders1, colors, 
					styles, legendEntriesSigBkg, weights1, setLinY, drawSame, 
					sigBkg);

  //compare data to MC in control region
  string dataVsMCOutputFile2p5InvFb(analysisFilePath + "results/dataVsMC_muHadNonIsoAnalysis" + 
				    tag2p5InvFb + outputVTag + fileExt);
  string dataVsMCOutputDiff2p5InvFb(analysisFilePath + "results/dataVsMC_muHadNonIsoDifference" + 
				    tag2p5InvFb + outputVTag + fileExt);
  vector<string> dataVsMCInputFiles;
  vector<string> dataVsMCDifferenceInputFiles;
  dataVsMCInputFiles.push_back(dataNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(QCDNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(QCDBNonIsoHaddOutputFile);
  //dataVsMCInputFiles.push_back(QCDBMuNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(DYJetsToLLNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(analysisFilePath + "TTJets/analysis/muHadNonIsoAnalysis_TTJets" + 
			       TTJetsVTag + fileExt);
  dataVsMCInputFiles.push_back(TNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(WNJetsToLNuNonIsoHaddOutputFile);
  //  dataVsMCInputFiles.push_back(analysisFilePath + "WJetsToLNu/analysis/muHadNonIsoAnalysis_WJetsToLNu" + WJetsVTag + 
  //			       fileExt);
  dataVsMCInputFiles.push_back(analysisFilePath + "WZ/analysis/muHadNonIsoAnalysis_WZ" + WZVTag + 
			       fileExt);
  dataVsMCInputFiles.push_back(analysisFilePath + "ZZ/analysis/muHadNonIsoAnalysis_ZZ" + ZZVTag + 
			       fileExt);
//   dataVsMCInputFiles.push_back(analysisFilePath + "WW/analysis/muHadNonIsoAnalysis_WW" + WWVTag + 
// 			       fileExt);
  std::reverse(dataVsMCInputFiles.begin() + 1, dataVsMCInputFiles.end());
  drawMultipleEfficiencyGraphsOn1Canvas(dataVsMCOutputFile2p5InvFb, dataVsMCInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders2p5InvFb, 
					colors, styles, legendEntriesMCData, 
					weightsMCData, setLogY, drawStack, dataMC);

  dataVsMCDifferenceInputFiles.push_back(dataNonIsoHaddOutputFile);
  dataVsMCDifferenceInputFiles.push_back(DYJetsToLLNonIsoHaddOutputFile);
  dataVsMCDifferenceInputFiles.push_back(analysisFilePath + "TTJets/analysis/muHadNonIsoAnalysis_TTJets" + 
			       TTJetsVTag + fileExt);
  dataVsMCDifferenceInputFiles.push_back(TNonIsoHaddOutputFile);
  dataVsMCDifferenceInputFiles.push_back(WNJetsToLNuNonIsoHaddOutputFile);
  //  dataVsMCDifferenceInputFiles.push_back(analysisFilePath + "WJetsToLNu/analysis/muHadNonIsoAnalysis_WJetsToLNu" + WJetsVTag + 
  //			       fileExt);
  dataVsMCDifferenceInputFiles.push_back(analysisFilePath + "WZ/analysis/muHadNonIsoAnalysis_WZ" + WZVTag + 
			       fileExt);
  dataVsMCDifferenceInputFiles.push_back(analysisFilePath + "ZZ/analysis/muHadNonIsoAnalysis_ZZ" + ZZVTag + 
			       fileExt);
  std::reverse(dataVsMCDifferenceInputFiles.begin() + 1, dataVsMCDifferenceInputFiles.end());

  
  drawDifferenceGraphsOn1Canvas(dataVsMCOutputDiff2p5InvFb,dataVsMCInputFiles,
				canvasNames1D, graphNames1D, legendHeaders2p5InvFb,
				colors, styles, legendEntriesMCData,
				weightsMCData, setLogY, dataMC);


  string outputFileNameA = "/data1/friccita/results/dataVsMC_RegionAQCDEstimate";
  outputFileNameA+=dataVTag;
  outputFileNameA+=".root";
  string inputFileNameB = "/data1/friccita/results/dataVsMC_muHadNonIsoDifference_2p5fb-1_v12.root";
  string inputFileNameC = "/data1/friccita/data/analysis/muHadIsoAnalysis_SingleMu_v9.root";
  string inputFileNameD = "/data1/friccita/data/analysis/muHadNonIsoAnalysis_SingleMu_v9.root";
  drawQCDRegionAHistograms(outputFileNameA,inputFileNameB,inputFileNameC,
			   inputFileNameD,canvasNames1D, graphNames1D,
			   legendHeaders2p5InvFb,colors, styles, legendEntriesMCData,
			   weightsMCData, setLogY, dataMC);
  
  //compare QCD search sample to control sample
  string QCDSearchVsControlOutputFile(analysisFilePath + 
					     "QCD/analysis/isoVsNonIsoTaus" + tag1 + 
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
  string QCDBSearchVsControlOutputFile(analysisFilePath + 
					     "QCDB/analysis/isoVsNonIsoTaus" + tag1 + 
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
  string QCDBMuSearchVsControlOutputFile(analysisFilePath + 
					     "QCDBMu/analysis/isoVsNonIsoTaus" + tag1 + 
					     outputVTag + fileExt);
  /* vector<string> QCDBMuSearchVsControlInputFiles;
  QCDBMuSearchVsControlInputFiles.push_back(QCDBMuIsoHaddOutputFile);
  QCDBMuSearchVsControlInputFiles.push_back(QCDBMuNonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(QCDBMuSearchVsControlOutputFile, 
					QCDBMuSearchVsControlInputFiles, canvasNames1D, 
					graphNames1D, legendHeaders1QCDBMu, colors, styles, 
					legendEntriesSearchVsControl, weights1, setLinY, drawSame, 
					sigBkg);
  */
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
//   drawMultipleEfficiencyGraphsOn1Canvas(WWSearchVsControlOutputFile, WWSearchVsControlInputFiles, 
// 					canvasNames1D, graphNames1D, legendHeaders1WW, colors, 
// 					styles, legendEntriesSearchVsControl, weights1, setLinY, 
// 					drawSame, sigBkg);

//   //unit strings
//   string unitPTTau("Reco #tau p_{T} (GeV)");
//   string unitPTMu("Reco #mu p_{T} (GeV)");
//   string unitEtaTau("Reco #tau #eta");
//   string unitEtaMu("Reco #mu #eta");
//   string unitDR("#DeltaR(visible gen #tau, gen #mu)");
//   string noUnit("");

//   //map of bin labels for certain efficiency plots
//   vector<string> binLabels;
//   binLabels.push_back("kNull");
//   binLabels.push_back("kOneProng0PiZero");
//   binLabels.push_back("kOneProng1PiZero");
//   binLabels.push_back("kOneProng2PiZero");
//   binLabels.push_back("kOneProng3PiZero");
//   binLabels.push_back("kOneProngNPiZero");
//   binLabels.push_back("kTwoProng0PiZero");
//   binLabels.push_back("kTwoProng1PiZero");
//   binLabels.push_back("kTwoProng2PiZero");
//   binLabels.push_back("kTwoProng3PiZero");
//   binLabels.push_back("kTwoProngNPiZero");
//   binLabels.push_back("kThreeProng0PiZero");
//   binLabels.push_back("kThreeProng1PiZero");
//   binLabels.push_back("kThreeProng2PiZero");
//   binLabels.push_back("kThreeProng3PiZero");
//   binLabels.push_back("kThreeProngNPiZero");
//   binLabels.push_back("kRareDecayMode");
//   map<string, vector<string> > binLabelMap;
//   binLabelMap["muHadGenDecayMode"] = binLabels;
//   binLabelMap["muHadCorrectRecoDecayModeGenDecayMode"] = binLabels;
//   binLabelMap["muHadRecoDecayMode"] = binLabels;
//   binLabelMap["muHadGen1ProngRecoDecayMode"] = binLabels;
//   binLabelMap["muHadGen1Prong1Pi0RecoDecayMode"] = binLabels;
//   binLabelMap["muHadGen3ProngRecoDecayMode"] = binLabels;

//   //map of inputs to efficiency histograms
//   map<string, pair<string, string> > effHistMap;
//   effHistMap["numeratorPT"] = make_pair(string("denominatorPT"), unitPTMu);
//   map<string, pair<string, string> > effHistMapMu;
//   effHistMapMu["numeratorPT"] = make_pair(string("denominatorPT"), unitPTMu);
//   effHistMapMu["numeratorEta"] = make_pair(string("denominatorEta"), unitEtaMu);
//   map<string, pair<string, string> > effHistMapTau;
//   effHistMapTau["numeratorPT"] = make_pair(string("denominatorPT"), unitPTTau);
//   effHistMapTau["numeratorEta"] = make_pair(string("denominatorEta"), unitEtaTau);

//   //map of inputs to 1D histograms
//   map<string, string> hist1DMap;
//   hist1DMap["numeratorPT"] = unitPTMu;
//   hist1DMap["denominatorPT"] = unitPTMu;
//   map<string, string> hist1DMapMu;
//   hist1DMapMu["numeratorPT"] = unitPTMu;
//   hist1DMapMu["denominatorPT"] = unitPTMu;
//   hist1DMapMu["numeratorEta"] = unitEtaMu;
//   hist1DMapMu["denominatorEta"] = unitEtaMu;
//   map<string, string> hist1DMapTau;
//   hist1DMapTau["numeratorPT"] = unitPTTau;
//   hist1DMapTau["denominatorPT"] = unitPTTau;
//   hist1DMapTau["numeratorEta"] = unitEtaTau;
//   hist1DMapTau["denominatorEta"] = unitEtaTau;

//   vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > > histMap;
//   histMap.push_back(make_pair(make_pair(), make_pair()));

//   map<pair<string, string>, vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > > > 
//     canvasMap;
//   canvasMap[make_pair(string("muHadGen1ProngRecoDecayMode"), noUnit)] = ;
//   canvasMap[make_pair(string("muHadGen1Prong1Pi0RecoDecayMode"), noUnit)] = ;
//   canvasMap[make_pair(string("muHadGen3ProngRecoDecayMode"), noUnit)] = ;

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
