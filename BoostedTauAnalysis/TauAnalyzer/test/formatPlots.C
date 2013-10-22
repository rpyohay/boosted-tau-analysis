{
  //load
  gROOT->Reset();
  gROOT->ProcessLine("#include <utility>");
  string macroPath("/afs/cern.ch/user/y/yohay/CMSSW_5_3_3_Git/src/BoostedTauAnalysis/");
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
  canvasNames1D.push_back("jet_pt_etacutCanvas");
  canvasNames1D.push_back("jet_etaCanvas");
  canvasNames1D.push_back("jet_phiCanvas");
  canvasNames1D.push_back("jet_mass_etacutCanvas");
  canvasNames1D.push_back("jet_ptmj_etacutCanvas");
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
  graphNames1D.push_back("jet_pt_etacut");
  graphNames1D.push_back("jet_eta");
  graphNames1D.push_back("jet_phi");
  graphNames1D.push_back("jet_mass_etacut");
  graphNames1D.push_back("jet_ptmj_etacut");
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

  //set up plot style options
  vector<string> legendHeaders20InvFb(canvasNames1D.size(), "Normalized to 20 fb^{-1}");
  vector<string> legendHeaders2p5InvFb(canvasNames1D.size(), "Normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbDYJetsToLL(canvasNames1D.size(), 
						 "Drell-Yan + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbTTJets(canvasNames1D.size(), 
					     "t#bar{t} + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbT(canvasNames1D.size(), 
					"t/#bar{t} normalized to 2.5 fb^{-1}");
  vector<string> 
    legendHeaders2p5InvFbWNJetsToLNu(canvasNames1D.size(), "W + jets normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbWZ(canvasNames1D.size(), "WZ normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbZZ(canvasNames1D.size(), "ZZ normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders2p5InvFbWW(canvasNames1D.size(), "WW normalized to 2.5 fb^{-1}");
  vector<string> legendHeaders1(canvasNames1D.size(), "Normalized to 1");
  vector<string> legendHeaders1DYJetsToLL(canvasNames1D.size(), 
					  "Drell-Yan + jets normalized to 1");
  vector<string> legendHeaders1TTJets(canvasNames1D.size(), "t#bar{t} + jets normalized to 1");
  vector<string> legendHeaders1T(canvasNames1D.size(), "t/#bar{t} normalized to 1");
  vector<string> legendHeaders1WNJetsToLNu(canvasNames1D.size(), "W + jets normalized to 1");
  vector<string> legendHeaders1WZ(canvasNames1D.size(), "WZ normalized to 1");
  vector<string> legendHeaders1ZZ(canvasNames1D.size(), "ZZ normalized to 1");
  vector<string> legendHeaders1WW(canvasNames1D.size(), "WW normalized to 1");
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
  vector<string> legendEntriesSigBkgInd;
  legendEntriesSigBkgInd.push_back("Wh_{1}");
  legendEntriesSigBkgInd.push_back("Drell-Yan + jets (10 < m_{l^{+}l^{-}} < 50) GeV");
  legendEntriesSigBkgInd.push_back("Drell-Yan + jets m_{l^{+}l^{-}} > 50 GeV");
  legendEntriesSigBkgInd.push_back("t#bar{t} + jets");
  legendEntriesSigBkgInd.push_back("t s-channel");
  legendEntriesSigBkgInd.push_back("#bar{t} s-channel");
  legendEntriesSigBkgInd.push_back("t t-channel");
  legendEntriesSigBkgInd.push_back("#bar{t} t-channel");
  legendEntriesSigBkgInd.push_back("W + 1 jet");
  legendEntriesSigBkgInd.push_back("W + 2 jets");
  legendEntriesSigBkgInd.push_back("W + 3 jets");
  legendEntriesSigBkgInd.push_back("W + 4 jets");
  legendEntriesSigBkgInd.push_back("WZ");
  legendEntriesSigBkgInd.push_back("ZZ");
  legendEntriesSigBkgInd.push_back("WW");
  std::reverse(legendEntriesSigBkgInd.begin() + 1, legendEntriesSigBkgInd.end());
  vector<string> legendEntriesMCDataInd(legendEntriesSigBkgInd);
  legendEntriesMCDataInd[0] = "Data 2.5 fb^{-1}";
  vector<string> legendEntriesSigBkg;
  legendEntriesSigBkg.push_back("Wh_{1}");
  legendEntriesSigBkg.push_back("Drell-Yan + jets");
  legendEntriesSigBkg.push_back("t#bar{t} + jets");
  legendEntriesSigBkg.push_back("t/#bar{t}");
  legendEntriesSigBkg.push_back("W + jets");
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

  //weights (sig. figs are probably wrong)
  //first number in parentheses is the PREP weight
  //second number in parentheses is the best available weight
  const float Wh1Weight20InvFb = 0.07208; /*(0.3604 pb(Pythia LO xs) * 20000 pb^-1)/
					   100000(no. events processed)*/
  const float Wh1Weight2p5InvFb = 0.00901; //20 fb^-1 weight * (2.5/20)
  vector<float> weights1(15, 0.0);
  vector<float> weightsSigBkgInd;
  weightsSigBkgInd.push_back(Wh1Weight20InvFb);
  weightsSigBkgInd.push_back(/*5.841109917*/1.18489425542698); /*Drell-Yan + jets 
								 (10 < ml+l- < 50) GeV weighted to 
								 20 fb^-1*/
  weightsSigBkgInd.push_back(/*1.936998118*/2.30056938223844); /*Drell-Yan + jets 
								 ml+l- > 50 GeV weighted to 
								 20 fb^-1*/
  weightsSigBkgInd.push_back(/*1.99738713*/3.60203783312072); //tt+jets weighted to 20 fb^-1
  weightsSigBkgInd.push_back(/*0.2169556203*/0.291582198868292); /*t s-channel weighted to 20 fb^-1
								   (2.82 pb(PREP xs) * 
								   20000 pb^-1)/
								   259961(no. events processed)*/
  weightsSigBkgInd.push_back(/*0.314*/0.352); /*tbar s-channel weighted to 20 fb^-1
						(1.57 pb(PREP xs) * 20000 pb^-1)/
						100000(no. events processed)*/
  weightsSigBkgInd.push_back(/*0.2718155864*/0.326178703711468); /*t t-channel weighted to 20 fb^-1
								   (47 pb(PREP xs) * 20000 pb^-1)/
								   3458227(no. events processed)*/
  weightsSigBkgInd.push_back(/*0.2583883184*/0.317300854955268); /*tbar t-channel weighted to 20 
								   fb^-1
								   (25 pb(PREP xs) * 20000 pb^-1)/
								   1935072(no. events processed)*/
  weightsSigBkgInd.push_back(/*4.666920582*/3.38697439995285); /*W + 1 jet weighted to 20 fb^-1
								 (5400 pb(PREP xs) * 20000 pb^-1)/
								 23141598(no. events processed)*/
  weightsSigBkgInd.push_back(/*1.034128577*/0.461044066257386); /*W + 2 jets weighted to 20 fb^-1
								  (1750 pb(PREP xs) * 20000 pb^-1)/
								  33844921(no. events processed)*/
  weightsSigBkgInd.push_back(/*0.6811245747*/0.164178582464271); /*W + 3 jets weighted to 20 fb^-1
								   (519 pb(PREP xs) * 20000 pb^-1)/
								   15239503(no. events processed)*/
  weightsSigBkgInd.push_back(/*0.3198134203*/0.0314134490360502); /*W + 4 jets weighted to 20 fb^-1
								    (214 pb(PREP xs) * 
								    20000 pb^-1)/
								    13382803(no. events 
								    processed)*/
  weightsSigBkgInd.push_back(/*0.0257747659*/0.067773553069845); /*WZ weighted to 20 fb^-1
								   (12.63 pb(PREP xs) * 
								   20000 pb^-1)/
								   9800283(no. events processed)*/
  weightsSigBkgInd.push_back(/*0.0112802265*/0.0383258502693219); /*ZZ weighted to 20 fb^-1
								    (5.196 pb(PREP xs) * 
								    20000 pb^-1)/
								    9212581(no. events processed)*/
  weightsSigBkgInd.push_back(/*0.0772605403111639*/0.126058122867706); /*WW weighted to 20 fb^-1
									 (33.61 pb(PREP xs) * 
									 20000 pb^-1)/
									 (no. events processed)*/
  std::reverse(weightsSigBkgInd.begin() + 1, weightsSigBkgInd.end());
  vector<float> weightsMCDataInd;
  weightsMCDataInd.push_back(1.0); //data
  weightsMCDataInd.push_back(/*0.7301387396*/0.148111781928372); /*Drell-Yan + jets 
								   (10 < ml+l- < 50) GeV weighted 
								   to 2.5 fb^-1
								   11050 pb(PREP xs) * 2500 pb^-1/
								   37835275(no. events processed)*/
  weightsMCDataInd.push_back(/*0.2421247648*/0.287571172779805); /*Drell-Yan + jets 
								   ml+l- > 50 GeV weighted 
								   to 2.5 fb^-1
								   2950 pb(PREP xs) * 2500 pb^-1/
								   30459503(no. events processed)*/
  weightsMCDataInd.push_back(/*0.2496733913*/0.45025472914009); //tt+jets weighted to 2.5 fb^-1
  weightsMCDataInd.push_back(/*0.0271194525*/0.0364477748585365); /*t s-channel weighted to 
								    2.5 fb^-1*/
  weightsMCDataInd.push_back(/*0.03925*/0.044); //tbar s-channel weighted to 2.5 fb^-1
  weightsMCDataInd.push_back(/*0.0339769483*/0.0407723379639335); /*t t-channel weighted to 
								    2.5 fb^-1*/
  weightsMCDataInd.push_back(/*0.0322985398*/0.0396626068694085); /*tbar t-channel weighted to 
								    2.5 fb^-1*/
  weightsMCDataInd.push_back(/*0.5833650728*/0.423371799994106); //W + 1 jet weighted to 2.5 fb^-1
  weightsMCDataInd.push_back(/*0.1292660721*/0.0576305082821733); /*W + 2 jets weighted to 
								    2.5 fb^-1*/
  weightsMCDataInd.push_back(/*0.0851405718*/0.0205223228080338); /*W + 3 jets weighted to 
								    2.5 fb^-1*/
  weightsMCDataInd.push_back(/*0.0399766775*/0.00392668112950628); /*W + 4 jets weighted to 
								     2.5 fb^-1*/
  weightsMCDataInd.push_back(/*0.0032218457*/0.00847169413373063); //WZ weighted to 2.5 fb^-1
  weightsMCDataInd.push_back(/*0.0014100283*/0.00479073128366524); //ZZ weighted to 2.5 fb^-1
  weightsMCDataInd.push_back(/*0.00965756753889549*/0.0157572653584633); //WW weighted to 2.5 fb^-1
  std::reverse(weightsMCDataInd.begin() + 1, weightsMCDataInd.end());
  vector<float> weightsSigBkg;
  weightsSigBkg.push_back(Wh1Weight20InvFb);
  weightsSigBkg.push_back(1.0); /*DYJetsToLLRelXSecWeights already weighted to 20 fb^-1 ==> 
				  multiply by 1 to get overall weight for 20 fb^-1*/
  weightsSigBkg.push_back(/*1.99738713*/3.60203783312072); //tt+jets weighted to 20 fb^-1
  weightsSigBkg.push_back(1.0); //t already weighted to 20 fb^-1
  weightsSigBkg.push_back(1.0); //W + jets already weighted to 20 fb^-1
//   weightsSigBkg.push_back(33.055892185598); /*W + jets weighted to 20 fb^-1
// 					      (30400 pb(PREP xs) * 20000 pb^-1)/
// 					      18393090(no. events processed)*/
  weightsSigBkg.push_back(/*0.0257747659*/0.067773553069845); //WZ weighted to 20 fb^-1
  weightsSigBkg.push_back(/*0.0112802265*/0.0383258502693219); //ZZ weighted to 20 fb^-1
  weightsSigBkg.push_back(/*0.0772605403111639*/0.126058122867706); //WW weighted to 20 fb^-1
  std::reverse(weightsSigBkg.begin() + 1, weightsSigBkg.end());
  vector<float> weightsMCData;
  weightsMCData.push_back(1.0); //data
  weightsMCData.push_back(0.125); /*DYJetsToLLRelXSecWeights already weighted to 20 fb^-1 ==> 
				    multiply by 2.5/20 to get overall weight for 2.5 fb^-1*/
  weightsMCData.push_back(/*0.2496733913*/0.45025472914009); /*tt+jets weighted to 2.5 fb^-1
							       (136.3 pb(PREP xs) * 2500 pb^-1)/
							       1364783(no. events processed)*/
  weightsMCData.push_back(0.125); /*TRelXSecWeights already weighted to 20 fb^-1 ==> 
				    multiply by 2.5/20 to get overall weight for 2.5 fb^-1*/
  weightsMCData.push_back(0.125); /*WNJetsToLNuRelXSecWeights already weighted to 20 fb^-1 ==> 
				    multiply by 2.5/20 to get overall weight for 2.5 fb^-1*/
//   weightsMCData.push_back(4.13198652319975); /*W+jets weighted to 2.5 fb^-1
// 					       20 fb^-1 weight * (2.5/20)*/
  weightsMCData.push_back(/*0.0032218457*/0.00847169413373063); //WZ weighted to 2.5 fb^-1
  weightsMCData.push_back(/*0.0014100283*/0.00479073128366524); //ZZ weighted to 2.5 fb^-1
  weightsMCData.push_back(/*0.00965756753889549*/0.0157572653584633); //WW weighted to 2.5 fb^-1
  std::reverse(weightsMCData.begin() + 1, weightsMCData.end());
  vector<float> WNJetsToLNuRelXSecWeights;
  WNJetsToLNuRelXSecWeights.push_back(/*4.666920582*/3.38697439995285); /*W + 1 jet weighted to 
									  20 fb^-1*/
  WNJetsToLNuRelXSecWeights.push_back(/*1.034128577*/0.461044066257386); /*W + 2 jets weighted to 
									   20 fb^-1*/
  WNJetsToLNuRelXSecWeights.push_back(/*0.6811245747*/0.164178582464271); /*W + 3 jets weighted to 
									    20 fb^-1*/
  WNJetsToLNuRelXSecWeights.push_back(/*0.3198134203*/0.0314134490360502); /*W + 4 jets weighted 
									     to 20 fb^-1*/
  vector<float> DYJetsToLLRelXSecWeights;
  DYJetsToLLRelXSecWeights.push_back(/*5.841109917*//*1.18489425542698*/9.24494976069056); /*(10 < m < 50) GeV 
									 weighted to 20 fb^-1*/
  DYJetsToLLRelXSecWeights.push_back(/*1.936998118*/2.30056938223844); /*m > 50 GeV 
									 weighted to 20 fb^-1*/
  vector<float> TRelXSecWeights;
  TRelXSecWeights.push_back(0.2169556203); //t s-channel weighted to 20 fb^-1
  TRelXSecWeights.push_back(0.314); //tbar s-channel weighted to 20 fb^-1
  TRelXSecWeights.push_back(0.2718155864); //t t-channel weighted to 20 fb^-1
  TRelXSecWeights.push_back(0.2583883184); //tbar t-channel weighted to 20 fb^-1

  //space-saving constant definitions
  const string analysisFilePath("/data1/yohay/");
  const string fileExt(".root");
  const string tag20InvFb("_20fb-1");
  const string tag2p5InvFb("_2p5fb-1");
  const string tag1("_normalizedTo1");

  //version tags
  const string outputVTag("_v57");
  const string dataVTag("_v57");
  const string sigVTag("_v57");
  const string WNJetsToLNuVTag("_v57");
  const string TTJetsVTag("_v57");
  const string TVTag("_v57");
  const string DYJetsToLLVTag("_v57");
  const string WZVTag("_v57");
  const string ZZVTag("_v57");
  const string WWVTag("_v57");

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
  haddCanvases(DYJetsToLLIsoHaddOutputFile, DYJetsToLLIsoHaddInputFiles, DYJetsToLLRelXSecWeights, 
	       canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(DYJetsToLLNonIsoHaddOutputFile, DYJetsToLLNonIsoHaddInputFiles, 
	       DYJetsToLLRelXSecWeights, canvasNames1D, graphNames1D, canvasNames2D, graphNames2D);

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

  //"hadd" ttbar sample just to get the formatting of the 2D plots the same
  string TTJetsIsoPrefix(analysisFilePath + "TTJets/analysis/muHadIsoAnalysis_TTJets");
  string TTJetsNonIsoPrefix(analysisFilePath + "TTJets/analysis/muHadNonIsoAnalysis_TTJets");
  string TTJetsSuffix(sigVTag + fileExt);
  string TTJetsIsoHaddOutputFile(TTJetsIsoPrefix + "_hadd" + TTJetsSuffix);
  string TTJetsNonIsoHaddOutputFile(TTJetsNonIsoPrefix + "_hadd" + TTJetsSuffix);
  vector<string> TTJetsIsoHaddInputFiles;
  vector<string> TTJetsNonIsoHaddInputFiles;
  stringstream TTJetsIsoName;
  TTJetsIsoName << TTJetsIsoPrefix << TTJetsSuffix;
  TTJetsIsoHaddInputFiles.push_back(TTJetsIsoName.str());
  stringstream TTJetsNonIsoName;
  TTJetsNonIsoName << TTJetsNonIsoPrefix << TTJetsSuffix;
  TTJetsNonIsoHaddInputFiles.push_back(TTJetsNonIsoName.str());
  haddCanvases(TTJetsIsoHaddOutputFile, TTJetsIsoHaddInputFiles, 
	       vector<float>(1, /*1.99738713*/3.60203783312072), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(TTJetsNonIsoHaddOutputFile, TTJetsNonIsoHaddInputFiles, 
	       vector<float>(1, /*1.99738713*/3.60203783312072), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);

  //"hadd" WZ sample just to get the formatting of the 2D plots the same
  string WZIsoPrefix(analysisFilePath + "WZ/analysis/muHadIsoAnalysis_WZ");
  string WZNonIsoPrefix(analysisFilePath + "WZ/analysis/muHadNonIsoAnalysis_WZ");
  string WZSuffix(sigVTag + fileExt);
  string WZIsoHaddOutputFile(WZIsoPrefix + "_hadd" + WZSuffix);
  string WZNonIsoHaddOutputFile(WZNonIsoPrefix + "_hadd" + WZSuffix);
  vector<string> WZIsoHaddInputFiles;
  vector<string> WZNonIsoHaddInputFiles;
  stringstream WZIsoName;
  WZIsoName << WZIsoPrefix << WZSuffix;
  WZIsoHaddInputFiles.push_back(WZIsoName.str());
  stringstream WZNonIsoName;
  WZNonIsoName << WZNonIsoPrefix << WZSuffix;
  WZNonIsoHaddInputFiles.push_back(WZNonIsoName.str());
  haddCanvases(WZIsoHaddOutputFile, WZIsoHaddInputFiles, 
	       vector<float>(1, /*0.0257747659*/0.067773553069845), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(WZNonIsoHaddOutputFile, WZNonIsoHaddInputFiles, 
	       vector<float>(1, /*0.0257747659*/0.067773553069845), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);

  //"hadd" ZZ sample just to get the formatting of the 2D plots the same
  string ZZIsoPrefix(analysisFilePath + "ZZ/analysis/muHadIsoAnalysis_ZZ");
  string ZZNonIsoPrefix(analysisFilePath + "ZZ/analysis/muHadNonIsoAnalysis_ZZ");
  string ZZSuffix(sigVTag + fileExt);
  string ZZIsoHaddOutputFile(ZZIsoPrefix + "_hadd" + ZZSuffix);
  string ZZNonIsoHaddOutputFile(ZZNonIsoPrefix + "_hadd" + ZZSuffix);
  vector<string> ZZIsoHaddInputFiles;
  vector<string> ZZNonIsoHaddInputFiles;
  stringstream ZZIsoName;
  ZZIsoName << ZZIsoPrefix << ZZSuffix;
  ZZIsoHaddInputFiles.push_back(ZZIsoName.str());
  stringstream ZZNonIsoName;
  ZZNonIsoName << ZZNonIsoPrefix << ZZSuffix;
  ZZNonIsoHaddInputFiles.push_back(ZZNonIsoName.str());
  haddCanvases(ZZIsoHaddOutputFile, ZZIsoHaddInputFiles, 
	       vector<float>(1, /*0.0112802265*/0.0383258502693219), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);
  haddCanvases(ZZNonIsoHaddOutputFile, ZZNonIsoHaddInputFiles, 
	       vector<float>(1, /*0.0112802265*/0.0383258502693219), canvasNames1D, graphNames1D, 
	       canvasNames2D, graphNames2D);

  //"hadd" WW sample just to get the formatting of the 2D plots the same
  string WWIsoPrefix(analysisFilePath + "WW/analysis/muHadIsoAnalysis_WW");
  string WWNonIsoPrefix(analysisFilePath + "WW/analysis/muHadNonIsoAnalysis_WW");
  string WWSuffix(sigVTag + fileExt);
  string WWIsoHaddOutputFile(WWIsoPrefix + "_hadd" + WWSuffix);
  string WWNonIsoHaddOutputFile(WWNonIsoPrefix + "_hadd" + WWSuffix);
  vector<string> WWIsoHaddInputFiles;
  vector<string> WWNonIsoHaddInputFiles;
  stringstream WWIsoName;
  WWIsoName << WWIsoPrefix << WWSuffix;
  WWIsoHaddInputFiles.push_back(WWIsoName.str());
  stringstream WWNonIsoName;
  WWNonIsoName << WWNonIsoPrefix << WWSuffix;
  WWNonIsoHaddInputFiles.push_back(WWNonIsoName.str());
  haddCanvases(WWIsoHaddOutputFile, WWIsoHaddInputFiles, 
	       vector<float>(1, /*0.0772605403111639*/0.126058122867706), canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D);
  haddCanvases(WWNonIsoHaddOutputFile, WWNonIsoHaddInputFiles, 
	       vector<float>(1, /*0.0772605403111639*/0.126058122867706), canvasNames1D, 
	       graphNames1D, canvasNames2D, graphNames2D);

  //compare MC signal to background
  string sigVsBkgOutputFile20InvFb(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + 
				   tag20InvFb + outputVTag + fileExt);
  string sigVsBkgOutputFile1(analysisFilePath + "results/sigVsBkg_muHadIsoAnalysis" + tag1 + 
			     outputVTag + fileExt);
  vector<string> sigVsBkgInputFiles;
  sigVsBkgInputFiles.push_back(analysisFilePath + "Wh1_Medium/muHadIsoAnalysis_Wh1" + sigVTag + 
			       fileExt);
  sigVsBkgInputFiles.push_back(DYJetsToLLIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(analysisFilePath + "TTJets/analysis/muHadIsoAnalysis_TTJets" + 
			       TTJetsVTag + fileExt);
  sigVsBkgInputFiles.push_back(TIsoHaddOutputFile);
  sigVsBkgInputFiles.push_back(WNJetsToLNuIsoHaddOutputFile);
//   sigVsBkgInputFiles.push_back(analysisFilePath + 
// 			       "WJetsToLNu/analysis/muHadIsoAnalysis_WJetsToLNu" + 
// 			       WNJetsToLNuVTag + fileExt);
  sigVsBkgInputFiles.push_back(analysisFilePath + "WZ/analysis/muHadIsoAnalysis_WZ" + WZVTag + 
			       fileExt);
  sigVsBkgInputFiles.push_back(analysisFilePath + "ZZ/analysis/muHadIsoAnalysis_ZZ" + ZZVTag + 
			       fileExt);
  sigVsBkgInputFiles.push_back(analysisFilePath + "WW/analysis/muHadIsoAnalysis_WW" + WWVTag + 
			       fileExt);
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
  vector<string> dataVsMCInputFiles;
  dataVsMCInputFiles.push_back(dataNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(DYJetsToLLNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(analysisFilePath + "TTJets/analysis/muHadNonIsoAnalysis_TTJets" + 
			       TTJetsVTag + fileExt);
  dataVsMCInputFiles.push_back(TNonIsoHaddOutputFile);
  dataVsMCInputFiles.push_back(WNJetsToLNuNonIsoHaddOutputFile);
//   dataVsMCInputFiles.push_back(analysisFilePath + 
// 			       "WJetsToLNu/analysis/muHadNonIsoAnalysis_WJetsToLNu" + 
// 			       WNJetsToLNuVTag + fileExt);
  dataVsMCInputFiles.push_back(analysisFilePath + "WZ/analysis/muHadNonIsoAnalysis_WZ" + WZVTag + 
			       fileExt);
  dataVsMCInputFiles.push_back(analysisFilePath + "ZZ/analysis/muHadNonIsoAnalysis_ZZ" + ZZVTag + 
			       fileExt);
  dataVsMCInputFiles.push_back(analysisFilePath + "WW/analysis/muHadNonIsoAnalysis_WW" + WWVTag + 
			       fileExt);
  std::reverse(dataVsMCInputFiles.begin() + 1, dataVsMCInputFiles.end());
  drawMultipleEfficiencyGraphsOn1Canvas(dataVsMCOutputFile2p5InvFb, dataVsMCInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders2p5InvFb, 
					colors, styles, legendEntriesMCData, 
					weightsMCData, setLogY, drawStack, dataMC);

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
  drawMultipleEfficiencyGraphsOn1Canvas(WWSearchVsControlOutputFile, WWSearchVsControlInputFiles, 
					canvasNames1D, graphNames1D, legendHeaders1WW, colors, 
					styles, legendEntriesSearchVsControl, weights1, setLinY, 
					drawSame, sigBkg);

//   //MET > 0 GeV MC closure with mu+had mass
//   vector<string> vars;
//   vars.push_back("muHadMass");
//   vars.push_back("muHadPTOverMuHadMass");
//   vector<string> units;
//   units.push_back("m_{#mu+had} (GeV)");
//   units.push_back("p_{T}^{#mu+had}/m^{#mu+had}");
//   vector<int> normRegionLowerBins;
//   normRegionLowerBins.push_back(1);
//   normRegionLowerBins.push_back(1);
//   vector<int> normRegionUpperBins;
//   normRegionUpperBins.push_back(1);
//   normRegionUpperBins.push_back(7);
//   makeMCClosurePlots("/data1/yohay/results/sigVsBkg_muHadIsoAnalysis_20fb-1_v46.root", vars, 
// 		     units, "/data1/yohay/results/dataVsMC_muHadNonIsoAnalysis_2p5fb-1_v46.root", 
// 		     20.0/2.5, normRegionLowerBins, normRegionUpperBins, 
// 		     "/data1/yohay/results/MC_closure_METGeq0_v46.root");

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
