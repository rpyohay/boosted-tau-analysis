void formatEfficiencyPlots(const bool compile)
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

  //unit strings
  string unitPT("p_{T} (GeV)");
  string unitPTTau("Reco #tau p_{T} (GeV)");
  string unitPTMu("Reco #mu p_{T} (GeV)");
  string unitEta("#eta");
  string unitAbsEta("|#eta|");
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

  //map of inputs to 1D efficiency histograms
  map<string, pair<string, string> > effHistMap1D;
  effHistMap1D["numeratorPT"] = make_pair(string("denominatorPT"), unitPT);
  effHistMap1D["numeratorEta"] = make_pair(string("denominatorEta"), unitEta);
  map<string, pair<string, string> > effHistMap1DMu;
  effHistMap1DMu["numeratorPT"] = make_pair(string("denominatorPT"), unitPTMu);
  effHistMap1DMu["numeratorEta"] = make_pair(string("denominatorEta"), unitEtaMu);
  map<string, pair<string, string> > effHistMap1DTau;
  effHistMap1DTau["numeratorPT"] = make_pair(string("denominatorPT"), unitPTTau);
  effHistMap1DTau["numeratorEta"] = make_pair(string("denominatorEta"), unitEtaTau);

  //map of inputs to 2D efficiency histograms
  map<pair<string, string>, pair<string, string> > effHistMap2D;
  effHistMap2D[pair<string, string>("numeratorPTAbsEta", "denominatorPTAbsEta")] = 
    make_pair(unitPT, unitAbsEta);

  //map of inputs to 1D histograms
  map<string, string> hist1DMap;
  hist1DMap["numeratorPT"] = unitPT;
  hist1DMap["denominatorPT"] = unitPT;
  hist1DMap["numeratorEta"] = unitEta;
  hist1DMap["denominatorEta"] = unitEta;
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

  //make individual efficiency plots for signal
//   vector<string> effInputFiles;
//   effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/HPSTauEff.root");
//   effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/HPSTauPTEff.root");
//   effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/HPSTauPTEtaEff.root");
//   effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/HPSTauPTEtaDMFEff.root");
//   effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/HPSTauPTEtaDMFIsoEff.root");
//   effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/softMuEtaEff.root");
//   effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/softMuPTEff.root");
//   for (vector<string>::const_iterator iFile = effInputFiles.begin(); iFile != effInputFiles.end(); 
//        ++iFile) {
//     const unsigned int strLen = iFile->find(".root");
//     const string outputFileName(iFile->substr(0, strLen) + "_final.root");
//     double tempWeight = 1.;
//     plotNice(*iFile, effHistMap1DTau, effHistMap2D, binLabelMap, hist1DMapTau, outputFileName, "noPDF", tempWeight);
//   }

//   //make HLT efficiency plots for signal
//   vector<string> effInputFiles;
//   effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/genWRecoMuonEtaHLTEff.root");
//   effInputFiles.
//     push_back("../../GenMatchedRecoObjectProducer/test/genWRecoMuonEtaPTIDHLTEff.root");
//   effInputFiles.
//     push_back("../../GenMatchedRecoObjectProducer/test/genA1TauRecoMuonEtaHLTEff.root");
//   effInputFiles.
//     push_back("../../GenMatchedRecoObjectProducer/test/genA1TauRecoMuonEtaPTIDHLTEff.root");
//   for (vector<string>::const_iterator iFile = effInputFiles.begin(); iFile != effInputFiles.end(); 
//        ++iFile) {
//     const unsigned int strLen = iFile->find(".root");
//     const string outputFileName(iFile->substr(0, strLen) + "_final.root");
//     plotNice(*iFile, effHistMap, binLabelMap, hist1DMap, outputFileName, "noPDF");
//   }

//   //make HLT efficiency plots for Z-->tautau
//   vector<string> effInputFiles;
//   effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/testoutput_mi_had_muHLT.root");
//   effInputFiles.
//     push_back("../../GenMatchedRecoObjectProducer/test/all_dmf_rerun_had_muHLT.root");
//   for (vector<string>::const_iterator iFile = effInputFiles.begin(); iFile != effInputFiles.end(); 
//        ++iFile) {
//     const unsigned int strLen = iFile->find(".root");
//     const string outputFileName(iFile->substr(0, strLen) + "_final.root");
//     plotNice(*iFile, effHistMap, binLabelMap, hist1DMap, outputFileName, "noPDF");
//   }

  //make signal b veto efficiency plots
  vector<string> effInputFiles;
  //low MT files
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_gg_a5.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_gg_a7.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_gg_a9.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_gg_a11.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_gg_a13.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_gg_a15.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_Wh1_a5.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_Wh1_a7.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_Wh1_a9.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_Wh1_a11.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_Wh1_a13.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_Wh1_a15.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_ZH_a9.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_lowMT_VBF_a9.root");
  //high MT files
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_gg_a5.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_gg_a7.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_gg_a9.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_gg_a11.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_gg_a13.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_gg_a15.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_Wh1_a5.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_Wh1_a7.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_Wh1_a9.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_Wh1_a11.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_Wh1_a13.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_Wh1_a15.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_ZH_a9.root");
  effInputFiles.push_back("../../GenMatchedRecoObjectProducer/test/v2/b_veto_eff_highMT_VBF_a9.root");
  for (vector<string>::const_iterator iFile = effInputFiles.begin(); iFile != effInputFiles.end(); 
       ++iFile) {
    const unsigned int strLen = iFile->find(".root");
    const string outputFileName(iFile->substr(0, strLen) + "_final.root");
    plotNice(*iFile, effHistMap1D, effHistMap2D, binLabelMap, hist1DMap, outputFileName, "noPDF", 1.0);
  }

  //put trigger efficiency and pT distribution on same plot
//   effInputFiles.push_back("/afs/cern.ch/user/y/yohay/CMSSW_5_3_3_Git/src/BoostedTauAnalysis/GenMatchedRecoObjectProducer/test/genA1TauRecoMuonEtaHLTEff.root");
//   effInputFiles.push_back("/afs/cern.ch/user/y/yohay/CMSSW_5_3_3_Git/src/BoostedTauAnalysis/GenMatchedRecoObjectProducer/test/genWRecoMuonEtaHLTEff.root");
//   effInputFiles.push_back("ggHWMuHLTEff.root");
//   effInputFiles.push_back("WHWMuHLTEff.root");
//   for (vector<string>::const_iterator iFile = effInputFiles.begin(); iFile != effInputFiles.end(); 
//        ++iFile) {
//     const unsigned int strLen = iFile->find(".root");
//     const string outputFileName(iFile->substr(0, strLen) + "_final.root");
//     // float weight = 1.95295217378794; //Pythia WH cross section scaled by SM ggH:WH ratio
//     float weight = 3.79619; //SM ggH cross section
//     // if ((iFile - effInputFiles.begin()) == 1) weight = 0.0709988; //Pythia WH cross section
//     if ((iFile - effInputFiles.begin()) == 1) weight = 0.138806; /*SM WH cross section * 
// 								   BR(W-->leptons)*/
//     plotNice(*iFile, effHistMap1D, effHistMap2D, binLabelMap, hist1DMap, outputFileName, "noPDF", 
// 	     weight);
//   }
}
