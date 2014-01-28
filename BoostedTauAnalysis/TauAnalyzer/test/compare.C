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
  gROOT->LoadMacro((macroPath + "Miscellaneous.C++").c_str());

  //space-saving constant definitions
  string user(gSystem->GetFromPipe("whoami").Data());
  const string analysisFilePath("/data1/" + user + "/");
  const string fileExt(".root");

  //versions to compare
  const vector<string> versions;
  versions.push_back("_v13");
  versions.push_back("_v14");

  //compare data samples
  vector<vector<string> > dataIsoHaddInputFilesBothVersions; //BLINDED!!!
  vector<vector<string> > dataNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > dataNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > dataAllHaddInputFilesBothVersions; //BLINDED!!!
  vector<string> runEras;
  runEras.push_back("_Run2012A");
  runEras.push_back("_Run2012B");
  runEras.push_back("_Run2012C");
  runEras.push_back("_Run2012D");
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string dataSuffix(*iVersion + fileExt);
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
    for (vector<string>::const_iterator iRunEra = runEras.begin(); iRunEra != runEras.end(); 
	 ++iRunEra) {
      stringstream dataIsoName;
      dataIsoName << dataIsoPrefix << *iRunEra << dataSuffix; //BLINDED!!!
      dataIsoHaddInputFiles.push_back(dataIsoName.str());
      stringstream dataNonIsoName;
      dataNonIsoName << dataNonIsoPrefix << *iRunEra << dataSuffix;
      dataNonIsoHaddInputFiles.push_back(dataNonIsoName.str());
      stringstream dataNonIsoReweightName;
      dataNonIsoReweightName << dataNonIsoReweightPrefix << *iRunEra << dataSuffix;
      dataNonIsoReweightHaddInputFiles.push_back(dataNonIsoReweightName.str());
      stringstream dataAllName; //BLINDED!!!
      dataAllName << dataAllPrefix << *iRunEra << dataSuffix;
      dataAllHaddInputFiles.push_back(dataAllName.str());
    }
    dataIsoHaddInputFilesBothVersions.push_back(dataIsoHaddInputFiles); //BLINDED!!!
    dataNonIsoHaddInputFilesBothVersions.push_back(dataNonIsoHaddInputFiles);
    dataNonIsoReweightHaddInputFilesBothVersions.push_back(dataNonIsoReweightHaddInputFiles);
    dataAllHaddInputFilesBothVersions.push_back(dataAllHaddInputFiles); //BLINDED!!!
  }
//   compareVersions(dataIsoHaddInputFilesBothVersions); //BLINDED!!!
  compareVersions(dataNonIsoHaddInputFilesBothVersions);
  compareVersions(dataNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(dataAllHaddInputFilesBothVersions); //BLINDED!!!
/*
  //compare non-isolated W data samples
  vector<vector<string> > nonIsoWDataIsoHaddInputFilesBothVersions;
  vector<vector<string> > nonIsoWDataNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > nonIsoWDataAllHaddInputFilesBothVersions;
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string nonIsoWDataSuffix(*iVersion + fileExt);
    string nonIsoWDataIsoPrefix(analysisFilePath + 
				"nonIsoWData/analysis/nonIsoW_muHadIsoAnalysis_SingleMu");
    string nonIsoWDataIsoHaddOutputFile(nonIsoWDataIsoPrefix + nonIsoWDataSuffix);
    string nonIsoWDataNonIsoPrefix(analysisFilePath + 
				   "nonIsoWData/analysis/nonIsoW_muHadNonIsoAnalysis_SingleMu");
    string nonIsoWDataNonIsoHaddOutputFile(nonIsoWDataNonIsoPrefix + nonIsoWDataSuffix);
    string nonIsoWDataAllPrefix(analysisFilePath + 
				"nonIsoWData/analysis/nonIsoW_muHadAnalysis_SingleMu");
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
    nonIsoWDataIsoHaddInputFilesBothVersions.push_back(nonIsoWDataIsoHaddInputFiles);
    nonIsoWDataNonIsoHaddInputFilesBothVersions.push_back(nonIsoWDataNonIsoHaddInputFiles);
    nonIsoWDataAllHaddInputFilesBothVersions.push_back(nonIsoWDataAllHaddInputFiles);
  }
  compareVersions(nonIsoWDataIsoHaddInputFilesBothVersions);
  compareVersions(nonIsoWDataNonIsoHaddInputFilesBothVersions);
//   compareVersions(nonIsoWDataAllHaddInputFilesBothVersions);

  //compare Wh1 sample
  vector<vector<string> > Wh1IsoHaddInputFilesBothVersions;
  vector<vector<string> > Wh1AllHaddInputFilesBothVersions;
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string Wh1Suffix(*iVersion + fileExt);
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
    Wh1IsoHaddInputFilesBothVersions.push_back(Wh1IsoHaddInputFiles);
    Wh1AllHaddInputFilesBothVersions.push_back(Wh1AllHaddInputFiles);
  }
  compareVersions(Wh1IsoHaddInputFilesBothVersions);
//   compareVersions(Wh1AllHaddInputFilesBothVersions);

  //compare gg sample
  vector<vector<string> > ggIsoHaddInputFilesBothVersions;
  vector<vector<string> > ggAllHaddInputFilesBothVersions;
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string ggSuffix(*iVersion + fileExt);
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
    ggIsoHaddInputFilesBothVersions.push_back(ggIsoHaddInputFiles);
    ggAllHaddInputFilesBothVersions.push_back(ggAllHaddInputFiles);
  }
  compareVersions(ggIsoHaddInputFilesBothVersions);
//   compareVersions(ggAllHaddInputFilesBothVersions);

  //compare QCD Mu-enriched Pt-binned samples
  vector<vector<string> > QCDIsoHaddInputFilesBothVersions;
  vector<vector<string> > QCDNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > QCDNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > QCDAllHaddInputFilesBothVersions;
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
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string QCDSuffix(*iVersion + fileExt);
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
    QCDIsoHaddInputFilesBothVersions.push_back(QCDIsoHaddInputFiles);
    QCDNonIsoHaddInputFilesBothVersions.push_back(QCDNonIsoHaddInputFiles);
    QCDNonIsoReweightHaddInputFilesBothVersions.push_back(QCDNonIsoReweightHaddInputFiles);
    QCDAllHaddInputFilesBothVersions.push_back(QCDAllHaddInputFiles);
  }
  compareVersions(QCDIsoHaddInputFilesBothVersions);
  compareVersions(QCDNonIsoHaddInputFilesBothVersions);
  compareVersions(QCDNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(QCDAllHaddInputFilesBothVersions);

  //compare QCD b-enriched Pt-binned samples
  vector<vector<string> > QCDBIsoHaddInputFilesBothVersions;
  vector<vector<string> > QCDBNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > QCDBNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > QCDBAllHaddInputFilesBothVersions;
  vector<string> ptBinsB;
  ptBinsB.push_back("_Pt-15To30");
  ptBinsB.push_back("_Pt-30To50");
  ptBinsB.push_back("_Pt-50To150");
  ptBinsB.push_back("_Pt-150");
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string QCDBSuffix(*iVersion + fileExt);
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
    QCDBIsoHaddInputFilesBothVersions.push_back(QCDBIsoHaddInputFiles);
    QCDBNonIsoHaddInputFilesBothVersions.push_back(QCDBNonIsoHaddInputFiles);
    QCDBNonIsoReweightHaddInputFilesBothVersions.push_back(QCDBNonIsoReweightHaddInputFiles);
    QCDBAllHaddInputFilesBothVersions.push_back(QCDBAllHaddInputFiles);
  }
  compareVersions(QCDBIsoHaddInputFilesBothVersions);
  compareVersions(QCDBNonIsoHaddInputFilesBothVersions);
  compareVersions(QCDBNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(QCDBAllHaddInputFilesBothVersions);

  //compare QCD bToMu-enriched Pt-binned samples
  vector<vector<string> > QCDBMuIsoHaddInputFilesBothVersions;
  vector<vector<string> > QCDBMuNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > QCDBMuNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > QCDBMuAllHaddInputFilesBothVersions;
  vector<string> ptBinsBMu;
  ptBinsBMu.push_back("_pt15to30");
  ptBinsBMu.push_back("_pt30to50");
  ptBinsBMu.push_back("_pt50to150");
  ptBinsBMu.push_back("_pt150");
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string QCDBMuSuffix(*iVersion + fileExt);
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
    QCDBMuIsoHaddInputFilesBothVersions.push_back(QCDBMuIsoHaddInputFiles);
    QCDBMuNonIsoHaddInputFilesBothVersions.push_back(QCDBMuNonIsoHaddInputFiles);
    QCDBMuNonIsoReweightHaddInputFilesBothVersions.push_back(QCDBMuNonIsoReweightHaddInputFiles);
    QCDBMuAllHaddInputFilesBothVersions.push_back(QCDBMuAllHaddInputFiles);
  }
  compareVersions(QCDBMuIsoHaddInputFilesBothVersions);
  compareVersions(QCDBMuNonIsoHaddInputFilesBothVersions);
  compareVersions(QCDBMuNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(QCDBMuAllHaddInputFilesBothVersions);

  //compare Drell-Yan+jets ml+l- binned samples
  vector<vector<string> > DYJetsToLLIsoHaddInputFilesBothVersions;
  vector<vector<string> > DYJetsToLLNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > DYJetsToLLNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > DYJetsToLLAllHaddInputFilesBothVersions;
  vector<string> massBins;
  massBins.push_back("_M-10To50");
  massBins.push_back("_M-50");
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string DYJetsToLLSuffix(*iVersion + fileExt);
    string 
      DYJetsToLLIsoPrefix(analysisFilePath + "DYJetsToLL/analysis/muHadIsoAnalysis_DYJetsToLL");
    string DYJetsToLLIsoHaddOutputFile(DYJetsToLLIsoPrefix + DYJetsToLLSuffix);
    string DYJetsToLLNonIsoPrefix(analysisFilePath + 
				  "DYJetsToLL/analysis/muHadNonIsoAnalysis_DYJetsToLL");
    string DYJetsToLLNonIsoHaddOutputFile(DYJetsToLLNonIsoPrefix + DYJetsToLLSuffix);
    string 
      DYJetsToLLNonIsoReweightPrefix(analysisFilePath + 
				     "DYJetsToLL/analysis/muHadNonIsoReweightAnalysis_DYJetsToLL");
    string 
      DYJetsToLLNonIsoReweightHaddOutputFile(DYJetsToLLNonIsoReweightPrefix + DYJetsToLLSuffix);
    string DYJetsToLLAllPrefix(analysisFilePath + "DYJetsToLL/analysis/muHadAnalysis_DYJetsToLL");
    string DYJetsToLLAllHaddOutputFile(DYJetsToLLAllPrefix + DYJetsToLLSuffix);
    vector<string> DYJetsToLLIsoHaddInputFiles;
    vector<string> DYJetsToLLNonIsoHaddInputFiles;
    vector<string> DYJetsToLLNonIsoReweightHaddInputFiles;
    vector<string> DYJetsToLLAllHaddInputFiles;
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
    DYJetsToLLIsoHaddInputFilesBothVersions.push_back(DYJetsToLLIsoHaddInputFiles);
    DYJetsToLLNonIsoHaddInputFilesBothVersions.push_back(DYJetsToLLNonIsoHaddInputFiles);
    DYJetsToLLNonIsoReweightHaddInputFilesBothVersions.
      push_back(DYJetsToLLNonIsoReweightHaddInputFiles);
    DYJetsToLLAllHaddInputFilesBothVersions.push_back(DYJetsToLLAllHaddInputFiles);
  }
  compareVersions(DYJetsToLLIsoHaddInputFilesBothVersions);
  compareVersions(DYJetsToLLNonIsoHaddInputFilesBothVersions);
  compareVersions(DYJetsToLLNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(DYJetsToLLAllHaddInputFilesBothVersions);

  //compare ttbar sample
  vector<vector<string> > TTJetsIsoHaddInputFilesBothVersions;
  vector<vector<string> > TTJetsNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > TTJetsNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > TTJetsAllHaddInputFilesBothVersions;
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string TTJetsSuffix(*iVersion + fileExt);
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
    TTJetsIsoHaddInputFilesBothVersions.push_back(TTJetsIsoHaddInputFiles);
    TTJetsNonIsoHaddInputFilesBothVersions.push_back(TTJetsNonIsoHaddInputFiles);
    TTJetsNonIsoReweightHaddInputFilesBothVersions.push_back(TTJetsNonIsoReweightHaddInputFiles);
    TTJetsAllHaddInputFilesBothVersions.push_back(TTJetsAllHaddInputFiles);
  }
  compareVersions(TTJetsIsoHaddInputFilesBothVersions);
  compareVersions(TTJetsNonIsoHaddInputFilesBothVersions);
  compareVersions(TTJetsNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(TTJetsAllHaddInputFilesBothVersions);

  //compare single top samples
  vector<vector<string> > TIsoHaddInputFilesBothVersions;
  vector<vector<string> > TNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > TNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > TAllHaddInputFilesBothVersions;
  vector<string> singleTopSamples;
  singleTopSamples.push_back("_s-channel");
  singleTopSamples.push_back("bar_s-channel");
  singleTopSamples.push_back("_t-channel");
  singleTopSamples.push_back("bar_t-channel");
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string TSuffix(*iVersion + fileExt);
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
    TIsoHaddInputFilesBothVersions.push_back(TIsoHaddInputFiles);
    TNonIsoHaddInputFilesBothVersions.push_back(TNonIsoHaddInputFiles);
    TNonIsoReweightHaddInputFilesBothVersions.push_back(TNonIsoReweightHaddInputFiles);
    TAllHaddInputFilesBothVersions.push_back(TAllHaddInputFiles);
  }
  compareVersions(TIsoHaddInputFilesBothVersions);
  compareVersions(TNonIsoHaddInputFilesBothVersions);
  compareVersions(TNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(TAllHaddInputFilesBothVersions);

  //compare W+>=1 jet samples
  vector<vector<string> > WNJetsToLNuIsoHaddInputFilesBothVersions;
  vector<vector<string> > WNJetsToLNuNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > WNJetsToLNuNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > WNJetsToLNuAllTauHaddInputFilesBothVersions;
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string WNJetsToLNuSuffix("JetsToLNu" + *iVersion + fileExt);
    string WNJetsToLNuIsoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadIsoAnalysis_W");
    string WNJetsToLNuIsoHaddOutputFile(WNJetsToLNuIsoPrefix + "N" + WNJetsToLNuSuffix);
    string 
      WNJetsToLNuNonIsoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadNonIsoAnalysis_W");
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
    WNJetsToLNuIsoHaddInputFilesBothVersions.push_back(WNJetsToLNuIsoHaddInputFiles);
    WNJetsToLNuNonIsoHaddInputFilesBothVersions.push_back(WNJetsToLNuNonIsoHaddInputFiles);
    WNJetsToLNuNonIsoReweightHaddInputFilesBothVersions.
      push_back(WNJetsToLNuNonIsoReweightHaddInputFiles);
    WNJetsToLNuAllTauHaddInputFilesBothVersions.push_back(WNJetsToLNuAllTauHaddInputFiles);
  }
  compareVersions(WNJetsToLNuIsoHaddInputFilesBothVersions);
  compareVersions(WNJetsToLNuNonIsoHaddInputFilesBothVersions);
  compareVersions(WNJetsToLNuNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(WNJetsToLNuAllTauHaddInputFilesBothVersions);

  //compare W+jets sample
  vector<vector<string> > WJetsToLNuIsoHaddInputFilesBothVersions;
  vector<vector<string> > WJetsToLNuNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > WJetsToLNuNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > WJetsToLNuAllHaddInputFilesBothVersions;
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string WJetsToLNuSuffix(*iVersion + fileExt);
    string 
      WJetsToLNuIsoPrefix(analysisFilePath + "WJetsToLNu/analysis/muHadIsoAnalysis_WJetsToLNu");
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
    WJetsToLNuIsoHaddInputFilesBothVersions.push_back(WJetsToLNuIsoHaddInputFiles);
    WJetsToLNuNonIsoHaddInputFilesBothVersions.push_back(WJetsToLNuNonIsoHaddInputFiles);
    WJetsToLNuNonIsoReweightHaddInputFilesBothVersions.
      push_back(WJetsToLNuNonIsoReweightHaddInputFiles);
    WJetsToLNuAllHaddInputFilesBothVersions.push_back(WJetsToLNuAllHaddInputFiles);
  }
  compareVersions(WJetsToLNuIsoHaddInputFilesBothVersions);
  compareVersions(WJetsToLNuNonIsoHaddInputFilesBothVersions);
  compareVersions(WJetsToLNuNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(WJetsToLNuAllHaddInputFilesBothVersions);

  //compare WZ sample
  vector<vector<string> > WZIsoHaddInputFilesBothVersions;
  vector<vector<string> > WZNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > WZNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > WZAllHaddInputFilesBothVersions;
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string WZSuffix(*iVersion + fileExt);
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
    WZIsoHaddInputFilesBothVersions.push_back(WZIsoHaddInputFiles);
    WZNonIsoHaddInputFilesBothVersions.push_back(WZNonIsoHaddInputFiles);
    WZNonIsoReweightHaddInputFilesBothVersions.push_back(WZNonIsoReweightHaddInputFiles);
    WZAllHaddInputFilesBothVersions.push_back(WZAllHaddInputFiles);
  }
  compareVersions(WZIsoHaddInputFilesBothVersions);
  compareVersions(WZNonIsoHaddInputFilesBothVersions);
  compareVersions(WZNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(WZAllHaddInputFilesBothVersions);

  //compare ZZ sample
  vector<vector<string> > ZZIsoHaddInputFilesBothVersions;
  vector<vector<string> > ZZNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > ZZNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > ZZAllHaddInputFilesBothVersions;
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string ZZSuffix(*iVersion + fileExt);
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
    ZZIsoHaddInputFilesBothVersions.push_back(ZZIsoHaddInputFiles);
    ZZNonIsoHaddInputFilesBothVersions.push_back(ZZNonIsoHaddInputFiles);
    ZZNonIsoReweightHaddInputFilesBothVersions.push_back(ZZNonIsoReweightHaddInputFiles);
    ZZAllHaddInputFilesBothVersions.push_back(ZZAllHaddInputFiles);
  }
  compareVersions(ZZIsoHaddInputFilesBothVersions);
  compareVersions(ZZNonIsoHaddInputFilesBothVersions);
  compareVersions(ZZNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(ZZAllHaddInputFilesBothVersions);

  //compare WW sample
  vector<vector<string> > WWIsoHaddInputFilesBothVersions;
  vector<vector<string> > WWNonIsoHaddInputFilesBothVersions;
  vector<vector<string> > WWNonIsoReweightHaddInputFilesBothVersions;
  vector<vector<string> > WWAllHaddInputFilesBothVersions;
  for (vector<string>::const_iterator iVersion = versions.begin(); iVersion != versions.end(); 
       ++iVersion) {
    string WWSuffix(*iVersion + fileExt);
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
    WWIsoHaddInputFilesBothVersions.push_back(WWIsoHaddInputFiles);
    WWNonIsoHaddInputFilesBothVersions.push_back(WWNonIsoHaddInputFiles);
    WWNonIsoReweightHaddInputFilesBothVersions.push_back(WWNonIsoReweightHaddInputFiles);
    WWAllHaddInputFilesBothVersions.push_back(WWAllHaddInputFiles);
  }
  compareVersions(WWIsoHaddInputFilesBothVersions);
  compareVersions(WWNonIsoHaddInputFilesBothVersions);
  compareVersions(WWNonIsoReweightHaddInputFilesBothVersions);
//   compareVersions(WWAllHaddInputFilesBothVersions);
*/
}
