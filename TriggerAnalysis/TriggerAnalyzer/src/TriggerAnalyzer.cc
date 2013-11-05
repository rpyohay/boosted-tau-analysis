// -*- C++ -*-
//
// Package:    TriggerAnalyzer
// Class:      TriggerAnalyzer
// 
/**\class TriggerAnalyzer TriggerAnalyzer.cc TriggerAnalysis/TriggerAnalyzer/src/TriggerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Tue Mar 22 17:13:40 CET 2011
// $Id: TriggerAnalyzer.cc,v 1.2 2012/11/30 10:14:57 yohay Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "TH1F.h"
#include "TFile.h"
//
// class declaration
//

class TriggerAnalyzer : public edm::EDAnalyzer {

public:
  explicit TriggerAnalyzer(const edm::ParameterSet&);
  ~TriggerAnalyzer();


private:
  virtual void beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup);
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  //format the trigger menu
  template<typename T>
  void format(const char* HLTPathName, T prescale)
  {
    std::cout.width(62);
    std::cout.fill(' ');
    std::cout << std::left << HLTPathName;
    std::cout.width(13);
    std::cout << std::right << prescale;
    std::cout.width(10);
    std::cout << std::endl;
  }

  //print HLT menu with prescale values
  void printTriggerMenu(const edm::Event&, edm::EventSetup const&);

  //retrieve collection from the event
  template<typename T>
  const bool getCollection_(T& pCollection, const edm::InputTag& tag, const edm::Event& iEvent)
  {
    bool collectionFound = false;
    try { collectionFound = iEvent.getByLabel(tag, pCollection); }
    catch (cms::Exception& ex) {}
    if (!collectionFound) {
      std::cerr << "No collection of type " << tag << " found in run " << iEvent.run();
      std::cerr << ", event " << iEvent.id().event() << ", lumi section ";
      std::cerr << iEvent.getLuminosityBlock().luminosityBlock() << ".\n";
    }
    return collectionFound;
  }

      // ----------member data ---------------------------

  //input
  std::string HLTProcessName_;
  edm::InputTag triggerResultsTag_;
  std::vector<std::string> unprescaledHLTPaths_;
  std::string outputFile_;

  //HLT stuff
  HLTConfigProvider HLTCfg_;
  TH1F* firedTriggers_;
  unsigned int evt_;

  //output
  TFile* out_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerAnalyzer::TriggerAnalyzer(const edm::ParameterSet& iConfig) :
  HLTProcessName_(iConfig.getUntrackedParameter<std::string>("HLTProcessName", "HLT")),
  triggerResultsTag_(iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag", 
								  edm::InputTag("TriggerResults", 
										"", 
										HLTProcessName_))),
  outputFile_(iConfig.getUntrackedParameter<std::string>("outputFile", "trigger_analysis.root"))

{
   //now do what ever initialization is needed
  std::vector<std::string> unprescaledHLTPaths;
  unprescaledHLTPaths.push_back("HLT_Photon75_CaloIdVL_IsoL_v1");
  unprescaledHLTPaths.push_back("HLT_DoublePhoton5_IsoVL_CEP_v1");
  unprescaledHLTPaths.push_back("HLT_DoublePhoton33_v1");
  unprescaledHLTPaths.push_back("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1");
  unprescaledHLTPaths.push_back("HLT_Photon26_IsoVL_Photon18_v1");
  unprescaledHLTPaths.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_v1");
  unprescaledHLTPaths.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1");
  unprescaledHLTPaths.push_back("HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1");
  unprescaledHLTPaths.push_back("HLT_Photon60_CaloIdL_HT200_v1");
  unprescaledHLTPaths.push_back("HLT_Photon70_CaloIdL_HT200_v1");
  unprescaledHLTPaths.push_back("HLT_Photon70_CaloIdL_HT300_v1");
  unprescaledHLTPaths.push_back("HLT_Photon70_CaloIdL_MHT30_v1");
  unprescaledHLTPaths.push_back("HLT_Photon70_CaloIdL_MHT50_v1");
  unprescaledHLTPaths_ = iConfig.getUntrackedParameter<std::vector<std::string> >("unprescaledHLTPaths", unprescaledHLTPaths);
}


TriggerAnalyzer::~TriggerAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ++evt_;
  printTriggerMenu(iEvent, iSetup);

  //get the trigger collections
  edm::Handle<edm::TriggerResults> pTrgResults;
  if (getCollection_(pTrgResults, triggerResultsTag_, iEvent)) {

    //loop over the HLT paths in this event
    const edm::TriggerNames& trgNames = iEvent.triggerNames(*pTrgResults);
    const std::vector<std::string> HLTPathNames = HLTCfg_.triggerNames();
    for (std::vector<std::string>::const_iterator iName = HLTPathNames.begin(); 
	 iName != HLTPathNames.end(); ++iName) {

      //if it's a photon path, fill the histogram
      if ((*iName).find("Photon") != std::string::npos) {
	const unsigned int trgIndex = trgNames.triggerIndex(*iName);
	if ((trgIndex < trgNames.size()) && (pTrgResults->accept(trgIndex))) {
	  for (std::vector<std::string>::const_iterator iPath = unprescaledHLTPaths_.begin(); 
	       iPath != unprescaledHLTPaths_.end(); ++iPath) {

	    /*the input unprescaled HLT paths don't have to be exactly the same as the actual HLT 
	      path names in the menu*/
	    if ((*iName).find(*iPath) != std::string::npos) {
	      firedTriggers_->Fill((*iPath).c_str(), 1);
	    }
	  }
	}
      }
    }
  }
}

// ---- method called once each run  ---
void TriggerAnalyzer::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
  bool changed = false;
  if (!HLTCfg_.init(iRun, iSetup, HLTProcessName_, changed) ) {
    edm::LogError("TriggerAnalyzer") << "Error: can't initialize HLTConfigProvider.\n";
    throw cms::Exception("HLTConfigProvider::init() returned non-0.\n");
  }
  if (changed) std::cout << "HLT configuration changed!\n";
}

// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyzer::beginJob()
{
  evt_ = 0;
  firedTriggers_ = new TH1F("firedTriggers", "", unprescaledHLTPaths_.size(), -0.5, 
			    unprescaledHLTPaths_.size() - 0.5);
  for (std::vector<std::string>::const_iterator iPath = unprescaledHLTPaths_.begin(); 
       iPath != unprescaledHLTPaths_.end(); ++iPath) {
    firedTriggers_->GetXaxis()->SetBinLabel(iPath - unprescaledHLTPaths_.begin() + 1, 
					    (*iPath).c_str());
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyzer::endJob() {

  //normalize trigger histogram to total number of events to get a fraction firing for each trigger
  firedTriggers_->Scale(1.0/(float)evt_);

  //write histograms to file and deallocate memory
  try { out_ = new TFile(outputFile_.c_str(), "RECREATE"); }
  catch (cms::Exception& ex) {}
  if (out_->IsOpen()) {
    out_->cd();
    firedTriggers_->Write();
    out_->Write();
    out_->Close();
  }
  else std::cerr << "Error: unable to open file " << outputFile_ << ".\n";
  delete out_;
  delete firedTriggers_;
}

void TriggerAnalyzer::printTriggerMenu(const edm::Event& iEvent, 
				       edm::EventSetup const& iSetup)
{
  std::cout << "---------------------------HLT Menu----------------------------------\n\n";
  std::cout << "Run " << iEvent.run() << ", event " << iEvent.id().event() << ", lumi section ";
  std::cout << iEvent.getLuminosityBlock().luminosityBlock() << std::endl << std::endl;
  format("HLT Path Name", "Prescale");
  const std::vector<std::string> HLTPathNames = HLTCfg_.triggerNames();
  for (std::vector<std::string>::const_iterator iName = HLTPathNames.begin(); 
       iName != HLTPathNames.end(); ++iName) {
    format((*iName).c_str(), HLTCfg_.prescaleValue(iEvent, iSetup, *iName));
  }
  std::cout << std::endl;
  std::cout << "---------------------------------------------------------------------\n\n";
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzer);
