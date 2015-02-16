// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      IsolationAnalyzer
// 
/**\class IsolationAnalyzer IsolationAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/IsolationAnalyzer.cc

   Description: analyze isolation properties of objects firing triggers

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: IsolationAnalyzer.cc,v 1.1 2012/09/19 10:57:14 yohay Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;

//
// class declaration
//

class IsolationAnalyzer : public edm::EDAnalyzer {
public:
  explicit IsolationAnalyzer(const edm::ParameterSet&);
  ~IsolationAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  //delete memory
  void reset(const bool);

  // ----------member data ---------------------------

  //hltConfig for trigger matching
  HLTConfigProvider hltConfig_;

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //gen particle tag
  edm::InputTag genParticleTag_;

  //vertex tag
  edm::InputTag vtxTag_;

  //muon tag
  edm::InputTag muonTag_;

  //PU subtraction coefficient for muon PF isolation
  double muonPFIsoPUSubtractionCoeff_;

  //muon PFRelIso cut
  double PFIsoMax_;

  //trigger stuff
  edm::InputTag triggerEventTag_;
  edm::InputTag triggerResultsTag_;
  double delRMatchingCut_;
  std::vector<edm::InputTag> hltTags_;
  edm::InputTag theRightHLTTag_;
  edm::InputTag theRightHLTSubFilter_;
  std::vector<edm::InputTag> HLTSubFilters_;

  //histogram of gen-matched reco mu PFRelIso
  TH1F* recoMuPFRelIso_;

  //histogram of combined particle isolation vs. muon pT
  TH2F* combParticleIsoVsMuonPT_;
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
IsolationAnalyzer::IsolationAnalyzer(const edm::ParameterSet& iConfig):hltConfig_(),
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  vtxTag_(iConfig.getParameter<edm::InputTag>("vtxTag")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  muonPFIsoPUSubtractionCoeff_(iConfig.getParameter<double>("muonPFIsoPUSubtractionCoeff")),
  PFIsoMax_(iConfig.getParameter<double>("PFIsoMax"))
{
  //now do what ever initialization is needed
  reset(false);
  const edm::InputTag dTriggerEventTag("hltTriggerSummaryAOD","","HLT");
  triggerEventTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag",dTriggerEventTag);
  const edm::InputTag dTriggerResults("TriggerResults","","HLT");
  // By default, trigger results are labeled "TriggerResults" with process name "HLT" in the event.
  triggerResultsTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag",dTriggerResults);
  delRMatchingCut_ = iConfig.getUntrackedParameter<double>("triggerDelRMatch", 0.30);
  hltTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("hltTags");
  //  hltConfig_ = iConfig.getParameter<HLTConfigProvider>("hltConfig");
  theRightHLTTag_ = iConfig.getParameter<edm::InputTag>("theRightHLTTag");
  theRightHLTSubFilter_ = iConfig.getParameter<edm::InputTag>("theRightHLTSubFilter");
  //Whether using HLT trigger path name or the actual trigger filter name. Trigger path is default.
  HLTSubFilters_ = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("HLTSubFilters",std::vector<edm::InputTag>());
}

IsolationAnalyzer::~IsolationAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void IsolationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //get gen particle collection
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //get muon collection
  edm::Handle<edm::View<reco::Muon> > pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  //get vertices
  edm::Handle<reco::VertexCollection> pVertices;
  iEvent.getByLabel(vtxTag_, pVertices);

  //identify the first good vertex (the "primary" (?))
  reco::Vertex* pPV = Common::getPrimaryVertex(pVertices);

   // Trigger Info
  edm::Handle<trigger::TriggerEvent> trgEvent;
  iEvent.getByLabel(triggerEventTag_,trgEvent);
  edm::Handle<edm::TriggerResults> pTrgResults;
  iEvent.getByLabel(triggerResultsTag_, pTrgResults);
  std::map<std::string, bool> triggerInMenu;
  std::string myHLTFilter = "";
  int index = 9999;

  // get names of active HLT paths in this event
  std::vector<std::string> activeHLTPathsInThisEvent = hltConfig_.triggerNames();
  std::cout << "no. of active HLT paths = " << hltConfig_.triggerNames().size() << std::endl;
  // loop over active HLT paths to search for desired path
  for (std::vector<std::string>::const_iterator iHLT = activeHLTPathsInThisEvent.begin(); 
       iHLT != activeHLTPathsInThisEvent.end(); ++iHLT) { // active paths loop
    for (std::vector<edm::InputTag>::const_iterator iMyHLT = hltTags_.begin(); 
	 iMyHLT != hltTags_.end(); ++iMyHLT) {
      if ((*iMyHLT).label() == *iHLT) {
 	cout << "######## " << *iHLT << endl;
	myHLTFilter = (*iMyHLT).label();
	triggerInMenu[(*iMyHLT).label()] = true;
 	std::cout << "(*iMyHLT).label() = " << (*iMyHLT).label() << std::endl;
 	std::cout << "hltConfig_.prescaleValue(iEvent, iSetup, *iHLT) = ";
 	std::cout << hltConfig_.prescaleValue(iEvent, iSetup, *iHLT) << std::endl;
      }
    }
  } // active paths loop

   edm::InputTag filterTag;
   // loop over these objects to see whether they match
   const trigger::TriggerObjectCollection& TOC( trgEvent->getObjects() );

   //choose the right sub-filter depending on the HLT path name
   std::vector<std::string> filters;
   try { filters = hltConfig_.moduleLabels( theRightHLTTag_.label() ); }
   catch (std::exception ex) { std::cout << "bad trigger\n"; }

   for(int i=0; i != trgEvent->sizeFilters(); ++i) {
     std::string label(trgEvent->filterTag(i).label());
     //if( label == theRightHLTSubFilter_.label() ) index = i;
     if( label.find(theRightHLTSubFilter_.label()) != std::string::npos )
       {
 	 std::cout << "filterTag label: " << label << std::endl;
	 std::cout << "found subfilter!" << std::endl;
	 index = i;
       }
   }
   //    cout << "index = " << index << endl;
   // find how many objects there are
   if (index == 9999)
     index = 0;
   const trigger::Keys& KEYS(trgEvent->filterKeys(index));
   const size_type nK(KEYS.size());
   std::cout << "nK = " << nK << std::endl;

   //did this event fire the HLT?
   const edm::TriggerNames &trgNames = iEvent.triggerNames(*pTrgResults);
   const unsigned int trgIndex = trgNames.triggerIndex(myHLTFilter);
   std::cout << "trgIndex = " << trgIndex << " and trgNames.size() = " << trgNames.size() << std::endl;
   bool firedHLT = (trgIndex < trgNames.size()) && (pTrgResults->accept(trgIndex));


   std::vector<reco::GenParticle*> genZMuPtrs;
   for(reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); iGenParticle != pGenParticles->end(); ++iGenParticle)
     {
       if ((iGenParticle->status() == 1) && (fabs(iGenParticle->pdgId()) == 13))
	 { // if status-1 muon
	   if ((iGenParticle->mother()->mother()->pdgId() == 23)/* || (iGenParticle->mother()->mother()->pdgId() == 22)*/)
	     { // if it came from a Z or gamma
	       genZMuPtrs.push_back(const_cast<reco::GenParticle*>(&*iGenParticle));
	     } // if it came from a Z or gamma
	 } // if status-1 muon
     }


   //plot combined particle isolation vs. muon pT
   std::vector<unsigned int> muonsToIgnore;
   for (unsigned int iMuon = 0; iMuon < pMuons->size(); ++iMuon) {
     edm::RefToBase<reco::Muon> muon(pMuons->refAt(iMuon));
     double dRmin = 999.;
     unsigned int pos = 1000;
     for (unsigned int i = 0; i < genZMuPtrs.size(); ++i)
       { // loop over DY gen muons and look for match
	 double dEta = muon->eta() - genZMuPtrs.at(i)->eta();
	 double dPhi = muon->phi() - genZMuPtrs.at(i)->phi();
	 double delR = sqrt(dPhi*dPhi + dEta*dEta);
	 if (delR < dRmin)
	   {
	     dRmin = delR;
	     pos = i;
	   }
       } // loop over DY gen muons and look for match
     if ((dRmin < 0.1) && (pos < genZMuPtrs.size()) && (std::find(muonsToIgnore.begin(), muonsToIgnore.end(), pos) == muonsToIgnore.end()))
       {
	 
	 bool trigger_matched = false;
	 if (firedHLT)
	   { // firedHLT
	     for (int ipart = 0; ipart != nK; ++ipart)
	       {
		 const trigger::TriggerObject& TO = TOC[KEYS[ipart]];
		 double dEta_t = TO.eta() - muon->eta();
		 double dPhi_t = TO.phi() - muon->phi();
		 double delR_t = sqrt(dPhi_t*dPhi_t + dEta_t*dEta_t);
		 if (delR_t < 0.1)
		   trigger_matched = true;
	       }
	   } // firedHLT
	 
	 double etaMax = 2.1;
	 bool isTightMu = false;
	 isTightMu = Common::isTightIsolatedRecoMuon(muon, pPV, true, muonPFIsoPUSubtractionCoeff_,
						     PFIsoMax_, etaMax, true);
	 
	 if ((muon->pt() > 25) && (fabs(muon->eta()) < etaMax)  && trigger_matched && isTightMu)
	   {
	     recoMuPFRelIso_->Fill(Common::getMuonCombPFIso(*muon, muonPFIsoPUSubtractionCoeff_)/muon->pt());
	     muonsToIgnore.push_back(pos);
	   }
       }
     combParticleIsoVsMuonPT_->Fill(muon->pt(), 
				    Common::getMuonCombPFIso(*muon, muonPFIsoPUSubtractionCoeff_));
   }
}


// ------------ method called once each job just before starting event loop  ------------
void IsolationAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book muon isolation histograms
  recoMuPFRelIso_ = new TH1F("recoMuPFRelIso", "", 2000, 0.0, 40.0);
  combParticleIsoVsMuonPT_ = 
    new TH2F("combParticleIsoVsMuonPT", "", 20, 0.0, 100.0, 40, 0.0, 40.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void IsolationAnalyzer::endJob() 
{
  //make the muon isolation canvases
  TCanvas recoMuPFRelIsoCanvas("recoMuPFRelIsoCanvas", "", 600, 600);
  Common::setCanvasOptions(recoMuPFRelIsoCanvas, 1, 0, 0);
  Common::setHistogramOptions(recoMuPFRelIso_, kBlack, 0.7, 20, 1.0, "reco #mu PFRelIso", "", 0.05);
  recoMuPFRelIso_->SetLineWidth(2);

  TCanvas combParticleIsoVsMuonPTCanvas("combParticleIsoVsMuonPTCanvas", "", 600, 600);
  Common::setCanvasOptions(combParticleIsoVsMuonPTCanvas, 0, 0, 0);
  Common::setCanvasMargins(combParticleIsoVsMuonPTCanvas, 0.2, 0.2, 0.2, 0.2);

  //format the muon isolation plots
  Common::setHistogramOptions(combParticleIsoVsMuonPT_, kBlack, 0.7, 20, 1.6, 1.0, "p_{T} (GeV)", 
			      "Combined particle isolation (GeV)");

  //draw muon isolation plots
//  recoMuPFRelIsoCanvas.cd();
//  recoMuPFRelIso_->Draw();
//  combParticleIsoVsMuonPTCanvas.cd();
//  combParticleIsoVsMuonPT_->Draw("COLZ");

  //write output file
  out_->cd();
//  recoMuPFRelIsoCanvas.Write();
  recoMuPFRelIso_->Write();
//  combParticleIsoVsMuonPTCanvas.Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void IsolationAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  std::cout << "beginRun has run" << std::endl;
  bool changed_ = true;
  if ( !hltConfig_.init(iRun,iSetup,hltTags_[0].process(),changed_) ){
    edm::LogError("TriggerObjectFilter") << 
      "Error! Can't initialize HLTConfigProvider";
    throw cms::Exception("HLTConfigProvider::init() returned non 0");
  }
}

// ------------ method called when ending the processing of a run  ------------
void IsolationAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void IsolationAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
						    edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void IsolationAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
						  edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void IsolationAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void IsolationAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (out_ != NULL)) delete out_;
  out_ = NULL;
  if ((doDelete) && (recoMuPFRelIso_ != NULL)) delete recoMuPFRelIso_;
  recoMuPFRelIso_ = NULL;
  if ((doDelete) && (combParticleIsoVsMuonPT_ != NULL)) delete combParticleIsoVsMuonPT_;
  combParticleIsoVsMuonPT_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsolationAnalyzer);
