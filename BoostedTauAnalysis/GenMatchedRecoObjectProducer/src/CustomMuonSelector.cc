// -*- C++ -*-
//
// Package:    CustomMuonSelector
// Class:      CustomMuonSelector
// 
/**\class CustomMuonSelector CustomMuonSelector.cc 
   BoostedTauAnalysis/CustomMuonSelector/src/CustomMuonSelector.cc

 Description: create a collection of custom selected muons to put in the event, and stop 
 processing if no muons are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: CustomMuonSelector.cc,v 1.4 2012/11/08 16:40:22 yohay Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"

//
// class declaration
//

class CustomMuonSelector : public edm::EDFilter {
public:
  explicit CustomMuonSelector(const edm::ParameterSet&);
  ~CustomMuonSelector();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------

  //input tag for base reco muon collection
  edm::InputTag baseMuonTag_;

  //input tag for reco muon collection
  edm::InputTag muonTag_;

  //input tag for reco vertex collection
  edm::InputTag vtxTag_;

  //input tag for muons that should not be allowed to pass (i.e. they passed some other selection)
  edm::InputTag vetoMuonTag_;

  //muon ID to apply
  std::string muonID_;

  //PF isolation cut
  double PFIsoMax_;

  //detector isolation cut
  double detectorIsoMax_;

  //PU subtraction coefficient for PF isolation
  double PUSubtractionCoeff_;

  //flag indicating whether PF isolation or detector isolation should be used
  bool usePFIso_;

  //flag indicating whether the selected muons should pass the isolation cut
  bool passIso_;

  //|eta| cut
  double etaMax_;

  //minimum number of objects that must be found to pass the filter
  unsigned int minNumObjsToPassFilter_;
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
CustomMuonSelector::CustomMuonSelector(const edm::ParameterSet& iConfig) :
  baseMuonTag_(iConfig.getParameter<edm::InputTag>("baseMuonTag")),
  muonTag_(iConfig.existsAs<edm::InputTag>("muonTag") ? 
	   iConfig.getParameter<edm::InputTag>("muonTag") : edm::InputTag()),
  vtxTag_(iConfig.getParameter<edm::InputTag>("vtxTag")),
  vetoMuonTag_(iConfig.existsAs<edm::InputTag>("vetoMuonTag") ? 
	   iConfig.getParameter<edm::InputTag>("vetoMuonTag") : edm::InputTag()),
  muonID_(iConfig.getParameter<std::string>("muonID")),
  PFIsoMax_(iConfig.getParameter<double>("PFIsoMax")),
  detectorIsoMax_(iConfig.getParameter<double>("detectorIsoMax")),
  PUSubtractionCoeff_(iConfig.getParameter<double>("PUSubtractionCoeff")),
  usePFIso_(iConfig.getParameter<bool>("usePFIso")),
  passIso_(iConfig.getParameter<bool>("passIso")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{
  produces<reco::MuonRefVector>();
}


CustomMuonSelector::~CustomMuonSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool CustomMuonSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<reco::MuonRefVector> muonColl(new reco::MuonRefVector);

  //get base muons
  edm::Handle<reco::MuonCollection> pBaseMuons;
  iEvent.getByLabel(baseMuonTag_, pBaseMuons);

  //get muons
  edm::Handle<reco::MuonRefVector> pMuons;
  if (muonTag_ == edm::InputTag()) {}
  else iEvent.getByLabel(muonTag_, pMuons);

  //get vertices
  edm::Handle<reco::VertexCollection> pVertices;
  iEvent.getByLabel(vtxTag_, pVertices);

  //get veto muons
  edm::Handle<reco::MuonRefVector> pVetoMuons;
  if (vetoMuonTag_ == edm::InputTag()) {}
  else iEvent.getByLabel(vetoMuonTag_, pVetoMuons);

  //identify the first good vertex (the "primary" (?))
  reco::Vertex* pPV = Common::getPrimaryVertex(pVertices);

  /*fill STL container with muons passing the 2012 tight selection, PF isolation, and |eta| 
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId and 
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1)*/
  std::vector<reco::MuonRef> muons;
  if (muonID_ == "tight") {
    if (usePFIso_) {
      muons = pMuons.isValid() ? 
	Common::getTightPFIsolatedRecoMuons(pMuons, pBaseMuons, pPV, PUSubtractionCoeff_, 
					    PFIsoMax_, etaMax_, passIso_) : 
	Common::getTightPFIsolatedRecoMuons(pBaseMuons, pPV, PUSubtractionCoeff_, PFIsoMax_, 
					    etaMax_, passIso_);
    }

    /*fill STL container with muons passing the 2012 tight selection, detector isolation, and 
      |eta| (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId and 
      https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1)*/
    else {
      muons = pMuons.isValid() ? 
	Common::getTightDetectorIsolatedRecoMuons(pMuons, pBaseMuons, pPV, detectorIsoMax_, 
						  etaMax_, passIso_) : 
	Common::getTightDetectorIsolatedRecoMuons(pBaseMuons, pPV, detectorIsoMax_, 
						  etaMax_, passIso_);
    }
  }

  /*fill STL container of soft muons (cf. 
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Soft_Muon)*/
  else if (muonID_ == "soft") {
    muons = pMuons.isValid() ? 
      Common::getSoftRecoMuons(pMuons, pBaseMuons, pPV, etaMax_) : 
      Common::getSoftRecoMuons(pBaseMuons, pPV, etaMax_);
  }

  //error: unsupported muon ID
  else throw cms::Exception("CustomMuonSelector") << "Error: unsupported muon ID.\n";

  //make an STL container of the veto muon ref keys
  std::vector<int> vetoMuonRefKeys;
  if (pVetoMuons.isValid()) {
    for (reco::MuonRefVector::const_iterator iVetoMuon = pVetoMuons->begin(); 
	 iVetoMuon != pVetoMuons->end(); ++iVetoMuon) {
      vetoMuonRefKeys.push_back(iVetoMuon->key());
    }
  }

  //fill output collection
  unsigned int nPassingMuons = 0;
  for (std::vector<reco::MuonRef>::const_iterator iMuon = muons.begin(); iMuon != muons.end(); 
       ++iMuon) {
    if (std::find(vetoMuonRefKeys.begin(), vetoMuonRefKeys.end(), 
		  iMuon->key()) == vetoMuonRefKeys.end()) {
      muonColl->push_back(*iMuon);
      ++nPassingMuons;

//       //debug
//       std::cerr << "Selected muon pT: " << (*iMuon)->pt() << " GeV\n";
//       std::cerr << "Selected muon eta: " << (*iMuon)->eta() << std::endl;
//       std::cerr << "Selected muon phi: " << (*iMuon)->phi() << std::endl;
//       std::cerr << "Selected muon ref key: " << iMuon->key() << std::endl;
    }
  }
  iEvent.put(muonColl);

  //if not enough muons passing cuts were found in this event, stop processing
  return (nPassingMuons >= minNumObjsToPassFilter_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
CustomMuonSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CustomMuonSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
CustomMuonSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
CustomMuonSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
CustomMuonSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
CustomMuonSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CustomMuonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(CustomMuonSelector);
