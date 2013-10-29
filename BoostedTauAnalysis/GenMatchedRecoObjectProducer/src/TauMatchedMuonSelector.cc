// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      TauMatchedMuonSelector
// 
/**\class TauMatchedMuonSelector TauMatchedMuonSelector.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/TauMatchedMuonSelector.cc

   Description: create a collection of muons with associated taus (e.g. together they form a 
   mu+had object) to put in the event, and stop processing if no muons are selected

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: TauMatchedMuonSelector.cc,v 1.6 2012/12/12 16:02:05 yohay Exp $
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
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"

//
// class declaration
//

class TauMatchedMuonSelector : public edm::EDFilter {
public:
  explicit TauMatchedMuonSelector(const edm::ParameterSet&);
  ~TauMatchedMuonSelector();
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

  //input tag for reco tau collection
  edm::InputTag tauTag_;

  //input tag for muon collection
  edm::InputTag muonTag_;

  //jet-muon map tag
  edm::InputTag jetMuonMapTag_;

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
TauMatchedMuonSelector::TauMatchedMuonSelector(const edm::ParameterSet& iConfig) :
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  jetMuonMapTag_(iConfig.getParameter<edm::InputTag>("jetMuonMapTag")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{
  produces<reco::MuonRefVector>();
}


TauMatchedMuonSelector::~TauMatchedMuonSelector()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool TauMatchedMuonSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<reco::MuonRefVector> muonColl(new reco::MuonRefVector);

  //get taus
  edm::Handle<edm::View<reco::PFTau> > pTaus;
  iEvent.getByLabel(tauTag_, pTaus);

  //get muons
  edm::Handle<reco::MuonRefVector> pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  iEvent.getByLabel(jetMuonMapTag_, pMuonJetMap);

  //loop over muons
  unsigned int nPassingMuons = 0;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) {

    //loop over taus
    for (unsigned int iTau = 0; iTau < pTaus->size(); ++iTau) {

      //get removed muon collection associated to tau
      const reco::PFJetRef& tauJetRef = pTaus->refAt(iTau)->jetRef();
      const reco::MuonRefVector& removedMuons = (*pMuonJetMap)[tauJetRef];

      //fill an STL container of removed muon ref keys
      std::vector<unsigned int> removedMuonRefKeys;
      for (reco::MuonRefVector::const_iterator iRemovedMuon = removedMuons.begin(); 
	   iRemovedMuon != removedMuons.end(); ++iRemovedMuon) {
	removedMuonRefKeys.push_back(iRemovedMuon->key());
      }

      //if muon is in removed muon collection for this tau...
      if (std::find(removedMuonRefKeys.begin(), removedMuonRefKeys.end(), 
		    iMuon->key()) != removedMuonRefKeys.end()) {

	//...this muon passes
	muonColl->push_back(*iMuon);
	++nPassingMuons;
      }
    }
  }
  iEvent.put(muonColl);

  //if not enough muons passing cuts were found in this event, stop processing
  return (nPassingMuons >= minNumObjsToPassFilter_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
TauMatchedMuonSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauMatchedMuonSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
TauMatchedMuonSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
TauMatchedMuonSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
TauMatchedMuonSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
TauMatchedMuonSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void
TauMatchedMuonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TauMatchedMuonSelector);
