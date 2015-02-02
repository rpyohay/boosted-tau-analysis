// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      TauMuTauXMassSelector
// 
/**\class TauMuTauXMassSelector TauMuTauXMassSelector.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/TauMuTauXMassSelector.cc

 Description: create a collection consisting of the highest mu+X mass tau to put in the event

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: TauMuTauXMassSelector.cc,v 1.6 2012/12/12 16:02:05 yohay Exp $
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
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

//
// class declaration
//

class TauMuTauXMassSelector : public edm::EDFilter {
public:
  explicit TauMuTauXMassSelector(const edm::ParameterSet&);
  ~TauMuTauXMassSelector();
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

  //input tag for jet muon map
  edm::InputTag jetMuonMapTag_;
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
TauMuTauXMassSelector::TauMuTauXMassSelector(const edm::ParameterSet& iConfig) :
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  jetMuonMapTag_(iConfig.getParameter<edm::InputTag>("jetMuonMapTag"))
{
  produces<reco::PFTauRefVector>();
  produces<reco::MuonRefVector>();
}

TauMuTauXMassSelector::~TauMuTauXMassSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool TauMuTauXMassSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointers to output collections
  std::auto_ptr<reco::PFTauRefVector> tauColl(new reco::PFTauRefVector);
  std::auto_ptr<reco::MuonRefVector> muonColl(new reco::MuonRefVector);

  //get taus
  edm::Handle<reco::PFTauRefVector> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  iEvent.getByLabel(jetMuonMapTag_, pMuonJetMap);

  //sort selected taus by descending order in mu+X mass
  std::vector<reco::PFTauRef> muXMassSortedTaus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) { muXMassSortedTaus.push_back(*iTau); }
  Common::sortByMass(pMuonJetMap, muXMassSortedTaus);

  //put the highest mu+X mass tau and muon in the event
  bool tauExists = (muXMassSortedTaus.size() > 0);
  if (tauExists) {
    maxM comp;
    comp.setMuonJetMap(&pMuonJetMap);
    reco::MuonRef highestPTTauMu = comp.highestPTMu(muXMassSortedTaus[0]);
    tauColl->push_back(muXMassSortedTaus[0]);
    if (highestPTTauMu.isNonnull()) muonColl->push_back(highestPTTauMu);
    iEvent.put(tauColl);
    iEvent.put(muonColl);
    comp.deleteMuonJetMap();
  }
  return tauExists;
}

// ------------ method called once each job just before starting event loop  ------------
void TauMuTauXMassSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void TauMuTauXMassSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool TauMuTauXMassSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool TauMuTauXMassSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool TauMuTauXMassSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool TauMuTauXMassSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TauMuTauXMassSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TauMuTauXMassSelector);
