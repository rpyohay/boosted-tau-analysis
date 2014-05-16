// -*- C++ -*-
//
// Package:    CustomJetSelector
// Class:      CustomJetSelector
// 
/**\class CustomJetSelector CustomJetSelector.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/CustomJetSelector.cc

 Description: create a collection of custom selected jets to put in the event, and stop 
 processing if not enough jets are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: CustomJetSelector.cc,v 1.4 2012/11/08 16:40:22 yohay Exp $
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
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"

//
// class declaration
//

class CustomJetSelector : public edm::EDFilter {
public:
  explicit CustomJetSelector(const edm::ParameterSet&);
  ~CustomJetSelector();
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


  //input tag for selected tau collection
  edm::InputTag tauTag_;

  //input tag for overlap candidate collection
  edm::InputTag overlapCandTag_;

  //input tag for old jet collection
  edm::InputTag oldJetTag_;

  //input tag for the jet-soft-muon map
  edm::InputTag jetMuonMapTag_;

  //pT cut
  double pTMin_;

  //|eta| cut
  double absEtaMax_;

  //dR matching distance
  double dR_;

  //minimum number of objects that must be found to pass the filter
  unsigned int minNumObjsToPassFilter_;

  //maximum number of objects that must be found to pass the filter
  int maxNumObjsToPassFilter_;
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
CustomJetSelector::CustomJetSelector(const edm::ParameterSet& iConfig) :
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  overlapCandTag_(iConfig.getParameter<edm::InputTag>("overlapCandTag")),
  oldJetTag_(iConfig.getParameter<edm::InputTag>("oldJetTag")),
  jetMuonMapTag_(iConfig.getParameter<edm::InputTag>("jetMuonMapTag")),
  pTMin_(iConfig.getParameter<double>("pTMin")),
  absEtaMax_(iConfig.getParameter<double>("absEtaMax")),
  dR_(iConfig.getParameter<double>("dR")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter")),
  maxNumObjsToPassFilter_(iConfig.getParameter<int>("maxNumObjsToPassFilter"))
{
  produces<reco::PFJetCollection>();
}


CustomJetSelector::~CustomJetSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool CustomJetSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<reco::PFJetCollection> jetColl(new reco::PFJetCollection);

  //get taus
  edm::Handle<reco::PFTauRefVector> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);

  //get muons
  edm::Handle<edm::View<reco::Candidate> > pOverlapCands;
  iEvent.getByLabel(overlapCandTag_, pOverlapCands);

  //get old jets
  edm::Handle<reco::PFJetCollection> pOldJets;
  iEvent.getByLabel(oldJetTag_, pOldJets);

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  iEvent.getByLabel(jetMuonMapTag_, pMuonJetMap);

  //get AK5 PF L1FastL2L3 jet correction service
  const JetCorrector* corrector = JetCorrector::getJetCorrector("ak5PFL1FastL2L3", iSetup);

  //sort the overlap candidates in ascending order of pT
  std::vector<reco::Candidate*> overlapCandPtrs;
  for (unsigned int iOverlapCand = 0; iOverlapCand < pOverlapCands->size(); 
       ++iOverlapCand) {
    overlapCandPtrs.push_back(const_cast<reco::Candidate*>
			      (pOverlapCands->refAt(iOverlapCand).get()));
  }
  Common::sortByPT(overlapCandPtrs);

  //sort selected taus by descending order in mu+had mass
  std::vector<reco::PFTauRef> muHadMassSortedTaus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) { muHadMassSortedTaus.push_back(*iTau); }
  Common::sortByMass(pMuonJetMap, muHadMassSortedTaus);

  /*fill collection of new corrected jets
    - excluding the highest pT overlap candidate
    - excluding the highest mu+had mass tau
    - passing the pT cut
    - passing the eta cut*/
  unsigned int nPassingJets = 0;
  for (reco::PFJetCollection::const_iterator iJet = pOldJets->begin(); iJet != pOldJets->end(); 
       ++iJet) {
    reco::PFJet correctedJet = *iJet;
    double JEC = corrector->correction(*iJet, iEvent, iSetup);
    correctedJet.scaleEnergy(JEC);
    if ((reco::deltaR(*overlapCandPtrs[overlapCandPtrs.size() - 1], correctedJet) >= dR_) && 
	(reco::deltaR(**muHadMassSortedTaus.begin(), correctedJet) >= dR_) && 
	(correctedJet.pt() > pTMin_) && (fabs(correctedJet.eta()) < absEtaMax_)) {
      jetColl->push_back(correctedJet);
      ++nPassingJets;
    }
  }
  iEvent.put(jetColl);

  //if not enough jets passing cuts were found in this event, stop processing
  return ((nPassingJets >= minNumObjsToPassFilter_) && 
	  ((maxNumObjsToPassFilter_ == -1) || ((int)nPassingJets <= maxNumObjsToPassFilter_)));
}

// ------------ method called once each job just before starting event loop  ------------
void 
CustomJetSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CustomJetSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
CustomJetSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
CustomJetSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
CustomJetSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
CustomJetSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CustomJetSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(CustomJetSelector);
