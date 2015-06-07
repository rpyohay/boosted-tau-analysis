// -*- C++ -*-
//
// Package:    TauSkimmer
// Class:      BVetoFilter
// 
/**\class BVetoFilter BVetoFilter.cc BoostedTauAnalysis/TauSkimmer/plugins/BVetoFilter.cc

 Description: veto events with a mu+had object built from a b jet and produce a collection of 
 passing mu+had HPS taus

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,6 R-025,+41227672274,
//         Created:  Wed May 28 18:04:55 CEST 2014
// $Id$
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
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

//
// class declaration
//

class BVetoFilter : public edm::EDFilter {
   public:
      explicit BVetoFilter(const edm::ParameterSet&);
      ~BVetoFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  edm::InputTag tauTag_;
  edm::InputTag oldJetTag_;
  edm::InputTag jetMuonMapTag_;
  edm::InputTag bTagInfoTag_;
  double CSVMax_;
  bool passFilter_;
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
BVetoFilter::BVetoFilter(const edm::ParameterSet& iConfig) :
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  oldJetTag_(iConfig.getParameter<edm::InputTag>("oldJetTag")),
  jetMuonMapTag_(iConfig.getParameter<edm::InputTag>("jetMuonMapTag")),
  bTagInfoTag_(iConfig.getParameter<edm::InputTag>("bTagInfoTag")),
  CSVMax_(iConfig.getParameter<double>("CSVMax")),
  passFilter_(iConfig.getParameter<bool>("passFilter")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{
   //now do what ever initialization is needed
  produces<reco::PFTauRefVector>();
}


BVetoFilter::~BVetoFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
BVetoFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //pointer to output collection
  std::auto_ptr<reco::PFTauRefVector> tauColl(new reco::PFTauRefVector);

  //get taus
  edm::Handle<reco::PFTauRefVector> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);

  //get old jets
  edm::Handle<reco::PFJetCollection> pOldJets;
  iEvent.getByLabel(oldJetTag_, pOldJets);

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  iEvent.getByLabel(jetMuonMapTag_, pMuonJetMap);

  //get b tag information
  edm::Handle<reco::JetTagCollection> pBTagInfo;
  iEvent.getByLabel(bTagInfoTag_, pBTagInfo);
  const reco::JetTagCollection& bTags = *(pBTagInfo.product());

  //sort selected taus by descending order in mu+had mass
  std::vector<reco::PFTauRef> muHadMassSortedTaus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) { muHadMassSortedTaus.push_back(*iTau); }
  Common::sortByMass(pMuonJetMap, muHadMassSortedTaus);

  //loop over selected taus
  std::vector<reco::PFTauRef>::const_iterator iTau = muHadMassSortedTaus.begin();
  std::vector<reco::PFTauRef>::const_iterator endTau = iTau + 1;
  while (iTau != endTau) {

    //loop over b tag information
    for (unsigned int iBTagInfo = 0; iBTagInfo != pBTagInfo->size(); ++iBTagInfo) {

      //CSV score for uncleaned tau parent jet is right, save the tau
      if (reco::PFJetRef(pOldJets, iBTagInfo).key() == ((*iTau)->jetRef()).key()) {
	if ((passFilter_ && (bTags[iBTagInfo].second < CSVMax_)) || 
	    (!passFilter_ && (bTags[iBTagInfo].second >= CSVMax_))) {
	  tauColl->push_back(*iTau);
	}
      }
    }
    ++iTau;
  }

  //count number of passing taus
  const unsigned int nPassingTaus = tauColl->size();

  //put the selected tau collection into the event
  iEvent.put(tauColl);

  //if not enough taus passing cuts were found in this event, stop processing
  return (nPassingTaus >= minNumObjsToPassFilter_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
BVetoFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BVetoFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
BVetoFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
BVetoFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
BVetoFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
BVetoFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BVetoFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(BVetoFilter);
