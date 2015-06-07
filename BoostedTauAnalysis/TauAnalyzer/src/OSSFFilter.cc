// -*- C++ -*-
//
// Package:    OSSFFilter
// Class:      OSSFFilter
// 
/**\class OSSFFilter OSSFFilter.cc BoostedTauAnalysis/OSSFFilter/src/OSSFFilter.cc

 Description: OS/SF veto for charges of W decay muon and tau decay muon

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//         Created:  Mon Aug 19 11:19:05 CEST 2013
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
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// class declaration
//

class OSSFFilter : public edm::EDFilter {
   public:
      explicit OSSFFilter(const edm::ParameterSet&);
      ~OSSFFilter();

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

  edm::InputTag WMuonTag_;
  edm::InputTag tauTag_;
  edm::InputTag jetMuonMapTag_;
  bool passFilter_;
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
OSSFFilter::OSSFFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  WMuonTag_ = iConfig.getParameter<edm::InputTag>("WMuonTag");
  tauTag_ = iConfig.getParameter<edm::InputTag>("tauTag");
  jetMuonMapTag_ = iConfig.getParameter<edm::InputTag>("jetMuonMapTag");
  passFilter_ = iConfig.getParameter<bool>("passFilter");
}


OSSFFilter::~OSSFFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
OSSFFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool SignSelector = false; 
  
  //get W muons
  edm::Handle<reco::MuonRefVector> pMuons;
  iEvent.getByLabel(WMuonTag_, pMuons);
  
  //find the highest pT W muon
  std::vector<reco::MuonRef> WMuonRefs;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) { WMuonRefs.push_back(*iMuon); }
  Common::sortByPT(WMuonRefs);  
  double chargeWMuon = WMuonRefs[WMuonRefs.size() - 1]->charge();

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  iEvent.getByLabel(jetMuonMapTag_, pMuonJetMap);

  //get taus
  edm::Handle<reco::PFTauRefVector> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);

  //sort selected taus by descending order in mu+had mass
  std::vector<reco::PFTauRef> muHadMassSortedTaus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) { muHadMassSortedTaus.push_back(*iTau); }
  Common::sortByMass(pMuonJetMap, muHadMassSortedTaus);

  //loop over selected taus
  std::vector<reco::PFTauRef>::const_iterator iTau = muHadMassSortedTaus.begin();
  std::vector<reco::PFTauRef>::const_iterator endTau = iTau + 1;

  while (iTau != endTau) {
    const reco::PFJetRef& tauJetRef = (*iTau)->jetRef();
    const reco::MuonRefVector& removedMuons = (*pMuonJetMap)[tauJetRef];
    
    //find the highest pT associated muon
    std::vector<reco::MuonRef> removedMuonRefs;
    for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
	 iMuon != removedMuons.end(); ++iMuon) { removedMuonRefs.push_back(*iMuon); }
    Common::sortByPT(removedMuonRefs);
    double chargeTauMuon = removedMuonRefs[removedMuonRefs.size() - 1]->charge();
    if ((passFilter_ && (chargeTauMuon*chargeWMuon > 0)) || //select same charge muons
	(!passFilter_ && (chargeTauMuon*chargeWMuon <= 0))) { //select opposite charge muons
      SignSelector = true;
    }
    
    ++iTau; 
  }

  return SignSelector;

}

// ------------ method called once each job just before starting event loop  ------------
void 
OSSFFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
OSSFFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
OSSFFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
OSSFFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
OSSFFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
OSSFFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
OSSFFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(OSSFFilter);
