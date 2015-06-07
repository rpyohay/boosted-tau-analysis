// -*- C++ -*-
//
// Package:    TriggerMuonNeighborFilter
// Class:      TriggerMuonNeighborFilter
// 
/**\class TriggerMuonNeighborFilter TriggerMuonNeighborFilter.cc BoostedTauAnalysis/TriggerMuonNeighborFilter/src/TriggerMuonNeighborFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//         Created:  Sun Apr 19 18:10:02 CEST 2015
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
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
#include "FWCore/Common/interface/TriggerNames.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;

//
// class declaration
//

template<class T>
class TriggerMuonNeighborFilter : public edm::EDFilter {
   public:
      explicit TriggerMuonNeighborFilter(const edm::ParameterSet&);
      ~TriggerMuonNeighborFilter();

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
  //input tag for reco muon collection
  edm::InputTag muonTag_;

  //input tag for reco object collection
  edm::InputTag recoObjTag_;

  //deltaR cut
  double delRMin_;

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
template<class T>
TriggerMuonNeighborFilter<T>::TriggerMuonNeighborFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  recoObjTag_ = iConfig.getParameter<edm::InputTag>("recoObjTag");
  muonTag_ = (iConfig.existsAs<edm::InputTag>("muonTag") ? 
	   iConfig.getParameter<edm::InputTag>("muonTag") : edm::InputTag());
  delRMin_ = iConfig.getParameter<double>("delRMin");

}

template<class T>
TriggerMuonNeighborFilter<T>::~TriggerMuonNeighborFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
template<class T>
bool TriggerMuonNeighborFilter<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool noNeighbor = true;

  //get muons
  edm::Handle<reco::MuonRefVector> pMuons;
  if (muonTag_ == edm::InputTag()) {}
  else iEvent.getByLabel(muonTag_, pMuons);

  //get reco objects to veto from muon isolation cone
  edm::Handle<edm::RefVector<std::vector<T> > > recoObjs;
  iEvent.getByLabel(recoObjTag_, recoObjs);

  //find the highest pT W muon
  std::vector<reco::MuonRef> trigMuonRefs;
  if (pMuons.isValid()) {
    for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
	 ++iMuon) { trigMuonRefs.push_back(*iMuon); }
  }
  Common::sortByPT(trigMuonRefs);

  //loop over reco objects
  //see which ones overlap with muon iso cone
  if (recoObjs->size() != 0)
    {
      for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj = 
	     recoObjs->begin(); iRecoObj != recoObjs->end(); 
	   ++iRecoObj)
	{
	  if (deltaR(**iRecoObj, *trigMuonRefs[trigMuonRefs.size() - 1]) < delRMin_)
	    noNeighbor = false;
	}
    }
  return noNeighbor;
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void TriggerMuonNeighborFilter<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void TriggerMuonNeighborFilter<T>::endJob() {
}

// ------------ method called when starting to processes a run  ------------
template<class T>
bool TriggerMuonNeighborFilter<T>::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
bool TriggerMuonNeighborFilter<T>::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
bool TriggerMuonNeighborFilter<T>::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
bool TriggerMuonNeighborFilter<T>::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void TriggerMuonNeighborFilter<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef TriggerMuonNeighborFilter<reco::Muon> TriggerMuonMuonFilter;
typedef TriggerMuonNeighborFilter<reco::GsfElectron> TriggerMuonElectronFilter;
typedef TriggerMuonNeighborFilter<reco::PFTau> TriggerMuonTauFilter;
DEFINE_FWK_MODULE(TriggerMuonMuonFilter);
DEFINE_FWK_MODULE(TriggerMuonElectronFilter);
DEFINE_FWK_MODULE(TriggerMuonTauFilter);
