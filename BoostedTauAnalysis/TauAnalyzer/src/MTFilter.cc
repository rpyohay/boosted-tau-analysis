// -*- C++ -*-
//
// Package:    MTFilter
// Class:      MTFilter
// 
/**\class MTFilter MTFilter.cc BoostedTauAnalysis/MTFilter/src/MTFilter.cc

 Description: [one line class summary]

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
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"

using namespace edm;
using namespace reco;
using namespace std;

//
// class declaration
//

class MTFilter : public edm::EDFilter {
   public:
      explicit MTFilter(const edm::ParameterSet&);
      ~MTFilter();

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

  double minMT_;
  edm::InputTag METTag_;
  edm::InputTag muonTag_;

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
MTFilter::MTFilter(const edm::ParameterSet& iConfig) :
  minMT_(iConfig.getParameter<double>("minMT")),
  METTag_(iConfig.getParameter<edm::InputTag>("METTag")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag"))
{
  //now do what ever initialization is needed
}


MTFilter::~MTFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MTFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{   
  //get MET
  edm::Handle<edm::View<reco::PFMET> > pMET;
  iEvent.getByLabel(METTag_, pMET);
  edm::RefToBase<reco::PFMET> METRefToBase = pMET->refAt(0);
  //get W muons
  edm::Handle<reco::MuonRefVector> pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  //find the highest pT W muon
  std::vector<reco::MuonRef> WMuonRefs;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) { WMuonRefs.push_back(*iMuon); }
  Common::sortByPT(WMuonRefs);


  double MT = sqrt(2*WMuonRefs[WMuonRefs.size() - 1]->pt()*METRefToBase->et()*
		   (1.0 - cos(reco::deltaPhi(WMuonRefs[WMuonRefs.size() - 1]->phi(), METRefToBase->phi()))));

  /*  if (MT > minMT_)
    return true;
  else
    return false;
  */

  /* if (MT < minMT_ && METRefToBase->et() < 30.)
    return false;
  else
    return true;
  */
  if (MT >= (60. - ((60./55.)*METRefToBase->et())))
    return true;
  else
    return false;

}

// ------------ method called once each job just before starting event loop  ------------
void 
MTFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MTFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
MTFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
MTFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
MTFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
MTFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MTFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MTFilter);
