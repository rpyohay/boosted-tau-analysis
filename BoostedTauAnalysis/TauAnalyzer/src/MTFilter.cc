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
#include "DataFormats/PatCandidates/interface/MET.h"

using namespace edm;
using namespace reco;
using namespace std;

//
// class declaration
//

template<class T, class U>
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
  bool passFilter_;
  edm::InputTag METTag_;
  edm::InputTag objTag_;

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
template<class T, class U>
MTFilter<T, U>::MTFilter(const edm::ParameterSet& iConfig) :
  minMT_(iConfig.getParameter<double>("minMT")),
  passFilter_(iConfig.getParameter<bool>("passFilter")),
  METTag_(iConfig.getParameter<edm::InputTag>("METTag")),
  objTag_(iConfig.getParameter<edm::InputTag>("objTag"))
{
  //now do what ever initialization is needed
}

template<class T, class U>
MTFilter<T, U>::~MTFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
template<class T, class U>
bool
MTFilter<T, U>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{   
  //get MET
  edm::Handle<edm::View<T> > pMET;
  iEvent.getByLabel(METTag_, pMET);
  edm::RefToBase<T> METRefToBase = pMET->refAt(0);
  //get objects
  edm::Handle<edm::View<U> > pObjs;
  iEvent.getByLabel(objTag_, pObjs);
  edm::RefToBase<U> objRefToBase = pObjs->refAt(0);

  bool objsExist = (METRefToBase.isNonnull() && objRefToBase.isNonnull());

//   //find the highest pT object
//   std::vector<reco::MuonRef> WMuonRefs;
//   for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
//        ++iMuon) { WMuonRefs.push_back(*iMuon); }
//   Common::sortByPT(WMuonRefs);

  if (objsExist) {
    double MT = sqrt(2*objRefToBase->pt()*METRefToBase->et()*
		     (1.0 - cos(reco::deltaPhi(objRefToBase->phi(), METRefToBase->phi()))));

    if ((passFilter_ && (MT > minMT_)) || (!passFilter_ && (MT <= minMT_)))
      return true;
    else
      return false;
  }
  else return false;

  /* if (MT < minMT_ && METRefToBase->et() < 30.)
    return false;
  else
    return true;
  */

  /*  if (MT >= (60. - ((60./55.)*METRefToBase->et())))
    return true;
  else
    return false;
  */
}

// ------------ method called once each job just before starting event loop  ------------
template<class T, class U>
void 
MTFilter<T, U>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T, class U>
void 
MTFilter<T, U>::endJob() {
}

// ------------ method called when starting to processes a run  ------------
template<class T, class U>
bool 
MTFilter<T, U>::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
template<class T, class U>
bool 
MTFilter<T, U>::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T, class U>
bool 
MTFilter<T, U>::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T, class U>
bool 
MTFilter<T, U>::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T, class U>
void
MTFilter<T, U>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef MTFilter<pat::MET, reco::Muon> PATMTFilter;
typedef MTFilter<pat::MET, reco::PFTau> HPSPATMTFilter;
DEFINE_FWK_MODULE(PATMTFilter);
DEFINE_FWK_MODULE(HPSPATMTFilter);
