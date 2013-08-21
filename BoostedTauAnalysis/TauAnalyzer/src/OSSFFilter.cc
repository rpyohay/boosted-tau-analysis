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
  edm::InputTag tauMuonTag_;

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
  tauMuonTag_ = iConfig.getParameter<edm::InputTag>("tauMuonTag");
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
  //get W muons
  edm::Handle<reco::MuonRefVector> WMuons;
  iEvent.getByLabel(WMuonTag_, WMuons);

  //get tau muons
  edm::Handle<reco::MuonRefVector> tauMuons;
  iEvent.getByLabel(tauMuonTag_, tauMuons);

  double Wmu_pt = -9999.;
  double taumu_pt = -9999.;
  double chargeW = 0.;
  double chargetau = 0.;

  if (WMuons->size() == 1)
    chargeW = WMuons->at(0)->charge();
  else
    { // if >1 W muons
      for (reco::MuonRefVector::const_iterator iMuon = WMuons->begin(); iMuon != WMuons->end(); ++iMuon)
	{ // loop to find highest-pT
	  if ((*iMuon)->pt() > Wmu_pt)
	    {
	      Wmu_pt = (*iMuon)->pt();
	      chargeW = (*iMuon)->charge();
	    }
	} // loop to find highest-pT
    } // if >1 W muons

  if (tauMuons->size() == 1)
    chargetau = tauMuons->at(0)->charge();
  else
    { // if >1 tau muons
      for (reco::MuonRefVector::const_iterator iMuon = tauMuons->begin(); iMuon != tauMuons->end(); ++iMuon)
	{ // loop to find highest-pT
	  if ((*iMuon)->pt() > Wmu_pt)
	    {
	      taumu_pt = (*iMuon)->pt();
	      chargetau = (*iMuon)->charge();
	    }
	} // loop to find highest-pT
    } // if >1 tau muons

  double charge_product = chargetau*chargeW;
  if (charge_product > 0)
    return true;
  else
    return false;
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
