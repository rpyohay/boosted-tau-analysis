// -*- C++ -*-
//
// Package:    DelPhiFilter
// Class:      DelPhiFilter
// 
/**\class DelPhiFilter DelPhiFilter.cc BoostedTauAnalysis/DelPhiFilter/src/DelPhiFilter.cc

 Description: filter based on |delPhi| between W decay muon and tau decay muon

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
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// class declaration
//

class DelPhiFilter : public edm::EDFilter {
   public:
      explicit DelPhiFilter(const edm::ParameterSet&);
      ~DelPhiFilter();

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

  double delPhiCut_;
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
DelPhiFilter::DelPhiFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  delPhiCut_ = iConfig.getParameter<double>("delPhiCut");
  WMuonTag_ = iConfig.getParameter<edm::InputTag>("WMuonTag");
  tauMuonTag_ = iConfig.getParameter<edm::InputTag>("tauMuonTag");
}


DelPhiFilter::~DelPhiFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
DelPhiFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get W muons
  edm::Handle<reco::MuonRefVector> WMuons;
  iEvent.getByLabel(WMuonTag_, WMuons);

  //get tau muons
  edm::Handle<reco::MuonRefVector> tauMuons;
  iEvent.getByLabel(tauMuonTag_, tauMuons);

  double Wmu_pt = -9999.;
  double Wmu_phi = 0.;
  double taumu_pt = -9999.;
  double taumu_phi = 0.;

  if (WMuons->size() == 1)
    Wmu_phi = WMuons->at(0)->phi();
  else
    { // if >1 W muons
      for (reco::MuonRefVector::const_iterator iMuon = WMuons->begin(); iMuon != WMuons->end(); ++iMuon)
	{ // loop to find highest-pT
	  if ((*iMuon)->pt() > Wmu_pt)
	    {
	      Wmu_pt = (*iMuon)->pt();
	      Wmu_phi = (*iMuon)->phi();
	    }
	} // loop to find highest-pT
    } // if >1 W muons

  if (tauMuons->size() == 1)
    taumu_phi = tauMuons->at(0)->phi();
  else
    { // if >1 tau muons
      for (reco::MuonRefVector::const_iterator iMuon = tauMuons->begin(); iMuon != tauMuons->end(); ++iMuon)
	{ // loop to find highest-pT
	  if ((*iMuon)->pt() > taumu_pt)
	    {
	      taumu_pt = (*iMuon)->pt();
	      taumu_phi = (*iMuon)->phi();
	    }
	} // loop to find highest-pT
    } // if >1 tau muons

  double Pi = 3.14159265359;
  double deltaphi = Wmu_phi - taumu_phi;
  while (deltaphi > Pi)
    deltaphi -= 2.0*Pi;
  while (deltaphi < -1.0*Pi)
    deltaphi += 2.0*Pi;
  
  if (fabs(deltaphi) > delPhiCut_)
    return true;
  else
    return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DelPhiFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DelPhiFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
DelPhiFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
DelPhiFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
DelPhiFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
DelPhiFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DelPhiFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(DelPhiFilter);
