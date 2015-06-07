// -*- C++ -*-
//
// Package:    CustomElectronSelector
// Class:      CustomElectronSelector
// 
/**\class CustomElectronSelector CustomElectronSelector.cc 
   BoostedTauAnalysis/CustomElectronSelector/src/CustomElectronSelector.cc

 Description: create a collection of custom selected electrons to put in the event, and stop 
 processing if no electrons are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: CustomElectronSelector.cc,v 1.4 2012/11/08 16:40:22 yohay Exp $
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
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"

#include <cmath>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>

using namespace std;
using namespace reco;
using namespace edm;

//
// class declaration
//

class CustomElectronSelector : public edm::EDFilter {
public:
  explicit CustomElectronSelector(const edm::ParameterSet&);
  ~CustomElectronSelector();
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

  //input tag for base reco candidate collection
  edm::InputTag baseCandidateTag_;

  //pt cut
  double pTMin_;

  //|eta| cut
  double etaMax_;

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
CustomElectronSelector::CustomElectronSelector(const edm::ParameterSet& iConfig) :
  baseCandidateTag_(iConfig.getParameter<edm::InputTag>("baseCandidateTag")),
  pTMin_(iConfig.getParameter<double>("pTMin")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{
  produces<reco::GsfElectronRefVector>();
}


CustomElectronSelector::~CustomElectronSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool CustomElectronSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<reco::GsfElectronRefVector> eleColl(new reco::GsfElectronRefVector);
  bool letPass = true;

  //get base PF candidates
  edm::Handle<reco::PFCandidateCollection> pBaseCandidates;
  iEvent.getByLabel(baseCandidateTag_, pBaseCandidates);

  //fill STL container with ref to candidates that are electrons
  std::vector<reco::GsfElectronRef> electrons;
  for (reco::PFCandidateCollection::const_iterator cand = pBaseCandidates->begin(); cand != pBaseCandidates->end();
       ++cand) { // loop over pElectrons

    if (cand->particleId() == PFCandidate::e)
      {
	cout << "found a PFElectron!" << endl;
	reco::GsfElectronRef eleref = cand->gsfElectronRef();
	electrons.push_back(eleref);
      }

  } // loop over pElectrons

  //fill output collection
  unsigned int nPassingElectrons = 0;
  for (std::vector<reco::GsfElectronRef>::const_iterator iElec = electrons.begin(); iElec != electrons.end(); 
       ++iElec) {
    if (((*iElec)->pt() > pTMin_) && (fabs((*iElec)->eta()) < etaMax_))
      {
	eleColl->push_back(*iElec);
	++nPassingElectrons;
	
	//       //debug
	//       std::cerr << "Selected electron pT: " << (*iElec)->pt() << " GeV\n";
	//       std::cerr << "Selected electron eta: " << (*iElec)->eta() << std::endl;
	//       std::cerr << "Selected electron phi: " << (*iElec)->phi() << std::endl;
	//       std::cerr << "Selected electron ref key: " << iElec->key() << std::endl;
      }
  }
  iEvent.put(eleColl);

  //if not enough electrons passing cuts were found in this event, stop processing
  //return (nPassingElectrons >= minNumObjsToPassFilter_);
  return letPass;
}

// ------------ method called once each job just before starting event loop  ------------
void 
CustomElectronSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CustomElectronSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
CustomElectronSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
CustomElectronSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
CustomElectronSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
CustomElectronSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CustomElectronSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(CustomElectronSelector);
