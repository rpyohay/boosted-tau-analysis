// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      GenPartonProducer
// 
/**\class GenPartonProducer GenPartonProducer.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/GenPartonProducer.cc

Description: produce a collection of gen partons (status 3 partons with no status 3 daughters)

Implementation:

*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Thu Aug 23 11:23:58 CEST 2012
// $Id: GenPartonProducer.cc,v 1.2 2012/11/08 16:41:48 yohay Exp $
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
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"

//
// class declaration
//

class GenPartonProducer : public edm::EDFilter {
public:
  explicit GenPartonProducer(const edm::ParameterSet&);
  ~GenPartonProducer();

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

  //input tag for gen particle collection
  edm::InputTag genParticleTag_;

  //max |eta| of the partons
  double partonAbsEtaMax_;

  //min pT of the partons
  double partonPTMin_;

  //minimum number of gen objects passing cuts that must be found for event to pass filter
  unsigned int minNumGenObjectsToPassFilter_;
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
GenPartonProducer::GenPartonProducer(const edm::ParameterSet& iConfig) :
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  partonAbsEtaMax_(iConfig.getParameter<double>("partonAbsEtaMax")),
  partonPTMin_(iConfig.getParameter<double>("partonPTMin")),
  minNumGenObjectsToPassFilter_
  (iConfig.getParameter<unsigned int>("minNumGenObjectsToPassFilter"))
{
  //register your products
  produces<reco::GenParticleRefVector>();

  //now do what ever other initialization is needed
  
}


GenPartonProducer::~GenPartonProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
bool GenPartonProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get GEN particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //declare pointer to output collection to produce
  std::auto_ptr<reco::GenParticleRefVector> genObjs(new reco::GenParticleRefVector);

  //loop over gen particles
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {

    //save any parton passing pT and |eta| cuts with status 3 or with status 3 mother...
    const unsigned int absPDGID = fabs(iGenParticle->pdgId());
    if ((((absPDGID >= GenTauDecayID::D) && (absPDGID <= GenTauDecayID::T)) || 
	 (absPDGID == GenTauDecayID::G)) && 
	((partonAbsEtaMax_ == -1.0) || (fabs(iGenParticle->eta()) < partonAbsEtaMax_)) && 
	((partonPTMin_ == -1.0) || (iGenParticle->pt() > partonPTMin_)) && 
	((iGenParticle->status() == 3) || 
	 ((iGenParticle->numberOfMothers() > 0) && (iGenParticle->mother()->status() == 3)))) {

      //...with no status 3 daughters
      bool foundStatus3Daughter = false;
//       reco::Candidate::const_iterator iDaughter = iGenParticle->begin();
//       while ((iDaughter != iGenParticle->end()) && (!foundStatus3Daughter)) {
// 	if (iDaughter->status() == 3) foundStatus3Daughter = true;
// 	++iDaughter;
//       }
      if (!foundStatus3Daughter) genObjs->
	push_back(reco::GenParticleRef(pGenParticles, iGenParticle - pGenParticles->begin()));
    }
  }

  //flag indicating whether enough gen-matched reco objects were found
  bool foundGenObject = genObjs->size() >= minNumGenObjectsToPassFilter_;

  //put output collection into event
  iEvent.put(genObjs); //this function frees the auto_ptr argument

  //stop processing if not enough gen-matched objects were found
  return foundGenObject;
}

// ------------ method called once each job just before starting event loop  ------------
void GenPartonProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GenPartonProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool GenPartonProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool GenPartonProducer::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool GenPartonProducer::beginLuminosityBlock(edm::LuminosityBlock&, 
							   edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool GenPartonProducer::endLuminosityBlock(edm::LuminosityBlock&, 
							 edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void GenPartonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenPartonProducer);
