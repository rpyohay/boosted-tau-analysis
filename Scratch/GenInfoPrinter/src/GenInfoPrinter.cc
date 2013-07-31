// -*- C++ -*-
//
// Package:    GenInfoPrinter
// Class:      GenInfoPrinter
// 
/**\class GenInfoPrinter GenInfoPrinter.cc Scratch/GenInfoPrinter/src/GenInfoPrinter.cc

Description: print information about gen particles

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Sep 28 11:34:48 CEST 2012
// $Id: GenInfoPrinter.cc,v 1.1 2012/10/04 15:24:27 yohay Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//
// class declaration
//

class GenInfoPrinter : public edm::EDAnalyzer {
public:
  explicit GenInfoPrinter(const edm::ParameterSet&);
  ~GenInfoPrinter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  //recursively print all ancestors of the given gen particle
  void printAncestors(const reco::GenParticleRef&) const;

  // ----------member data ---------------------------
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
GenInfoPrinter::GenInfoPrinter(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed

}


GenInfoPrinter::~GenInfoPrinter()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenInfoPrinter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get gen particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel("genParticles", pGenParticles);

  //loop over gen particles
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {

    //identify status 3 partons
    const unsigned int PDGID = fabs(iGenParticle->pdgId());
    if /*(*/(iGenParticle->status() == 3)/* && (((PDGID >= 1) && (PDGID <= 6)) || (PDGID == 21)))*/ {

      //print all ancestors
      printAncestors(reco::GenParticleRef(pGenParticles, iGenParticle - pGenParticles->begin()));
    }
  }
  std::cout << std::endl << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
GenInfoPrinter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenInfoPrinter::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
GenInfoPrinter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenInfoPrinter::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenInfoPrinter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenInfoPrinter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void
GenInfoPrinter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//recursively print all ancestors of the given gen particle
void GenInfoPrinter::printAncestors(const reco::GenParticleRef& genParticleRef) const
{
  for (unsigned int iMother = 0; iMother < genParticleRef->numberOfMothers(); ++iMother) {
    reco::GenParticleRef motherRef = genParticleRef->motherRef(iMother);
    if (motherRef->status() == 3) {
      std::cout << "Mother " << iMother << " of gen particle " << genParticleRef.key();
      std::cout << " (with PDG ID " << genParticleRef->pdgId() << " and pT ";
      std::cout << genParticleRef->pt() << " GeV) is gen particle " << motherRef.key();
      std::cout << " with PDG ID " << motherRef->pdgId() << " and pT " << motherRef->pt();
      std::cout << " GeV.\n";
    }
    printAncestors(motherRef);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenInfoPrinter);
