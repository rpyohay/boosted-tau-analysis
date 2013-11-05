// -*- C++ -*-
//
// Package:    HadronicTauDecayFinder
// Class:      HadronicTauDecayFinder
// 
/**\class HadronicTauDecayFinder HadronicTauDecayFinder.cc BoostedTauAnalysis/TauSkimmer/plugins/HadronicTauDecayFinder.cc

 Description: return true if at least 1 hadronic tau decay was found at GEN level

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 11:29:40 CEST 2012
// $Id: HadronicTauDecayFinder.cc,v 1.4 2012/08/27 15:07:45 yohay Exp $
//
//


// system include files
#include <memory>
#include <algorithm>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"

//
// class declaration
//

class HadronicTauDecayFinder : public edm::EDFilter {
   public:
      explicit HadronicTauDecayFinder(const edm::ParameterSet&);
      ~HadronicTauDecayFinder();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      void reset();

      // ----------member data ---------------------------

      //input
      edm::InputTag genParticleTag_;
      int momPDGID_;
      double maxAbsEta_;
      edm::ParameterSet* cfg_;

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
HadronicTauDecayFinder::HadronicTauDecayFinder(const edm::ParameterSet& iConfig) :
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  momPDGID_(iConfig.getParameter<int>("momPDGID")),
  maxAbsEta_(iConfig.getParameter<double>("maxAbsEta")),
  cfg_(const_cast<edm::ParameterSet*>(&iConfig))
{
   //now do what ever initialization is needed

}


HadronicTauDecayFinder::~HadronicTauDecayFinder()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   reset();

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HadronicTauDecayFinder::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //get gen particles
   Handle<reco::GenParticleCollection> pGenParticles;
   iEvent.getByLabel(genParticleTag_, pGenParticles);

   //loop over gen particles
   bool foundHadronicDecay = false;
   bool foundGenTauOutsideEtaAcceptance = false;
   reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin();
   while ((iGenParticle != pGenParticles->end()) && !foundGenTauOutsideEtaAcceptance) {

     //instantiate GenTauDecayID object, turning off warnings for missing pT cut parameters
     try {
       GenTauDecayID tauDecay(*cfg_, pGenParticles, iGenParticle - pGenParticles->begin(), false);

       //look for a status 3 tau from boson decay
       if (tauDecay.tauIsStatus3DecayProduct()) {

	 //found a gen tau outside eta acceptance, so quit
	 if (fabs(iGenParticle->eta()) >= maxAbsEta_) foundGenTauOutsideEtaAcceptance = true;

	 //search for a leptonic decay (stop searching when an e-nu or mu-nu pair is found)
	 else {
	   if (!foundHadronicDecay) {

	     //leptonic decay not found ==> found a hadronic decay ==> stop searching, event passes
	     if (tauDecay.tauDecayType() == GenTauDecayID::HAD) foundHadronicDecay = true;
	   }
	 }
       }
     }
     catch (std::string& ex) { throw cms::Exception("HadronicTauDecayFinder") << ex; }

     //advance to the next gen particle
     ++iGenParticle;
   }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   return (foundHadronicDecay && !foundGenTauOutsideEtaAcceptance);
}

// ------------ method called once each job just before starting event loop  ------------
void 
HadronicTauDecayFinder::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HadronicTauDecayFinder::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
HadronicTauDecayFinder::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
HadronicTauDecayFinder::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
HadronicTauDecayFinder::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
HadronicTauDecayFinder::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HadronicTauDecayFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void HadronicTauDecayFinder::reset()
{
  momPDGID_ = 0;
  maxAbsEta_ = 0.0;
  cfg_ = NULL;
}
//define this as a plug-in
DEFINE_FWK_MODULE(HadronicTauDecayFinder);
