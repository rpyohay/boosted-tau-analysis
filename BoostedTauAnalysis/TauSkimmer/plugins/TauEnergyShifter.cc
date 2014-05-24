// -*- C++ -*-
//
// Package:    TauEnergyShifter
// Class:      TauEnergyShifter
// 
/**\class TauEnergyShifter TauEnergyShifter.cc BoostedTauAnalysis/TauEnergyShifter/src/TauEnergyShifter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//         Created:  Tue May 20 13:01:48 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/LorentzVector.h"

using namespace edm;
using namespace reco;
using namespace std;

//
// class declaration
//

class TauEnergyShifter : public edm::EDProducer {
   public:
      explicit TauEnergyShifter(const edm::ParameterSet&);
      ~TauEnergyShifter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  // PFTau collection tag
  //  edm::InputTag baseTauTag_;

  // selected tau refs tag
  edm::InputTag tauTag_;

  // minimum tau pT
  double pTMin_;

  // percent pT shift
  double pTShift_;

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
TauEnergyShifter::TauEnergyShifter(const edm::ParameterSet& iConfig) :
  //  baseTauTag_(iConfig.getParameter<edm::InputTag>("baseTauTag")),
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  pTMin_(iConfig.getParameter<double>("pTMin")),
  pTShift_(iConfig.getParameter<double>("pTShift"))
{
   //register your products
  produces<reco::PFTauCollection>("hpsTausUpShifted");
  produces<reco::PFTauCollection>("hpsTausDownShifted");
  produces<reco::PFTauRefVector>("hpsTausUpShifted");
  produces<reco::PFTauRefVector>("hpsTausDownShifted");
  produces<edm::ValueMap<reco::PFTauRef> >("hpsTausUpShifted");
  produces<edm::ValueMap<reco::PFTauRef> >("hpsTausDownShifted");
   //now do what ever other initialization is needed
  
}


TauEnergyShifter::~TauEnergyShifter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TauEnergyShifter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::auto_ptr<reco::PFTauCollection> hpsTausUpShifted(new reco::PFTauCollection);
  std::auto_ptr<reco::PFTauCollection> hpsTausDownShifted(new reco::PFTauCollection);
  std::auto_ptr<reco::PFTauRefVector> hpsTausUpShiftedRefVector(new reco::PFTauRefVector);
  std::auto_ptr<reco::PFTauRefVector> hpsTausDownShiftedRefVector(new reco::PFTauRefVector);

  //get base tau collection
  //  edm::Handle<reco::PFTauCollection> pBaseTaus;
  //  iEvent.getByLabel(baseTauTag_, pBaseTaus);

  //get taus
  edm::Handle<reco::PFTauRefVector> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);

  std::vector<reco::PFTauRef> taus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
    { taus.push_back(*iTau); }

  //create vectors of tau refs to the original tau that was up- or down-shifted
  std::vector<reco::PFTauRef> originalUpShiftedTauRefs;
  std::vector<reco::PFTauRef> originalDownShiftedTauRefs;

  for (std::vector<reco::PFTauRef>::const_iterator iTau = taus.begin(); iTau != taus.end(); ++iTau)
    {
      //cerr << "Original pT = " << (*iTau)->pT() << endl;
      //cerr << "Original eta = " << (*iTau)->eta() << endl;
      //cerr << "Original phi = " << (*iTau)->phi() << endl;

      // shift up --> still greater than pTMin_?
      if (((*iTau)->pt())*(1.+pTShift_) > pTMin_)
	{
	  reco::PFTau upShiftedTau = *((**iTau).clone());
	  math::PtEtaPhiMLorentzVector tauP4 = upShiftedTau.polarP4();
	  double newPT = (1.+pTShift_)*tauP4.pt();
	  tauP4.SetPt(newPT);
	  upShiftedTau.setP4(tauP4);
	  //cerr << "Upshifted pT = " << upShiftedTau.pt() << endl;
	  //cerr << "Upshifted eta = " << upShiftedTau.eta() << endl;
	  //cerr << "Upshifted phi = " << upShiftedTau.phi() << endl;
	  hpsTausUpShifted->push_back(upShiftedTau);
	  originalUpShiftedTauRefs.push_back(*iTau);
	}
      // shift down --> still greater than pTMin_?
      if (((*iTau)->pt())*(1.-pTShift_) > pTMin_)
	{
	  reco::PFTau downShiftedTau = *((**iTau).clone());
	  math::PtEtaPhiMLorentzVector tauP4 = downShiftedTau.polarP4();
	  double newPT = (1.-pTShift_)*tauP4.pt();
	  tauP4.SetPt(newPT);
	  downShiftedTau.setP4(tauP4);
	  //cerr << "Downshifted pT = " << downShiftedTau.pt() << endl;
	  //cerr << "Downshifted eta = " << downShiftedTau.eta() << endl;
	  //cerr << "Downshifted phi = " << downShiftedTau.phi() << endl;
	  hpsTausDownShifted->push_back(downShiftedTau);
	  originalDownShiftedTauRefs.push_back(*iTau);
	}
    }

  const unsigned int nUpShifted = hpsTausUpShifted->size();
  const unsigned int nDownShifted = hpsTausDownShifted->size();
  const edm::OrphanHandle<reco::PFTauCollection> upShiftedRefProd = 
    iEvent.put(hpsTausUpShifted, "hpsTausUpShifted");
  const edm::OrphanHandle<reco::PFTauCollection> downShiftedRefProd = 
    iEvent.put(hpsTausDownShifted, "hpsTausDownShifted");
  for (unsigned int iTau = 0; iTau < nUpShifted; ++iTau) {
    hpsTausUpShiftedRefVector->
      push_back(reco::PFTauRef(upShiftedRefProd, iTau));
  }
  for (unsigned int iTau = 0; iTau < nDownShifted; ++iTau) {
    hpsTausDownShiftedRefVector->
      push_back(reco::PFTauRef(downShiftedRefProd, iTau));
  }
//   const edm::OrphanHandle<reco::PFTauRefVector> originalUpShiftedRefProd = 
    iEvent.put(hpsTausUpShiftedRefVector, "hpsTausUpShifted");
//   const edm::OrphanHandle<reco::PFTauRefVector> originalDownShiftedRefProd = 
  iEvent.put(hpsTausDownShiftedRefVector, "hpsTausDownShifted");

  //fill the value map of original tau refs for each energy up-shifted tau
  std::auto_ptr<edm::ValueMap<reco::PFTauRef> > 
    tauUpShiftedValMap(new edm::ValueMap<reco::PFTauRef>());
  edm::ValueMap<reco::PFTauRef>::Filler tauUpShiftedFiller(*tauUpShiftedValMap);
  tauUpShiftedFiller.insert(upShiftedRefProd, 
			    originalUpShiftedTauRefs.begin(), 
			    originalUpShiftedTauRefs.end());
  tauUpShiftedFiller.fill();
  iEvent.put(tauUpShiftedValMap, "hpsTausUpShifted");

  //fill the value map of original tau refs for each energy down-shifted tau
  std::auto_ptr<edm::ValueMap<reco::PFTauRef> > 
    tauDownShiftedValMap(new edm::ValueMap<reco::PFTauRef>());
  edm::ValueMap<reco::PFTauRef>::Filler tauDownShiftedFiller(*tauDownShiftedValMap);
  tauDownShiftedFiller.insert(downShiftedRefProd, 
			      originalDownShiftedTauRefs.begin(), 
			      originalDownShiftedTauRefs.end());
  tauDownShiftedFiller.fill();
  iEvent.put(tauDownShiftedValMap, "hpsTausDownShifted");
}

// ------------ method called once each job just before starting event loop  ------------
void 
TauEnergyShifter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauEnergyShifter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
TauEnergyShifter::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TauEnergyShifter::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TauEnergyShifter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TauEnergyShifter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauEnergyShifter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauEnergyShifter);
