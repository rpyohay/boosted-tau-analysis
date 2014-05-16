// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      CustomPhotonSelector
// 
/**\class CustomPhotonSelector CustomPhotonSelector.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/CustomPhotonSelector.cc

 Description: create a collection of custom selected photons to put in the event, and stop 
 processing if no photons are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: CustomPhotonSelector.cc,v 1.6 2012/12/12 16:02:05 yohay Exp $
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
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"

//
// class declaration
//

class CustomPhotonSelector : public edm::EDFilter {
public:
  explicit CustomPhotonSelector(const edm::ParameterSet&);
  ~CustomPhotonSelector();
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

  //input tag for reco photon collection
  edm::InputTag photonTag_;

  //input tag for base photon collection
  edm::InputTag basePhotonTag_;

  //input tag for beam spot collection
  edm::InputTag beamSpotTag_;

  //input tag for conversion collection
  edm::InputTag conversionTag_;

  //input tag for electron collection
  edm::InputTag electronTag_;

  //input tag for charged hadron isolation
  edm::InputTag chargedHadronIsoTag_;

  //input tag for charged hadron isolation
  edm::InputTag neutralHadronIsoTag_;

  //input tag for charged hadron isolation
  edm::InputTag photonIsoTag_;

  //input tag for rho (PU energy per unit area)
  edm::InputTag rhoTag_;

  //pT cut
  double pTMin_;

  //|eta| cut
  double etaMin_;

  //|eta| cut
  double etaMax_;

  //flag indicating whether the selected photons should pass or fail the CSEV
  bool passCSEV_;

  //maximum H/E cut
  double HOverEMax_;

  //flag indicating whether the selected photons should pass or fail the H/E cut
  bool passHOverE_;

  //maximum sigmaIetaIeta cut
  double sigmaIetaIetaMax_;

  //flag indicating whether the selected photons should pass or fail the sigmaIetaIeta cut
  bool passSigmaIetaIeta_;

  //maximum charged hadron isolation cut
  double chargedHadronIsoMax_;

  /*flag indicating whether the selected photons should pass or fail the charged hadron isolation 
    cut*/
  bool passChargedHadronIso_;

  //neutral hadron isolation cut constant
  double neutralHadronIsoConst_;

  //neutral hadron isolation cut pT multiplier
  double neutralHadronIsoPTMultiplier_;

  /*flag indicating whether the selected photons should pass or fail the neutral hadron isolation 
    cut*/
  bool passNeutralHadronIso_;

  //photon isolation cut constant
  double photonIsoConst_;

  //photon isolation cut pT multiplier
  double photonIsoPTMultiplier_;

  //flag indicating whether the selected photons should pass or fail the photon isolation cut
  bool passPhotonIso_;

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
CustomPhotonSelector::CustomPhotonSelector(const edm::ParameterSet& iConfig) :
  photonTag_(iConfig.existsAs<edm::InputTag>("photonTag") ? 
	     iConfig.getParameter<edm::InputTag>("photonTag") : edm::InputTag()),
  basePhotonTag_(iConfig.getParameter<edm::InputTag>("basePhotonTag")),
  beamSpotTag_(iConfig.getParameter<edm::InputTag>("beamSpotTag")),
  conversionTag_(iConfig.getParameter<edm::InputTag>("conversionTag")),
  electronTag_(iConfig.getParameter<edm::InputTag>("electronTag")),
  chargedHadronIsoTag_(iConfig.getParameter<edm::InputTag>("chargedHadronIsoTag")),
  neutralHadronIsoTag_(iConfig.getParameter<edm::InputTag>("neutralHadronIsoTag")),
  photonIsoTag_(iConfig.getParameter<edm::InputTag>("photonIsoTag")),
  rhoTag_(iConfig.getParameter<edm::InputTag>("rhoTag")),
  pTMin_(iConfig.getParameter<double>("pTMin")),
  etaMin_(iConfig.getParameter<double>("etaMin")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  passCSEV_(iConfig.getParameter<bool>("passCSEV")),
  HOverEMax_(iConfig.getParameter<double>("HOverEMax")),
  passHOverE_(iConfig.getParameter<bool>("passHOverE")),
  sigmaIetaIetaMax_(iConfig.getParameter<double>("sigmaIetaIetaMax")),
  passSigmaIetaIeta_(iConfig.getParameter<bool>("passSigmaIetaIeta")),
  chargedHadronIsoMax_(iConfig.getParameter<double>("chargedHadronIsoMax")),
  passChargedHadronIso_(iConfig.getParameter<bool>("passChargedHadronIso")),
  neutralHadronIsoConst_(iConfig.getParameter<double>("neutralHadronIsoConst")),
  neutralHadronIsoPTMultiplier_(iConfig.getParameter<double>("neutralHadronIsoPTMultiplier")),
  passNeutralHadronIso_(iConfig.getParameter<bool>("passNeutralHadronIso")),
  photonIsoConst_(iConfig.getParameter<double>("photonIsoConst")),
  photonIsoPTMultiplier_(iConfig.getParameter<double>("photonIsoPTMultiplier")),
  passPhotonIso_(iConfig.getParameter<bool>("passPhotonIso")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{
  produces<reco::PhotonRefVector>();
}


CustomPhotonSelector::~CustomPhotonSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool CustomPhotonSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<reco::PhotonRefVector> photonColl(new reco::PhotonRefVector);

  //get photons
  edm::Handle<reco::PhotonRefVector> pPhotons;
  if (photonTag_ == edm::InputTag()) {}
  else iEvent.getByLabel(photonTag_, pPhotons);

  //get base photon collection
  edm::Handle<reco::PhotonCollection> pBasePhotons;
  iEvent.getByLabel(basePhotonTag_, pBasePhotons);

  //get beam spot collection
  edm::Handle<reco::BeamSpot> pBeamSpot;
  iEvent.getByLabel(beamSpotTag_, pBeamSpot);

  //get conversion collection
  edm::Handle<reco::ConversionCollection> pConversions;
  iEvent.getByLabel(conversionTag_, pConversions);

  //get electron collection
  edm::Handle<reco::GsfElectronCollection> pElectrons;
  iEvent.getByLabel(electronTag_, pElectrons);

  //get charged hadron isolation
  edm::Handle<edm::ValueMap<double> > pChargedHadronIso;
  iEvent.getByLabel(chargedHadronIsoTag_, pChargedHadronIso);

  //get neutral hadron isolation
  edm::Handle<edm::ValueMap<double> > pNeutralHadronIso;
  iEvent.getByLabel(neutralHadronIsoTag_, pNeutralHadronIso);

  //get photon isolation
  edm::Handle<edm::ValueMap<double> > pPhotonIso;
  iEvent.getByLabel(photonIsoTag_, pPhotonIso);

  //get rho (PU energy per unit area)
  edm::Handle<double> pRho;
  iEvent.getByLabel(rhoTag_, pRho);

  //fill STL container with photons passing specified discriminators in specified eta and pT range
  std::vector<reco::PhotonRef> photons = pPhotons.isValid() ? 
    Common::getRecoPhotons(pPhotons, pBasePhotons, pBeamSpot, pConversions, pElectrons, 
			   pChargedHadronIso, pNeutralHadronIso, pPhotonIso, pRho, pTMin_, 
			   etaMin_, etaMax_, passCSEV_, HOverEMax_, passHOverE_, 
			   sigmaIetaIetaMax_, passSigmaIetaIeta_, chargedHadronIsoMax_, 
			   passChargedHadronIso_, neutralHadronIsoConst_, 
			   neutralHadronIsoPTMultiplier_, passNeutralHadronIso_, photonIsoConst_, 
			   photonIsoPTMultiplier_, passPhotonIso_) : 
    Common::getRecoPhotons(pBasePhotons, pBeamSpot, pConversions, pElectrons, 
			   pChargedHadronIso, pNeutralHadronIso, pPhotonIso, pRho, pTMin_, 
			   etaMin_, etaMax_, passCSEV_, HOverEMax_, passHOverE_, 
			   sigmaIetaIetaMax_, passSigmaIetaIeta_, chargedHadronIsoMax_, 
			   passChargedHadronIso_, neutralHadronIsoConst_, 
			   neutralHadronIsoPTMultiplier_, passNeutralHadronIso_, photonIsoConst_, 
			   photonIsoPTMultiplier_, passPhotonIso_);

  //put the photon collection into the event
  unsigned int nPassingPhotons = photons.size();
  for (std::vector<reco::PhotonRef>::const_iterator iPhoton = photons.begin(); 
       iPhoton != photons.end(); ++iPhoton) { photonColl->push_back(*iPhoton); }
  iEvent.put(photonColl);

  //if not enough photons passing cuts were found in this event, stop processing
  return (nPassingPhotons >= minNumObjsToPassFilter_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
CustomPhotonSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CustomPhotonSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
CustomPhotonSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
CustomPhotonSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
CustomPhotonSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
CustomPhotonSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void
CustomPhotonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(CustomPhotonSelector);
