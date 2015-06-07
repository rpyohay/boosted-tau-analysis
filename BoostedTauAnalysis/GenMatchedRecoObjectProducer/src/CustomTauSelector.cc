// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      CustomTauSelector
// 
/**\class CustomTauSelector CustomTauSelector.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/CustomTauSelector.cc

 Description: create a collection of custom selected hadronic taus to put in the event, and stop 
 processing if no taus are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: CustomTauSelector.cc,v 1.6 2012/12/12 16:02:05 yohay Exp $
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
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

//
// class declaration
//

template<class T>
class CustomTauSelector : public edm::EDFilter {
public:
  explicit CustomTauSelector(const edm::ParameterSet&);
  ~CustomTauSelector();
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

  //input tag for reco tau collection
  edm::InputTag tauTag_;

  //input tag for base tau collection
  edm::InputTag baseTauTag_;

  //input tag for tau isolation energy
  edm::InputTag tauHadIsoTag_;

  //input tag for clean jet collection
  edm::InputTag jetTag_;

  //input tag for map jet muon removal decisions
  edm::InputTag muonRemovalDecisionTag_;

  //input tag for overlap candidate collection
  edm::InputTag overlapCandTag_;

  //vector of input tags, 1 for each discriminator the tau should pass
  std::vector<edm::InputTag> tauDiscriminatorTags_;

  //flag indicating whether the selected taus should pass or fail the discriminator
  bool passDiscriminator_;

  //pT cut
  double pTMin_;

  //|eta| cut
  double etaMax_;

  //maximum isolation cut
  double isoMax_;

  //overlap candidate matching cut
  double dR_;

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
template<class T>
CustomTauSelector<T>::CustomTauSelector(const edm::ParameterSet& iConfig) :
  tauTag_(iConfig.existsAs<edm::InputTag>("tauTag") ? 
	  iConfig.getParameter<edm::InputTag>("tauTag") : edm::InputTag()),
  baseTauTag_(iConfig.getParameter<edm::InputTag>("baseTauTag")),
  tauHadIsoTag_(iConfig.getParameter<edm::InputTag>("tauHadIsoTag")),
  jetTag_(iConfig.existsAs<edm::InputTag>("jetTag") ? 
	  iConfig.getParameter<edm::InputTag>("jetTag") : edm::InputTag()),
  muonRemovalDecisionTag_(iConfig.existsAs<edm::InputTag>("muonRemovalDecisionTag") ? 
			  iConfig.getParameter<edm::InputTag>("muonRemovalDecisionTag") : 
			  edm::InputTag()),
  overlapCandTag_(iConfig.existsAs<edm::InputTag>("overlapCandTag") ? 
		  iConfig.getParameter<edm::InputTag>("overlapCandTag") : edm::InputTag()),
  tauDiscriminatorTags_(iConfig.getParameter<std::vector<edm::InputTag> >("tauDiscriminatorTags")),
  passDiscriminator_(iConfig.getParameter<bool>("passDiscriminator")),
  pTMin_(iConfig.getParameter<double>("pTMin")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  isoMax_(iConfig.getParameter<double>("isoMax")),
  dR_(iConfig.getParameter<double>("dR")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{
  if (((jetTag_ == edm::InputTag()) && !(muonRemovalDecisionTag_ == edm::InputTag())) || 
      (!(jetTag_ == edm::InputTag()) && (muonRemovalDecisionTag_ == edm::InputTag()))) {
    std::cerr << "Warning: only one of jetTag or muonRemovalDecisionTag was supplied.  No ";
    std::cerr << "decision on tau seed jet will be made.\n";
  }
  if ((overlapCandTag_ == edm::InputTag()) && !(jetTag_ == edm::InputTag()) && 
      !(muonRemovalDecisionTag_ == edm::InputTag())) {
    std::cerr << "Warning: both jetTag and muonRemovalDecisionTag were supplied, but not ";
    std::cerr << "overlapCandTag.  Overlap with overlap candidates will not be checked.\n";
  }
  produces<reco::PFTauRefVector>();
}

template<class T>
CustomTauSelector<T>::~CustomTauSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
template<class T>
bool CustomTauSelector<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<reco::PFTauRefVector> tauColl(new reco::PFTauRefVector);

  //get taus
  edm::Handle<reco::PFTauRefVector> pTaus;
  if (tauTag_ == edm::InputTag()) {}
  else iEvent.getByLabel(tauTag_, pTaus);

  //get base tau collection
  edm::Handle<reco::PFTauCollection> pBaseTaus;
  iEvent.getByLabel(baseTauTag_, pBaseTaus);

  //get hadronic tau deltaBeta-corrected isolation
  edm::Handle<reco::PFTauDiscriminator> pTauHadIso;
  iEvent.getByLabel(tauHadIsoTag_, pTauHadIso);

  //get tau discriminators
  std::vector<edm::Handle<reco::PFTauDiscriminator> > 
    pTauDiscriminators(tauDiscriminatorTags_.size(), edm::Handle<reco::PFTauDiscriminator>());
  for (std::vector<edm::InputTag>::const_iterator iTag = tauDiscriminatorTags_.begin(); 
       iTag != tauDiscriminatorTags_.end(); ++iTag) {
    iEvent.getByLabel(*iTag, pTauDiscriminators[iTag - tauDiscriminatorTags_.begin()]);
  }

  //get jet collection
  edm::Handle<reco::PFJetCollection> pJets;
  if (jetTag_ == edm::InputTag()) {}
  else iEvent.getByLabel(jetTag_, pJets);

  //get map of jet muon removal decisions
  edm::Handle<edm::ValueMap<bool> > pMuonRemovalDecisions;
  if (muonRemovalDecisionTag_ == edm::InputTag()) {}
  else iEvent.getByLabel(muonRemovalDecisionTag_, pMuonRemovalDecisions);

  //get overlap candidates
  edm::Handle<edm::RefVector<std::vector<T> > > pOverlapCands;
  if (overlapCandTag_ == edm::InputTag()) {}
  else iEvent.getByLabel(overlapCandTag_, pOverlapCands);

  //fill STL container of pointers to overlap candidates
  std::vector<T*> overlapCandPtrs;
  if (pOverlapCands.isValid()) {
    for (typename edm::RefVector<std::vector<T> >::const_iterator iOverlapCand = 
	   pOverlapCands->begin(); iOverlapCand != pOverlapCands->end(); 
	 ++iOverlapCand) { overlapCandPtrs.push_back(const_cast<T*>(iOverlapCand->get())); }
  }

//   //debug
//   std::cerr << "Jets " << pJets.isValid() << std::endl;
//   std::cerr << "Value map " << pMuonRemovalDecisions.isValid() << std::endl;
//   std::cerr << "Taus " << pTaus.isValid() << std::endl;
//   std::cerr << "Base taus " << pBaseTaus.isValid() << std::endl;

  //fill STL container with taus passing specified discriminators in specified eta and pT range
  std::vector<reco::PFTauRef> taus = pTaus.isValid() ? 
    Common::getRecoTaus(pTaus, pBaseTaus, pTauDiscriminators, pTauHadIso, pTMin_, etaMax_, 
			passDiscriminator_, isoMax_) : 
    Common::getRecoTaus(pBaseTaus, pTauDiscriminators, pTauHadIso, pTMin_, etaMax_, 
			passDiscriminator_, isoMax_);

  //loop over selected taus
  unsigned int nPassingTaus = 0;
  for (std::vector<reco::PFTauRef>::const_iterator iTau = taus.begin(); iTau != taus.end(); 
       ++iTau) {

    //find the nearest overlap candidate to the tau
    int nearestMuonIndex = -1;
    const reco::Candidate* nearestMuon = 
      Common::nearestObject(*iTau, overlapCandPtrs, nearestMuonIndex);

    //if tau doesn't overlap with overlap candidate (or no overlap checking requested)...
    if (!(pOverlapCands.isValid()) || 
	((nearestMuon != NULL) && (reco::deltaR(**iTau, *nearestMuon) > dR_))) {

      /*...if jet collection and muon removal decision map exist, fill output collection if tau is 
	matched to jet tagged for muon removal*/
      if (pJets.isValid() && pMuonRemovalDecisions.isValid()) {
	if ((*pMuonRemovalDecisions)[(*iTau)->jetRef()]) {
	  tauColl->push_back(*iTau);
	  ++nPassingTaus;
	}
      }

      /*...if jet collection and muon removal decision map do not exist, assume no selection on 
	tau seed jet is desired and fill output collection*/
      else {
	tauColl->push_back(*iTau);
	++nPassingTaus;
      }
    }
  }
  iEvent.put(tauColl);

  //if not enough taus passing cuts were found in this event, stop processing
  return (nPassingTaus >= minNumObjsToPassFilter_);
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void CustomTauSelector<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void CustomTauSelector<T>::endJob() {
}

// ------------ method called when starting to processes a run  ------------
template<class T>
bool CustomTauSelector<T>::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
bool CustomTauSelector<T>::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
bool CustomTauSelector<T>::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
bool CustomTauSelector<T>::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void CustomTauSelector<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef CustomTauSelector<reco::Muon> CustomTauSepFromMuonSelector;
typedef CustomTauSelector<reco::Photon> CustomTauSepFromPhotonSelector;
DEFINE_FWK_MODULE(CustomTauSepFromMuonSelector);
DEFINE_FWK_MODULE(CustomTauSepFromPhotonSelector);
