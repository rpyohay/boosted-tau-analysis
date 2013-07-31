// -*- C++ -*-
//
// Package:    GenObjectProducer
// Class:      GenObjectProducer
// 
/**\class GenObjectProducer GenObjectProducer.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/GenObjectProducer.cc

Description: produce a collection of gen objects from boosted di-tau objects

Implementation:

*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Thu Aug 23 11:23:58 CEST 2012
// $Id: GenObjectProducer.cc,v 1.4 2012/10/10 09:11:43 yohay Exp $
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

//code for any tau decay
#define TAU_ALL 3

//
// class declaration
//

class GenObjectProducer : public edm::EDFilter {
public:
  explicit GenObjectProducer(const edm::ParameterSet&);
  ~GenObjectProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  //fill vector of decay mode names so it can be used by all functions
  void fillDecayModes();

  //find correct sister type if applicable
  bool goodSister(GenTauDecayID&, const std::vector<unsigned int>&) const;

  /*save this gen object
    - if primaryTauDecayType_ == GenTauDecayID::HAD, the reco::GenParticle object saved is the 
    status 3 tau from the decay of momPDGID_
    - if primaryTauDecayType_ == GenTauDecayID::MU or primaryTauDecayType_ == GenTauDecayID::E, 
    the reco::GenParticle object saved is the charged lepton daughter of the status 2 tau that is 
    itself the sole daughter of the status 3 tau from the decay of momPDGID_*/
  int keyToStore(const unsigned int, edm::Handle<reco::GenParticleCollection>&, 
		 GenTauDecayID::DecayType) const;

  // ----------member data ---------------------------

  //input tag for gen particle collection
  edm::InputTag genParticleTag_;

  //list of fabs(PDG ID) to match
  std::vector<unsigned int> absMatchPDGIDs_;

  //sister fabs(PDG ID) to match
  unsigned int sisterAbsMatchPDGID_;

  //set of parameters for GenTauDecayID class
  edm::ParameterSet genTauDecayIDPSet_;

  //decay type of one member of the tau pair, arbitrarily designated the primary
  GenTauDecayID::DecayType primaryTauDecayType_;

  //decay type of the second member of the tau pair, arbitrarily designated the sister
  GenTauDecayID::DecayType sisterTauDecayType_;

  //pT rank of the primary tau
  int primaryTauPTRank_;

  //hadronic decay type of the primary tau
  int primaryTauHadronicDecayType_;

  //hadronic decay type of the sister
  int sisterHadronicDecayType_;

  //max |eta| of the primary tau
  double primaryTauAbsEtaMax_;

  //min pT of the primary tau
  double primaryTauPTMin_;

  //flag indicating whether sister should be counted
  bool countSister_;

  //flag indicating whether pT cuts should be applied in determining valid gen objects
  bool applyPTCuts_;

  //flag indicating whether KShorts should be counted as neutral hadrons
  double countKShort_;

  //minimum number of gen objects passing cuts that must be found for event to pass filter
  unsigned int minNumGenObjectsToPassFilter_;

  //flag indicating whether all collections should be produced or just 1
  bool makeAllCollections_;

  //vector of decay mode names
  std::vector<std::string> decayModes_;
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
GenObjectProducer::GenObjectProducer(const edm::ParameterSet& iConfig) :
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  absMatchPDGIDs_(iConfig.getParameter<std::vector<unsigned int> >("absMatchPDGIDs")),
  sisterAbsMatchPDGID_(iConfig.getParameter<unsigned int>("sisterAbsMatchPDGID")),
  genTauDecayIDPSet_(iConfig.getParameter<edm::ParameterSet>("genTauDecayIDPSet")),
  primaryTauDecayType_(static_cast<GenTauDecayID::DecayType>
		       (iConfig.getParameter<unsigned int>("primaryTauDecayType"))),
  sisterTauDecayType_(static_cast<GenTauDecayID::DecayType>
		      (iConfig.getParameter<unsigned int>("sisterTauDecayType"))),
  primaryTauPTRank_(iConfig.getParameter<int>("primaryTauPTRank")),
  primaryTauHadronicDecayType_(static_cast<reco::PFTau::hadronicDecayMode>
			       (iConfig.getParameter<int>("primaryTauHadronicDecayType"))),
  sisterHadronicDecayType_(static_cast<reco::PFTau::hadronicDecayMode>
			   (iConfig.getParameter<int>("sisterHadronicDecayType"))),
  primaryTauAbsEtaMax_(iConfig.getParameter<double>("primaryTauAbsEtaMax")),
  primaryTauPTMin_(iConfig.getParameter<double>("primaryTauPTMin")),
  countSister_(iConfig.getParameter<bool>("countSister")),
  applyPTCuts_(iConfig.getParameter<bool>("applyPTCuts")),
  countKShort_(iConfig.getParameter<bool>("countKShort")),
  minNumGenObjectsToPassFilter_
  (iConfig.getParameter<unsigned int>("minNumGenObjectsToPassFilter")),
  makeAllCollections_(iConfig.getParameter<bool>("makeAllCollections"))
{
  //fill vector of decay mode names so it can be used by all functions
  fillDecayModes();

  //register your products
  if (makeAllCollections_) {
    for (std::vector<std::string>::const_iterator iMode = decayModes_.begin(); 
	 iMode != decayModes_.end(); ++iMode) {
      for (unsigned int iPTRank = 0; iPTRank <= 3; ++iPTRank) {
	std::stringstream instance;
	instance << "decayMode" << *iMode << "PTRank" << iPTRank;
	produces<reco::GenParticleRefVector>(instance.str());
      }
    }
  }
  else produces<reco::GenParticleRefVector>();

  //now do what ever other initialization is needed
  
}


GenObjectProducer::~GenObjectProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
bool GenObjectProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get GEN particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //fill STL container of all decay products of pseudoscalar Higgses
  std::vector<GenTauDecayID> aDecayProducts;
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {
    try {
      GenTauDecayID tauDecay(genTauDecayIDPSet_, pGenParticles, 
			     iGenParticle - pGenParticles->begin());
      if (tauDecay.isStatus3DecayProduct()) aDecayProducts.push_back(tauDecay);
    }
    catch (std::string& ex) { throw cms::Exception("GenObjectProducer") << ex; }
  }

  /*find the decay type of each a decay product if it's a tau (needed to get the visible 
    4-vector), otherwise use the status 3 4-vector*/
  for (std::vector<GenTauDecayID>::iterator iDecayProduct = aDecayProducts.begin(); 
       iDecayProduct != aDecayProducts.end(); ++iDecayProduct) {
    try { iDecayProduct->tauDecayType(applyPTCuts_, countKShort_); }
    catch (std::string& ex) { throw cms::Exception("GenObjectProducer") << ex; }
  }

  //sort the a decay products in ascending order by visible pT
  if (primaryTauPTRank_ != GenTauDecayID::ANY_PT_RANK) {
    Common::sortByPT(aDecayProducts);

    //set the pT rank of the a decay product (highest rank is 0, next highest is 1, etc.)
    for (std::vector<GenTauDecayID>::iterator iDecayProduct = aDecayProducts.begin(); 
	 iDecayProduct != aDecayProducts.end(); ++iDecayProduct) {
      iDecayProduct->setPTRank(aDecayProducts.end() - iDecayProduct - 1);
    }
  }

  //declare pointer to output collection to produce
  std::vector<std::auto_ptr<reco::GenParticleRefVector> > genObjs;
  for (std::vector<std::string>::const_iterator iMode = decayModes_.begin(); 
       iMode != decayModes_.end(); ++iMode) {
    for (unsigned int iPTRank = 0; iPTRank <= 3; ++iPTRank) {
      genObjs.push_back(std::auto_ptr<reco::GenParticleRefVector>(new reco::GenParticleRefVector));
    }
  }

  //loop over gen tau decays
  std::vector<unsigned int> keysToIgnore;
  for (std::vector<GenTauDecayID>::iterator iTau = aDecayProducts.begin(); 
       iTau != aDecayProducts.end(); ++iTau) {
    try {

      //select gen tau decays or light leptons...
      bool isStatus3 = false;
      std::vector<unsigned int>::const_iterator iPDGID = absMatchPDGIDs_.begin();
      while ((iPDGID != absMatchPDGIDs_.end()) && !isStatus3) {
	isStatus3 = iTau->isStatus3DecayProduct(*iPDGID);
	++iPDGID;
      }
      if (isStatus3) {
	const unsigned int tauKey = iTau->getTauIndex();

	//...with the right sister
	bool rightSister = goodSister(*iTau, keysToIgnore);
	if (rightSister) {

	  //get pT rank and decay type
	  unsigned int pTRank = 0;
	  if (primaryTauPTRank_ != GenTauDecayID::ANY_PT_RANK) pTRank = iTau->getPTRank();
	  const std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> decayType = 
	    iTau->tauDecayType(applyPTCuts_, countKShort_);

	  //choose taus in the allowed fiducial region
	  bool passAcceptanceCut = 
	    ((primaryTauAbsEtaMax_ == -1.0) || 
	     (fabs(iTau->getVisibleTauP4().Eta()) < primaryTauAbsEtaMax_)) && 
	    ((primaryTauPTMin_ == -1.0) || 
	     (fabs(iTau->getVisibleTauP4().Pt()) > primaryTauPTMin_));

	  //if making all collections, get hashed index to access proper collection
	  unsigned int hashedIndex = 0;
	  if (makeAllCollections_) {

	    //ignore unsupported hadronic decay modes
	    const bool ignore = 
	      (decayType.first != reco::PFTau::kOneProng0PiZero) && 
	      (decayType.first != reco::PFTau::kOneProng1PiZero) && 
	      (decayType.first != reco::PFTau::kOneProng2PiZero) && 
	      (decayType.first != reco::PFTau::kThreeProng0PiZero);

	    //count tau-->mu as decay mode index 0
	    int decayModeIndex = -1;
	    if (decayType.second == GenTauDecayID::MU) decayModeIndex = 0;

	    /*count tau-->had decays as indices 1-4, in order of appearance in 
	      reco::PFTau::hadronicDecayMode*/
	    if (decayType.second == GenTauDecayID::HAD) {
	      decayModeIndex = 1;
	      if (!ignore) decayModeIndex+=static_cast<int>(decayType.first);
	      const int maxDecayModeIndex = decayModes_.size() - 1;
	      if (decayModeIndex > maxDecayModeIndex) decayModeIndex = maxDecayModeIndex;
	    }

	    //make hashed index out of pT rank and decay type
	    hashedIndex = 
	      ((decayModeIndex == -1) ? genObjs.size() : (4*decayModeIndex + pTRank));
	  }

	  /*in the case where only 1 collection is made, verify that this tau has the right decay 
	    mode and pT rank*/
	  bool right = makeAllCollections_ || 
	    (((primaryTauPTRank_ == GenTauDecayID::ANY_PT_RANK) || 
	      ((int)pTRank == primaryTauPTRank_)) && 
	     ((primaryTauDecayType_ == TAU_ALL) || (decayType.second == primaryTauDecayType_)) && 
	     ((primaryTauHadronicDecayType_ == reco::PFTau::kNull) || /*using kNull to mean don't 
								       select on hadronic tau 
								       decay type*/
	      (decayType.first == primaryTauHadronicDecayType_)));

	  //get the key of the gen object to store
	  if (right && passAcceptanceCut) {
	    int genObjKey = keyToStore(tauKey, pGenParticles, decayType.second);

	    //store
	    if (genObjKey != -1) {
	      if (hashedIndex < genObjs.size()) {
		genObjs[hashedIndex]->push_back(reco::GenParticleRef(pGenParticles, genObjKey));
	      }
	    }
	    else throw cms::Exception("GenObjectProducer") << "Invalid gen object key.\n";

	    /*add this tau's key to the list of keys to ignore when we get to its sister if 
	      applicable*/
	    keysToIgnore.push_back(tauKey);
	  }
	}
      }
    }
    catch (std::string& ex) { throw cms::Exception("GenObjectProducer") << ex; }
  }

  //flag indicating whether enough gen-matched reco objects were found
  bool foundGenObject = genObjs[0]->size() >= minNumGenObjectsToPassFilter_;
  if (makeAllCollections_) {
    unsigned int iColl = 0;
    while ((iColl < genObjs.size()) && foundGenObject) {
      foundGenObject = genObjs[iColl]->size() >= minNumGenObjectsToPassFilter_;
      ++iColl;
    }
  }

  //put output collection into event
  if (makeAllCollections_) {
    for (std::vector<std::string>::const_iterator iMode = decayModes_.begin(); 
	 iMode != decayModes_.end(); ++iMode) {
      for (unsigned int iPTRank = 0; iPTRank <= 3; ++iPTRank) {
	std::stringstream instance;
	instance << "decayMode" << *iMode << "PTRank" << iPTRank;
	iEvent.put(genObjs[4*(iMode - decayModes_.begin()) + iPTRank], instance.str());
      }
    }
  }
  else iEvent.put(genObjs[0]); //this function frees the auto_ptr argument

  //stop processing if not enough gen-matched objects were found
  return foundGenObject;
}

// ------------ method called once each job just before starting event loop  ------------
void GenObjectProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GenObjectProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool GenObjectProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool GenObjectProducer::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool GenObjectProducer::beginLuminosityBlock(edm::LuminosityBlock&, 
							   edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool GenObjectProducer::endLuminosityBlock(edm::LuminosityBlock&, 
							 edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void GenObjectProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void GenObjectProducer::fillDecayModes()
{
  decayModes_.push_back("Mu");
  decayModes_.push_back("1Prong");
  decayModes_.push_back("1Prong1Pi0");
  decayModes_.push_back("1Prong2Pi0");
  decayModes_.push_back("3Prong");
}

bool GenObjectProducer::goodSister(GenTauDecayID& tau, 
				   const std::vector<unsigned int>& keysToIgnore) const
{
  bool goodSister = false;
  try {
    tau.findSister();
    const unsigned int iSister = tau.getSisterIndex();
    if (countSister_ || (std::find(keysToIgnore.begin(), keysToIgnore.end(), iSister) == 
			 keysToIgnore.end())) { /*keysToIgnore keeps track of whether you've 
						  looped over the other half of the boosted di-tau 
						  pair before*/
      std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> decayType = 
	tau.sisterDecayType(applyPTCuts_, countKShort_);
      if (((sisterTauDecayType_ == TAU_ALL) || (decayType.second == sisterTauDecayType_)) && 
	  ((sisterHadronicDecayType_ == reco::PFTau::kNull) || /*using kNull to mean don't select 
								 on hadronic tau decay type*/
	   (decayType.first == sisterHadronicDecayType_)) && 
	  ((sisterAbsMatchPDGID_ == GenTauDecayID::ANY_PDGID) || 
	   (fabs(reco::GenParticleRef(tau.getGenParticleHandle(), iSister)->pdgId()) == 
	    sisterAbsMatchPDGID_))) goodSister = true;
    }
  }
  catch (std::string& ex) {
    if (sisterAbsMatchPDGID_ != GenTauDecayID::ANY_PDGID) throw ex; /*assuming you wanted a 
								      sister, so if one isn't 
								      found an exception should be 
								      thrown*/
    //else assume particle really didn't have sister and you knew that, so just move on
  }
  return goodSister;
}

int GenObjectProducer::keyToStore(const unsigned int tauKey, 
				  edm::Handle<reco::GenParticleCollection>& pGenParticles, 
				  GenTauDecayID::DecayType primaryTauDecayType) const
{
  int genObjKey = -1;
  reco::GenParticleRef tauRef(pGenParticles, tauKey);
  if ((primaryTauDecayType == GenTauDecayID::HAD) || 
      (fabs(tauRef->pdgId()) != GenTauDecayID::TAUPDGID)) genObjKey = tauKey;
  if ((primaryTauDecayType == GenTauDecayID::MU) || 
      (primaryTauDecayType == GenTauDecayID::E)) {
    int chargedLeptonDaughterKey = -1;
    reco::GenParticleRef status2DaughterRef = tauRef->daughterRef(0);
    unsigned int iDaughter = 0;
    while ((iDaughter < status2DaughterRef->numberOfDaughters()) && 
	   (chargedLeptonDaughterKey == -1)) {
      reco::GenParticleRef decayProductRef = status2DaughterRef->daughterRef(iDaughter);
      const unsigned int daughterPDGID = fabs(decayProductRef->pdgId());
      if ((daughterPDGID == GenTauDecayID::EPDGID) || 
	  (daughterPDGID == GenTauDecayID::MUPDGID)) {
	chargedLeptonDaughterKey = decayProductRef.key();
      }
      ++iDaughter;
    }
    genObjKey = chargedLeptonDaughterKey;
  }
  return genObjKey;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenObjectProducer);
