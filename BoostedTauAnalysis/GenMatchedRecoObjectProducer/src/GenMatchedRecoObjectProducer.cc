// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      GenMatchedRecoObjectProducer
// 
/**\class GenMatchedRecoObjectProducer GenMatchedRecoObjectProducer.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/GenMatchedRecoObjectProducer.cc

Description: produce a collection of reco objects matched to gen boosted di-tau objects

Implementation:

*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Thu Aug 23 11:23:58 CEST 2012
// $Id: GenMatchedRecoObjectProducer.cc,v 1.3 2012/09/25 11:44:34 yohay Exp $
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
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"

//
// class declaration
//

template<class T>
class GenMatchedRecoObjectProducer : public edm::EDFilter {
public:
  explicit GenMatchedRecoObjectProducer(const edm::ParameterSet&);
  ~GenMatchedRecoObjectProducer();

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

  //input tag for base gen particle collection
  edm::InputTag genParticleTag_;

  /*input tag for gen particle collection to match
    count on the user to pass in a collection that will not lead to the same reco object being 
    matched to multiple different gen objects
    for example, if the input object is a boosted di-tau pair, only 1 member of the pair should be 
    in the input collection*/
  edm::InputTag selectedGenParticleTag_;

  //input tag for reco object collection
  edm::InputTag recoObjTag_;

  //input tag for base reco object collection
  edm::InputTag baseRecoObjTag_;

  //set of parameters for GenTauDecayID class
  edm::ParameterSet genTauDecayIDPSet_;

  //flag indicating whether pT cuts should be applied in determining valid gen objects
  bool applyPTCuts_;

  //flag indicating whether KShorts should be counted as neutral hadrons
  bool countKShort_;

  //pT rank of the matched reco object in the event
  int pTRank_;

  //flag indicating whether all collections should be produced or just 1
  bool makeAllCollections_;

  //flag indicating whether pT rank should be assessed against reco object or matching gen object
  bool useGenObjPTRank_;

  //number of output collections in the makeAllCollections_ = true case
  unsigned int nOutputColls_;

  //dR matching cut
  double dR_;

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
template<class T>
GenMatchedRecoObjectProducer<T>::GenMatchedRecoObjectProducer(const edm::ParameterSet& iConfig) :
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  selectedGenParticleTag_(iConfig.getParameter<edm::InputTag>("selectedGenParticleTag")),
  recoObjTag_(iConfig.getParameter<edm::InputTag>("recoObjTag")),
  baseRecoObjTag_(iConfig.getParameter<edm::InputTag>("baseRecoObjTag")),
  genTauDecayIDPSet_(iConfig.getParameter<edm::ParameterSet>("genTauDecayIDPSet")),
  applyPTCuts_(iConfig.getParameter<bool>("applyPTCuts")),
  countKShort_(iConfig.getParameter<bool>("countKShort")),
  pTRank_(iConfig.getParameter<int>("pTRank")),
  makeAllCollections_(iConfig.getParameter<bool>("makeAllCollections")),
  useGenObjPTRank_(iConfig.getParameter<bool>("useGenObjPTRank")),
  nOutputColls_(iConfig.getParameter<unsigned int>("nOutputColls")),
  dR_(iConfig.getParameter<double>("dR")),
  minNumGenObjectsToPassFilter_
  (iConfig.getParameter<unsigned int>("minNumGenObjectsToPassFilter"))
{
  //register your products
  if (makeAllCollections_) {
    for (unsigned int i = 0; i < nOutputColls_; ++i) {
      std::stringstream instance;
      instance << "coll" << i;
      produces<edm::RefVector<std::vector<T> > >(instance.str());
    }
  }
  else produces<edm::RefVector<std::vector<T> > >();

  //now do what ever other initialization is needed
  
}


template<class T>
GenMatchedRecoObjectProducer<T>::~GenMatchedRecoObjectProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get base gen particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //get selected gen particles
  edm::Handle<reco::GenParticleRefVector> pSelectedGenParticles;
  iEvent.getByLabel(selectedGenParticleTag_, pSelectedGenParticles);

  //get reco object collection
  edm::Handle<edm::RefVector<std::vector<T> > > pRecoObjs;
  iEvent.getByLabel(recoObjTag_, pRecoObjs);

  //get base reco object collection
  edm::Handle<std::vector<T> > pBaseRecoObjs;
  iEvent.getByLabel(baseRecoObjTag_, pBaseRecoObjs);

  //fill STL container of pointers to reco objects
  std::vector<T*> recoObjPtrs;
  for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj = pRecoObjs->begin(); 
       iRecoObj != pRecoObjs->end(); ++iRecoObj) {
    recoObjPtrs.push_back(const_cast<T*>(iRecoObj->get()));
  }

  //make a copy of the reco object vector
  std::vector<T*> recoObjPtrsCopy = recoObjPtrs;

  //sort the reco objects in the copied vector by ascending order by pT in the copied vector
  Common::sortByPT(recoObjPtrsCopy);

  //fill STL container of selected gen objects
  std::vector<GenTauDecayID> selectedGenObjs;
  for (unsigned int iGenParticle = 0; iGenParticle < pSelectedGenParticles->size(); 
       ++iGenParticle) {
    try {

      //status 2 gen particles in the selected gen particle collection are not supported
      GenTauDecayID tauDecay(genTauDecayIDPSet_, pGenParticles, 
			     Common::getStatus3Key(pSelectedGenParticles, pGenParticles, 
						   iGenParticle));
      if (tauDecay.isStatus3DecayProduct()) selectedGenObjs.push_back(tauDecay);
    }
    catch (std::string& ex) { throw cms::Exception("GenMatchedRecoObjectProducer") << ex; }
  }

  /*find the decay type of each gen object if it's a tau (needed to get the visible 
    4-vector), otherwise use the status 3 4-vector*/
  for (std::vector<GenTauDecayID>::iterator iGenObj = selectedGenObjs.begin(); 
       iGenObj != selectedGenObjs.end(); ++iGenObj) {
    try { iGenObj->tauDecayType(applyPTCuts_, countKShort_); }
    catch (std::string& ex) { throw cms::Exception("GenObjectProducer") << ex; }
  }

  //sort the gen objects in ascending order by visible pT
  Common::sortByPT(selectedGenObjs);

  //set the pT rank of the a decay product (highest rank is 0, next highest is 1, etc.)
  for (std::vector<GenTauDecayID>::iterator iGenObj = selectedGenObjs.begin(); 
       iGenObj != selectedGenObjs.end(); ++iGenObj) {
    iGenObj->setPTRank(selectedGenObjs.end() - iGenObj - 1);
  }

  //declare pointers to output collection to produce
  std::vector<std::auto_ptr<edm::RefVector<std::vector<T> > > > genMatchedRecoObjs;
  for (unsigned int i = 0; i < nOutputColls_; ++i) {
    genMatchedRecoObjs.push_back(std::auto_ptr<edm::RefVector<std::vector<T> > >
				 (new edm::RefVector<std::vector<T> >));
  }

  //debug
  std::vector<edm::Ref<std::vector<T> > > recoObjsToSave;

  //loop over selected gen particles
  for (std::vector<GenTauDecayID>::iterator iGenObj = selectedGenObjs.begin(); 
       iGenObj != selectedGenObjs.end(); ++iGenObj) {

    //make a dummy LeafCandidate and ref out of the visible 4-vector
    reco::LeafCandidate::LorentzVector visibleGenP4 = iGenObj->getVisibleTauP4();
    std::vector<reco::LeafCandidate> 
      visibleGenParticle(1, reco::LeafCandidate(0.0, visibleGenP4));
    edm::Ref<std::vector<reco::LeafCandidate> > visibleGenParticleRef(&visibleGenParticle, 0);

    //find the nearest reco object to the gen particle
    int nearestRecoObjPTRank = -1; /*this is the index into recoObjPtrsCopy of the nearest object
				     since recoObjPtrsCopy is sorted in ascending order by pT, the 
				     pT rank of the nearest object is recoObjPtrsCopy.size() - 
				     nearestRecoObjPTRank - 1*/
    const T* nearestRecoObj = 
      Common::nearestObject(visibleGenParticleRef, recoObjPtrsCopy, nearestRecoObjPTRank);
    if (nearestRecoObjPTRank != -1) {
      nearestRecoObjPTRank = recoObjPtrsCopy.size() - nearestRecoObjPTRank - 1;
    }

    /*still need the index into the original collection pRecoObjs of the nearest object, so repeat 
      call to nearestObject, but with original recoObjPtrs vector
      ideally when sorting the vector in the first place we'd save the original object keys*/
    int nearestRecoObjKey = -1;
    nearestRecoObj = Common::nearestObject(visibleGenParticleRef, recoObjPtrs, nearestRecoObjKey);

    //if nearest reco object is within dR_ of the gen object... 
//     bool save = false;
    if ((nearestRecoObj != NULL) && (nearestRecoObjPTRank >= 0) && 
	(reco::deltaR(*nearestRecoObj, *visibleGenParticleRef) < dR_)) {
      int matchedGenObjPTRank = (int)iGenObj->getPTRank();

      /*...and in the case of makeAllCollections_ = true has a pT rank higher than or equal to the 
	max supported, ...*/
      if (makeAllCollections_) {

// 	if ((useGenObjPTRank_ && /*debug*//*(matchedGenObjPTRank < nOutputColls_)*/true) || 
// 	    (!useGenObjPTRank_ && 
// 	     /*debug*//*(nearestRecoObjPTRank < (int)nOutputColls_)*/true)) save = true;

	//debug
	/*nearestRecoObjKey is the index into recoObjPtrs of the matched object
	  recoObjPtrs is in the same order as pRecoObjs
	  the element of pRecoObjs with index nearestRecoObjKey is the ref of the matched object
	  its key MUST BE the index into the original collection pBaseRecoObjs (or this fails)
	  so, a ref to the original collection is saved*/
	recoObjsToSave.
	  push_back(edm::Ref<std::vector<T> >(pBaseRecoObjs, 
					      pRecoObjs->at(nearestRecoObjKey).key()));
      }

      //or in the case of makeAllCollections_ = false, the pTRank is the one wanted, ...
      else if ((pTRank_ == GenTauDecayID::ANY_PT_RANK) || 
	       (useGenObjPTRank_ && (matchedGenObjPTRank == pTRank_)) || 
	       (!useGenObjPTRank_ && (nearestRecoObjPTRank == pTRank_))) {
// 	save = true;
	nearestRecoObjPTRank = 0; /*since only 1 output collection is produced in this case, the 
				    index into genMatchedRecoObjs should be 0*/

	//debug
	genMatchedRecoObjs[nearestRecoObjPTRank]->
	  push_back(edm::Ref<std::vector<T> >(pBaseRecoObjs, 
					      pRecoObjs->at(nearestRecoObjKey).key()));
      }
    }

//     //comment out following when debugging
//     //...save this reco object
//     if (save) {
//       genMatchedRecoObjs[nearestRecoObjPTRank]->
// 	push_back(edm::Ref<std::vector<T> >(pBaseRecoObjs, 
// 					    pRecoObjs->at(nearestRecoObjKey).key()));
//     }
  }

  /*debug - here we pay no attention to the pT rank of the reco object among other reco object or 
    the gen object among other gen objects, and just sort the passing reco objects by pT to 
    determine their pT rank
    for instance, if 2 reco objects passed and the leading one had pT rank 0 by the normal measure 
    and the trailing one had pT rank 10 by the normal measure, the trailing one would be saved 
    here as pT rank 1, instead of normally being discarded*/
  Common::sortByPT(recoObjsToSave);
  if (makeAllCollections_) {
    for (unsigned int i = 0; i < nOutputColls_; ++i) {
      if ((recoObjsToSave.size() - i - 1) < genMatchedRecoObjs.size()) {
	genMatchedRecoObjs[i]->push_back(edm::Ref<std::vector<T> >
					 (pBaseRecoObjs, 
					  recoObjsToSave[recoObjsToSave.size() - i - 1].key()));
      }
    }
  }

  //flag indicating whether right number of gen-matched reco objects were found
  bool foundGenMatchedRecoObject = genMatchedRecoObjs[0]->size() >= minNumGenObjectsToPassFilter_;
  if (makeAllCollections_) {
    unsigned int iColl = 0;
    while ((iColl < nOutputColls_) && foundGenMatchedRecoObject) {
      foundGenMatchedRecoObject = 
	genMatchedRecoObjs[iColl]->size() >= minNumGenObjectsToPassFilter_;
      ++iColl;
    }
  }

  //put output collection into event
  if (makeAllCollections_) {
    for (unsigned int i = 0; i < nOutputColls_; ++i) {
      std::stringstream instance;
      instance << "coll" << i;
      iEvent.put(genMatchedRecoObjs[i], instance.str());
    }
  }
  else iEvent.put(genMatchedRecoObjs[0]); //this function frees the auto_ptr argument

  //stop processing if no gen-matched objects were found
  return foundGenMatchedRecoObject;
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void GenMatchedRecoObjectProducer<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void GenMatchedRecoObjectProducer<T>::endJob() {
}

// ------------ method called when starting to processes a run  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::beginLuminosityBlock(edm::LuminosityBlock&, 
							   edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::endLuminosityBlock(edm::LuminosityBlock&, 
							 edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
template<class T>
void 
GenMatchedRecoObjectProducer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
typedef GenMatchedRecoObjectProducer<reco::Muon> GenMatchedMuonProducer;
typedef GenMatchedRecoObjectProducer<reco::PFJet> GenMatchedJetProducer;
typedef GenMatchedRecoObjectProducer<reco::PFTau> GenMatchedTauProducer;
DEFINE_FWK_MODULE(GenMatchedMuonProducer);
DEFINE_FWK_MODULE(GenMatchedJetProducer);
DEFINE_FWK_MODULE(GenMatchedTauProducer);
