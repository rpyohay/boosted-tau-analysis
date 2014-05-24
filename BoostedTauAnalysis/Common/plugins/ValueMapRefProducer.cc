// -*- C++ -*-
//
// Package:    Common
// Class:      ValueMapRefProducer
// 
/**\class ValueMapRefProducer ValueMapRefProducer.cc 
   BoostedTauAnalysis/Common/plugins/ValueMapRefProducer.cc

 Description: create a collection of tau refs from a value map

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: ValueMapRefProducer.cc,v 1.6 2012/12/12 16:02:05 yohay Exp $
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
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TauReco/interface/PFTau.h"

//
// class declaration
//

template<class T, class U>
class ValueMapRefProducer : public edm::EDProducer {
public:
  explicit ValueMapRefProducer(const edm::ParameterSet&);
  ~ValueMapRefProducer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------

  //input tag for value map
  edm::InputTag valMapTag_;

  //input tag for the keys with mapped values
  edm::InputTag keyTag_;
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
template<class T, class U>
ValueMapRefProducer<T, U>::ValueMapRefProducer(const edm::ParameterSet& iConfig) :
  valMapTag_(iConfig.getParameter<edm::InputTag>("valMapTag")),
  keyTag_(iConfig.getParameter<edm::InputTag>("keyTag"))
{
  produces<edm::RefVector<std::vector<T> > >();
}

template<class T, class U>
ValueMapRefProducer<T, U>::~ValueMapRefProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
template<class T, class U>
void ValueMapRefProducer<T, U>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<edm::RefVector<std::vector<T> > > coll(new edm::RefVector<std::vector<T> >);

  //get value map
  edm::Handle<edm::ValueMap<edm::Ref<std::vector<T> > > > pValMap;
  iEvent.getByLabel(valMapTag_, pValMap);

  //get keys with mapped values
  edm::Handle<edm::RefVector<std::vector<U> > > pKeys;
  iEvent.getByLabel(keyTag_, pKeys);

  //loop over keys
  for (typename edm::RefVector<std::vector<U> >::const_iterator iKey = pKeys->begin(); 
       iKey != pKeys->end(); ++iKey) {

    //store mapped value, which is an object of type T
    coll->push_back((*pValMap)[*iKey]);
  }

  //put the collection in the event
  iEvent.put(coll);
}

// ------------ method called once each job just before starting event loop  ------------
template<class T, class U>
void ValueMapRefProducer<T, U>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T, class U>
void ValueMapRefProducer<T, U>::endJob() {
}

// ------------ method called when starting to processes a run  ------------
template<class T, class U>
void ValueMapRefProducer<T, U>::beginRun(edm::Run&, edm::EventSetup const&)
{ 
}

// ------------ method called when ending the processing of a run  ------------
template<class T, class U>
void ValueMapRefProducer<T, U>::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T, class U>
void ValueMapRefProducer<T, U>::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T, class U>
void ValueMapRefProducer<T, U>::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T, class U>
void ValueMapRefProducer<T, U>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef ValueMapRefProducer<reco::PFTau, reco::PFTau> ValueMapTauRefProducer;
DEFINE_FWK_MODULE(ValueMapTauRefProducer);
