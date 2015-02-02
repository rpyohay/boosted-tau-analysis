// -*- C++ -*-
//
// Package:    Common
// Class:      HighestPTRefSelector
// 
/**\class HighestPTRefSelector HighestPTRefSelector.cc 
   BoostedTauAnalysis/Common/plugins/HighestPTRefSelector.cc

 Description: create a collection consisting of a ref to the highest pT object to put in the event

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: HighestPTRefSelector.cc,v 1.6 2012/12/12 16:02:05 yohay Exp $
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

//
// class declaration
//

template<class T>
class HighestPTRefSelector : public edm::EDFilter {
public:
  explicit HighestPTRefSelector(const edm::ParameterSet&);
  ~HighestPTRefSelector();
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

  //input tag for object ref collection
  edm::InputTag objRefTag_;
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
HighestPTRefSelector<T>::HighestPTRefSelector(const edm::ParameterSet& iConfig) :
  objRefTag_(iConfig.getParameter<edm::InputTag>("objRefTag"))
{
  produces<edm::RefVector<std::vector<T> > >();
}

template<class T>
HighestPTRefSelector<T>::~HighestPTRefSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
template<class T>
bool HighestPTRefSelector<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<edm::RefVector<std::vector<T> > > objColl(new edm::RefVector<std::vector<T> >);

  //get objects
  edm::Handle<edm::RefVector<std::vector<T> > > pObjs;
  iEvent.getByLabel(objRefTag_, pObjs);

  //sort objects by descending order in pT
  std::vector<edm::Ref<std::vector<T> > > pTSortedObjs;
  for (typename edm::RefVector<std::vector<T> >::const_iterator iObj = pObjs->begin(); 
       iObj != pObjs->end(); ++iObj) { pTSortedObjs.push_back(*iObj); }
  Common::sortByPT(pTSortedObjs);
  std::reverse(pTSortedObjs.begin(), pTSortedObjs.end());

  //put the highest pT object in the event
  bool objExists = (pTSortedObjs.size() > 0);
  if (objExists) objColl->push_back(pTSortedObjs[0]);
  iEvent.put(objColl);
  return objExists;
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void HighestPTRefSelector<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void HighestPTRefSelector<T>::endJob() {
}

// ------------ method called when starting to processes a run  ------------
template<class T>
bool HighestPTRefSelector<T>::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
bool HighestPTRefSelector<T>::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
bool HighestPTRefSelector<T>::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
bool HighestPTRefSelector<T>::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
template<class T>
void HighestPTRefSelector<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef HighestPTRefSelector<reco::Muon> HighestPTMuonRefSelector;
DEFINE_FWK_MODULE(HighestPTMuonRefSelector);
