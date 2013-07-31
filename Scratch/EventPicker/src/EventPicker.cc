// -*- C++ -*-
//
// Package:    EventPicker
// Class:      EventPicker
// 
/**\class EventPicker EventPicker.cc Scratch/EventPicker/src/EventPicker.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Oct 10 13:29:07 CEST 2012
// $Id: EventPicker.cc,v 1.1 2012/10/24 14:27:32 yohay Exp $
//
//


// system include files
#include <memory>
// #include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class EventPicker : public edm::EDFilter {
public:
  explicit EventPicker(const edm::ParameterSet&);
  ~EventPicker();

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

  //vector of events to pick
  std::vector<edm::EventRange> evts_;

  //vector of run numbers for wanted events
  std::vector<edm::RunNumber_t> runNumbers_;

  //vector of event numbers for wanted events
  std::vector<edm::EventNumber_t> evtNumbers_;
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
EventPicker::EventPicker(const edm::ParameterSet& iConfig) :
  evts_(iConfig.getParameter<std::vector<edm::EventRange> >("evts"))
{
  //now do what ever initialization is needed

  /*fill vectors of run, lumi section, and event numbers, assuming each range contains exactly 1 
    event*/
  for (std::vector<edm::EventRange>::const_iterator iEvt = evts_.begin(); iEvt != evts_.end(); 
       ++iEvt) {
    runNumbers_.push_back(iEvt->startRun());
    evtNumbers_.push_back(iEvt->startEvent());
  }
}


EventPicker::~EventPicker()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
EventPicker::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //if event is in the specified list, return true
  bool ret = false;
  std::vector<edm::RunNumber_t>::const_iterator iRun = 
    std::find(runNumbers_.begin(), runNumbers_.end(), iEvent.run());
  if ((iRun != runNumbers_.end()) && 
      (evtNumbers_[iRun - runNumbers_.begin()] == iEvent.id().event())) {
    ret = true;

    //debug
    std::cerr << "Saving event " << iEvent.run() << ":" << iEvent.id().event() << std::endl;

    //speed this up by erasing the found elements?
  }
  return ret;
}

// ------------ method called once each job just before starting event loop  ------------
void 
EventPicker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventPicker::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
EventPicker::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
EventPicker::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
EventPicker::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
EventPicker::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void
EventPicker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(EventPicker);
