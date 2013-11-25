// -*- C++ -*-
//
// Package:    TriggerObjectFilter
// Class:      TriggerObjectFilter
// 
/**\class TriggerObjectFilter TriggerObjectFilter.cc BoostedTauAnalysis/TriggerObjectFilter/src/TriggerObjectFilter.cc

 Description: matching the W_mu to the trigger muon

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//         Created:  Mon Aug 19 11:19:05 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <cmath>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
#include "FWCore/Common/interface/TriggerNames.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;

//
// class declaration
//

class TriggerObjectFilter : public edm::EDFilter {
   public:
      explicit TriggerObjectFilter(const edm::ParameterSet&);
      ~TriggerObjectFilter();

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

  edm::InputTag WMuonTag_;
  edm::InputTag triggerEventTag_;
  edm::InputTag triggerResultsTag_;
  double delRMatchingCut_;
  std::vector<edm::InputTag> hltTags_;
  HLTConfigProvider hltConfig_;
  edm::InputTag theRightHLTTag_;
  edm::InputTag theRightHLTSubFilter_;
  std::vector<edm::InputTag> HLTSubFilters_;
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
TriggerObjectFilter::TriggerObjectFilter(const edm::ParameterSet& iConfig):hltConfig_()
{
   //now do what ever initialization is needed
  WMuonTag_ = iConfig.getParameter<edm::InputTag>("WMuonTag");
  const edm::InputTag dTriggerEventTag("hltTriggerSummaryAOD","","HLT");
  triggerEventTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag",dTriggerEventTag);
  const edm::InputTag dTriggerResults("TriggerResults","","HLT");
  // By default, trigger results are labeled "TriggerResults" with process name "HLT" in the event.
  triggerResultsTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag",dTriggerResults);
  delRMatchingCut_ = iConfig.getUntrackedParameter<double>("triggerDelRMatch", 0.30);
  hltTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("hltTags");
  //  hltConfig_ = iConfig.getParameter<HLTConfigProvider>("hltConfig");
  theRightHLTTag_ = iConfig.getParameter<edm::InputTag>("theRightHLTTag");
  theRightHLTSubFilter_ = iConfig.getParameter<edm::InputTag>("theRightHLTSubFilter");
  //Whether using HLT trigger path name or the actual trigger filter name. Trigger path is default.
  HLTSubFilters_ = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("HLTSubFilters",std::vector<edm::InputTag>());

}


TriggerObjectFilter::~TriggerObjectFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TriggerObjectFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool trigger_matched = false;
  int index = 9999;

  //get W muons
  edm::Handle<reco::MuonRefVector> WMuons;
  iEvent.getByLabel(WMuonTag_, WMuons);

   // Trigger Info
  edm::Handle<trigger::TriggerEvent> trgEvent;
  iEvent.getByLabel(triggerEventTag_,trgEvent);
  edm::Handle<edm::TriggerResults> pTrgResults;
  iEvent.getByLabel(triggerResultsTag_, pTrgResults);
  std::map<std::string, bool> triggerInMenu;
  std::string myHLTFilter = "";

  double Wmu_pt = -9999.;
  double Wmu_phi = 0.;
  double Wmu_eta = -9999.;

  if (WMuons->size() == 1)
    {
      Wmu_eta = WMuons->at(0)->eta();
      Wmu_phi = WMuons->at(0)->phi();
    }
  else
    { // if >1 W muons
      for (reco::MuonRefVector::const_iterator iMuon = WMuons->begin(); iMuon != WMuons->end(); ++iMuon)
        { // loop to find highest-pT
          if ((*iMuon)->pt() > Wmu_pt)
            {
              Wmu_pt = (*iMuon)->pt();
	      Wmu_eta = (*iMuon)->eta();
              Wmu_phi = (*iMuon)->phi();
            }
        } // loop to find highest-pT
    } // if >1 W muons


  // get names of active HLT paths in this event
  std::vector<std::string> activeHLTPathsInThisEvent = hltConfig_.triggerNames();
  // loop over active HLT paths to search for desired path
  for (std::vector<std::string>::const_iterator iHLT = activeHLTPathsInThisEvent.begin(); 
       iHLT != activeHLTPathsInThisEvent.end(); ++iHLT) { // active paths loop
    for (std::vector<edm::InputTag>::const_iterator iMyHLT = hltTags_.begin(); 
	 iMyHLT != hltTags_.end(); ++iMyHLT) {
      if ((*iMyHLT).label() == *iHLT) {
	cout << "######## " << *iHLT << endl;
	myHLTFilter = (*iMyHLT).label();
	triggerInMenu[(*iMyHLT).label()] = true;
	std::cout << "(*iMyHLT).label() = " << (*iMyHLT).label() << std::endl;
	std::cout << "hltConfig_.prescaleValue(iEvent, iSetup, *iHLT) = ";
	std::cout << hltConfig_.prescaleValue(iEvent, iSetup, *iHLT) << std::endl;
      }
    }
  } // active paths loop
  
   edm::InputTag filterTag;
   // loop over these objects to see whether they match
   const trigger::TriggerObjectCollection& TOC( trgEvent->getObjects() );

   //choose the right sub-filter depending on the HLT path name
   std::vector<std::string> filters;
   try { filters = hltConfig_.moduleLabels( theRightHLTTag_.label() ); }
   catch (std::exception ex) { cout << "bad trigger\n"; }
   //try { filters = hltConfig_.moduleLabels( myHLTFilter ); }
   //catch (std::exception ex) { cout << "bad trigger 2\n"; }

   //loop over moduleLabels of hltConfig
   //loop over sub-filters and see if they match any moduleLabel
   std::vector<std::string>::const_iterator filter = filters.begin();
   bool foundMatch = false;
   /*   while ((filter != filters.end()) && !foundMatch) {
     std::vector<edm::InputTag>::const_iterator iHLTSubFilter = HLTSubFilters_.begin();
     while ((iHLTSubFilter != HLTSubFilters_.end()) && !foundMatch) {
       if (*filter == iHLTSubFilter->label()) {
	 foundMatch = true;
	 theRightHLTSubFilter_ = *iHLTSubFilter;
       }
       ++iHLTSubFilter;
     }
     ++filter;
     }*/

   //loop over filterTags of the trgEvent
   //store the position of the one that matches the right sub-filter
   cout << "theRightSubFilter label = " << theRightHLTSubFilter_.label() << endl;
   for(int i=0; i != trgEvent->sizeFilters(); ++i) {
     std::string label(trgEvent->filterTag(i).label());
     //if( label == theRightHLTSubFilter_.label() ) index = i;
     if( label.find(theRightHLTSubFilter_.label()) != std::string::npos )
       {
	 std::cout << "filterTag label: " << label << std::endl;
	 cout << "found subfilter!" << endl;
	 index = i;
       }
   }
   cout << "index = " << index << endl;
   // find how many objects there are
   if (index == 9999)
     index = 0;
   const trigger::Keys& KEYS(trgEvent->filterKeys(index));
   const size_type nK(KEYS.size());

   //did this event fire the HLT?
   const edm::TriggerNames &trgNames = iEvent.triggerNames(*pTrgResults);
   /* for (unsigned int i = 0; i < trgNames.size(); i++)
      { cout << "trgName = " << trgNames.triggerName(i) << endl; } */
   //const unsigned int trgIndex = trgNames.triggerIndex(theRightHLTTag_.label());
   const unsigned int trgIndex = trgNames.triggerIndex(myHLTFilter);
   bool firedHLT = (trgIndex < trgNames.size()) && (pTrgResults->accept(trgIndex));
   cout << "trgIndex = " << trgIndex << " and trgNames size = " << trgNames.size() << endl;
   cout << "firedHLT = " << firedHLT << endl;

   // If the event fired the HLT,
   // loop over the trigger object collection 
   if (firedHLT)
     { // firedHLT
       // Get cut decision for each candidate
       // Did this candidate cause an HLT trigger?
       bool hltTrigger = false;
       for(int ipart = 0; ipart != nK; ++ipart) { 
	 
	 const trigger::TriggerObject& TO = TOC[KEYS[ipart]];	
	 double dRval = deltaR(Wmu_eta, Wmu_phi, TO.eta(), TO.phi());
	 cout << "dRval = " << dRval << endl;
	 hltTrigger = dRval < delRMatchingCut_;
	 if( hltTrigger )
	   {
	     trigger_matched = true;
	     break;
	   }
       }        
     } //firedHLT

  return trigger_matched;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TriggerObjectFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerObjectFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
TriggerObjectFilter::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup)
{ 
  bool changed_ = true;
  if ( !hltConfig_.init(iRun,iSetup,hltTags_[0].process(),changed_) ){
    edm::LogError("TriggerObjectFilter") << 
      "Error! Can't initialize HLTConfigProvider";
    throw cms::Exception("HLTConfigProvider::init() returned non 0");
  }
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
TriggerObjectFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
TriggerObjectFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
TriggerObjectFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerObjectFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TriggerObjectFilter);
