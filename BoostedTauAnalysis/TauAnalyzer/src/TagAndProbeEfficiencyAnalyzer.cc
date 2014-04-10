// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      TagAndProbeEfficiencyAnalyzer
// 
/**\class TagAndProbeEfficiencyAnalyzer TagAndProbeEfficiencyAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/TagAndProbeEfficiencyAnalyzer.cc

   Description: make efficiency plots for tag and probe analysis

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: TagAndProbeEfficiencyAnalyzer.cc,v 1.2 2012/09/25 11:49:57 yohay Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

//
// class declaration
//

template<class T>
class TagAndProbeEfficiencyAnalyzer : public edm::EDAnalyzer {
public:
  explicit TagAndProbeEfficiencyAnalyzer(const edm::ParameterSet&);
  ~TagAndProbeEfficiencyAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  //delete memory
  void reset(const bool);

  //fill pT histogram with object of same type as template parameter
  void fillPTHistogram(edm::Handle<edm::View<T> >&, TH1F*);

  //fill eta histogram
  void fillEtaHistogram(edm::Handle<edm::View<T> >&, TH1F*);

  // ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //tag-tag input tag
  edm::InputTag tagTagTag_;

  //tag-pass input tag
  edm::InputTag tagPassTag_;

  //tag-fail input tag
  edm::InputTag tagFailTag_;

  //marker colors for histograms with different pT rank
  std::vector<unsigned int> pTRankColors_;

  //marker colors for histograms with different decay mode
  std::vector<unsigned int> decayModeColors_;

  //marker styles for histograms with different pT rank
  std::vector<unsigned int> pTRankStyles_;

  //marker styles for histograms with different decay mode
  std::vector<unsigned int> decayModeStyles_;

  //legend entries for histograms with different pT rank
  std::vector<std::string> pTRankEntries_;

  //legend entries for histograms with different decay mode
  std::vector<std::string> decayModeEntries_;

  //histogram of denominator pT
  TH1F* denominatorPT_;

  //histogram of numerator pT
  TH1F* numeratorPT_;

  //histogram of denominator eta
  TH1F* denominatorEta_;

  //histogram of numerator eta
  TH1F* numeratorEta_;
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
TagAndProbeEfficiencyAnalyzer<T>::TagAndProbeEfficiencyAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  tagTagTag_(iConfig.getParameter<edm::InputTag>("tagTagTag")),
  tagFailTag_(iConfig.getParameter<edm::InputTag>("tagFailTag")),
  pTRankColors_(iConfig.getParameter<std::vector<unsigned int> >("pTRankColors")),
  decayModeColors_(iConfig.getParameter<std::vector<unsigned int> >("decayModeColors")),
  pTRankStyles_(iConfig.getParameter<std::vector<unsigned int> >("pTRankStyles")),
  decayModeStyles_(iConfig.getParameter<std::vector<unsigned int> >("decayModeStyles")),
  pTRankEntries_(iConfig.getParameter<std::vector<std::string> >("pTRankEntries")),
  decayModeEntries_(iConfig.getParameter<std::vector<std::string> >("decayModeEntries"))
{
  //now do what ever initialization is needed
  if (iConfig.existsAs<edm::InputTag>("tagPassTag")) {
    tagPassTag_ = iConfig.getParameter<edm::InputTag>("tagPassTag");
  }
  reset(false);
}

template<class T>
TagAndProbeEfficiencyAnalyzer<T>::~TagAndProbeEfficiencyAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
template<class T>
void TagAndProbeEfficiencyAnalyzer<T>::analyze(const edm::Event& iEvent, 
					       const edm::EventSetup& iSetup)
{
  //get tag-tag collection
  edm::Handle<edm::View<T> > pTagTagView;
  iEvent.getByLabel(tagTagTag_, pTagTagView);

  //get tag-pass collection
  edm::Handle<edm::View<T> > pTagPassView;
  if (tagPassTag_ == edm::InputTag()) {}
  else iEvent.getByLabel(tagPassTag_, pTagPassView);

  //get tag-fail collection
  edm::Handle<edm::View<T> > pTagFailView;
  iEvent.getByLabel(tagFailTag_, pTagFailView);

  /*fill numerator and denominator distributions with tag-tag events (fill twice: once for first 
    object, once for second object ==> pTagTagView->size() == 2)*/
  const unsigned int tagTagCollSize = pTagTagView->size();
  if (tagTagCollSize == 2) {
    fillPTHistogram(pTagTagView, denominatorPT_);
    fillPTHistogram(pTagTagView, numeratorPT_);
    fillEtaHistogram(pTagTagView, denominatorEta_);
    fillEtaHistogram(pTagTagView, numeratorEta_);
  }

  /*fill numerator and denominator distributions with tag-pass events (fill once for probe object 
    ==> pTagTagView->size() == 1, pTagPassView->size() == 2)*/
  if ((tagTagCollSize == 1) && (pTagPassView.isValid()) && (pTagPassView->size() == 2)) {

    //fill only the probe that is not the tag
    for (unsigned int i = 0; i < pTagPassView->size(); ++i) {
      if (pTagPassView->refAt(i).key() != pTagTagView->refAt(0).key()) {
	const double pT = pTagPassView->refAt(i)->pt();
	const double eta = pTagPassView->refAt(i)->eta();
	denominatorPT_->Fill(pT);
	numeratorPT_->Fill(pT);
	denominatorEta_->Fill(eta);
	numeratorEta_->Fill(eta);
      }
    }
  }

  /*fill denominator distribution with tag-fail events (fill once for probe object 
    ==> pTagTagView->size() == 1, pTagFailView->size() == 2)*/
  if ((tagTagCollSize == 1) && (pTagFailView->size() == 2)) {

    //fill only the probe that is not the tag
    for (unsigned int i = 0; i < pTagFailView->size(); ++i) {
      if (pTagFailView->refAt(i).key() != pTagTagView->refAt(0).key()) {
	const double pT = pTagFailView->refAt(i)->pt();
	const double eta = pTagFailView->refAt(i)->eta();
	denominatorPT_->Fill(pT);
	denominatorEta_->Fill(eta);
      }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
template<class T>
void TagAndProbeEfficiencyAnalyzer<T>::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book pT histograms
//   const Double_t bins[11] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 100.0};
  denominatorPT_ = new TH1F("denominatorPT", "", 20, 0.0, 100.0/*10, bins*/);
  numeratorPT_ = new TH1F("numeratorPT", "", 20, 0.0, 100.0/*10, bins*/);

  //book eta histograms
  denominatorEta_ = new TH1F("denominatorEta", "", 20, -5.0, 5.0);
  numeratorEta_ = new TH1F("numeratorEta", "", 20, -5.0, 5.0);
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void TagAndProbeEfficiencyAnalyzer<T>::endJob() 
{
  //write output file
  out_->cd();
  denominatorPT_->Write();
  numeratorPT_->Write();
  denominatorEta_->Write();
  numeratorEta_->Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
template<class T>
void TagAndProbeEfficiencyAnalyzer<T>::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
void TagAndProbeEfficiencyAnalyzer<T>::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
void TagAndProbeEfficiencyAnalyzer<T>::beginLuminosityBlock(edm::LuminosityBlock const&, 
						    edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
void TagAndProbeEfficiencyAnalyzer<T>::endLuminosityBlock(edm::LuminosityBlock const&, 
						  edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
template<class T>
void
TagAndProbeEfficiencyAnalyzer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template<class T>
void TagAndProbeEfficiencyAnalyzer<T>::fillPTHistogram(edm::Handle<edm::View<T> >& pView, 
						       TH1F* hist)
{
  for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->pt());
}

template<class T>
void TagAndProbeEfficiencyAnalyzer<T>::fillEtaHistogram(edm::Handle<edm::View<T> >& pView, 
							TH1F* hist)
{
  for (unsigned int i = 0; i < pView->size(); ++i) {
    if (pView->refAt(i)->pt() > 26.0) hist->Fill(pView->refAt(i)->eta());
  }
}

template<class T>
void TagAndProbeEfficiencyAnalyzer<T>::reset(const bool doDelete)
{
  if ((doDelete) && (out_ != NULL)) delete out_;
  out_ = NULL;
  if ((doDelete) && (denominatorPT_ != NULL)) delete denominatorPT_;
  denominatorPT_ = NULL;
  if ((doDelete) && (numeratorPT_ != NULL)) delete numeratorPT_;
  numeratorPT_ = NULL;
  if ((doDelete) && (denominatorEta_ != NULL)) delete denominatorEta_;
  denominatorEta_ = NULL;
  if ((doDelete) && (numeratorEta_ != NULL)) delete numeratorEta_;
  numeratorEta_ = NULL;
}

//define this as a plug-in
typedef TagAndProbeEfficiencyAnalyzer<reco::GenParticle> GenParticleTagAndProbeEfficiencyAnalyzer;
typedef TagAndProbeEfficiencyAnalyzer<reco::Muon> MuonTagAndProbeEfficiencyAnalyzer;
typedef TagAndProbeEfficiencyAnalyzer<reco::PFTau> TauTagAndProbeEfficiencyAnalyzer;
DEFINE_FWK_MODULE(GenParticleTagAndProbeEfficiencyAnalyzer);
DEFINE_FWK_MODULE(MuonTagAndProbeEfficiencyAnalyzer);
DEFINE_FWK_MODULE(TauTagAndProbeEfficiencyAnalyzer);
