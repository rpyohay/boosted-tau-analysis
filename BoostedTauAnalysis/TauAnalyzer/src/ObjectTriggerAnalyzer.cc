// -*- C++ -*-
//
// Package:    ObjectTriggerAnalyzer
// Class:      ObjectTriggerAnalyzer
// 
/**\class ObjectTriggerAnalyzer ObjectTriggerAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/ObjectTriggerAnalyzer.cc

   Description: analyze objects firing triggers

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: ObjectTriggerAnalyzer.cc,v 1.2 2012/09/19 10:57:14 yohay Exp $
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
class ObjectTriggerAnalyzer : public edm::EDAnalyzer {
public:
  explicit ObjectTriggerAnalyzer(const edm::ParameterSet&);
  ~ObjectTriggerAnalyzer();

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

  //fill pT histogram with object of arbitrary type
  template<typename U>
  void fillPTHistogramArbitrary(edm::Handle<edm::View<U> >& pView, TH1F* hist)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->pt());
  }

  //fill eta histogram
  void fillEtaHistogram(edm::Handle<edm::View<T> >&, TH1F*);

  //fill ET histogram
  template<typename U>
  void fillETHistogram(edm::Handle<edm::View<U> >& pView, TH1F* hist)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->et());
  }

  //fill visible pT histogram
  void fillPTHistogram(const edm::Handle<reco::GenParticleRefVector>&, TH1F*, 
		       edm::Handle<reco::GenParticleCollection>&);

  //fill HT histogram
  template<typename U>
  void fillHTHistogram(std::map<TH1F*, edm::Handle<edm::View<U> > >& map, 
		       const unsigned int numObjs)
  {
    double HT = 0.0;
    unsigned int nObjs = 0;
    for (typename std::map<TH1F*, edm::Handle<edm::View<U> > >::iterator i = map.begin(); 
	 i != map.end(); ++i) {
      nObjs+=i->second->size();
      for (unsigned int j = 0; j < i->second->size(); ++j) HT+=(i->second->refAt(j)->et());
    }
    if (nObjs == 2) HT_->Fill(HT);
  }

  //fill multiple pT histograms
  void fillPTHistograms(std::map<TH1F*, edm::Handle<reco::GenParticleRefVector> >&, 
			edm::Handle<reco::GenParticleCollection>&);

  //make visible gen tau pT canvas for 1 decay mode, multiple pT ranks
  void makePTRankCanvas(TCanvas&, TLegend&, const std::string&, std::vector<TH1F*>&);

  //make visible gen tau pT canvas for 1 pT rank, multiple decay modes
  void makeDecayModeCanvas(TCanvas&, TLegend&, const std::string&, const unsigned int);

  //format and draw multiple pT histograms on one canvas
  void drawMultiplePTHistograms(TCanvas&, std::vector<TH1F*>&, const std::vector<unsigned int>&, 
				const std::vector<unsigned int>&, TLegend&, 
				const std::vector<std::string>&, const std::string&);

  // ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //denominator input tag
  edm::InputTag denominatorTag_;

  //numerator input tag
  edm::InputTag numeratorTag_;

  //vector of tau-->mu input tags, highest pT rank first
  std::vector<edm::InputTag> tauMuInputTags_;

  //vector of tau-->1-prong input tags, highest pT rank first
  std::vector<edm::InputTag> tau1ProngInputTags_;

  //vector of tau-->1-prong+1-pi0 input tags, highest pT rank first
  std::vector<edm::InputTag> tau1Prong1Pi0InputTags_;

  //vector of tau-->1-prong+2-pi0 input tags, highest pT rank first
  std::vector<edm::InputTag> tau1Prong2Pi0InputTags_;

  //vector of tau-->3-prong input tags, highest pT rank first
  std::vector<edm::InputTag> tau3ProngInputTags_;

  //base gen particle tag
  edm::InputTag genParticleTag_;

  //tag for single input collection (i.e. for plots not split by pT rank and/or decay mode)
  edm::InputTag singleCollTag_;

  //vector of jet tags
  std::vector<edm::InputTag> jetTags_;

  //MET tag
  edm::InputTag METTag_;

  //set of parameters for GenTauDecayID class
  edm::ParameterSet genTauDecayIDPSet_;

  //flag indicating whether pT cuts should be applied in determining valid gen objects
  bool applyPTCuts_;

  //flag indicating whether KShorts should be counted as neutral hadrons
  double countKShort_;

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

  //vector of tau-->mu visible pT histograms, highest pT rank first
  std::vector<TH1F*> tauMuPTHists_;

  //vector of tau-->1-prong visible pT histograms, highest pT rank first
  std::vector<TH1F*> tau1ProngPTHists_;

  //vector of tau-->1-prong+1-pi0 visible pT histograms, highest pT rank first
  std::vector<TH1F*> tau1Prong1Pi0PTHists_;

  //vector of tau-->1-prong+2-pi0 visible pT histograms, highest pT rank first
  std::vector<TH1F*> tau1Prong2Pi0PTHists_;

  //vector of tau-->3-prong visible pT histograms, highest pT rank first
  std::vector<TH1F*> tau3ProngPTHists_;

  //histogram of eta for the single collection
  TH1F* eta_;

  //vector of jet pT histograms, highest pT rank first
  std::vector<TH1F*> jetPTHists_;

  //histogram of MET
  TH1F* MET_;

  //histogram of HT
  TH1F* HT_;
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
ObjectTriggerAnalyzer<T>::ObjectTriggerAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  denominatorTag_(iConfig.getParameter<edm::InputTag>("denominatorTag")),
  numeratorTag_(iConfig.getParameter<edm::InputTag>("numeratorTag")),
//   tauMuInputTags_(iConfig.getParameter<std::vector<edm::InputTag> >("tauMuInputTags")),
//   tau1ProngInputTags_(iConfig.getParameter<std::vector<edm::InputTag> >("tau1ProngInputTags")),
//   tau1Prong1Pi0InputTags_(iConfig.getParameter<std::vector<edm::InputTag> >
// 			  ("tau1Prong1Pi0InputTags")),
//   tau1Prong2Pi0InputTags_(iConfig.getParameter<std::vector<edm::InputTag> >
// 			  ("tau1Prong2Pi0InputTags")),
//   tau3ProngInputTags_(iConfig.getParameter<std::vector<edm::InputTag> >("tau3ProngInputTags")),
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
//   singleCollTag_(iConfig.getParameter<edm::InputTag>("singleCollTag")),
//   jetTags_(iConfig.getParameter<std::vector<edm::InputTag> >("jetTags")),
  METTag_(iConfig.getParameter<edm::InputTag>("METTag")),
  genTauDecayIDPSet_(iConfig.getParameter<edm::ParameterSet>("genTauDecayIDPSet")),
  applyPTCuts_(iConfig.getParameter<bool>("applyPTCuts")),
  countKShort_(iConfig.getParameter<bool>("countKShort")),
  pTRankColors_(iConfig.getParameter<std::vector<unsigned int> >("pTRankColors")),
  decayModeColors_(iConfig.getParameter<std::vector<unsigned int> >("decayModeColors")),
  pTRankStyles_(iConfig.getParameter<std::vector<unsigned int> >("pTRankStyles")),
  decayModeStyles_(iConfig.getParameter<std::vector<unsigned int> >("decayModeStyles")),
  pTRankEntries_(iConfig.getParameter<std::vector<std::string> >("pTRankEntries")),
  decayModeEntries_(iConfig.getParameter<std::vector<std::string> >("decayModeEntries"))
{
  //now do what ever initialization is needed
  reset(false);
}

template<class T>
ObjectTriggerAnalyzer<T>::~ObjectTriggerAnalyzer()
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
void ObjectTriggerAnalyzer<T>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get denominator collection
  edm::Handle<edm::View<T> > pDenominatorView;
  iEvent.getByLabel(denominatorTag_, pDenominatorView);

  //get numerator collection
  edm::Handle<edm::View<T> > pNumeratorView;
  iEvent.getByLabel(numeratorTag_, pNumeratorView);

//   //get gen object collections split by decay type and pT rank and associate each to a pT histogram
//   std::map<TH1F*, edm::Handle<reco::GenParticleRefVector> > genParticleCollMap;
//   for (std::vector<edm::InputTag>::const_iterator iTag = tauMuInputTags_.begin(); 
//        iTag != tauMuInputTags_.end(); ++iTag) {
//     edm::Handle<reco::GenParticleRefVector> pRefVector;
//     iEvent.getByLabel(*iTag, pRefVector);
//     genParticleCollMap[tauMuPTHists_[iTag - tauMuInputTags_.begin()]] = pRefVector;
//   }
//   for (std::vector<edm::InputTag>::const_iterator iTag = tau1ProngInputTags_.begin(); 
//        iTag != tau1ProngInputTags_.end(); ++iTag) {
//     edm::Handle<reco::GenParticleRefVector> pRefVector;
//     iEvent.getByLabel(*iTag, pRefVector);
//     genParticleCollMap[tau1ProngPTHists_[iTag - tau1ProngInputTags_.begin()]] = pRefVector;
//   }
//   for (std::vector<edm::InputTag>::const_iterator iTag = tau1Prong1Pi0InputTags_.begin(); 
//        iTag != tau1Prong1Pi0InputTags_.end(); ++iTag) {
//     edm::Handle<reco::GenParticleRefVector> pRefVector;
//     iEvent.getByLabel(*iTag, pRefVector);
//     genParticleCollMap[tau1Prong1Pi0PTHists_[iTag - tau1Prong1Pi0InputTags_.begin()]] = pRefVector;
//   }
//   for (std::vector<edm::InputTag>::const_iterator iTag = tau1Prong2Pi0InputTags_.begin(); 
//        iTag != tau1Prong2Pi0InputTags_.end(); ++iTag) {
//     edm::Handle<reco::GenParticleRefVector> pRefVector;
//     iEvent.getByLabel(*iTag, pRefVector);
//     genParticleCollMap[tau1Prong2Pi0PTHists_[iTag - tau1Prong2Pi0InputTags_.begin()]] = pRefVector;
//   }
//   for (std::vector<edm::InputTag>::const_iterator iTag = tau3ProngInputTags_.begin(); 
//        iTag != tau3ProngInputTags_.end(); ++iTag) {
//     edm::Handle<reco::GenParticleRefVector> pRefVector;
//     iEvent.getByLabel(*iTag, pRefVector);
//     genParticleCollMap[tau3ProngPTHists_[iTag - tau3ProngInputTags_.begin()]] = pRefVector;
//   }

  //get base gen particle collection
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

//   //get single collection tag
//   edm::Handle<edm::View<T> > pSingleColl;
//   iEvent.getByLabel(singleCollTag_, pSingleColl);

//   //get jet tags
//   std::map<TH1F*, edm::Handle<edm::View<reco::PFJet> > > jetCollMap;
//   for (std::vector<edm::InputTag>::const_iterator iTag = jetTags_.begin(); 
//        iTag != jetTags_.end(); ++iTag) {
//     edm::Handle<edm::View<reco::PFJet> > pView;
//     iEvent.getByLabel(*iTag, pView);
//     jetCollMap[jetPTHists_[iTag - jetTags_.begin()]] = pView;
//   }

  //get MET tag
  edm::Handle<edm::View<reco::PFMET> > pMET;
  iEvent.getByLabel(METTag_, pMET);

//   //plot pT distributions
//   fillPTHistogram(pDenominatorView, denominatorPT_);
//   fillPTHistogram(pNumeratorView, numeratorPT_);

//plot tag/probe pT distributions (debug)
  if (pNumeratorView->size() == 2) {
    fillPTHistogram(pNumeratorView, denominatorPT_);
    fillPTHistogram(pNumeratorView, numeratorPT_);
    fillEtaHistogram(pNumeratorView, denominatorEta_);
    fillEtaHistogram(pNumeratorView, numeratorEta_);
  }
  if ((pNumeratorView->size() == 1) && (pDenominatorView->size() == 2)) {
    //fill only the probe that is not the tag
    for (unsigned int i = 0; i < pDenominatorView->size(); ++i) {
      if (pDenominatorView->refAt(i).key() != pNumeratorView->refAt(0).key()) {
	denominatorPT_->Fill(pDenominatorView->refAt(i)->pt());
	if (pDenominatorView->refAt(i)->pt() > 26.0) {
	  denominatorEta_->Fill(pDenominatorView->refAt(i)->eta());
	}
      }
    }
  }

  //debug
  fillEtaHistogram(pDenominatorView, eta_);

//   //fill pT histograms for gen objects, 1 per decay type and pT rank
//   fillPTHistograms(genParticleCollMap, pGenParticles);

//   //plot eta distribution
//   fillEtaHistogram(pSingleColl, eta_);

//   //fill pT histograms for jets, 1 per pT rank
//   for (std::map<TH1F*, edm::Handle<edm::View<reco::PFJet> > >::iterator i = jetCollMap.begin(); 
//        i != jetCollMap.end(); ++i) { fillPTHistogramArbitrary(i->second, i->first); }

  //plot MET distribution
  fillETHistogram(pMET, MET_);

//   /*plot HT distribution (sum ET of jets matched to di-tau objects in the case where there are 2 
//     di-tau objects with jet matches)*/
//   fillHTHistogram(jetCollMap, 2);
}


// ------------ method called once each job just before starting event loop  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book pT histograms
  const Double_t bins[11] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 100.0};
  denominatorPT_ = new TH1F("denominatorPT", "", /*20, 0.0, 100.0*/10, bins);
  numeratorPT_ = new TH1F("numeratorPT", "", /*20, 0.0, 100.0*/10, bins);

  //book eta histograms
  denominatorEta_ = new TH1F("denominatorEta", "", 20, -5.0, 5.0);
  numeratorEta_ = new TH1F("numeratorEta", "", 20, -5.0, 5.0);

//   //book gen particle pT histograms split by decay type and pT rank
//   for (std::vector<edm::InputTag>::const_iterator iTag = tauMuInputTags_.begin(); 
//        iTag != tauMuInputTags_.end(); ++iTag) {
//     tauMuPTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
// 				     .c_str(), "", 20, 0.0, 100.0));
//   }
//   for (std::vector<edm::InputTag>::const_iterator iTag = tau1ProngInputTags_.begin(); 
//        iTag != tau1ProngInputTags_.end(); ++iTag) {
//     tau1ProngPTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
// 					 .c_str(), "", 20, 0.0, 100.0));
//   }
//   for (std::vector<edm::InputTag>::const_iterator iTag = tau1Prong1Pi0InputTags_.begin(); 
//        iTag != tau1Prong1Pi0InputTags_.end(); ++iTag) {
//     tau1Prong1Pi0PTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
// 					     .c_str(), "", 20, 0.0, 100.0));
//   }
//   for (std::vector<edm::InputTag>::const_iterator iTag = tau1Prong2Pi0InputTags_.begin(); 
//        iTag != tau1Prong2Pi0InputTags_.end(); ++iTag) {
//     tau1Prong2Pi0PTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
// 					     .c_str(), "", 20, 0.0, 100.0));
//   }
//   for (std::vector<edm::InputTag>::const_iterator iTag = tau3ProngInputTags_.begin(); 
//        iTag != tau3ProngInputTags_.end(); ++iTag) {
//     tau3ProngPTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
// 					 .c_str(), "", 20, 0.0, 100.0));
//   }

  //book gen particle eta histogram
  eta_ = new TH1F("eta", "", 20, -5.0, 5.0);

//   //book jet pT histograms split by pT rank
//   for (std::vector<edm::InputTag>::const_iterator iTag = jetTags_.begin(); 
//        iTag != jetTags_.end(); ++iTag) {
//     jetPTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
// 				   .c_str(), "", 20, 0.0, 100.0));
//   }

  //book reco MET histogram
  MET_ = new TH1F("MET", "", 20, 0.0, 100.0);

  //book reco HT histogram
  HT_ = new TH1F("HT", "", 100, 0.0, 500.0);
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::endJob() 
{
//   //make the gen particle pT canvases, 1 per decay mode
//   TCanvas tauMuPTRankCanvas("tauMuPTRankCanvas", "", 600, 600);
//   TCanvas tau1ProngPTRankCanvas("tau1ProngPTRankCanvas", "", 600, 600);
//   TCanvas tau1Prong1Pi0PTRankCanvas("tau1Prong1Pi0PTRankCanvas", "", 600, 600);
//   TCanvas tau1Prong2Pi0PTRankCanvas("tau1Prong2Pi0PTRankCanvas", "", 600, 600);
//   TCanvas tau3ProngPTRankCanvas("tau3ProngPTRankCanvas", "", 600, 600);
//   TLegend tauMuPTRankLegend(0.4, 0.6, 0.8, 0.8);
//   TLegend tau1ProngPTRankLegend(0.4, 0.6, 0.8, 0.8);
//   TLegend tau1Prong1Pi0PTRankLegend(0.4, 0.6, 0.8, 0.8);
//   TLegend tau1Prong2Pi0PTRankLegend(0.4, 0.6, 0.8, 0.8);
//   TLegend tau3ProngPTRankLegend(0.4, 0.6, 0.8, 0.8);
//   makePTRankCanvas(tauMuPTRankCanvas, tauMuPTRankLegend, 
// 		   "gg fusion NMSSM Higgs #tau_{#mu}", tauMuPTHists_);
//   makePTRankCanvas(tau1ProngPTRankCanvas, tau1ProngPTRankLegend, 
// 		   "gg fusion NMSSM Higgs #tau_{had}, 1 prong", tau1ProngPTHists_);
//   makePTRankCanvas(tau1Prong1Pi0PTRankCanvas, tau1Prong1Pi0PTRankLegend, 
// 		   "gg fusion NMSSM Higgs #tau_{had}, 1 prong + 1 #pi^{0}", tau1Prong1Pi0PTHists_);
//   makePTRankCanvas(tau1Prong2Pi0PTRankCanvas, tau1Prong2Pi0PTRankLegend, 
// 		   "gg fusion NMSSM Higgs #tau_{had}, 1 prong + 2 #pi^{0}", tau1Prong2Pi0PTHists_);
//   makePTRankCanvas(tau3ProngPTRankCanvas, tau3ProngPTRankLegend, 
// 		   "gg fusion NMSSM Higgs #tau_{had}, 3 prong", tau3ProngPTHists_);

//   //make the gen particle pT canvases, 1 per pT rank
//   TCanvas pTRank0DecayModeCanvas("pTRank0DecayModeCanvas", "", 600, 600);
//   TCanvas pTRank1DecayModeCanvas("pTRank1DecayModeCanvas", "", 600, 600);
//   TCanvas pTRank2DecayModeCanvas("pTRank2DecayModeCanvas", "", 600, 600);
//   TCanvas pTRank3DecayModeCanvas("pTRank3DecayModeCanvas", "", 600, 600);
//   TLegend pTRank0DecayModeLegend(0.4, 0.6, 0.8, 0.8);
//   TLegend pTRank1DecayModeLegend(0.4, 0.6, 0.8, 0.8);
//   TLegend pTRank2DecayModeLegend(0.4, 0.6, 0.8, 0.8);
//   TLegend pTRank3DecayModeLegend(0.4, 0.6, 0.8, 0.8);
//   makeDecayModeCanvas(pTRank0DecayModeCanvas, pTRank0DecayModeLegend, 
// 		      "gg fusion NMSSM Higgs highest p_{T} object", 0);
//   makeDecayModeCanvas(pTRank1DecayModeCanvas, pTRank1DecayModeLegend, 
// 		      "gg fusion NMSSM Higgs second highest p_{T} object", 1);
//   makeDecayModeCanvas(pTRank2DecayModeCanvas, pTRank2DecayModeLegend, 
// 		      "gg fusion NMSSM Higgs third highest p_{T} object", 2);
//   makeDecayModeCanvas(pTRank3DecayModeCanvas, pTRank3DecayModeLegend, 
// 		      "gg fusion NMSSM Higgs lowest p_{T} object", 3);

//   //make the jet pT canvas
//   TCanvas jetPTRankCanvas("jetPTRankCanvas", "", 600, 600);
//   TLegend jetPTRankLegend(0.4, 0.6, 0.8, 0.8);
//   makePTRankCanvas(jetPTRankCanvas, jetPTRankLegend, 
// 		   "gg fusion NMSSM Higgs-matched AK5 jets", jetPTHists_);

  //write output file
  out_->cd();
  denominatorPT_->Write();
  numeratorPT_->Write();
  denominatorEta_->Write();
  numeratorEta_->Write();
//   for (std::vector<TH1F*>::iterator iHist = tauMuPTHists_.begin(); 
//        iHist != tauMuPTHists_.end(); ++iHist) { (*iHist)->Write(); }
//   for (std::vector<TH1F*>::iterator iHist = tau1ProngPTHists_.begin(); 
//        iHist != tau1ProngPTHists_.end(); ++iHist) { (*iHist)->Write(); }
//   for (std::vector<TH1F*>::iterator iHist = tau1Prong1Pi0PTHists_.begin(); 
//        iHist != tau1Prong1Pi0PTHists_.end(); ++iHist) { (*iHist)->Write(); }
//   for (std::vector<TH1F*>::iterator iHist = tau1Prong2Pi0PTHists_.begin(); 
//        iHist != tau1Prong2Pi0PTHists_.end(); ++iHist) { (*iHist)->Write(); }
//   for (std::vector<TH1F*>::iterator iHist = tau3ProngPTHists_.begin(); 
//        iHist != tau3ProngPTHists_.end(); ++iHist) { (*iHist)->Write(); }
  eta_->Write();
//   for (std::vector<TH1F*>::iterator iHist = jetPTHists_.begin(); 
//        iHist != jetPTHists_.end(); ++iHist) { (*iHist)->Write(); }
  MET_->Write();
  HT_->Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::beginLuminosityBlock(edm::LuminosityBlock const&, 
						    edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::endLuminosityBlock(edm::LuminosityBlock const&, 
						  edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
template<class T>
void ObjectTriggerAnalyzer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template<class T>
void ObjectTriggerAnalyzer<T>::fillPTHistogram(edm::Handle<edm::View<T> >& pView, TH1F* hist)
{
  for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->pt());
}

template<class T>
void ObjectTriggerAnalyzer<T>::fillEtaHistogram(edm::Handle<edm::View<T> >& pView, TH1F* hist)
{
  for (unsigned int i = 0; i < pView->size(); ++i) {
    if (pView->refAt(i)->pt() > 26.0) hist->Fill(pView->refAt(i)->eta());
  }
}

template<class T>
void ObjectTriggerAnalyzer<T>::fillPTHistogram(const edm::Handle<reco::GenParticleRefVector>& 
					       pRefVector, TH1F* hist, 
					       edm::Handle<reco::GenParticleCollection>& 
					       pGenParticles)
{
  for (unsigned int i = 0; i < pRefVector->size(); ++i) {
    try {
      GenTauDecayID tauDecay(genTauDecayIDPSet_, pGenParticles, 
			     Common::getStatus3Key(pRefVector, pGenParticles, i));
      tauDecay.tauDecayType(applyPTCuts_, countKShort_);
      hist->Fill(tauDecay.getVisibleTauP4().Pt());
    }
    catch (std::string& ex) { throw cms::Exception("ObjectTriggerAnalyzer<T>") << ex; }
  }
}

template<class T>
void ObjectTriggerAnalyzer<T>::fillPTHistograms(std::map<TH1F*, 
						edm::Handle<reco::GenParticleRefVector> >& map, 
						edm::Handle<reco::GenParticleCollection>& 
						pGenParticles)
{
  for (std::map<TH1F*, edm::Handle<reco::GenParticleRefVector> >::iterator i = map.begin(); 
       i != map.end(); ++i) { fillPTHistogram(i->second, i->first, pGenParticles); }
}

template<class T>
void ObjectTriggerAnalyzer<T>::makePTRankCanvas(TCanvas& canvas, TLegend& legend, 
						const std::string& header, 
						std::vector<TH1F*>& hists)
{
  drawMultiplePTHistograms(canvas, hists, pTRankColors_, pTRankStyles_, legend, 
			   pTRankEntries_, header);
}

template<class T>
void ObjectTriggerAnalyzer<T>::makeDecayModeCanvas(TCanvas& canvas, TLegend& legend, 
						   const std::string& header, 
						   const unsigned int pTRank)
{
  if ((pTRank >= tauMuPTHists_.size()) || (pTRank >= tau1ProngPTHists_.size()) || 
      (pTRank >= tau1Prong1Pi0PTHists_.size()) || (pTRank >= tau1Prong2Pi0PTHists_.size()) || 
      (pTRank >= tau3ProngPTHists_.size())) {
    std::cerr << "Warning: rank " << pTRank << " is too large.  Skipping this canvas.\n";
    return;
  }
  std::vector<TH1F*> hists;
  hists.push_back(tauMuPTHists_[pTRank]);
  hists.push_back(tau1ProngPTHists_[pTRank]);
  hists.push_back(tau1Prong1Pi0PTHists_[pTRank]);
  hists.push_back(tau1Prong2Pi0PTHists_[pTRank]);
  hists.push_back(tau3ProngPTHists_[pTRank]);
  drawMultiplePTHistograms(canvas, hists, decayModeColors_, decayModeStyles_, legend, 
			   decayModeEntries_, header);
}

template<class T>
void ObjectTriggerAnalyzer<T>::drawMultiplePTHistograms(TCanvas& canvas, 
							std::vector<TH1F*>& hists, 
							const std::vector<unsigned int>& colors, 
							const std::vector<unsigned int>& styles, 
							TLegend& legend, 
							const std::vector<std::string>& entries, 
							const std::string& header)
{
  Common::setLegendOptions(legend, header.c_str());
  canvas.cd();
  Common::setCanvasOptions(canvas, 1, 0, 0);
  TH1F* pHistWithMaxMaxBin = NULL;
  Double_t maxBinContent = 0.0;
  for (std::vector<TH1F*>::iterator iHist = hists.begin(); 
       iHist != hists.end(); ++iHist) {
    const unsigned int i = iHist - hists.begin();
    Common::setHistogramOptions(*iHist, colors[i], 0.7, styles[i], 1.0, 
				"p_{T} (GeV)", "", 0.04);
    legend.AddEntry(*iHist, entries[i].c_str(), "l");
    Int_t histMaxBin = (*iHist)->GetMaximumBin();
    Double_t histMaxBinContent = (*iHist)->GetBinContent(histMaxBin);
    if (histMaxBinContent > maxBinContent) {
      maxBinContent = histMaxBinContent;
      pHistWithMaxMaxBin = *iHist;
    }
  }
  pHistWithMaxMaxBin->Draw();
  for (std::vector<TH1F*>::iterator iHist = hists.begin(); 
       iHist != hists.end(); ++iHist) {
    (*iHist)->Draw("SAME");
  }
  legend.Draw();
  canvas.Write();
}

template<class T>
void ObjectTriggerAnalyzer<T>::reset(const bool doDelete)
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
  for (std::vector<TH1F*>::iterator iHist = tauMuPTHists_.begin(); 
       iHist != tauMuPTHists_.end(); ++iHist) {
    if ((doDelete) && (*iHist != NULL)) delete *iHist;
    *iHist = NULL;
  }
  for (std::vector<TH1F*>::iterator iHist = tau1ProngPTHists_.begin(); 
       iHist != tau1ProngPTHists_.end(); ++iHist) {
    if ((doDelete) && (*iHist != NULL)) delete *iHist;
    *iHist = NULL;
  }
  for (std::vector<TH1F*>::iterator iHist = tau1Prong1Pi0PTHists_.begin(); 
       iHist != tau1Prong1Pi0PTHists_.end(); ++iHist) {
    if ((doDelete) && (*iHist != NULL)) delete *iHist;
    *iHist = NULL;
  }
  for (std::vector<TH1F*>::iterator iHist = tau1Prong2Pi0PTHists_.begin(); 
       iHist != tau1Prong2Pi0PTHists_.end(); ++iHist) {
    if ((doDelete) && (*iHist != NULL)) delete *iHist;
    *iHist = NULL;
  }
  for (std::vector<TH1F*>::iterator iHist = tau3ProngPTHists_.begin(); 
       iHist != tau3ProngPTHists_.end(); ++iHist) {
    if ((doDelete) && (*iHist != NULL)) delete *iHist;
    *iHist = NULL;
  }
  if ((doDelete) && (eta_ != NULL)) delete eta_;
  eta_ = NULL;
  for (std::vector<TH1F*>::iterator iHist = jetPTHists_.begin(); 
       iHist != jetPTHists_.end(); ++iHist) {
    if ((doDelete) && (*iHist != NULL)) delete *iHist;
    *iHist = NULL;
  }
  if ((doDelete) && (MET_ != NULL)) delete MET_;
  MET_ = NULL;
  if ((doDelete) && (HT_ != NULL)) delete HT_;
  HT_ = NULL;
}

//define this as a plug-in
typedef ObjectTriggerAnalyzer<reco::GenParticle> GenParticleTriggerAnalyzer;
typedef ObjectTriggerAnalyzer<reco::Muon> MuonTriggerAnalyzer;
DEFINE_FWK_MODULE(GenParticleTriggerAnalyzer);
DEFINE_FWK_MODULE(MuonTriggerAnalyzer);
