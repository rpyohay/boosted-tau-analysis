// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      MuonIsolationAnalyzer
// 
/**\class MuonIsolationAnalyzer MuonIsolationAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/MuonIsolationAnalyzer.cc

   Description: analyze isolation properties of objects firing triggers

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: MuonIsolationAnalyzer.cc,v 1.1 2012/09/19 10:57:14 yohay Exp $
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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;
using namespace edm;
using namespace reco;
//
// class declaration
//

class MuonIsolationAnalyzer : public edm::EDAnalyzer {
public:
  explicit MuonIsolationAnalyzer(const edm::ParameterSet&);
  ~MuonIsolationAnalyzer();

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

  // ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //muon tag
  edm::InputTag muonTag_;

  //collection of associated HPS reco taus (for signal only!)
  edm::InputTag tauTag_;

  //map of reco muons to to associated cleaned jets (for signal only!)
  edm::InputTag muonJetMapTag_;

  //flag for signal
  bool isSignal_;

  //PU subtraction coefficient for muon PF isolation
  double muonPFIsoPUSubtractionCoeff_;

  //histogram of combined particle isolation vs. muon pT
  TH1F* recoMuPFIso_;
  
  //histogram of combined particle relative isolation vs. muon pT
  TH1F* recoMuPFRelIso_;
  
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
MuonIsolationAnalyzer::MuonIsolationAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  tauTag_(iConfig.existsAs<edm::InputTag>("tauTag") ?
  	  iConfig.getParameter<edm::InputTag>("tauTag") : edm::InputTag()),
  muonJetMapTag_(iConfig.existsAs<edm::InputTag>("muonJetMapTag") ?
  	  iConfig.getParameter<edm::InputTag>("muonJetMapTag") : edm::InputTag()),
  isSignal_(iConfig.getParameter<bool>("isSignal")),
  muonPFIsoPUSubtractionCoeff_(iConfig.getParameter<double>("muonPFIsoPUSubtractionCoeff"))
{
  //now do what ever initialization is needed
  reset(false);
}

MuonIsolationAnalyzer::~MuonIsolationAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void MuonIsolationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get muon collection
  // edm::Handle<edm::View<reco::Muon> > pMuons;
  //iEvent.getByLabel(muonTag_, pMuons);
  edm::Handle<reco::MuonRefVector> pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  edm::Handle<reco::PFTauRefVector> pTaus;
  if (tauTag_ == edm::InputTag()) {} // do nothing
  else iEvent.getByLabel(tauTag_, pTaus);

  std::vector<reco::PFTauRef> pTSortedTaus;
  std::vector<reco::PFTauRef> taus;
  if (isSignal_ && pTaus.isValid())
    {
      //sort selected taus by descending order in pT
      for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
	   ++iTau) { pTSortedTaus.push_back(*iTau); }
      taus = pTSortedTaus;
      Common::sortByPT(pTSortedTaus);
      std::reverse(pTSortedTaus.begin(), pTSortedTaus.end());
      cout << "number of taus = " << taus.size() << endl;
    }

  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  if (muonJetMapTag_ == edm::InputTag()) {} // do nothing
  else iEvent.getByLabel(muonJetMapTag_, pMuonJetMap);

  std::vector<reco::MuonRef> muonRefs;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); ++iMuon)
    {
	muonRefs.push_back(*iMuon);
    }
  Common::sortByPT(muonRefs);
  cout << "number of muon refs = " << muonRefs.size() << endl;
  //plot combined particle isolation vs. muon pT
  //  for (std::vector<reco::MuonRef>::const_iterator iMuon = muonRefs.begin(); iMuon != muonRefs.end(); ++iMuon) {
  //edm::RefToBase<reco::Muon> muon(pMuons->refAt(iMuon));

  double iso;
  double pT = (muonRefs[muonRefs.size() - 1])->pt();
  
  if (isSignal_ && pTaus.isValid())
    { // calculate modified tau isolation
      bool foundRemovedMu = false;
      std::vector<reco::PFTauRef>::const_iterator iTau = taus.begin();
      while ((iTau != taus.end()) && !foundRemovedMu)
	{ // loop over tau refs
	  const reco::MuonRefVector& removedMuons = (*pMuonJetMap)[(*iTau)->jetRef()]; // get removed muons for tau
	  reco::LeafCandidate::LorentzVector tauHadVisibleP4 = (*iTau)->p4();
	  for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
	       iMuon != removedMuons.end(); ++iMuon) { // loop over removed muons
	    if ((*iMuon).key() == muonRefs[muonRefs.size() - 1].key())
	      { // if this is the highest-pT mu
		//		iso = Common::getMuonCombPFIsoMinusTau(*(muonRefs[muonRefs.size() - 1]), tauHadVisibleP4, muonPFIsoPUSubtractionCoeff_); // this is wrong!!!!
		iso = Common::getMuonCombPFIsoModified(*(muonRefs[muonRefs.size() - 1]), *iTau, muonPFIsoPUSubtractionCoeff_);
		//iso = Common::getMuonCombPFIso(*(muonRefs[muonRefs.size() - 1]), muonPFIsoPUSubtractionCoeff_);
		recoMuPFRelIso_->Fill(iso/pT);
		recoMuPFIso_->Fill(iso);
		foundRemovedMu = true;
		break;
	      } // if this is the highest-pT mu
	  } // loop over removed muons
	  ++iTau;  
	} // loop over tau refs
    } // calculate modified tau isolation
  else
    {
      iso = Common::getMuonCombPFIso(*(muonRefs[muonRefs.size() - 1]), muonPFIsoPUSubtractionCoeff_);
      recoMuPFRelIso_->Fill(iso/pT);
      recoMuPFIso_->Fill(iso);
    }
  //  }
}


// ------------ method called once each job just before starting event loop  ------------
void MuonIsolationAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book muon isolation histograms
  recoMuPFIso_ = 
    new TH1F("recoMuPFIso", "", 1000, 0.0, 20.0);
  recoMuPFRelIso_ = 
    new TH1F("recoMuPFRelIso", "", 1000, 0.0, 20.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void MuonIsolationAnalyzer::endJob() 
{
  //make the muon isolation canvases
  //TCanvas recoMuPFRelIsoCanvas("recoMuPFRelIsoCanvas", "", 600, 600);
  //Common::setCanvasOptions(recoMuPFRelIsoCanvas, 0, 0, 0);
  //Common::setCanvasMargins(recoMuPFRelIsoCanvas, 0.2, 0.2, 0.2, 0.2);

  //format the muon isolation plots
  Common::setHistogramOptions(recoMuPFIso_, kBlack, 0.7, 20, 1.0, "Muon PFIso", "", 0.04);
  Common::setHistogramOptions(recoMuPFRelIso_, kBlack, 0.7, 20, 1.0, "Muon PFRelIso", "", 0.04);

  //draw muon isolation plots
  //recoMuPFRelIsoCanvas.cd();
  //recoMuPFRelIso_->Draw();

  //write output file
  out_->cd();
  recoMuPFIso_->Write();
  recoMuPFRelIso_->Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void MuonIsolationAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void MuonIsolationAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void MuonIsolationAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
						    edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void MuonIsolationAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
						  edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void MuonIsolationAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void MuonIsolationAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (out_ != NULL)) delete out_;
  out_ = NULL;
  if ((doDelete) && (recoMuPFIso_ != NULL)) delete recoMuPFIso_;
  recoMuPFIso_ = NULL;
  if ((doDelete) && (recoMuPFRelIso_ != NULL)) delete recoMuPFRelIso_;
  recoMuPFRelIso_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonIsolationAnalyzer);
