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

  //PU subtraction coefficient for muon PF isolation
  double muonPFIsoPUSubtractionCoeff_;

  //histogram of combined particle isolation vs. muon pT
  TH1F* combParticleIsoOverMuonPT_;
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
  edm::Handle<edm::View<reco::Muon> > pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  double highestpT = -99999.9;
  double iso;

  //plot combined particle isolation vs. muon pT
  for (unsigned int iMuon = 0; iMuon < pMuons->size(); ++iMuon) {
    edm::RefToBase<reco::Muon> muon(pMuons->refAt(iMuon));
//     if (muon->pt() > highestpT)
//       {
	highestpT = muon->pt();
	iso = Common::getMuonCombPFIso(*muon, muonPFIsoPUSubtractionCoeff_);
  combParticleIsoOverMuonPT_->Fill(iso/highestpT);
//       }
  }
//   combParticleIsoOverMuonPT_->Fill(iso/highestpT);
}


// ------------ method called once each job just before starting event loop  ------------
void MuonIsolationAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book muon isolation histograms
  combParticleIsoOverMuonPT_ = 
    new TH1F("combParticleIsoOverMuonPT", "", 100, 0.0, 10.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void MuonIsolationAnalyzer::endJob() 
{
  //make the muon isolation canvases
  TCanvas combParticleIsoOverMuonPTCanvas("combParticleIsoOverMuonPTCanvas", "", 600, 600);
  Common::setCanvasOptions(combParticleIsoOverMuonPTCanvas, 0, 0, 0);
  Common::setCanvasMargins(combParticleIsoOverMuonPTCanvas, 0.2, 0.2, 0.2, 0.2);

  //format the muon isolation plots
  Common::setHistogramOptions(combParticleIsoOverMuonPT_, kBlack, 0.7, 20, 1.0, "Muon isolation/pT", "", 0.04);

  //draw muon isolation plots
  combParticleIsoOverMuonPTCanvas.cd();
  combParticleIsoOverMuonPT_->Draw();

  //write output file
  out_->cd();
  combParticleIsoOverMuonPTCanvas.Write();
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
  if ((doDelete) && (combParticleIsoOverMuonPT_ != NULL)) delete combParticleIsoOverMuonPT_;
  combParticleIsoOverMuonPT_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonIsolationAnalyzer);
