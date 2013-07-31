// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      IsolationAnalyzer
// 
/**\class IsolationAnalyzer IsolationAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/IsolationAnalyzer.cc

   Description: analyze isolation properties of objects firing triggers

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: IsolationAnalyzer.cc,v 1.1 2012/09/19 10:57:14 yohay Exp $
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
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

//
// class declaration
//

class IsolationAnalyzer : public edm::EDAnalyzer {
public:
  explicit IsolationAnalyzer(const edm::ParameterSet&);
  ~IsolationAnalyzer();

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
  TH2F* combParticleIsoVsMuonPT_;
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
IsolationAnalyzer::IsolationAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  muonPFIsoPUSubtractionCoeff_(iConfig.getParameter<double>("muonPFIsoPUSubtractionCoeff"))
{
  //now do what ever initialization is needed
  reset(false);
}

IsolationAnalyzer::~IsolationAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void IsolationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get muon collection
  edm::Handle<edm::View<reco::Muon> > pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  //plot combined particle isolation vs. muon pT
  for (unsigned int iMuon = 0; iMuon < pMuons->size(); ++iMuon) {
    edm::RefToBase<reco::Muon> muon(pMuons->refAt(iMuon));
    combParticleIsoVsMuonPT_->Fill(muon->pt(), 
				   Common::getMuonCombPFIso(*muon, muonPFIsoPUSubtractionCoeff_));
  }
}


// ------------ method called once each job just before starting event loop  ------------
void IsolationAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book muon isolation histograms
  combParticleIsoVsMuonPT_ = 
    new TH2F("combParticleIsoVsMuonPT", "", 20, 0.0, 100.0, 40, 0.0, 40.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void IsolationAnalyzer::endJob() 
{
  //make the muon isolation canvases
  TCanvas combParticleIsoVsMuonPTCanvas("combParticleIsoVsMuonPTCanvas", "", 600, 600);
  Common::setCanvasOptions(combParticleIsoVsMuonPTCanvas, 0, 0, 0);
  Common::setCanvasMargins(combParticleIsoVsMuonPTCanvas, 0.2, 0.2, 0.2, 0.2);

  //format the muon isolation plots
  Common::setHistogramOptions(combParticleIsoVsMuonPT_, kBlack, 0.7, 20, 1.6, 1.0, "p_{T} (GeV)", 
			      "Combined particle isolation (GeV)");

  //draw muon isolation plots
  combParticleIsoVsMuonPTCanvas.cd();
  combParticleIsoVsMuonPT_->Draw("COLZ");

  //write output file
  out_->cd();
  combParticleIsoVsMuonPTCanvas.Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void IsolationAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void IsolationAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void IsolationAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
						    edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void IsolationAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
						  edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void IsolationAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void IsolationAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (out_ != NULL)) delete out_;
  out_ = NULL;
  if ((doDelete) && (combParticleIsoVsMuonPT_ != NULL)) delete combParticleIsoVsMuonPT_;
  combParticleIsoVsMuonPT_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsolationAnalyzer);
