// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      GenAnalyzer
// 
/**\class GenAnalyzer GenAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/GenAnalyzer.cc

   Description: plot gen quantities

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: GenAnalyzer.cc,v 1.1 2012/09/25 11:50:20 yohay Exp $
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
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

//
// class declaration
//

class GenAnalyzer : public edm::EDAnalyzer {
public:
  explicit GenAnalyzer(const edm::ParameterSet&);
  ~GenAnalyzer();

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

  //gen particle tag
  edm::InputTag genParticleTag_;

  //PU info tag
  edm::InputTag PUTag_;

  //set of parameters for GenTauDecayID class
  edm::ParameterSet genTauDecayIDPSet_;

  //histogram of dR between gen objects from a1 decay
  TH1F* dRA1TauDaughters_;

  //histogram of mu+had mu pT
  TH1F* tauMuPT_;

  //histogram of mu+had had pT
  TH1F* tauHadPT_;

  //histogram of true no. in-time interactions
  TH1D* trueNInt_;

  //a2 tau pair decay type vs. a1 tau pair decay type
  TH2F* a2TauPairDecayVsA1TauPairDecay_;
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
GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  PUTag_(iConfig.getParameter<edm::InputTag>("PUTag")),
  genTauDecayIDPSet_(iConfig.getParameter<edm::ParameterSet>("genTauDecayIDPSet"))
{
  //now do what ever initialization is needed
  reset(false);
}

GenAnalyzer::~GenAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get gen particle collection
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //get PU info
  edm::Handle<std::vector<PileupSummaryInfo> > pPU;
  iEvent.getByLabel(PUTag_, pPU);

//   //look at muon parentage
//   for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
//        iGenParticle != pGenParticles->end(); ++iGenParticle) {
//     if (fabs(iGenParticle->pdgId()) == 13)
//       { // if it's a muon
// 	std::cout << "Gen muon found with pT = " << iGenParticle->pt() << ", status " << iGenParticle->status() << std::endl;
// 	std::cout << "Its mother's pdgId was: " << iGenParticle->mother()->pdgId() << std::endl;
// 	std::cout << "Its grandmother's pdgId was: " << iGenParticle->mother()->mother()->pdgId() << std::endl;
// 	std::cout << "Its great-grandmothers' pdgId was: " << iGenParticle->mother()->mother()->mother()->pdgId() << std::endl;
//       } // if it's a muon
//   }

  //find a1 tau decay products
  std::vector<GenTauDecayID> aDecayProducts;
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {
    //    if (iGenParticle->pdgId() == 35)
    //     std::cout << "Mass of a = " << iGenParticle->mass() << std::endl;
    try {
      GenTauDecayID tauDecay(genTauDecayIDPSet_, pGenParticles, 
			     iGenParticle - pGenParticles->begin());
      if (tauDecay.isStatus3DecayProduct()) aDecayProducts.push_back(tauDecay);
    }
    catch (std::string& ex) { throw cms::Exception("GenAnalyzer") << ex; }
  }

  //containers for a tau decay types
  std::vector<std::pair<GenTauDecayID::DecayType, GenTauDecayID::DecayType> > aDecay;

  //loop over a1 tau daughters
  std::vector<unsigned int> keysToIgnore;
  for (std::vector<GenTauDecayID>::iterator iTau = aDecayProducts.begin(); 
       iTau != aDecayProducts.end(); ++iTau) {
    try {
      const unsigned int tauKey = iTau->getTauIndex();

      //find sister
      iTau->findSister();
      const unsigned int iSister = iTau->getSisterIndex();

      //if sister wasn't already looped over...
      if (std::find(keysToIgnore.begin(), keysToIgnore.end(), iSister) == keysToIgnore.end()) {

	//...plot dR(sisters)
	dRA1TauDaughters_->Fill(reco::deltaR(*reco::GenParticleRef(pGenParticles, tauKey), 
					     *reco::GenParticleRef(pGenParticles, iSister)));

	//...save pair decay mode
	aDecay.push_back(std::pair<GenTauDecayID::DecayType, 
			 GenTauDecayID::DecayType>(iTau->tauDecayType(false, true).second, 
						   iTau->sisterDecayType(false, true).second));

	//ignore this tau in the future
	keysToIgnore.push_back(tauKey);
      }

      //is this a mu+had decay?
      std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> thisDecay = 
	iTau->tauDecayType(false, true);
      std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> sisterDecay = 
	iTau->sisterDecayType(false, true);
      if (((thisDecay.second == GenTauDecayID::MU) && 
	   (sisterDecay.second == GenTauDecayID::HAD)) || 
	  ((thisDecay.second == GenTauDecayID::HAD) && 
	   (sisterDecay.second == GenTauDecayID::MU))) {

	//plot mu pT and had pT
	reco::LeafCandidate::LorentzVector visibleP4 = iTau->getVisibleTauP4();
	if (thisDecay.second == GenTauDecayID::MU) tauMuPT_->Fill(visibleP4.Pt());
	if (thisDecay.second == GenTauDecayID::HAD) tauHadPT_->Fill(visibleP4.Pt());
      }
    }
    catch (std::string& ex) { throw cms::Exception("GenAnalyzer") << ex; }
  }

  //plot a2 decay type vs. a1 decay type
  int val1 = -1;
  int val2 = -1;
  if (aDecay.size() == 2) {
    for (std::vector<std::pair<GenTauDecayID::DecayType, 
	   GenTauDecayID::DecayType> >::const_iterator iADecay = aDecay.begin(); 
	 iADecay != aDecay.end(); ++iADecay) {
      int val = -1;
      if ((iADecay->first == GenTauDecayID::MU) && (iADecay->second == GenTauDecayID::MU)) val = 0;
      if (((iADecay->first == GenTauDecayID::MU) && (iADecay->second == GenTauDecayID::E)) || 
	  ((iADecay->first == GenTauDecayID::E) && (iADecay->second == GenTauDecayID::MU))) val = 1;
      if (((iADecay->first == GenTauDecayID::MU) && (iADecay->second == GenTauDecayID::HAD)) || 
	  ((iADecay->first == GenTauDecayID::HAD) && (iADecay->second == GenTauDecayID::MU))) {
	val = 2;
      }
      if ((iADecay->first == GenTauDecayID::E) && (iADecay->second == GenTauDecayID::E)) val = 3;
      if (((iADecay->first == GenTauDecayID::E) && (iADecay->second == GenTauDecayID::HAD)) || 
	  ((iADecay->first == GenTauDecayID::HAD) && (iADecay->second == GenTauDecayID::E))) {
	val = 4;
      }
      if ((iADecay->first == GenTauDecayID::HAD) && (iADecay->second == GenTauDecayID::HAD)) {
	val = 5;
      }
      if ((iADecay - aDecay.begin()) == 0) val1 = val;
      if ((iADecay - aDecay.begin()) == 1) val2 = val;
    }
  }
  a2TauPairDecayVsA1TauPairDecay_->Fill(val1, val2);

  //plot distribution of true no. in-time interactions
  float trueNInt = -1;
  std::vector<PileupSummaryInfo>::const_iterator iPU = pPU->begin();
  int BX = 0;
  while ((iPU != pPU->end()) && (BX == 0)) {
    int BX = iPU->getBunchCrossing();
    if (BX == 0) { 
      trueNInt = iPU->getTrueNumInteractions();
      BX = -1;
    }
    ++iPU;
  }
  trueNInt_->Fill(trueNInt);
}


// ------------ method called once each job just before starting event loop  ------------
void GenAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book histograms
  dRA1TauDaughters_ = new TH1F("dRA1TauDaughters", "", 60, 0.0, 3.0);
  tauMuPT_ = new TH1F("tauMuPT", "", 50, 0.0, 100.0);
  tauHadPT_ = new TH1F("tauHadPT", "", 50, 0.0, 100.0);
  trueNInt_ = new TH1D("trueNInt", "", 60, 0.0, 60.0);
  a2TauPairDecayVsA1TauPairDecay_ = new TH2F("a2TauPairDecayVsA1TauPairDecay", 
					     ";a_{1} di-tau decay;a_{2} di-tau decay", 
					     6, -0.5, 5.5, 6, -0.5, 5.5);

  //set bin labels where appropriate
  a2TauPairDecayVsA1TauPairDecay_->GetXaxis()->SetBinLabel(1, "#tau_{#mu}#tau_{#mu}");
  a2TauPairDecayVsA1TauPairDecay_->GetXaxis()->SetBinLabel(2, "#tau_{#mu}#tau_{e}");
  a2TauPairDecayVsA1TauPairDecay_->GetXaxis()->SetBinLabel(3, "#tau_{#mu}#tau_{had}");
  a2TauPairDecayVsA1TauPairDecay_->GetXaxis()->SetBinLabel(4, "#tau_{e}#tau_{e}");
  a2TauPairDecayVsA1TauPairDecay_->GetXaxis()->SetBinLabel(5, "#tau_{e}#tau_{had}");
  a2TauPairDecayVsA1TauPairDecay_->GetXaxis()->SetBinLabel(6, "#tau_{had}#tau_{had}");
  a2TauPairDecayVsA1TauPairDecay_->GetYaxis()->SetBinLabel(1, "#tau_{#mu}#tau_{#mu}");
  a2TauPairDecayVsA1TauPairDecay_->GetYaxis()->SetBinLabel(2, "#tau_{#mu}#tau_{e}");
  a2TauPairDecayVsA1TauPairDecay_->GetYaxis()->SetBinLabel(3, "#tau_{#mu}#tau_{had}");
  a2TauPairDecayVsA1TauPairDecay_->GetYaxis()->SetBinLabel(4, "#tau_{e}#tau_{e}");
  a2TauPairDecayVsA1TauPairDecay_->GetYaxis()->SetBinLabel(5, "#tau_{e}#tau_{had}");
  a2TauPairDecayVsA1TauPairDecay_->GetYaxis()->SetBinLabel(6, "#tau_{had}#tau_{had}");
}

// ------------ method called once each job just after ending the event loop  ------------
void GenAnalyzer::endJob() 
{
  //make the canvases
  TCanvas dRA1TauDaughtersCanvas("dRA1TauDaughtersCanvas", "", 600, 600);
  Common::setCanvasOptions(dRA1TauDaughtersCanvas, 1, 0, 0);
  TCanvas tauMuPTCanvas("tauMuPTCanvas", "", 600, 600);
  Common::setCanvasOptions(tauMuPTCanvas, 1, 0, 0);
  TCanvas tauHadPTCanvas("tauHadPTCanvas", "", 600, 600);
  Common::setCanvasOptions(tauHadPTCanvas, 1, 0, 0);
  TCanvas trueNIntCanvas("trueNIntCanvas", "", 600, 600);
  Common::setCanvasOptions(trueNIntCanvas, 1, 0, 0);
  TCanvas a2TauPairDecayVsA1TauPairDecayCanvas("a2TauPairDecayVsA1TauPairDecayCanvas", "", 
					       600, 600);

  //format the plots
  Common::setHistogramOptions(dRA1TauDaughters_, kBlack, 0.7, 20, 1.0, "#DeltaR", "", 0.05);
  dRA1TauDaughters_->SetLineWidth(2);
  Common::setHistogramOptions(tauMuPT_, kBlack, 0.7, 20, 1.0, "No. interactions", "", 0.05);
  tauMuPT_->SetLineWidth(2);
  Common::setHistogramOptions(tauHadPT_, kBlack, 0.7, 20, 1.0, "No. interactions", "", 0.05);
  tauHadPT_->SetLineWidth(2);
  Common::setHistogramOptions(trueNInt_, kBlack, 0.7, 20, 1.0, "No. interactions", "", 0.05);
  trueNInt_->SetLineWidth(2);

  //draw plots
  dRA1TauDaughtersCanvas.cd();
  dRA1TauDaughters_->Draw();
  tauMuPTCanvas.cd();
  tauMuPT_->Draw();
  tauHadPTCanvas.cd();
  tauHadPT_->Draw();
  trueNIntCanvas.cd();
  trueNInt_->Draw();
  Common::draw2DHistograms(a2TauPairDecayVsA1TauPairDecayCanvas, a2TauPairDecayVsA1TauPairDecay_);

  //write output file
  out_->cd();
  dRA1TauDaughtersCanvas.Write();
  tauMuPTCanvas.Write();
  tauHadPTCanvas.Write();
  trueNIntCanvas.Write();
  a2TauPairDecayVsA1TauPairDecayCanvas.Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void GenAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void GenAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void GenAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
						    edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void GenAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
						  edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void GenAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (out_ != NULL)) delete out_;
  out_ = NULL;
  if ((doDelete) && (dRA1TauDaughters_ != NULL)) delete dRA1TauDaughters_;
  dRA1TauDaughters_ = NULL;
  if ((doDelete) && (tauMuPT_ != NULL)) delete tauMuPT_;
  tauMuPT_ = NULL;
  if ((doDelete) && (tauHadPT_ != NULL)) delete tauHadPT_;
  tauHadPT_ = NULL;
  if ((doDelete) && (trueNInt_ != NULL)) delete trueNInt_;
  trueNInt_ = NULL;
  if (doDelete && (a2TauPairDecayVsA1TauPairDecay_ != NULL)) {
    delete a2TauPairDecayVsA1TauPairDecay_;
  }
  a2TauPairDecayVsA1TauPairDecay_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
