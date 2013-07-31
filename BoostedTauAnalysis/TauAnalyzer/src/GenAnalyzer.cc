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
#include "TFile.h"
#include "TH1F.h"
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

  //set of parameters for GenTauDecayID class
  edm::ParameterSet genTauDecayIDPSet_;

  //histogram of dR between gen objects from a1 decay
  TH1F* dRA1TauDaughters_;
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

  //find a1 tau decay products
  std::vector<GenTauDecayID> aDecayProducts;
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {
    try {
      GenTauDecayID tauDecay(genTauDecayIDPSet_, pGenParticles, 
			     iGenParticle - pGenParticles->begin());
      if (tauDecay.isStatus3DecayProduct()) aDecayProducts.push_back(tauDecay);
    }
    catch (std::string& ex) { throw cms::Exception("GenAnalyzer") << ex; }
  }

  //loop over a1 tau daughters
  std::vector<unsigned int> keysToIgnore;
  for (std::vector<GenTauDecayID>::iterator iTau = aDecayProducts.begin(); 
       iTau != aDecayProducts.end(); ++iTau) {
    try {
      const unsigned int tauKey = iTau->getTauIndex();

      //find sister
      iTau->findSister();
      const unsigned int iSister = iTau->getSisterIndex();

      //if sister wasn't already looped over, plot dR(sisters)
      if (std::find(keysToIgnore.begin(), keysToIgnore.end(), iSister) == keysToIgnore.end()) {
	dRA1TauDaughters_->Fill(reco::deltaR(*reco::GenParticleRef(pGenParticles, tauKey), 
					     *reco::GenParticleRef(pGenParticles, iSister)));

	//ignore this tau in the future
	keysToIgnore.push_back(tauKey);
      }
    }
    catch (std::string& ex) { throw cms::Exception("GenAnalyzer") << ex; }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void GenAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book histograms
  dRA1TauDaughters_ = new TH1F("dRA1TauDaughters", "", 60, 0.0, 3.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void GenAnalyzer::endJob() 
{
  //make the canvases
  TCanvas dRA1TauDaughtersCanvas("dRA1TauDaughtersCanvas", "", 600, 600);
  Common::setCanvasOptions(dRA1TauDaughtersCanvas, 1, 0, 0);

  //format the plots
  Common::setHistogramOptions(dRA1TauDaughters_, kBlack, 0.7, 20, 1.0, "#DeltaR", "", 0.05);
  dRA1TauDaughters_->SetLineWidth(2);

  //draw plots
  dRA1TauDaughtersCanvas.cd();
  dRA1TauDaughters_->Draw();

  //write output file
  out_->cd();
  dRA1TauDaughtersCanvas.Write();
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
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
