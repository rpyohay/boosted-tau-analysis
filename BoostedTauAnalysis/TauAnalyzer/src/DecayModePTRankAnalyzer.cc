// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      DecayModePTRankAnalyzer
// 
/**\class DecayModePTRankAnalyzer DecayModePTRankAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/DecayModePTRankAnalyzer.cc

   Description: analyze tau decays split by decay mode and pT rank

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: DecayModePTRankAnalyzer.cc,v 1.2 2012/09/25 11:49:24 yohay Exp $
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

class DecayModePTRankAnalyzer : public edm::EDAnalyzer {
public:
  explicit DecayModePTRankAnalyzer(const edm::ParameterSet&);
  ~DecayModePTRankAnalyzer();

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

  //fill visible pT histogram
  void fillPTHistogram(const edm::Handle<reco::GenParticleRefVector>&, TH1F*, 
		       edm::Handle<reco::GenParticleCollection>&);

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
DecayModePTRankAnalyzer::DecayModePTRankAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  tauMuInputTags_(iConfig.getParameter<std::vector<edm::InputTag> >("tauMuInputTags")),
  tau1ProngInputTags_(iConfig.getParameter<std::vector<edm::InputTag> >("tau1ProngInputTags")),
  tau1Prong1Pi0InputTags_(iConfig.getParameter<std::vector<edm::InputTag> >
			  ("tau1Prong1Pi0InputTags")),
  tau1Prong2Pi0InputTags_(iConfig.getParameter<std::vector<edm::InputTag> >
			  ("tau1Prong2Pi0InputTags")),
  tau3ProngInputTags_(iConfig.getParameter<std::vector<edm::InputTag> >("tau3ProngInputTags")),
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
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

DecayModePTRankAnalyzer::~DecayModePTRankAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void DecayModePTRankAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get gen object collections split by decay type and pT rank and associate each to a pT histogram
  std::map<TH1F*, edm::Handle<reco::GenParticleRefVector> > genParticleCollMap;
  for (std::vector<edm::InputTag>::const_iterator iTag = tauMuInputTags_.begin(); 
       iTag != tauMuInputTags_.end(); ++iTag) {
    edm::Handle<reco::GenParticleRefVector> pRefVector;
    iEvent.getByLabel(*iTag, pRefVector);
    genParticleCollMap[tauMuPTHists_[iTag - tauMuInputTags_.begin()]] = pRefVector;
  }
  for (std::vector<edm::InputTag>::const_iterator iTag = tau1ProngInputTags_.begin(); 
       iTag != tau1ProngInputTags_.end(); ++iTag) {
    edm::Handle<reco::GenParticleRefVector> pRefVector;
    iEvent.getByLabel(*iTag, pRefVector);
    genParticleCollMap[tau1ProngPTHists_[iTag - tau1ProngInputTags_.begin()]] = pRefVector;
  }
  for (std::vector<edm::InputTag>::const_iterator iTag = tau1Prong1Pi0InputTags_.begin(); 
       iTag != tau1Prong1Pi0InputTags_.end(); ++iTag) {
    edm::Handle<reco::GenParticleRefVector> pRefVector;
    iEvent.getByLabel(*iTag, pRefVector);
    genParticleCollMap[tau1Prong1Pi0PTHists_[iTag - tau1Prong1Pi0InputTags_.begin()]] = pRefVector;
  }
  for (std::vector<edm::InputTag>::const_iterator iTag = tau1Prong2Pi0InputTags_.begin(); 
       iTag != tau1Prong2Pi0InputTags_.end(); ++iTag) {
    edm::Handle<reco::GenParticleRefVector> pRefVector;
    iEvent.getByLabel(*iTag, pRefVector);
    genParticleCollMap[tau1Prong2Pi0PTHists_[iTag - tau1Prong2Pi0InputTags_.begin()]] = pRefVector;
  }
  for (std::vector<edm::InputTag>::const_iterator iTag = tau3ProngInputTags_.begin(); 
       iTag != tau3ProngInputTags_.end(); ++iTag) {
    edm::Handle<reco::GenParticleRefVector> pRefVector;
    iEvent.getByLabel(*iTag, pRefVector);
    genParticleCollMap[tau3ProngPTHists_[iTag - tau3ProngInputTags_.begin()]] = pRefVector;
  }

  //get base gen particle collection
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //fill pT histograms for gen objects, 1 per decay type and pT rank
  fillPTHistograms(genParticleCollMap, pGenParticles);
}


// ------------ method called once each job just before starting event loop  ------------
void DecayModePTRankAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book gen particle pT histograms split by decay type and pT rank
  for (std::vector<edm::InputTag>::const_iterator iTag = tauMuInputTags_.begin(); 
       iTag != tauMuInputTags_.end(); ++iTag) {
    tauMuPTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
				     .c_str(), "", 20, 0.0, 100.0));
  }
  for (std::vector<edm::InputTag>::const_iterator iTag = tau1ProngInputTags_.begin(); 
       iTag != tau1ProngInputTags_.end(); ++iTag) {
    tau1ProngPTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
					 .c_str(), "", 20, 0.0, 100.0));
  }
  for (std::vector<edm::InputTag>::const_iterator iTag = tau1Prong1Pi0InputTags_.begin(); 
       iTag != tau1Prong1Pi0InputTags_.end(); ++iTag) {
    tau1Prong1Pi0PTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
					     .c_str(), "", 20, 0.0, 100.0));
  }
  for (std::vector<edm::InputTag>::const_iterator iTag = tau1Prong2Pi0InputTags_.begin(); 
       iTag != tau1Prong2Pi0InputTags_.end(); ++iTag) {
    tau1Prong2Pi0PTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
					     .c_str(), "", 20, 0.0, 100.0));
  }
  for (std::vector<edm::InputTag>::const_iterator iTag = tau3ProngInputTags_.begin(); 
       iTag != tau3ProngInputTags_.end(); ++iTag) {
    tau3ProngPTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
					 .c_str(), "", 20, 0.0, 100.0));
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void DecayModePTRankAnalyzer::endJob() 
{
  //make the gen particle pT canvases, 1 per decay mode
  out_->cd();
  TCanvas tauMuPTRankCanvas("tauMuPTRankCanvas", "", 600, 600);
  TCanvas tau1ProngPTRankCanvas("tau1ProngPTRankCanvas", "", 600, 600);
  TCanvas tau1Prong1Pi0PTRankCanvas("tau1Prong1Pi0PTRankCanvas", "", 600, 600);
  TCanvas tau1Prong2Pi0PTRankCanvas("tau1Prong2Pi0PTRankCanvas", "", 600, 600);
  TCanvas tau3ProngPTRankCanvas("tau3ProngPTRankCanvas", "", 600, 600);
  TLegend tauMuPTRankLegend(0.4, 0.6, 0.8, 0.8);
  TLegend tau1ProngPTRankLegend(0.4, 0.6, 0.8, 0.8);
  TLegend tau1Prong1Pi0PTRankLegend(0.4, 0.6, 0.8, 0.8);
  TLegend tau1Prong2Pi0PTRankLegend(0.4, 0.6, 0.8, 0.8);
  TLegend tau3ProngPTRankLegend(0.4, 0.6, 0.8, 0.8);
  makePTRankCanvas(tauMuPTRankCanvas, tauMuPTRankLegend, 
		   "gg fusion NMSSM Higgs #tau_{#mu}", tauMuPTHists_);
  makePTRankCanvas(tau1ProngPTRankCanvas, tau1ProngPTRankLegend, 
		   "gg fusion NMSSM Higgs #tau_{had}, 1 prong", tau1ProngPTHists_);
  makePTRankCanvas(tau1Prong1Pi0PTRankCanvas, tau1Prong1Pi0PTRankLegend, 
		   "gg fusion NMSSM Higgs #tau_{had}, 1 prong + 1 #pi^{0}", tau1Prong1Pi0PTHists_);
  makePTRankCanvas(tau1Prong2Pi0PTRankCanvas, tau1Prong2Pi0PTRankLegend, 
		   "gg fusion NMSSM Higgs #tau_{had}, 1 prong + 2 #pi^{0}", tau1Prong2Pi0PTHists_);
  makePTRankCanvas(tau3ProngPTRankCanvas, tau3ProngPTRankLegend, 
		   "gg fusion NMSSM Higgs #tau_{had}, 3 prong", tau3ProngPTHists_);

  //make the gen particle pT canvases, 1 per pT rank
  TCanvas pTRank0DecayModeCanvas("pTRank0DecayModeCanvas", "", 600, 600);
  TCanvas pTRank1DecayModeCanvas("pTRank1DecayModeCanvas", "", 600, 600);
  TCanvas pTRank2DecayModeCanvas("pTRank2DecayModeCanvas", "", 600, 600);
  TCanvas pTRank3DecayModeCanvas("pTRank3DecayModeCanvas", "", 600, 600);
  TLegend pTRank0DecayModeLegend(0.4, 0.6, 0.8, 0.8);
  TLegend pTRank1DecayModeLegend(0.4, 0.6, 0.8, 0.8);
  TLegend pTRank2DecayModeLegend(0.4, 0.6, 0.8, 0.8);
  TLegend pTRank3DecayModeLegend(0.4, 0.6, 0.8, 0.8);
  makeDecayModeCanvas(pTRank0DecayModeCanvas, pTRank0DecayModeLegend, 
		      "gg fusion NMSSM Higgs highest p_{T} object", 0);
  makeDecayModeCanvas(pTRank1DecayModeCanvas, pTRank1DecayModeLegend, 
		      "gg fusion NMSSM Higgs second highest p_{T} object", 1);
  makeDecayModeCanvas(pTRank2DecayModeCanvas, pTRank2DecayModeLegend, 
		      "gg fusion NMSSM Higgs third highest p_{T} object", 2);
  makeDecayModeCanvas(pTRank3DecayModeCanvas, pTRank3DecayModeLegend, 
		      "gg fusion NMSSM Higgs lowest p_{T} object", 3);

  //write output file
  out_->cd();
  for (std::vector<TH1F*>::iterator iHist = tauMuPTHists_.begin(); 
       iHist != tauMuPTHists_.end(); ++iHist) { (*iHist)->Write(); }
  for (std::vector<TH1F*>::iterator iHist = tau1ProngPTHists_.begin(); 
       iHist != tau1ProngPTHists_.end(); ++iHist) { (*iHist)->Write(); }
  for (std::vector<TH1F*>::iterator iHist = tau1Prong1Pi0PTHists_.begin(); 
       iHist != tau1Prong1Pi0PTHists_.end(); ++iHist) { (*iHist)->Write(); }
  for (std::vector<TH1F*>::iterator iHist = tau1Prong2Pi0PTHists_.begin(); 
       iHist != tau1Prong2Pi0PTHists_.end(); ++iHist) { (*iHist)->Write(); }
  for (std::vector<TH1F*>::iterator iHist = tau3ProngPTHists_.begin(); 
       iHist != tau3ProngPTHists_.end(); ++iHist) { (*iHist)->Write(); }
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void DecayModePTRankAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void DecayModePTRankAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void DecayModePTRankAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
						   edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void DecayModePTRankAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
						 edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void DecayModePTRankAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void DecayModePTRankAnalyzer::fillPTHistogram(const edm::Handle<reco::GenParticleRefVector>& 
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
    catch (std::string& ex) { throw cms::Exception("DecayModePTRankAnalyzer") << ex; }
  }
}

void DecayModePTRankAnalyzer::fillPTHistograms(std::map<TH1F*, 
					       edm::Handle<reco::GenParticleRefVector> >& map, 
					       edm::Handle<reco::GenParticleCollection>& 
					       pGenParticles)
{
  for (std::map<TH1F*, edm::Handle<reco::GenParticleRefVector> >::iterator i = map.begin(); 
       i != map.end(); ++i) { fillPTHistogram(i->second, i->first, pGenParticles); }
}

void DecayModePTRankAnalyzer::makePTRankCanvas(TCanvas& canvas, TLegend& legend, 
					       const std::string& header, 
					       std::vector<TH1F*>& hists)
{
  drawMultiplePTHistograms(canvas, hists, pTRankColors_, pTRankStyles_, legend, 
			   pTRankEntries_, header);
}

void DecayModePTRankAnalyzer::makeDecayModeCanvas(TCanvas& canvas, TLegend& legend, 
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

void DecayModePTRankAnalyzer::drawMultiplePTHistograms(TCanvas& canvas, 
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
    (*iHist)->SetLineWidth(2);
    legend.AddEntry(*iHist, entries[i].c_str(), "l");
    Int_t histMaxBin = (*iHist)->GetMaximumBin();
    Double_t histMaxBinContent = (*iHist)->GetBinContent(histMaxBin);
    if (histMaxBinContent >= maxBinContent) {
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

void DecayModePTRankAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (out_ != NULL)) delete out_;
  out_ = NULL;
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
}

//define this as a plug-in
DEFINE_FWK_MODULE(DecayModePTRankAnalyzer);
