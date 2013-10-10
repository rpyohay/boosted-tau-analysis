// -*- C++ -*-
//
// Package:    DrellYanAnalyzer
// Class:      DrellYanAnalyzer
// 
/**\class DrellYanAnalyzer DrellYanAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/DrellYanAnalyzer.cc

   Description: analyze Drell-Yan enriched sample to understand MC normalization to data

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: DrellYanAnalyzer.cc,v 1.12 2013/07/31 08:10:29 yohay Exp $
//
//

// system include files
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// class declaration
//

class DrellYanAnalyzer : public edm::EDAnalyzer {
public:
  explicit DrellYanAnalyzer(const edm::ParameterSet&);
  ~DrellYanAnalyzer();

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

  //fill pT histogram with object of arbitrary type
  template<typename U>
  void fillPTHistogramArbitrary(edm::Handle<edm::View<U> >& pView, TH1F* hist)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->pt());
  }

  //fill ET histogram
  template<typename U>
  void fillETHistogram(edm::Handle<edm::View<U> >& pView, TH1F* hist, const double weight)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->et(), weight);
  }

  //calculate MT for 2 vectors
  double MT(const reco::LeafCandidate::LorentzVector&, 
	    const reco::LeafCandidate::LorentzVector&) const;

  //plot histogram of HT for user's choice of input candidates
  void plotHT(const std::vector<reco::Candidate*>&, TH1F*, const double);

  // ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //MET tag
  edm::InputTag METTag_;

  //muon tag
  edm::InputTag muonTag_;

  //PU info tag
  edm::InputTag PUTag_;

  //vertex tag
  edm::InputTag vtxTag_;

  //muon PU subtraction coefficient
  double PUSubtractionCoeff_;

  //PU reweighting scenario
  std::string PUScenario_;

  //MC flag
  bool MC_;

  //minimum dimuon invariant mass
  double minMDimuon_;

  //minimum muon pT
  double minPTMuon_;

  //histogram of first muon pT
  TH1F* mu1PT_;

  //histogram of second muon pT
  TH1F* mu2PT_;

  //histogram of dimuon system pT
  TH1F* dimuonPT_;

  //histogram of first muon eta
  TH1F* mu1Eta_;

  //histogram of second muon eta
  TH1F* mu2Eta_;

  //histogram of dimuon system eta
  TH1F* dimuonEta_;

  //histogram of first muon phi
  TH1F* mu1Phi_;

  //histogram of second muon phi
  TH1F* mu2Phi_;

  //histogram of dimuon system phi
  TH1F* dimuonPhi_;

  //histogram of transverse mass of first muon and MET
  TH1F* mu1METMT_;

  //histogram of transverse mass of second muon and MET
  TH1F* mu2METMT_;

  //histogram of transverse mass of dimuon system and MET
  TH1F* dimuonMETMT_;

  //histogram of dPhi(first muon, MET)
  TH1F* dPhiMu1MET_;

  //histogram of dPhi(second muon, MET)
  TH1F* dPhiMu2MET_;

  //histogram of dPhi(dimuon system, MET)
  TH1F* dPhiDimuonMET_;

  //histogram of first muon isolation (PF, combined, PU-subtracted)
  TH1F* mu1Iso_;

  //histogram of second muon isolation (PF, combined, PU-subtracted)
  TH1F* mu2Iso_;

  //histogram of dEta(first muon, second muon)
  TH1F* dEtaMu1Mu2_;

  //histogram of dPhi(first muon, second muon)
  TH1F* dPhiMu1Mu2_;

  //histogram of dR(first muon, second muon)
  TH1F* dRMu1Mu2_;

  //histogram of MET
  TH1F* MET_;

  //histogram of number of good vertices (to check PU reweighting)
  TH1F* nGoodVtx_;

  //histogram of dimuon invariant mass
  TH1F* mDimuon_;

  //PU reweighting object
  edm::LumiReWeighting PUReweight_;
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
DrellYanAnalyzer::DrellYanAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  METTag_(iConfig.getParameter<edm::InputTag>("METTag")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  PUTag_(iConfig.getParameter<edm::InputTag>("PUTag")),
  vtxTag_(iConfig.getParameter<edm::InputTag>("vtxTag")),
  PUSubtractionCoeff_(iConfig.getParameter<double>("PUSubtractionCoeff")),
  PUScenario_(iConfig.getParameter<std::string>("PUScenario")),
  MC_(iConfig.getParameter<bool>("MC")),
  minMDimuon_(iConfig.getParameter<double>("minMDimuon")),
  minPTMuon_(iConfig.getParameter<double>("minPTMuon"))
{
  //now do what ever initialization is needed
  reset(false);

  //instantiate the PU reweighting object
  float Data2012PUDistArray[60] = {1.226e+04, 3.285e+04, 9.233e+04, 3.395e+05, 6.185e+05, 
				   3.050e+06, 1.772e+07, 5.414e+07, 1.305e+08, 2.590e+08, 
				   4.463e+08, 6.856e+08, 8.816e+08, 9.991e+08, 1.079e+09, 
				   1.138e+09, 1.172e+09, 1.182e+09, 1.177e+09, 1.161e+09, 
				   1.136e+09, 1.105e+09, 1.068e+09, 1.021e+09, 9.556e+08, 
				   8.671e+08, 7.587e+08, 6.389e+08, 5.164e+08, 3.999e+08, 
				   2.963e+08, 2.101e+08, 1.424e+08, 9.205e+07, 5.654e+07, 
				   3.291e+07, 1.815e+07, 9.512e+06, 4.764e+06, 2.300e+06, 
				   1.081e+06, 5.020e+05, 2.337e+05, 1.111e+05, 5.483e+04, 
				   2.840e+04, 1.549e+04, 8.845e+03, 5.236e+03, 3.180e+03, 
				   1.964e+03, 1.225e+03, 7.678e+02, 4.813e+02, 3.006e+02, 
				   1.866e+02, 1.147e+02, 6.969e+01, 4.179e+01, 2.470e+01};
  float Data201219p7InvFbPUDistArray[60] = {1.353e+04, 5.161e+04, 2.116e+06, 3.323e+05, 
					    5.569e+05, 3.645e+06, 2.024e+07, 6.093e+07, 
					    1.430e+08, 2.762e+08, 4.717e+08, 7.054e+08, 
					    8.867e+08, 9.951e+08, 1.068e+09, 1.120e+09, 
					    1.148e+09, 1.153e+09, 1.149e+09, 1.135e+09, 
					    1.111e+09, 1.082e+09, 1.050e+09, 1.009e+09, 
					    9.496e+08, 8.655e+08, 7.598e+08, 6.420e+08, 
					    5.221e+08, 4.084e+08, 3.072e+08, 2.223e+08, 
					    1.547e+08, 1.033e+08, 6.614e+07, 4.073e+07, 
					    2.438e+07, 1.450e+07, 8.859e+06, 5.779e+06, 
					    4.131e+06, 3.233e+06, 2.706e+06, 2.356e+06, 
					    2.092e+06, 1.871e+06, 1.674e+06, 1.492e+06, 
					    1.324e+06, 1.168e+06, 1.022e+06, 8.887e+05, 
					    7.663e+05, 6.553e+05, 5.556e+05, 4.668e+05, 
					    3.887e+05, 3.207e+05, 2.623e+05, 2.129e+05};
  std::vector<float> Data2012PUDist(Data2012PUDistArray, Data2012PUDistArray + 
				    sizeof(Data2012PUDistArray)/sizeof(float));
  std::vector<float> 
    Data201219p7InvFbPUDist(Data201219p7InvFbPUDistArray, Data201219p7InvFbPUDistArray + 
			   sizeof(Data201219p7InvFbPUDistArray)/sizeof(float));
  if (PUScenario_ == "S7") {
//     float S7PUDistArray[60] = {2.344E-05, 2.344E-05, 2.344E-05, 2.344E-05, 4.687E-04, 4.687E-04, 
// 			       7.032E-04, 9.414E-04, 1.234E-03, 1.603E-03, 2.464E-03, 3.250E-03, 
// 			       5.021E-03, 6.644E-03, 8.502E-03, 1.121E-02, 1.518E-02, 2.033E-02, 
// 			       2.608E-02, 3.171E-02, 3.667E-02, 4.060E-02, 4.338E-02, 4.520E-02, 
// 			       4.641E-02, 4.735E-02, 4.816E-02, 4.881E-02, 4.917E-02, 4.909E-02, 
// 			       4.842E-02, 4.707E-02, 4.501E-02, 4.228E-02, 3.896E-02, 3.521E-02, 
// 			       3.118E-02, 2.702E-02, 2.287E-02, 1.885E-02, 1.508E-02, 1.166E-02, 
// 			       8.673E-03, 6.190E-03, 4.222E-03, 2.746E-03, 1.698E-03, 9.971E-04, 
// 			       5.549E-04, 2.924E-04, 1.457E-04, 6.864E-05, 3.054E-05, 1.282E-05, 
// 			       5.081E-06, 1.898E-06, 6.688E-07, 2.221E-07, 6.947E-08, 2.047E-08};
    float S7PUDistArray[60] = {3.200e+02, 1.089e+03, 2.438e+03, 4.308e+03, 5.883e+03, 6.872e+03, 
			       7.273e+03, 7.123e+03, 6.709e+03, 6.038e+03, 5.470e+03, 5.052e+03, 
			       4.431e+03, 4.115e+03, 3.816e+03, 3.601e+03, 3.415e+03, 3.026e+03, 
			       2.860e+03, 2.550e+03, 2.369e+03, 1.997e+03, 1.776e+03, 1.522e+03, 
			       1.322e+03, 1.100e+03, 8.970e+02, 7.040e+02, 5.190e+02, 3.840e+02, 
			       3.290e+02, 2.130e+02, 1.730e+02, 9.600e+01, 7.700e+01, 4.800e+01, 
			       3.100e+01, 2.300e+01, 1.000e+01, 1.300e+01, 4.000e+00, 3.000e+00, 
			       0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 1.000e+00, 0.000e+00, 
			       0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 
			       0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00};
    std::vector<float> S7PUDist(S7PUDistArray, S7PUDistArray + 
				sizeof(S7PUDistArray)/sizeof(float));
//     PUReweight_ = edm::LumiReWeighting(S7PUDist, Data2012PUDist);
    PUReweight_ = edm::LumiReWeighting(S7PUDist, Data201219p7InvFbPUDist);
  }
  else if (PUScenario_ == "S10") {
    float S10PUDistArray[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 
				1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 
				2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 
				5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 
				4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 
				2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 
				1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 
				2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 
				5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 
				7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06};
    std::vector<float> S10PUDist(S10PUDistArray, S10PUDistArray + 
				 sizeof(S10PUDistArray)/sizeof(float));
//     PUReweight_ = edm::LumiReWeighting(S10PUDist, Data2012PUDist);
    PUReweight_ = edm::LumiReWeighting(S10PUDist, Data201219p7InvFbPUDist);
  }
}

DrellYanAnalyzer::~DrellYanAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void DrellYanAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get MET
  edm::Handle<edm::View<reco::PFMET> > pMET;
  iEvent.getByLabel(METTag_, pMET);

  //get muons
  edm::Handle<reco::MuonRefVector> pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  //get PU info
  edm::Handle<std::vector<PileupSummaryInfo> > pPU;
  if (MC_) iEvent.getByLabel(PUTag_, pPU);

  //get vertices
  edm::Handle<reco::VertexCollection> pVtx;
  iEvent.getByLabel(vtxTag_, pVtx);

  //get PU weight
  double PUWeight = 1.0;
  float trueNInt = -1;
  if (MC_ && ((PUScenario_ == "S7") || (PUScenario_ == "S10"))) {
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
    PUWeight = PUReweight_.weight(trueNInt);
  }

  //find the 2 highest pT muons
  std::vector<reco::MuonRef> muonRefs;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) { muonRefs.push_back(*iMuon); }
  Common::sortByPT(muonRefs);
  std::reverse(muonRefs.begin(), muonRefs.end());
  reco::LeafCandidate::LorentzVector dimuonSystem = muonRefs[0]->p4() + muonRefs[1]->p4();

  //apply cuts
  if ((dimuonSystem.M() > minMDimuon_) && 
      (muonRefs[0]->pt() > minPTMuon_) && (muonRefs[1]->pt() > minPTMuon_)) {

    //plot first muon pT
    mu1PT_->Fill(muonRefs[0]->pt(), PUWeight);

    //plot second muon pT
    mu2PT_->Fill(muonRefs[1]->pt(), PUWeight);

    //plot dimuon system pT
    dimuonPT_->Fill(dimuonSystem.Pt(), PUWeight);

    //plot first muon eta
    mu1Eta_->Fill(muonRefs[0]->eta(), PUWeight);

    //plot second muon eta
    mu2Eta_->Fill(muonRefs[1]->eta(), PUWeight);

    //plot dimuon system eta
    dimuonEta_->Fill(dimuonSystem.Eta(), PUWeight);

    //plot first muon phi
    mu1Phi_->Fill(muonRefs[0]->phi(), PUWeight);

    //plot second muon phi
    mu2Phi_->Fill(muonRefs[1]->phi(), PUWeight);

    //plot dimuon system phi
    dimuonPhi_->Fill(dimuonSystem.Phi(), PUWeight);

    //plot transverse mass of first muon and MET
    reco::LeafCandidate::LorentzVector MET = pMET->refAt(0)->p4();
    mu1METMT_->Fill(MT(muonRefs[0]->p4(), MET), PUWeight);

    //plot transverse mass of second muon and MET
    mu2METMT_->Fill(MT(muonRefs[1]->p4(), MET), PUWeight);

    //plot transverse mass of dimuon system and MET
    dimuonMETMT_->Fill(MT(dimuonSystem, MET), PUWeight);

    //plot dPhi(first muon, MET)
    dPhiMu1MET_->Fill(fabs(reco::deltaPhi(muonRefs[0]->phi(), MET.Phi())), PUWeight);

    //plot dPhi(second muon, MET)
    dPhiMu2MET_->Fill(fabs(reco::deltaPhi(muonRefs[1]->phi(), MET.Phi())), PUWeight);

    //plot dPhi(dimuon system, MET)
    dPhiDimuonMET_->Fill(fabs(reco::deltaPhi(dimuonSystem.Phi(), MET.Phi())), PUWeight);

    //plot first muon isolation (PF, combined, PU-subtracted)
    mu1Iso_->Fill(Common::getMuonCombPFIso(*muonRefs[0], PUSubtractionCoeff_), PUWeight);

    //plot second muon isolation (PF, combined, PU-subtracted)
    mu2Iso_->Fill(Common::getMuonCombPFIso(*muonRefs[1], PUSubtractionCoeff_), PUWeight);

    //plot dEta(first muon, second muon)
    dEtaMu1Mu2_->Fill(muonRefs[0]->eta() - muonRefs[1]->eta(), PUWeight);

    //plot dPhi(first muon, second muon)
    dPhiMu1Mu2_->Fill(fabs(reco::deltaPhi(muonRefs[0]->phi(), muonRefs[1]->phi())), PUWeight);

    //plot dR(first muon, second muon)
    dRMu1Mu2_->Fill(reco::deltaR(*muonRefs[0], *muonRefs[1]), PUWeight);

    //plot MET distribution
    fillETHistogram(pMET, MET_, PUWeight);  

    //plot the number of good vertices
    nGoodVtx_->Fill(Common::numGoodVertices(pVtx), PUWeight);

    //plot dimuon invariant mass
    mDimuon_->Fill(dimuonSystem.M());
  }
}


// ------------ method called once each job just before starting event loop  ------------
void DrellYanAnalyzer::beginJob()
{
  //open output files
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book histograms
  mu1PT_ = new TH1F("mu1PT", ";p_{T} (GeV);", 40, 0.0, 200.0);
  mu2PT_ = new TH1F("mu2PT", ";p_{T} (GeV);", 40, 0.0, 200.0);
  dimuonPT_ = new TH1F("dimuonPT", ";p_{T} (GeV);", 40, 0.0, 200.0);
  mu1Eta_ = new TH1F("mu1Eta", ";#eta;", 42, -2.1, 2.1);
  mu2Eta_ = new TH1F("mu2Eta", ";#eta;", 42, -2.1, 2.1);
  dimuonEta_ = new TH1F("dimuonEta", ";#eta;", 42, -2.1, 2.1);
  mu1Phi_ = new TH1F("mu1Phi", ";#phi;", 60, -3.15, 3.15);
  mu2Phi_ = new TH1F("mu2Phi", ";#phi;", 60, -3.15, 3.15);
  dimuonPhi_ = new TH1F("dimuonPhi", ";#phi;", 60, -3.15, 3.15);
  mu1METMT_ = new TH1F("mu1METMT", ";M_{T} (GeV);", 50, 0.0, 200.0);
  mu2METMT_ = new TH1F("mu2METMT", ";M_{T} (GeV);", 50, 0.0, 200.0);
  dimuonMETMT_ = new TH1F("dimuonMETMT", ";M_{T} (GeV);", 50, 0.0, 200.0);
  dPhiMu1MET_ = new TH1F("dPhiMu1MET", ";#Delta#phi;", 30, 0.0, 3.15);
  dPhiMu2MET_ = new TH1F("dPhiMu2MET", ";#Delta#phi;", 30, 0.0, 3.15);
  dPhiDimuonMET_ = new TH1F("dPhiDimuonMET", ";#Delta#phi;", 30, 0.0, 3.15);
  mu1Iso_ = new TH1F("mu1Iso", ";Isolation energy (GeV);", 20, 0.0, 20.0);
  mu2Iso_ = new TH1F("mu2Iso", ";Isolation energy (GeV);", 20, 0.0, 20.0);
  dEtaMu1Mu2_ = new TH1F("dEtaMu1Mu2", ";#Delta#eta;", 42, -4.2, 4.2);
  dPhiMu1Mu2_ = new TH1F("dPhiMu1Mu2", ";#Delta#phi;", 30, 0.0, 3.15);
  dRMu1Mu2_ = new TH1F("dRMu1Mu2", ";#DeltaR;", 35, 0.0, 5.25);
  MET_ = new TH1F("MET", ";#slash{E}_{T} (GeV);", 40, 0.0, 200.0);
  nGoodVtx_ = new TH1F("nGoodVtx", ";N_{PV};", 60, -0.5, 59.5);
  mDimuon_ = new TH1F("mDimuon", ";m_{#mu#mu} (GeV);", 36, 0.0, 180.0);

  //set sumw2
  mu1PT_->Sumw2();
  mu2PT_->Sumw2();
  dimuonPT_->Sumw2();
  mu1Eta_->Sumw2();
  mu2Eta_->Sumw2();
  dimuonEta_->Sumw2();
  mu1Phi_->Sumw2();
  mu2Phi_->Sumw2();
  dimuonPhi_->Sumw2();
  mu1METMT_->Sumw2();
  mu2METMT_->Sumw2();
  dimuonMETMT_->Sumw2();
  dPhiMu1MET_->Sumw2();
  dPhiMu2MET_->Sumw2();
  dPhiDimuonMET_->Sumw2();
  mu1Iso_->Sumw2();
  mu2Iso_->Sumw2();
  dEtaMu1Mu2_->Sumw2();
  dPhiMu1Mu2_->Sumw2();
  dRMu1Mu2_->Sumw2();
  MET_->Sumw2();
  nGoodVtx_->Sumw2();
  mDimuon_->Sumw2();
}

// ------------ method called once each job just after ending the event loop  ------------
void DrellYanAnalyzer::endJob() 
{
  //make canvases
  out_->cd();
  TCanvas mu1PTCanvas("mu1PTCanvas", "", 600, 600);
  TCanvas mu2PTCanvas("mu2PTCanvas", "", 600, 600);
  TCanvas dimuonPTCanvas("dimuonPTCanvas", "", 600, 600);
  TCanvas mu1EtaCanvas("mu1EtaCanvas", "", 600, 600);
  TCanvas mu2EtaCanvas("mu2EtaCanvas", "", 600, 600);
  TCanvas dimuonEtaCanvas("dimuonEtaCanvas", "", 600, 600);
  TCanvas mu1PhiCanvas("mu1PhiCanvas", "", 600, 600);
  TCanvas mu2PhiCanvas("mu2PhiCanvas", "", 600, 600);
  TCanvas dimuonPhiCanvas("dimuonPhiCanvas", "", 600, 600);
  TCanvas mu1METMTCanvas("mu1METMTCanvas", "", 600, 600);
  TCanvas mu2METMTCanvas("mu2METMTCanvas", "", 600, 600);
  TCanvas dimuonMETMTCanvas("dimuonMETMTCanvas", "", 600, 600);
  TCanvas dPhiMu1METCanvas("dPhiMu1METCanvas", "", 600, 600);
  TCanvas dPhiMu2METCanvas("dPhiMu2METCanvas", "", 600, 600);
  TCanvas dPhiDimuonMETCanvas("dPhiDimuonMETCanvas", "", 600, 600);
  TCanvas mu1IsoCanvas("mu1IsoCanvas", "", 600, 600);
  TCanvas mu2IsoCanvas("mu2IsoCanvas", "", 600, 600);
  TCanvas dEtaMu1Mu2Canvas("dEtaMu1Mu2Canvas", "", 600, 600);
  TCanvas dPhiMu1Mu2Canvas("dPhiMu1Mu2Canvas", "", 600, 600);
  TCanvas dRMu1Mu2Canvas("dRMu1Mu2Canvas", "", 600, 600);
  TCanvas METCanvas("METCanvas", "", 600, 600);
  TCanvas nGoodVtxCanvas("nGoodVtxCanvas", "", 600, 600);
  TCanvas mDimuonCanvas("mDimuonCanvas", "", 600, 600);

  //format and draw 1D plots
  Common::draw1DHistograms(mu1PTCanvas, mu1PT_);
  Common::draw1DHistograms(mu2PTCanvas, mu2PT_);
  Common::draw1DHistograms(dimuonPTCanvas, dimuonPT_);
  Common::draw1DHistograms(mu1EtaCanvas, mu1Eta_);
  Common::draw1DHistograms(mu2EtaCanvas, mu2Eta_);
  Common::draw1DHistograms(dimuonEtaCanvas, dimuonEta_);
  Common::draw1DHistograms(mu1PhiCanvas, mu1Phi_);
  Common::draw1DHistograms(mu2PhiCanvas, mu2Phi_);
  Common::draw1DHistograms(dimuonPhiCanvas, dimuonPhi_);
  Common::draw1DHistograms(mu1METMTCanvas, mu1METMT_);
  Common::draw1DHistograms(mu2METMTCanvas, mu2METMT_);
  Common::draw1DHistograms(dimuonMETMTCanvas, dimuonMETMT_);
  Common::draw1DHistograms(dPhiMu1METCanvas, dPhiMu1MET_);
  Common::draw1DHistograms(dPhiMu2METCanvas, dPhiMu2MET_);
  Common::draw1DHistograms(dPhiDimuonMETCanvas, dPhiDimuonMET_);
  Common::draw1DHistograms(mu1IsoCanvas, mu1Iso_);
  Common::draw1DHistograms(mu2IsoCanvas, mu2Iso_);
  Common::draw1DHistograms(dEtaMu1Mu2Canvas, dEtaMu1Mu2_);
  Common::draw1DHistograms(dPhiMu1Mu2Canvas, dPhiMu1Mu2_);
  Common::draw1DHistograms(dRMu1Mu2Canvas, dRMu1Mu2_);
  Common::draw1DHistograms(METCanvas, MET_);
  Common::draw1DHistograms(nGoodVtxCanvas, nGoodVtx_);
  Common::draw1DHistograms(mDimuonCanvas, mDimuon_);

  //write output files
  out_->cd();
  mu1PTCanvas.Write();
  mu2PTCanvas.Write();
  dimuonPTCanvas.Write();
  mu1EtaCanvas.Write();
  mu2EtaCanvas.Write();
  dimuonEtaCanvas.Write();
  mu1PhiCanvas.Write();
  mu2PhiCanvas.Write();
  dimuonPhiCanvas.Write();
  mu1METMTCanvas.Write();
  mu2METMTCanvas.Write();
  dimuonMETMTCanvas.Write();
  dPhiMu1METCanvas.Write();
  dPhiMu2METCanvas.Write();
  dPhiDimuonMETCanvas.Write();
  mu1IsoCanvas.Write();
  mu2IsoCanvas.Write();
  dEtaMu1Mu2Canvas.Write();
  dPhiMu1Mu2Canvas.Write();
  dRMu1Mu2Canvas.Write();
  METCanvas.Write();
  nGoodVtxCanvas.Write();
  mDimuonCanvas.Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void DrellYanAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void DrellYanAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void DrellYanAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
				       edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void DrellYanAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
				     edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void DrellYanAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

double DrellYanAnalyzer::MT(const reco::LeafCandidate::LorentzVector& p41, 
			    const reco::LeafCandidate::LorentzVector& p42) const
{
  return sqrt(2*p41.Pt()*p42.Et()*(1.0 - cos(reco::deltaPhi(p41.Phi(), p42.Phi()))));
}
void DrellYanAnalyzer::plotHT(const std::vector<reco::Candidate*>& cands, TH1F* hist, 
			 const double weight)
{
  double HT = 0.0;
  for (std::vector<reco::Candidate*>::const_iterator iCand = cands.begin(); iCand != cands.end(); 
       ++iCand) {
    HT+=(*iCand)->pt();
  }
  hist->Fill(HT, weight);
}

void DrellYanAnalyzer::reset(const bool doDelete)
{
  if (doDelete && (out_ != NULL)) delete out_;
  out_ = NULL;
  if (doDelete && (mu1PT_ != NULL)) delete mu1PT_;
  mu1PT_ = NULL;
  if (doDelete && (mu2PT_ != NULL)) delete mu2PT_;
  mu2PT_ = NULL;
  if (doDelete && (dimuonPT_ != NULL)) delete dimuonPT_;
  dimuonPT_ = NULL;
  if (doDelete && (mu1Eta_ != NULL)) delete mu1Eta_;
  mu1Eta_ = NULL;
  if (doDelete && (mu2Eta_ != NULL)) delete mu2Eta_;
  mu2Eta_ = NULL;
  if (doDelete && (dimuonEta_ != NULL)) delete dimuonEta_;
  dimuonEta_ = NULL;
  if (doDelete && (mu1Phi_ != NULL)) delete mu1Phi_;
  mu1Phi_ = NULL;
  if (doDelete && (mu2Phi_ != NULL)) delete mu2Phi_;
  mu2Phi_ = NULL;
  if (doDelete && (dimuonPhi_ != NULL)) delete dimuonPhi_;
  dimuonPhi_ = NULL;
  if (doDelete && (mu1METMT_ != NULL)) delete mu1METMT_;
  mu1METMT_ = NULL;
  if (doDelete && (mu2METMT_ != NULL)) delete mu2METMT_;
  mu2METMT_ = NULL;
  if (doDelete && (dimuonMETMT_ != NULL)) delete dimuonMETMT_;
  dimuonMETMT_ = NULL;
  if (doDelete && (dPhiMu1MET_ != NULL)) delete dPhiMu1MET_;
  dPhiMu1MET_ = NULL;
  if (doDelete && (dPhiMu2MET_ != NULL)) delete dPhiMu2MET_;
  dPhiMu2MET_ = NULL;
  if (doDelete && (dPhiDimuonMET_ != NULL)) delete dPhiDimuonMET_;
  dPhiDimuonMET_ = NULL;
  if (doDelete && (mu1Iso_ != NULL)) delete mu1Iso_;
  mu1Iso_ = NULL;
  if (doDelete && (mu2Iso_ != NULL)) delete mu2Iso_;
  mu2Iso_ = NULL;
  if (doDelete && (dEtaMu1Mu2_ != NULL)) delete dEtaMu1Mu2_;
  dEtaMu1Mu2_ = NULL;
  if (doDelete && (dPhiMu1Mu2_ != NULL)) delete dPhiMu1Mu2_;
  dPhiMu1Mu2_ = NULL;
  if (doDelete && (dRMu1Mu2_ != NULL)) delete dRMu1Mu2_;
  dRMu1Mu2_ = NULL;
  if (doDelete && (MET_ != NULL)) delete MET_;
  MET_ = NULL;
  if (doDelete && (nGoodVtx_ != NULL)) delete nGoodVtx_;
  nGoodVtx_ = NULL;
  if (doDelete && (mDimuon_ != NULL)) delete mDimuon_;
  mDimuon_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(DrellYanAnalyzer);
