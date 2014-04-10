// -*- C++ -*-

// Package:    TauAnalyzer
// Class:      TauAnalyzer

// \class TauAnalyzer TauAnalyzer.cc 
//    BoostedTauAnalysis/TauAnalyzer/src/TauAnalyzer.cc

//    Description: analyze tau variables that can discriminate signal to background

//    Implementation:
//    [Notes on implementation]

// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: TauAnalyzer.cc,v 1.12 2013/07/31 08:10:29 yohay Exp $



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
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Pruner.hh"
#include "BoostedTauAnalysis/TauAnalyzer/interface/Nsubjettiness.h"
#include "BoostedTauAnalysis/TauAnalyzer/interface/Njettiness.hh"

using namespace std;
using namespace edm;
using namespace reco;
using namespace fastjet;


// class declaration


class TauAnalyzer : public edm::EDAnalyzer {
public:
  explicit TauAnalyzer(const edm::ParameterSet&);
  ~TauAnalyzer();

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

  //calculate HT for a given set of candidates
  double HT(const std::vector<reco::Candidate*>&) const;

  //make visible gen tau pT canvas for 1 decay mode, multiple pT ranks
  void makePTRankCanvas(TCanvas&, TLegend&, const std::string&, std::vector<TH1F*>&);

  //format and draw multiple pT histograms on one canvas
  void drawMultiplePTHistograms(TCanvas&, std::vector<TH1F*>&, const std::vector<unsigned int>&, 
				const std::vector<unsigned int>&, TLegend&, 
				const std::vector<std::string>&, const std::string&);

  //plot histogram of dPhi(muon, MET) for user's choice of muon
  void plotDPhiMuMet(const reco::MuonRef&, const edm::Handle<edm::View<reco::PFMET> >&, TH1F*, 
		     const double);

  //plot histogram of MT for user's choice of muon
  void plotMT(const reco::MuonRef&, const edm::Handle<edm::View<reco::PFMET> >&, TH1F*, 
	      const double);

  //plot histogram of HT for user's choice of input candidates
  void plotHT(const std::vector<reco::Candidate*>&, TH1F*, const double);

//   ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //tau tag
  edm::InputTag tauTag_;

  //MET tag
  edm::InputTag METTag_;

  //muon tag
  edm::InputTag muonTag_;

  //gen-matched muon tag
  edm::InputTag genMatchedMuonTag_;

  //old jet tag
  edm::InputTag oldJetTag_;

  //new jet tag
  edm::InputTag newJetTag_;

  //jet-muon map tag
  edm::InputTag jetMuonMapTag_;

  //old-new jet map tag
  edm::InputTag oldNewJetMapTag_;

  //gen particle tag
  edm::InputTag genParticleTag_;

  //gen tau-->mu tag
  edm::InputTag genTauMuTag_;

  //gen W-->tau-->mu tag
  edm::InputTag genWTauMuTag_;

  //hadronic tau deltaBeta-corrected isolation energy tag
  edm::InputTag tauHadIsoTag_;

  //all muons tag
  edm::InputTag allMuonTag_;

  //muon gen particle tag
  edm::InputTag muonGenParticleTag_;

  //PU subtraction coefficient for muon PF isolation
  double muonPFIsoPUSubtractionCoeff_;

  //PU info tag
  edm::InputTag PUTag_;

  //vertex tag
  edm::InputTag vtxTag_;

  //all gen particles tag
  edm::InputTag allGenParticleTag_;

  //selected corrected jet tag
  edm::InputTag corrJetTag_;

  //b jet tag
  edm::InputTag bJetTag_;

  //dR matching distance
  double dR_;

  //minimum tau pT
  double tauPTMin_;

  //tau decay mode
  reco::PFTau::hadronicDecayMode tauDecayMode_;

  //minimum uncorrected jet pT
  double uncorrJetPTMin_;

  //tau arbitration method when there are multiple selected taus per event
  std::string tauArbitrationMethod_;

  //PU reweighting scenario
  std::string PUScenario_;

  //zCut for jet pruning
  double zCut_;

  //RcutFactor for jet pruning
  double RcutFactor_;

  //MC flag
  bool MC_;

  //reweight flag
  bool reweight_;

  //marker colors for histograms with different pT rank
  std::vector<unsigned int> pTRankColors_;

  //marker styles for histograms with different pT rank
  std::vector<unsigned int> pTRankStyles_;

  //legend entries for histograms with different pT rank
  std::vector<std::string> pTRankEntries_;

  //histogram of MET
  TH1F* MET_;

  //histograms of selected tau associated muon multiplicity
  TH1F* hadTauAssociatedMuMultiplicity_;

  //histogram of mu+had mass
  TH1F* muHadMass_;

  /*histogram of mu+had mass weighted by the stat. error on hadronic tau pT weights (i.e. each 
    event is filled with a weight of sigma_w, where sigma_w is the stat. error on the hadronic tau 
    pT weight for that event)*/
  TH1F* muHadMassReweightErrSq_;

  //histogram of mu charge + had charge
  TH1F* muHadCharge_;

  //histogram of b-tagging discriminators for the b-tagged jets used to reconstruct the tau_h
  TH1F* bTagDiscrim_;

  //histogram of W muon isolation
  TH1F* WMuIso_;

  //histogram of W muon transverse mass
  TH1F* WMuMT_;

  //histogram of tau muon transverse mass
  TH1F* tauMuMT_;

  //histogram of hadronic tau transverse mass
  TH1F* tauHadMT_;

  //histogram of dPhi(W muon, MET)
  TH1F* dPhiWMuMET_;

  //histogram of dPhi(tau muon, MET)
  TH1F* dPhiTauMuMET_;

  //histogram of HT (tau muon + hadronic tau + leading distinct corrected jet)
  TH1F* tauMuTauHadJetHT_;

  //histogram of HT (two leading corrected jets)
  TH1F* diJetHT_;

  //histogram of HT (corrected jet associated to hadronic tau + leading distinct corrected jet)
  TH1F* jetTauJetHT_;

  //histogram of HT (tau muon + hadronic tau + leading distinct corrected jet + W muon)
  TH1F* tauMuTauHadJetWMuHT_;

  //histogram of HT (tau muon + hadronic tau + leading distinct corrected jet + W muon + MET)
  TH1F* tauMuTauHadJetWMuMETHT_;

  //histogram of HT (two leading corrected jets + W muon)
  TH1F* diJetWMuHT_;

  /*histogram of HT (corrected jet associated to hadronic tau + leading distinct corrected jet + W 
    muon)*/
  TH1F* jetTauJetWMuHT_;

  //histogram of dR(W muon, leading soft muon per jet)
  TH1F* dRWMuSoftMu_;

  //histogram of dPhi(W muon, leading soft muon per jet)
  TH1F* dPhiWMuSoftMu_;

  //histogram of dPhi(W muon, leading soft muon per jet)
  TH1F* dPhiWMuSoftMu_withCut_;

  //histogram of dR(W muon, leading soft gen-matched muon per jet)
  TH1F* dRWMuSoftGenMatchedMu_;

  //histogram of dR(W muon, leading soft muon per jet) when mu+had mass > 2 GeV
  TH1F* dRWMuSoftMuMuHadMassGe2_;

  //histogram of tau muon pT
  TH1F* tauMuPT_;

  //histogram of hadronic tau pT
  TH1F* tauHadPT_;

  //histogram of hadronic tau deltaBeta corrected isolation
  TH1F* tauHadIso_;

  //histogram of W muon relative isolation computed from leptons only
  TH1F* WMuLeptonRelIso_;

  //histogram of hadronic tau eta
  TH1F* tauHadEta_;

  //histogram of pT(soft muon)/(mu+had mass)
  TH1F* softMuPTOverMuHadMass_;

  //histogram of pT(mu+had)/(mu+had mass)
  TH1F* muHadPTOverMuHadMass_;

  //histogram of dR(soft muon, nearest gen muon)
  TH1F* dRSoftMuNearestGenMuHist_;

  //histogram of mu+had pT
  TH1F* muHadPT_;

  //histogram of number of good vertices (to check PU reweighting)
  TH1F* nGoodVtx_;

  //histogram of W muon + tau muon invariant mass
  TH1F* mWMuTauMu_;

  //histogram of the PDG ID of the nearest status 1 particle to the soft muon
  TH1F* PDGIDNearestStatus1GenParticleToSoftMu_;

  //histogram of the PDG ID of the source particle of the gen muon nearest to the soft muon 
  TH1F* PDGIDMuSrc_;

  //histogram of the invariant mass of the mu+had object and the hardest distinct jet
  TH1F* mSecondJetMuHad_;

  //histogram of Dphi(mu+had, second jet)
  TH1F* dPhiMuHadSecondJet_;

  //histogram of pT rank of mu+had parent uncleaned jet
  TH1F* muHadUncleanedJetPTRank_;

  //histogram of the number of additional loose isolated tight muons with pT > 20 GeV in the event
  TH1F* nAddlHardMuons_;

  //histogram of the number of additional jets with pT > 20 GeV in the event
  TH1F* nAddlJetsPTGeq20_;

  //histogram of the number of additional jets with pT > 30 GeV in the event
  TH1F* nAddlJetsPTGeq30_;

  //histogram of the number of additional jets with pT > 40 GeV in the event
  TH1F* nAddlJetsPTGeq40_;

  //histogram of hadronic tau decay mode
  TH1F* tauHadDecayMode_;

  //histogram of dR(soft muon, nearest gen Z decay muon)
  TH1F* dRSoftMuNearestGenZOrTTMu_;

  //histogram of dR(soft muon, nearest gen Z decay muon that radiated)
  TH1F* dRSoftMuNearestGenZOrTTMuFSR_;

  //histogram of gen muon information
  TH1F* genMuInfo_;

  //histogram of charge product of W muon and third tight muon
  TH1F* WMu3rdTightMuChargeProduct_;

  //histogram of photon energy fraction of hadronic tau
  TH1F* tauHadPhotonEnergyFraction_;

  //histogram of the opening angle between the photon and the other hadronic tau constituents
  TH1F* dThetaPhotonOtherTauConstituents_;

  //histogram of cleaned jet pT vs. cleaned tau pT
  TH2F* cleanedJetPTVsCleanedTauPT_;

  //histogram of uncleaned jet pT vs. cleaned tau pT
  TH2F* uncleanedJetPTVsCleanedTauPT_;

  //histogram of mu+had mass vs. dR(tagged soft muon, tau axis)
  TH2F* muHadMassVsDRSoftMuTau_;

  //histogram of mu+had mass vs. tagged soft muon pT
  TH2F* muHadMassVsSoftMuPT_;

  //histogram of hadronic tau isolation vs. soft muon pT
  TH2F* tauHadIsoVsSoftMuPT_;

  //histogram of properties of nearest muon to soft muon
  TH2F* genMuExistsVsSoftMuNearestMuProperties_;

  //histogram of mu+had multiplicity
  TH1F* muHadMultiplicity_;

  //histogram of cleaned jet tau3/tau1
  TH1F* muHad_t3t1_;

  //histogram of cleaned jet tau2/tau1
  TH1F* muHad_t2t1_;

  //histograms of cleaned jet tau3/tau1, binned by pT
  TH1F* muHad_t3t1_pT1020_;
  TH1F* muHad_t3t1_pT2030_;
  TH1F* muHad_t3t1_pT3040_;
  TH1F* muHad_t3t1_pT4050_;
  TH1F* muHad_t3t1_pT50Up_;

  //histograms of cleaned jet tau3/tau1, binned by NJets
  TH1F* muHad_t3t1_0Jets_;
  TH1F* muHad_t3t1_1Jets_;
  TH1F* muHad_t3t1_2Jets_;
  TH1F* muHad_t3t1_3Jets_;
  TH1F* muHad_t3t1_4Jets_;
  TH1F* muHad_t3t1_5Jets_;
  TH1F* muHad_t3t1_MoreJets_;

  //histogram of no. of charged tracks in cleaned jet (pT > 0 GeV)
  TH1F* muHad_Nchtrk_0_;

  //histogram of no. of charged tracks in cleaned jet (pT > 1 GeV)
  TH1F* muHad_Nchtrk_1_;

  //histogram of no. of charged tracks in cleaned jet (pT > 10 GeV)
  TH1F* muHad_Nchtrk_10_;

  //histogram of no. of charged tracks in cleaned jet (pT > 30 GeV)
  TH1F* muHad_Nchtrk_30_;

  //histogram of no. of charged tracks in second jet (pT > 0 GeV)
  TH1F* second_Nchtrk_0_;

  //histogram of no. of charged tracks in second jet (pT > 1 GeV)
  TH1F* second_Nchtrk_1_;

  //histogram of no. of charged tracks in second jet (pT > 10 GeV)
  TH1F* second_Nchtrk_10_;

  //histogram of no. of charged tracks in second jet (pT > 30 GeV)
  TH1F* second_Nchtrk_30_;

  //histogram of mu+had mass vs. hadronic tau eta
  TH2F* muHadMassVsTauHadEta_;

  //histogram of mu+had mass vs. soft muon eta
  TH2F* muHadMassVsSoftMuEta_;

  //histogram of mu+had mass vs. hadronic tau isolation
  TH2F* muHadMassVsTauHadIso_;

  //histogram of mu+had mass vs. hadronic tau pT
  TH2F* muHadMassVsTauHadPT_;

  //histogram of hadronic tau isolation vs. eta
  TH2F* tauHadIsoVsEta_;

  //histogram of hadronic tau eta vs. soft muon eta
  TH2F* tauHadEtaVsSoftMuEta_;

  //histogram of dEta(hadronic tau, soft muon) vs. dPhi(hadronic tau, soft muon)
  TH2F* dEtaTauHadSoftMuVsDPhiTauHadSoftMu_;

  //histogram of pT(hadronic tau)/(mu+had mass) vs. hadronic tau isolation
  TH2F* tauHadPTOverMuHadMassVsTauHadIso_;

  //histogram of pT(soft muon)/(mu+had mass) vs. hadronic tau isolation
  TH2F* softMuPTOverMuHadMassVsTauHadIso_;

  //histogram of (pT(soft muon) + pT(hadronic tau)/(2(mu+had mass)) vs. hadronic tau isolation
  TH2F* avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_;

  //histogram of pT(mu+had)/(mu+had mass) vs. hadronic tau isolation
  TH2F* muHadPTOverMuHadMassVsTauHadIso_;

  //histogram of W muon isolation vs. hadronic tau isolation
  TH2F* WMuIsoVsTauHadIso_;

  //histogram of soft muon pT vs. hadronic tau pT
  TH2F* softMuPTVsTauHadPT_;

  //histogram of pT/m vs. dimuon invariant mass
  TH2F* muHadPTOverMuHadMassVsMWMuSoftMu_;

  //histogram of soft muon PF pT vs. RECO pT
  TH2F* softMuPFPTVsRECOPT_;

  //histogram of soft muon PF eta vs. RECO eta
  TH2F* softMuPFEtaVsRECOEta_;

  /*histogram of the invariant mass of the mu+had object and the hardest distinct jet vs. mu+had 
    mass*/
  TH2F* mSecondJetMuHadVsMuHadMass_;

  //histogram of mu+had mass vs. 2nd jet mass
  TH2F* muHadMassVsSecondJetMass_;

  //histogram of mu+had pT/m vs. 2nd jet pT/m
  TH2F* muHadPTOverMVsSecondJetPTOverM_;

  //histogram of mu+had mass vs. 2nd jet pT/m
  TH2F* muHadMassVsSecondJetPTOverM_;

  //histogram of HT vs. MET
  TH2F* tauMuTauHadJetWMuHTVsMET_;

  //histogram of dimuon invariant mass vs. mumutau invariant mass
  TH2F* mWMuTauMuVsMWMuTauMuTauHad_;

  /*histogram of dimuon invariant mass vs. mumutau invariant mass for muons which radiated at the 
    gen level*/
  TH2F* mWMuTauMuVsMWMuTauMuTauHadGenFSR_;

  //histogram of mu+had invariant mass vs. mumutau invariant mass
  TH2F* muHadMassVsMWMuTauMuTauHad_;

  //histogram of W_mu MT vs. MET
  TH2F* WMuMTVsMET_;

  //histogram of pruned, cleaned tau3/tau1 vs unpruned, uncleaned pt/m for muHad jet
  TH2F* muHad_t3t1Vsptmj_;

  //histogram of pruned, cleaned tau3/tau1 vs HPS tau decay mode
  TH2F* muHad_t3t1VsDecayMode_;

  //histogram of mu+had mass vs. Nj
  TH2F* muHadMassVsNAddlJets_;

  //PU reweighting object
  edm::LumiReWeighting PUReweight_;

  //hadronic tau pT weights
  std::vector<double> tauHadPTWeights_;

  //hadronic tau pT weight errors
  std::vector<double> tauHadPTWeightErrs_;

  //hadronic tau pT bins
  std::vector<double> tauHadPTBins_;

  //mu+had mass bins
  std::vector<double> muHadMassBins_;

  //hadronic tau pT bins for fake rate corrections
  std::vector<double> fakeRateTauHadPTBins_;

  //MC fake rate corrections
  std::vector<double> fakeRateCorrs_;

  //MC fake rate correction errors
  std::vector<double> fakeRateCorrErrs_;

  /*second jet pT for |eta| < 2.4, second jet is highest pT jet in the event excluding W muon and 
    mu+had*/
  TH1F *jet_pt_etacut;

  //second jet eta, second jet is highest pT jet in the event excluding W muon and mu+had
  TH1F *jet_eta;

  //second jet phi, second jet is highest pT jet in the event excluding W muon and mu+had
  TH1F *jet_phi;

  /*second jet mass for |eta| < 2.4, second jet is highest pT jet in the event excluding W muon 
    and mu+had*/
  TH1F *jet_mass_etacut;

  /*second jet pT/m for |eta| < 2.4, second jet is highest pT jet in the event excluding W muon 
    and mu+had*/
  TH1F *jet_ptmj_etacut;

  //dPhi between W muon and second jet (highest pT jet in the event excluding W muon)
  TH1F *dPhiWMuSecJet_;

  //maximum hadronic tau pT
  double maxTauHadPT_;

  //maximum mu+had mass
  double maxMuHadMass_;
};


// constants, enums and typedefs



// static data member definitions



// constructors and destructor

TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  METTag_(iConfig.getParameter<edm::InputTag>("METTag")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  genMatchedMuonTag_(iConfig.getParameter<edm::InputTag>("genMatchedMuonTag")),
  oldJetTag_(iConfig.getParameter<edm::InputTag>("oldJetTag")),
  newJetTag_(iConfig.getParameter<edm::InputTag>("newJetTag")),
  jetMuonMapTag_(iConfig.getParameter<edm::InputTag>("jetMuonMapTag")),
  oldNewJetMapTag_(iConfig.getParameter<edm::InputTag>("oldNewJetMapTag")),
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  genTauMuTag_(iConfig.getParameter<edm::InputTag>("genTauMuTag")),
  genWTauMuTag_(iConfig.getParameter<edm::InputTag>("genWTauMuTag")),
  tauHadIsoTag_(iConfig.getParameter<edm::InputTag>("tauHadIsoTag")),
  allMuonTag_(iConfig.getParameter<edm::InputTag>("allMuonTag")),
  muonGenParticleTag_(iConfig.getParameter<edm::InputTag>("muonGenParticleTag")),
  muonPFIsoPUSubtractionCoeff_(iConfig.getParameter<double>("muonPFIsoPUSubtractionCoeff")),
  PUTag_(iConfig.getParameter<edm::InputTag>("PUTag")),
  vtxTag_(iConfig.getParameter<edm::InputTag>("vtxTag")),
  allGenParticleTag_(iConfig.getParameter<edm::InputTag>("allGenParticleTag")),
  corrJetTag_(iConfig.getParameter<edm::InputTag>("corrJetTag")),
  bJetTag_(iConfig.getParameter<edm::InputTag>("bJetTag")),
  dR_(iConfig.getParameter<double>("dR")),
  tauPTMin_(iConfig.getParameter<double>("tauPTMin")),
  tauDecayMode_(static_cast<reco::PFTau::hadronicDecayMode>
		(iConfig.getParameter<int>("tauDecayMode"))),
  uncorrJetPTMin_(iConfig.getParameter<double>("uncorrJetPTMin")),
  tauArbitrationMethod_(iConfig.getParameter<std::string>("tauArbitrationMethod")),
  PUScenario_(iConfig.getParameter<std::string>("PUScenario")),
  zCut_(iConfig.getParameter<double>("zCut")),
  RcutFactor_(iConfig.getParameter<double>("RcutFactor")),
  MC_(iConfig.getParameter<bool>("MC")),
  reweight_(iConfig.getParameter<bool>("reweight")),
  pTRankColors_(iConfig.getParameter<std::vector<unsigned int> >("pTRankColors")),
  pTRankStyles_(iConfig.getParameter<std::vector<unsigned int> >("pTRankStyles")),
  pTRankEntries_(iConfig.getParameter<std::vector<std::string> >("pTRankEntries"))
{
  //now do what ever initialization is needed
  reset(false);

  //check that tau arbitration method is valid
  if ((tauArbitrationMethod_ != "pT") && (tauArbitrationMethod_ != "m") && 
      (tauArbitrationMethod_ != "none")) {
    std::stringstream err;
    err << "Error: unsupported tau arbitration method.  Valid choices are:\n";
    err << "---\"pT\": highest pT selected tau chosen per event\n";
    err << "---\"m\": highest mu+had mass selected tau chosen per event\n";
    err << "--\"none\": all selected taus are used\n";
    throw cms::Exception("TauAnalyzer") << err.str();
  }

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
  float Data20122p5InvFbPUDistArray[60] = {6.633e-03, 7.767e-02, 9.649e+02, 1.449e+05, 
					   3.185e+05, 5.658e+05, 4.224e+06, 1.779e+07, 
					   4.045e+07, 6.108e+07, 8.996e+07, 1.319e+08, 
					   1.670e+08, 1.925e+08, 2.075e+08, 2.081e+08, 
					   1.960e+08, 1.755e+08, 1.537e+08, 1.349e+08, 
					   1.186e+08, 1.054e+08, 9.523e+07, 8.614e+07, 
					   7.637e+07, 6.542e+07, 5.384e+07, 4.252e+07, 
					   3.212e+07, 2.310e+07, 1.580e+07, 1.035e+07, 
					   6.544e+06, 4.001e+06, 2.354e+06, 1.320e+06, 
					   6.986e+05, 3.473e+05, 1.620e+05, 7.126e+04, 
					   2.994e+04, 1.232e+04, 5.175e+03, 2.318e+03, 
					   1.136e+03, 6.016e+02, 3.322e+02, 1.849e+02, 
					   1.015e+02, 5.438e+01, 2.826e+01, 1.422e+01, 
					   6.923e+00, 3.258e+00, 1.482e+00, 6.517e-01, 
					   2.770e-01, 1.138e-01, 4.518e-02, 1.734e-02};
  std::vector<float> Data2012PUDist(Data2012PUDistArray, Data2012PUDistArray + 
				    sizeof(Data2012PUDistArray)/sizeof(float));
  std::vector<float> 
    Data20122p5InvFbPUDist(Data20122p5InvFbPUDistArray, Data20122p5InvFbPUDistArray + 
			   sizeof(Data20122p5InvFbPUDistArray)/sizeof(float));
  if (PUScenario_ == "S7") {
//     float S7PUDistArray[60] = {2.344E-05, 2.344E-05, 2.344E-05, 2.344E-05, 4.687E-04, 4.687E-04, 
//     			       7.032E-04, 9.414E-04, 1.234E-03, 1.603E-03, 2.464E-03, 3.250E-03, 
//     			       5.021E-03, 6.644E-03, 8.502E-03, 1.121E-02, 1.518E-02, 2.033E-02, 
//     			       2.608E-02, 3.171E-02, 3.667E-02, 4.060E-02, 4.338E-02, 4.520E-02, 
//     			       4.641E-02, 4.735E-02, 4.816E-02, 4.881E-02, 4.917E-02, 4.909E-02, 
//     			       4.842E-02, 4.707E-02, 4.501E-02, 4.228E-02, 3.896E-02, 3.521E-02, 
//     			       3.118E-02, 2.702E-02, 2.287E-02, 1.885E-02, 1.508E-02, 1.166E-02, 
//     			       8.673E-03, 6.190E-03, 4.222E-03, 2.746E-03, 1.698E-03, 9.971E-04, 
//     			       5.549E-04, 2.924E-04, 1.457E-04, 6.864E-05, 3.054E-05, 1.282E-05, 
//     			       5.081E-06, 1.898E-06, 6.688E-07, 2.221E-07, 6.947E-08, 2.047E-08};
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
    PUReweight_ = edm::LumiReWeighting(S7PUDist, Data2012PUDist);
//     PUReweight_ = edm::LumiReWeighting(S7PUDist, Data20122p5InvFbPUDist);
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
    PUReweight_ = edm::LumiReWeighting(S10PUDist, Data2012PUDist);
//     PUReweight_ = edm::LumiReWeighting(S10PUDist, Data20122p5InvFbPUDist);
  }

  //instantiate the vector of weights based on hadronic tau pT bin
//   tauHadPTWeights_.push_back(0.0); //derived from 2.5 fb^-1 data with trigger matching enabled
//   tauHadPTWeights_.push_back(0.0);
//   tauHadPTWeights_.push_back(1.52423);
//   tauHadPTWeights_.push_back(1.29072);
//   tauHadPTWeights_.push_back(0.927203);
//   tauHadPTWeights_.push_back(0.782077);
//   tauHadPTWeights_.push_back(0.626892);
//   tauHadPTWeights_.push_back(0.458541);
//   tauHadPTWeights_.push_back(0.397239);
//   tauHadPTWeights_.push_back(0.323495);
//   tauHadPTWeights_.push_back(0.34062);
//   tauHadPTWeights_.push_back(0.290588);
//   tauHadPTWeights_.push_back(0.137788);
//   tauHadPTWeights_.push_back(0.313786);
  tauHadPTWeights_.push_back(0.0); /*derived from full 2012 dataset with trigger matching enabled 
				     and second jet requirement (pT > 20 GeV)*/
  tauHadPTWeights_.push_back(0.0);
  tauHadPTWeights_.push_back(1.62967);
  tauHadPTWeights_.push_back(1.36438);
  tauHadPTWeights_.push_back(1.09932);
  tauHadPTWeights_.push_back(0.732214);
  tauHadPTWeights_.push_back(0.505928);
  tauHadPTWeights_.push_back(0.234365);
  tauHadPTWeights_.push_back(0.328963);
  tauHadPTWeights_.push_back(0.254085);
  tauHadPTWeights_.push_back(0.452579);
  tauHadPTWeights_.push_back(0.0864345);

  //instantiate the vector of hadronic tau pT weight errors
//   tauHadPTWeightErrs_.push_back(0.0); //derived from 2.5 fb^-1 data with trigger matching enabled
//   tauHadPTWeightErrs_.push_back(0.0);
//   tauHadPTWeightErrs_.push_back(0.11641);
//   tauHadPTWeightErrs_.push_back(0.114134);
//   tauHadPTWeightErrs_.push_back(0.105168);
//   tauHadPTWeightErrs_.push_back(0.112472);
//   tauHadPTWeightErrs_.push_back(0.11584);
//   tauHadPTWeightErrs_.push_back(0.12018);
//   tauHadPTWeightErrs_.push_back(0.100846);
//   tauHadPTWeightErrs_.push_back(0.123222);
//   tauHadPTWeightErrs_.push_back(0.171247);
//   tauHadPTWeightErrs_.push_back(0.206207);
//   tauHadPTWeightErrs_.push_back(0.138026);
//   tauHadPTWeightErrs_.push_back(0.222703);
  tauHadPTWeightErrs_.push_back(0.0); /*derived from full 2012 dataset with trigger matching 
					enabled and second jet requirement (pT > 20 GeV)*/
  tauHadPTWeightErrs_.push_back(0.0);
  tauHadPTWeightErrs_.push_back(0.177485);
  tauHadPTWeightErrs_.push_back(0.170437);
  tauHadPTWeightErrs_.push_back(0.165546);
  tauHadPTWeightErrs_.push_back(0.152739);
  tauHadPTWeightErrs_.push_back(0.143558);
  tauHadPTWeightErrs_.push_back(0.118059);
  tauHadPTWeightErrs_.push_back(0.125919);
  tauHadPTWeightErrs_.push_back(0.147566);
  tauHadPTWeightErrs_.push_back(0.263058);
  tauHadPTWeightErrs_.push_back(0.0866064);
 
//   //fill hadronic tau pT bins
//   tauHadPTBins_.push_back(0.0);
//   tauHadPTBins_.push_back(5.0);
//   tauHadPTBins_.push_back(10.0);
//   tauHadPTBins_.push_back(15.0);
//   tauHadPTBins_.push_back(20.0);
//   tauHadPTBins_.push_back(25.0);
//   tauHadPTBins_.push_back(30.0);
//   tauHadPTBins_.push_back(35.0);
//   tauHadPTBins_.push_back(40.0);
//   tauHadPTBins_.push_back(50.0);
//   tauHadPTBins_.push_back(60.0);
//   tauHadPTBins_.push_back(70.0);
// //   tauHadPTBins_.push_back(80.0);
// //   tauHadPTBins_.push_back(100.0);
//   tauHadPTBins_.push_back(600.0);

  //fill hadronic tau pT bins
  tauHadPTBins_.push_back(0.0);
  tauHadPTBins_.push_back(5.0);
  tauHadPTBins_.push_back(10.0);
  tauHadPTBins_.push_back(15.0);
  tauHadPTBins_.push_back(20.0);
  tauHadPTBins_.push_back(25.0);
  tauHadPTBins_.push_back(30.0);
  tauHadPTBins_.push_back(35.0);
  tauHadPTBins_.push_back(40.0);
  tauHadPTBins_.push_back(45.0);
  tauHadPTBins_.push_back(50.0);
  tauHadPTBins_.push_back(55.0);
  tauHadPTBins_.push_back(60.0);
  tauHadPTBins_.push_back(65.0);
  tauHadPTBins_.push_back(70.0);
  tauHadPTBins_.push_back(75.0);
  tauHadPTBins_.push_back(80.0);
  tauHadPTBins_.push_back(85.0);
  tauHadPTBins_.push_back(90.0);
  tauHadPTBins_.push_back(95.0);
  tauHadPTBins_.push_back(100.0);
  tauHadPTBins_.push_back(105.0);
  tauHadPTBins_.push_back(110.0);
  tauHadPTBins_.push_back(115.0);
  tauHadPTBins_.push_back(120.0);
  tauHadPTBins_.push_back(125.0);
  tauHadPTBins_.push_back(130.0);
  tauHadPTBins_.push_back(135.0);
  tauHadPTBins_.push_back(140.0);
  tauHadPTBins_.push_back(145.0);
  tauHadPTBins_.push_back(150.0);
  tauHadPTBins_.push_back(155.0);
  tauHadPTBins_.push_back(160.0);
  tauHadPTBins_.push_back(165.0);
  tauHadPTBins_.push_back(170.0);
  tauHadPTBins_.push_back(175.0);
  tauHadPTBins_.push_back(180.0);
  tauHadPTBins_.push_back(185.0);
  tauHadPTBins_.push_back(190.0);
  tauHadPTBins_.push_back(195.0);
  tauHadPTBins_.push_back(200.0);
  tauHadPTBins_.push_back(205.0);
  tauHadPTBins_.push_back(210.0);
  tauHadPTBins_.push_back(215.0);
  tauHadPTBins_.push_back(220.0);
  tauHadPTBins_.push_back(225.0);
  tauHadPTBins_.push_back(230.0);
  tauHadPTBins_.push_back(235.0);
  tauHadPTBins_.push_back(240.0);
  tauHadPTBins_.push_back(245.0);
  tauHadPTBins_.push_back(250.0);
  tauHadPTBins_.push_back(255.0);
  tauHadPTBins_.push_back(260.0);
  tauHadPTBins_.push_back(265.0);
  tauHadPTBins_.push_back(270.0);
  tauHadPTBins_.push_back(275.0);
  tauHadPTBins_.push_back(280.0);
  tauHadPTBins_.push_back(285.0);
  tauHadPTBins_.push_back(290.0);
  tauHadPTBins_.push_back(295.0);
  tauHadPTBins_.push_back(300.0);
  tauHadPTBins_.push_back(305.0);
  tauHadPTBins_.push_back(310.0);
  tauHadPTBins_.push_back(315.0);
  tauHadPTBins_.push_back(320.0);
  tauHadPTBins_.push_back(325.0);
  tauHadPTBins_.push_back(330.0);
  tauHadPTBins_.push_back(335.0);
  tauHadPTBins_.push_back(340.0);
  tauHadPTBins_.push_back(345.0);
  tauHadPTBins_.push_back(350.0);
  tauHadPTBins_.push_back(355.0);
  tauHadPTBins_.push_back(360.0);
  tauHadPTBins_.push_back(365.0);
  tauHadPTBins_.push_back(370.0);
  tauHadPTBins_.push_back(375.0);
  tauHadPTBins_.push_back(380.0);
  tauHadPTBins_.push_back(385.0);
  tauHadPTBins_.push_back(390.0);
  tauHadPTBins_.push_back(395.0);
  tauHadPTBins_.push_back(400.0);
  tauHadPTBins_.push_back(405.0);
  tauHadPTBins_.push_back(410.0);
  tauHadPTBins_.push_back(415.0);
  tauHadPTBins_.push_back(420.0);
  tauHadPTBins_.push_back(425.0);
  tauHadPTBins_.push_back(430.0);
  tauHadPTBins_.push_back(435.0);
  tauHadPTBins_.push_back(440.0);
  tauHadPTBins_.push_back(445.0);
  tauHadPTBins_.push_back(450.0);
  tauHadPTBins_.push_back(455.0);
  tauHadPTBins_.push_back(460.0);
  tauHadPTBins_.push_back(465.0);
  tauHadPTBins_.push_back(470.0);
  tauHadPTBins_.push_back(475.0);
  tauHadPTBins_.push_back(480.0);
  tauHadPTBins_.push_back(485.0);
  tauHadPTBins_.push_back(490.0);
  tauHadPTBins_.push_back(495.0);
  tauHadPTBins_.push_back(500.0);
  tauHadPTBins_.push_back(505.0);
  tauHadPTBins_.push_back(510.0);
  tauHadPTBins_.push_back(515.0);
  tauHadPTBins_.push_back(520.0);
  tauHadPTBins_.push_back(525.0);
  tauHadPTBins_.push_back(530.0);
  tauHadPTBins_.push_back(535.0);
  tauHadPTBins_.push_back(540.0);
  tauHadPTBins_.push_back(545.0);
  tauHadPTBins_.push_back(550.0);
  tauHadPTBins_.push_back(555.0);
  tauHadPTBins_.push_back(560.0);
  tauHadPTBins_.push_back(565.0);
  tauHadPTBins_.push_back(570.0);
  tauHadPTBins_.push_back(575.0);
  tauHadPTBins_.push_back(580.0);
  tauHadPTBins_.push_back(585.0);
  tauHadPTBins_.push_back(590.0);
  tauHadPTBins_.push_back(595.0);
  tauHadPTBins_.push_back(600.0);

  //fill mu+had mass bins
  muHadMassBins_.push_back(0.0);
  muHadMassBins_.push_back(1.0);
  muHadMassBins_.push_back(2.0);
  muHadMassBins_.push_back(3.0);
  muHadMassBins_.push_back(4.0);
  muHadMassBins_.push_back(5.0);
  muHadMassBins_.push_back(6.0);
  muHadMassBins_.push_back(8.0);
//   muHadMassBins_.push_back(10.0);
  muHadMassBins_.push_back(60.0); //maximum over all MC and full 2012 dataset

  //fill hadronic tau pT bins for isolated MC fake rate correction
  fakeRateTauHadPTBins_.push_back(0.0);
  fakeRateTauHadPTBins_.push_back(5.0);
  fakeRateTauHadPTBins_.push_back(10.0);
  fakeRateTauHadPTBins_.push_back(15.0);
  fakeRateTauHadPTBins_.push_back(20.0);
  fakeRateTauHadPTBins_.push_back(25.0);
  fakeRateTauHadPTBins_.push_back(30.0);
  fakeRateTauHadPTBins_.push_back(35.0);
  fakeRateTauHadPTBins_.push_back(40.0);
  fakeRateTauHadPTBins_.push_back(45.0);
  fakeRateTauHadPTBins_.push_back(50.0);
  fakeRateTauHadPTBins_.push_back(600.0);

  //fill fake rate corrections
  fakeRateCorrs_.push_back(0.0);
  fakeRateCorrs_.push_back(0.0);
  fakeRateCorrs_.push_back(1.22837);
  fakeRateCorrs_.push_back(1.1525);
  fakeRateCorrs_.push_back(1.15403);
  fakeRateCorrs_.push_back(0.861062);
  fakeRateCorrs_.push_back(0.572035);
  fakeRateCorrs_.push_back(0.597443);
  fakeRateCorrs_.push_back(0.643079);
  fakeRateCorrs_.push_back(0.907813);
  fakeRateCorrs_.push_back(1.87002);

  //fill fake rate correction errors
  fakeRateCorrErrs_.push_back(0.0);
  fakeRateCorrErrs_.push_back(0.0);
  fakeRateCorrErrs_.push_back(0.201149);
  fakeRateCorrErrs_.push_back(0.191108);
  fakeRateCorrErrs_.push_back(0.234955);
  fakeRateCorrErrs_.push_back(0.206023);
  fakeRateCorrErrs_.push_back(0.174063);
  fakeRateCorrErrs_.push_back(0.222943);
  fakeRateCorrErrs_.push_back(0.295332);
  fakeRateCorrErrs_.push_back(0.505194);
  fakeRateCorrErrs_.push_back(0.734409);
}

TauAnalyzer::~TauAnalyzer()
{
 
//   do anything here that needs to be done at desctruction time
//   (e.g. close files, deallocate resources etc.)
  reset(true);
}



// member functions


// ------------ method called for each event  ------------
void TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ //analyzer

//   //print event number
//   std::cerr << "Run " << iEvent.run() << ", event " << iEvent.id().event() << ", lumi ";
//   std::cerr << iEvent.luminosityBlock() << std::endl;

  //get taus
  edm::Handle<reco::PFTauRefVector> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);
  
  //get MET
  edm::Handle<edm::View<reco::PFMET> > pMET;
  iEvent.getByLabel(METTag_, pMET);

  //get W muons
  edm::Handle<reco::MuonRefVector> pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  //get gen-matched muons
  edm::Handle<reco::MuonRefVector> pGenMatchedMuons;
  if (MC_) iEvent.getByLabel(genMatchedMuonTag_, pGenMatchedMuons);

  //get old jets
  edm::Handle<reco::PFJetCollection> pOldJets;
  iEvent.getByLabel(oldJetTag_, pOldJets);

  //get new jets
  edm::Handle<reco::PFJetCollection> pNewJets;
  iEvent.getByLabel(newJetTag_, pNewJets);

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  iEvent.getByLabel(jetMuonMapTag_, pMuonJetMap);

  //get old-new jet map
  edm::Handle<edm::ValueMap<reco::PFJetRef> > pOldNewJetMap;
  iEvent.getByLabel(oldNewJetMapTag_, pOldNewJetMap);

  //get gen particles
  edm::Handle<reco::GenParticleRefVector> pGenParticles;
  if (MC_) iEvent.getByLabel(genParticleTag_, pGenParticles);

  //get gen tau-->mu muons
  edm::Handle<reco::GenParticleRefVector> pGenTauMus;
  if (MC_) iEvent.getByLabel(genTauMuTag_, pGenTauMus);

  //get gen W-->tau-->mu muons
  edm::Handle<reco::GenParticleRefVector> pGenWTauMus;
  if (MC_) iEvent.getByLabel(genWTauMuTag_, pGenWTauMus);

  //get hadronic tau deltaBeta-corrected isolation
  edm::Handle<reco::PFTauDiscriminator> pTauHadIso;
  iEvent.getByLabel(tauHadIsoTag_, pTauHadIso);

  //get muons
  edm::Handle<reco::MuonCollection> pAllMuons;
  iEvent.getByLabel(allMuonTag_, pAllMuons);

  //get muon gen particles
  edm::Handle<reco::CandidateView> pMuonGenParticles;
  if (MC_) iEvent.getByLabel(muonGenParticleTag_, pMuonGenParticles);

  //get AK5 PF L1FastL2L3 jet correction service
  const JetCorrector* corrector = JetCorrector::getJetCorrector("ak5PFL1FastL2L3", iSetup);

  //get PU info
  edm::Handle<std::vector<PileupSummaryInfo> > pPU;
  if (MC_) iEvent.getByLabel(PUTag_, pPU);

  //get vertices
  edm::Handle<reco::VertexCollection> pVtx;
  iEvent.getByLabel(vtxTag_, pVtx);

  //get all gen particles
  edm::Handle<reco::GenParticleCollection> pAllGenParticles;
  if (MC_) iEvent.getByLabel(allGenParticleTag_, pAllGenParticles);

  //get selected corrected jets
  edm::Handle<reco::PFJetCollection> pCorrJets;
  iEvent.getByLabel(corrJetTag_, pCorrJets);

//   //debug
//   for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
//        ++iTau) {
//     const reco::PFJetRef& tauJetRef = (*iTau)->jetRef();
//     std::cerr << "Tau jet ref key: " << tauJetRef.key() << std::endl;
//     const reco::MuonRefVector& removedMuons = (*pMuonJetMap)[tauJetRef];
//     for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
//   	 iMuon != removedMuons.end(); ++iMuon) {
//       std::cerr << "Muon ref key: " << iMuon->key() << std::endl;
//       std::cerr << "Muon pT: " << (*iMuon)->pt() << " GeV\n";
//       std::cerr << "Muon eta: " << (*iMuon)->eta() << std::endl;
//       std::cerr << "Muon phi: " << (*iMuon)->phi() << std::endl;
//     }
//   }

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

  //Get b tag information
  edm::Handle<reco::JetTagCollection> bTagHandle;
  iEvent.getByLabel(bJetTag_, bTagHandle);
  const reco::JetTagCollection & bTags = *(bTagHandle.product());

  //find the highest pT W muon
  std::vector<reco::MuonRef> WMuonRefs;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) { WMuonRefs.push_back(*iMuon); }
  Common::sortByPT(WMuonRefs);

  //fill STL containers of pointers to the gen particles
  std::vector<reco::GenParticle*> genParticlePtrs;
  std::vector<reco::GenParticle*> genTauMuPtrs;
  std::vector<reco::GenParticle*> genWTauMuPtrs;
  std::vector<reco::GenParticle*> status1GenParticlePtrs;
  std::vector<reco::GenParticle*> ZOrTTMuGenParticlePtrs;
  std::vector<reco::GenParticle*> ZOrTTMuFSRGenParticlePtrs;
  std::vector<int> DrellYanOrTTMothers;
  DrellYanOrTTMothers.push_back(GenTauDecayID::ZPDGID);
  DrellYanOrTTMothers.push_back(GenTauDecayID::GAMMAPDGID);
  DrellYanOrTTMothers.push_back(GenTauDecayID::TPDGID);
  if (MC_) {
    for (reco::GenParticleRefVector::const_iterator iGenParticle = pGenParticles->begin(); 
	 iGenParticle != pGenParticles->end(); ++iGenParticle) {
      genParticlePtrs.push_back(const_cast<reco::GenParticle*>((*iGenParticle).get()));
    }
    for (reco::GenParticleRefVector::const_iterator iGenParticle = pGenWTauMus->begin(); 
	 iGenParticle != pGenWTauMus->end(); ++iGenParticle) {
      genParticlePtrs.push_back(const_cast<reco::GenParticle*>((*iGenParticle).get()));
    }
    for (reco::GenParticleRefVector::const_iterator iGenParticle = pGenTauMus->begin(); 
	 iGenParticle != pGenTauMus->end(); ++iGenParticle) {
      genTauMuPtrs.push_back(const_cast<reco::GenParticle*>((*iGenParticle).get()));
    }
    for (reco::GenParticleCollection::const_iterator iGenParticle = pAllGenParticles->begin(); 
	 iGenParticle != pAllGenParticles->end(); ++iGenParticle) {
      if (iGenParticle->status() == 1) {
	status1GenParticlePtrs.push_back(const_cast<reco::GenParticle*>(&*iGenParticle));
// 	std::cerr << iGenParticle->pdgId() << std::endl;
      }
      if ((iGenParticle->status() == 1) && 
	  (fabs(iGenParticle->pdgId()) == GenTauDecayID::MUPDGID) && 
	  Common::hasAncestor(reco::GenParticleRef(pAllGenParticles, 
						   iGenParticle - pAllGenParticles->begin()), 
			      DrellYanOrTTMothers)) {
	ZOrTTMuGenParticlePtrs.push_back(const_cast<reco::GenParticle*>(&*iGenParticle));
	reco::GenParticleRef mom = iGenParticle->motherRef();
	if ((fabs(mom->pdgId()) == GenTauDecayID::MUPDGID) && (mom->status() == 2)) {
	  bool isFSR = false;
	  unsigned int iDaughter = 0;
	  while ((iDaughter < mom->numberOfDaughters()) && !isFSR) {
	    if (fabs(mom->daughterRef(iDaughter)->pdgId()) == GenTauDecayID::GAMMAPDGID) {
	      isFSR = true;
	    }
	    ++iDaughter;
	  }
	  if (isFSR) {
	    ZOrTTMuFSRGenParticlePtrs.push_back(const_cast<reco::GenParticle*>(&*iGenParticle));
	  }
	}
      }
    }
  }

  //fill collection of all corrected jets (preserve order of original jet collection)
  std::vector<reco::PFJet> correctedJets;
  for (reco::PFJetCollection::const_iterator iJet = pOldJets->begin(); iJet != pOldJets->end(); 
       ++iJet) {
    reco::PFJet correctedJet = *iJet;
    double JEC = corrector->correction(*iJet, iEvent, iSetup);
    correctedJet.scaleEnergy(JEC);
    correctedJets.push_back(correctedJet);
  }  

  //fill collection of all corrected jets excluding the jet associated to the W muon
  std::vector<reco::PFJet> correctedOldJets;
  std::vector<reco::PFJetRef> correctedOldJetRefs; //ref keys refer to original jet collection
  for (reco::PFJetCollection::const_iterator iJet = correctedJets.begin(); 
       iJet != correctedJets.end(); ++iJet) {
    if (reco::deltaR(*WMuonRefs[WMuonRefs.size() - 1], *iJet) >= dR_) { /*technically you should 
									  check that the W muon 
									  ref key is not 
									  identical to the ref 
									  key of any jet 
									  constituent, but that 
									  would take a much 
									  longer time*/
      correctedOldJets.push_back(*iJet);
      correctedOldJetRefs.push_back(reco::PFJetRef(&correctedJets, iJet - correctedJets.begin()));
    }
  }

    /*number of jets in event with |eta| < 2.4 and corrected pT > 20
      that do not overlap with the W muon or HPS tau*/
    int NJets = pCorrJets->size();

    /*make and sort by pT refs to AK5 PF jets in the selected corrected collection (distinct from 
      this tau and the W muon and passing |eta| < 2.4)*/
    std::vector<reco::PFJetRef> oldJetRefsExcludingTau;
    for (reco::PFJetCollection::const_iterator iCorrectedJetExcludingTau = pCorrJets->begin(); 
	 iCorrectedJetExcludingTau != pCorrJets->end(); ++iCorrectedJetExcludingTau) {
      oldJetRefsExcludingTau.
	push_back(reco::PFJetRef(pCorrJets, iCorrectedJetExcludingTau - pCorrJets->begin()));
    }
    Common::sortByPT(oldJetRefsExcludingTau);

    //make a collection of corrected old jets in |eta| < 2.4 not overlapping the W muon
    //make a collection of corrected old jets in |eta| < 2.4 not overlapping the W muon or tau
    std::vector<reco::PFJetRef> oldJetRefsInclTauNoPTCut;
    for (reco::PFJetCollection::const_iterator iCorrectedJet = correctedOldJets.begin(); 
	 iCorrectedJet != correctedOldJets.end(); ++iCorrectedJet) {
      if (fabs(iCorrectedJet->eta()) < 2.4) {
	oldJetRefsInclTauNoPTCut.
	  push_back(reco::PFJetRef(&correctedOldJets, iCorrectedJet - correctedOldJets.begin()));
      }
    }
    Common::sortByPT(oldJetRefsInclTauNoPTCut);


//   //debug
//   for (reco::MuonCollection::const_iterator iMu = pAllMuons->begin(); iMu != pAllMuons->end(); 
//        ++iMu) {
//     std::cerr << "Reco muon\n----------\n";
//     std::cerr << "pT = " << iMu->pt() << " GeV\n";
//     std::cerr << "eta = " << iMu->eta() << std::endl;
//     std::cerr << "phi = " << iMu->phi() << std::endl;
//   }

  //fill an STL container of muon gen particles
  std::vector<reco::Candidate*> muGenParticlePtrs;
  if (MC_) {
    for (unsigned int iMuGenParticle = 0; iMuGenParticle < pMuonGenParticles->size(); 
	 ++iMuGenParticle) {
      muGenParticlePtrs.
	push_back(const_cast<reco::Candidate*>(pMuonGenParticles->refAt(iMuGenParticle).get()));
    }
  }

  //sort selected taus by descending order in pT
  std::vector<reco::PFTauRef> pTSortedTaus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) { pTSortedTaus.push_back(*iTau); }
  std::vector<reco::PFTauRef> taus = pTSortedTaus;
  Common::sortByPT(pTSortedTaus);
  std::reverse(pTSortedTaus.begin(), pTSortedTaus.end());

  //sort selected taus by descending order in mu+had mass
  std::vector<reco::PFTauRef> muHadMassSortedTaus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) { muHadMassSortedTaus.push_back(*iTau); }
  Common::sortByMass(pMuonJetMap, muHadMassSortedTaus);

  //mu+had object counter
  unsigned int nMuHad = 0;

  //loop over selected taus
  std::vector<reco::PFTauRef>::const_iterator iTau = taus.begin();
  std::vector<reco::PFTauRef>::const_iterator endTau = taus.end();
  if (tauArbitrationMethod_ == "pT") {
    iTau = pTSortedTaus.begin();
    endTau = iTau + 1;
  }
  if (tauArbitrationMethod_ == "m") {
    iTau = muHadMassSortedTaus.begin();
    endTau = iTau + 1;
  }
  while (iTau != endTau) {
    const reco::PFJetRef& tauJetRef = (*iTau)->jetRef();
    const reco::PFJetRef& tauOldJetRef = (*pOldNewJetMap)[tauJetRef];
    const reco::PFJet correctedTauJet = correctedJets[tauOldJetRef.key()];
    const reco::MuonRefVector& removedMuons = (*pMuonJetMap)[tauJetRef];
    
    double b_discriminant = 999.;
    for (unsigned int i = 0; i != bTags.size(); ++i)
      { // loop over btags
	PFJetRef myBJetRef(pOldJets,i);
	if (myBJetRef.key() == tauJetRef.key())
	  {
	    bTagDiscrim_->Fill(bTags[i].second);
	    b_discriminant = bTags[i].second;
	    break;
	  }
      } // loop over btags

    //make a collection of corrected old jets in |eta| < 2.4 not overlapping the W muon or tau
    std::vector<reco::PFJetRef> oldJetRefsExclTauNoPTCut;
    for (reco::PFJetCollection::const_iterator iCorrectedJet = correctedOldJets.begin(); 
	 iCorrectedJet != correctedOldJets.end(); ++iCorrectedJet) {
      if ((fabs(iCorrectedJet->eta()) < 2.4) && (reco::deltaR(*iCorrectedJet, **iTau) >= dR_)) {
	oldJetRefsExclTauNoPTCut.
	  push_back(reco::PFJetRef(&correctedOldJets, iCorrectedJet - correctedOldJets.begin()));
      }
    }
    Common::sortByPT(oldJetRefsExclTauNoPTCut);

    //get weight based on hadronic tau pT
    bool foundBin = false;
    std::vector<double>::const_iterator iBinEdge = tauHadPTBins_.begin();
    while ((iBinEdge != tauHadPTBins_.end()) && !foundBin) {
      if (((*iTau)->pt() >= *iBinEdge) && ((*iTau)->pt() < *(iBinEdge + 1))) {
	foundBin = true;
      }
      else ++iBinEdge;
    }
//     const int tauHadPTBin = (int)(iBinEdge - fakeRateTauHadPTBins_.begin());
    double tauHadPTWeight = 1.0;
    double tauHadPTWeightErr = 0.0;
    if (reweight_) {
      tauHadPTWeight = exp(1.2 - 0.059*(*iTau)->pt());
      tauHadPTWeightErr = tauHadPTWeight*sqrt(0.1*0.1 + (*iTau)->pt()*(*iTau)->pt()*0.005*0.005);
//       if (foundBin) {
// 	tauHadPTWeight = tauHadPTWeights_[tauHadPTBin];
// 	tauHadPTWeightErr = tauHadPTWeightErrs_[tauHadPTBin];
//       }
    }

    //find the highest pT associated muon
    std::vector<reco::MuonRef> removedMuonRefs;
    for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
	 iMuon != removedMuons.end(); ++iMuon) { removedMuonRefs.push_back(*iMuon); }
    Common::sortByPT(removedMuonRefs);

    /*calculate invariant mass of W muon, tau muon, and hadronic tau (sensitive to Z-->mumugamma 
      FSR)*/
    const double mWMuTauMuTauHad = (WMuonRefs[WMuonRefs.size() - 1]->p4() + 
				    removedMuonRefs[removedMuonRefs.size() - 1]->p4() + 
				    (*iTau)->p4()).M();

    //calculate invariant mass of W muon and tau muon (sensitive to Z-->mumu)
    const double mWMuTauMu = (WMuonRefs[WMuonRefs.size() - 1]->p4() + 
			      removedMuonRefs[removedMuonRefs.size() - 1]->p4()).M();

    //calculate mu + tau invariant mass for the highest pT muon
    const double muHadMass = 
      (removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).M();

    /*get charge product of W muon and all third tight muons
      - 3rd muon not overlapping the tau muon or the W muon and passing tight ID*/
    std::vector<reco::MuonRef> thirdMuons = 
      Common::getTightPFIsolatedRecoMuons(pAllMuons, Common::getPrimaryVertex(pVtx), 
					  muonPFIsoPUSubtractionCoeff_, -1.0, -1.0, true);
    bool passThirdMuonCut = true;
    std::vector<double> chargeProduct;
    for (std::vector<reco::MuonRef>::const_iterator i3rdMu = thirdMuons.begin(); 
	 i3rdMu != thirdMuons.end(); ++i3rdMu) {
      if ((i3rdMu->key() != removedMuonRefs[removedMuonRefs.size() - 1].key()) && 
	  (i3rdMu->key() != WMuonRefs[WMuonRefs.size() - 1].key())) {
	chargeProduct.push_back(WMuonRefs[WMuonRefs.size() - 1]->charge()*(*i3rdMu)->charge());
	passThirdMuonCut = passThirdMuonCut && (chargeProduct[chargeProduct.size() - 1] == 1.0);
      }
    }

    //impose pT and decay mode cut on tau
    if (((*iTau)->pt() > tauPTMin_) && 
	((tauDecayMode_ == reco::PFTau::kNull) || ((*iTau)->decayMode() == tauDecayMode_)) && 
	(b_discriminant < 0.679)) {

      //plot the number of good vertices
      nGoodVtx_->Fill(Common::numGoodVertices(pVtx)/*trueNInt*/, PUWeight);

      //plot MET distribution
      fillETHistogram(pMET, MET_, PUWeight);

      //plot dPhi(highest pT W muon, MET)
      plotDPhiMuMet(WMuonRefs[WMuonRefs.size() - 1], pMET, dPhiWMuMET_, PUWeight);

      //plot the transverse mass for the highest pT W muon
      plotMT(WMuonRefs[WMuonRefs.size() - 1], pMET, WMuMT_, PUWeight);

      //plot W muon MT versus MET
      edm::RefToBase<reco::PFMET> METRefToBase = pMET->refAt(0);
      double MT = sqrt(2*WMuonRefs[WMuonRefs.size() - 1]->pt()*METRefToBase->et()*(1.0 - cos(reco::deltaPhi(WMuonRefs[WMuonRefs.size() - 1]->phi(), METRefToBase->phi()))));
      double met = METRefToBase->et();
      WMuMTVsMET_->Fill(met, MT, PUWeight);

      //plot multiplicity of muons associated to this tau
      hadTauAssociatedMuMultiplicity_->Fill(removedMuons.size(), PUWeight);

      //plot the mu + tau invariant mass for the highest pT muon
      muHadMass_->Fill(muHadMass, PUWeight*tauHadPTWeight);

      //track maximum mu+had mass
      if (muHadMass > maxMuHadMass_) maxMuHadMass_ = muHadMass;

      /*plot the mu + tau invariant mass for the highest pT muon weighted by the statistical error 
	on the hadronic tau pT weight for this event*/
      muHadMassReweightErrSq_->Fill(muHadMass, tauHadPTWeightErr);

      //plot the mu + tau charge for the highest pT muon
      muHadCharge_->
	Fill(removedMuonRefs[removedMuonRefs.size() - 1]->charge() + (*iTau)->charge(), PUWeight);

      //plot dPhi(highest pT muon, MET)
      plotDPhiMuMet(removedMuonRefs[removedMuonRefs.size() - 1], pMET, dPhiTauMuMET_, PUWeight);

      //plot the transverse mass for the highest pT muon
      plotMT(removedMuonRefs[removedMuonRefs.size() - 1], pMET, tauMuMT_, PUWeight);

      //plot the transverse mass for the hadronic tau
      double tauHad_TransverseMass = 
	sqrt(2*(*iTau)->pt()*METRefToBase->et()*
	     (1.0 - cos(reco::deltaPhi((*iTau)->phi(), METRefToBase->phi()))));
      tauHadMT_->Fill(tauHad_TransverseMass, PUWeight);

      //plot HT (tau muon + hadronic tau + leading distinct corrected jet)
      reco::Candidate* leadDistinctCorrJet = 
	dynamic_cast<reco::Candidate*>(const_cast<reco::PFJet*>
				       (oldJetRefsExclTauNoPTCut
					[oldJetRefsExclTauNoPTCut.size() - 1].get()));
      std::vector<reco::Candidate*> tauMuTauHadJet;
      tauMuTauHadJet.push_back(dynamic_cast<reco::Candidate*>
			       (const_cast<reco::Muon*>
				(removedMuonRefs[removedMuonRefs.size() - 1].get())));
      tauMuTauHadJet.push_back(dynamic_cast<reco::Candidate*>
			       (const_cast<reco::PFTau*>(((*iTau).get()))));
      tauMuTauHadJet.push_back(leadDistinctCorrJet);
      plotHT(tauMuTauHadJet, tauMuTauHadJetHT_, PUWeight);

      //plot HT (tau muon + hadronic tau + leading distinct corrected jet + W muon)
      tauMuTauHadJet.push_back(dynamic_cast<reco::Candidate*>
			       (const_cast<reco::Muon*>(WMuonRefs[WMuonRefs.size() - 1].get())));
      plotHT(tauMuTauHadJet, tauMuTauHadJetWMuHT_, PUWeight);

      //plot HT (tau muon + hadronic tau + leading distinct corrected jet + W muon + MET)
      tauMuTauHadJet.push_back(dynamic_cast<reco::Candidate*>
			       (const_cast<reco::PFMET*>(&*(pMET->refAt(0)))));
      plotHT(tauMuTauHadJet, tauMuTauHadJetWMuMETHT_, PUWeight);

      //plot HT (two leading corrected jets if both exist)
      std::vector<reco::Candidate*> diJet;
      if (oldJetRefsInclTauNoPTCut.size() > 0) {
	diJet.push_back(dynamic_cast<reco::Candidate*>
			(const_cast<reco::PFJet*>
			 (oldJetRefsInclTauNoPTCut[oldJetRefsInclTauNoPTCut.size() - 1].get())));
      }
      if (oldJetRefsInclTauNoPTCut.size() > 1) {
	diJet.push_back(dynamic_cast<reco::Candidate*>
			(const_cast<reco::PFJet*>
			 (oldJetRefsInclTauNoPTCut[oldJetRefsInclTauNoPTCut.size() - 2].get())));
      }
      plotHT(diJet, diJetHT_, PUWeight);

      //plot HT (two leading corrected jets + W muon)
      diJet.push_back(dynamic_cast<reco::Candidate*>
		      (const_cast<reco::Muon*>(WMuonRefs[WMuonRefs.size() - 1].get())));
      plotHT(diJet, diJetWMuHT_, PUWeight);

      //plot HT (corrected jet associated to hadronic tau + leading distinct corrected jet)
      std::vector<reco::Candidate*> jetTauJet;
      jetTauJet.push_back(dynamic_cast<reco::Candidate*>
			  (const_cast<reco::PFJet*>(&correctedTauJet)));
      jetTauJet.push_back(leadDistinctCorrJet);
      plotHT(jetTauJet, jetTauJetHT_, PUWeight);

      /*plot HT (corrected jet associated to hadronic tau + leading distinct corrected jet + W 
	muon)*/
      jetTauJet.push_back(dynamic_cast<reco::Candidate*>
			  (const_cast<reco::Muon*>(WMuonRefs[WMuonRefs.size() - 1].get())));
      plotHT(jetTauJet, jetTauJetWMuHT_, PUWeight);

      //get the nearest gen W decay muon to the reco W muon
      int nearestGenWMuIndex = -1;
      const reco::GenParticle* nearestGenWMu = NULL;
      if (MC_) {
	nearestGenWMu = Common::nearestObject(WMuonRefs[WMuonRefs.size() - 1], genParticlePtrs, 
					      nearestGenWMuIndex);
      }

      //get the nearest gen muon from a1-->tau(-->mu)tau(-->had) decay to the soft muon
      int iNearestGenTauMu = -1;
      const reco::GenParticle* nearestGenTauMu = NULL;
      if (MC_) {
	nearestGenTauMu = Common::nearestObject(removedMuonRefs[removedMuonRefs.size() - 1], 
						genTauMuPtrs, iNearestGenTauMu);
      }

      //categorize based on the gen matches to the W and tau muons
      int genMuInfoVal = -1;
      if ((iNearestGenTauMu != -1) && (nearestGenWMuIndex != -1)) {
	const double dRWMuNearestGenWMu = 
	  reco::deltaR(*nearestGenWMu, *WMuonRefs[WMuonRefs.size() - 1]);
	const double dRTauMuNearestGenTauMu = 
	  reco::deltaR(*nearestGenTauMu, *removedMuonRefs[removedMuonRefs.size() - 1]);
	if (dRWMuNearestGenWMu < dR_) {
	  if (dRTauMuNearestGenTauMu < dR_) genMuInfoVal = 0; //correct assignment
	  else genMuInfoVal = 1; //W muon correct, gen tau muon not reconstructed
	}
	else if (dRTauMuNearestGenTauMu < dR_) genMuInfoVal = 2; /*tau muon correct, gen W muon 
								   not reconstructed*/
	else { //incorrect assignment
	  int iNearestGenTauMuToWMu = -1;
	  int iNearestGenWMuToTauMu = -1;
	  const reco::GenParticle* nearestGenTauMuToWMu = 
	    Common::nearestObject(WMuonRefs[WMuonRefs.size() - 1], genTauMuPtrs, 
				  iNearestGenTauMuToWMu);
	  const reco::GenParticle* nearestGenWMuToTauMu = 
	    Common::nearestObject(removedMuonRefs[removedMuonRefs.size() - 1], genParticlePtrs, 
				  iNearestGenWMuToTauMu);
	  const double dRWMuNearestGenTauMu = 
	    reco::deltaR(*nearestGenTauMuToWMu, *WMuonRefs[WMuonRefs.size() - 1]);
	  const double dRTauMuNearestGenWMu = 
	    reco::deltaR(*nearestGenWMuToTauMu, *removedMuonRefs[removedMuonRefs.size() - 1]);
	  if (dRWMuNearestGenTauMu < dR_) {
	    if (dRTauMuNearestGenWMu < dR_) genMuInfoVal = 3; //switched assignment
	    else genMuInfoVal = 4; //tau-->mu reconstructed as W muon, W muon not reconstructed
	  }
	  else if (dRTauMuNearestGenWMu < dR_) genMuInfoVal = 5; /*W muon reconstructed as tau 
								   muon, tau muon not 
								   reconstructed*/
	  else genMuInfoVal = 6; //neither lepton reconstructed
	}
      }

      //plot information about muon gen matching
      genMuInfo_->Fill(genMuInfoVal, PUWeight);

      //plot dR(W muon, leading soft muon per jet) when mu+had mass > 2 GeV
      if (muHadMass > 2.0/*GeV*/) {
	dRWMuSoftMuMuHadMassGe2_->Fill(reco::deltaR(*WMuonRefs[WMuonRefs.size() - 1], 
						    *removedMuonRefs[removedMuonRefs.size() - 1]), 
				       PUWeight);
      }
    
      //plot tau muon pT
      const double tauMuPT = removedMuonRefs[removedMuonRefs.size() - 1]->pt();
      tauMuPT_->Fill(tauMuPT, PUWeight);

      //plot hadronic tau pT
      tauHadPT_->Fill((*iTau)->pt(), PUWeight);

      //track maximum hadronic tau pT
      if ((*iTau)->pt() > maxTauHadPT_) maxTauHadPT_ = (*iTau)->pt();

      //plot hadronic tau deltaBeta-corrected isolation energy
      tauHadIso_->Fill((*pTauHadIso)[*iTau], PUWeight);

      //plot W mu PF isolation vs hadronic tau isolation energy
      double WMuIso = 
	Common::getMuonCombPFIso(*(WMuonRefs[WMuonRefs.size() - 1]), muonPFIsoPUSubtractionCoeff_);
      double WMuRelIso = WMuIso/WMuonRefs[WMuonRefs.size() - 1]->pt();
      WMuIso_->Fill(WMuRelIso, PUWeight);
      WMuIsoVsTauHadIso_->Fill((*pTauHadIso)[*iTau], WMuRelIso, PUWeight);

      //plot W muon PF relative isolation, counting leptons in the isolation sum
      double WMuLeptonRelIso = Common::getMuonLeptonPFIso(*(WMuonRefs[WMuonRefs.size() - 1]));
      WMuLeptonRelIso_->Fill(WMuLeptonRelIso/WMuIso, PUWeight);

      //fill vector of reco::Muon pointers of muons not overlapping the tau muon
      std::vector<reco::Muon*> muPtrs;
      for (reco::MuonCollection::const_iterator iMu = pAllMuons->begin(); iMu != pAllMuons->end(); 
	   ++iMu) {
	if ((iMu - pAllMuons->begin()) != removedMuonRefs[removedMuonRefs.size() - 1].key()) {
	  muPtrs.push_back(const_cast<reco::Muon*>(&*iMu));
	}
      }

      //find nearest muon to soft muon
      int nearestMuIndex = -1;
      const reco::Muon* nearestMu = 
	Common::nearestObject(removedMuonRefs[removedMuonRefs.size() - 1], muPtrs, nearestMuIndex);

      //find nearest gen muon to soft muon
      int nearestGenMuIndex = -1;
      const reco::Candidate* nearestGenMu = NULL;
      if (MC_) {
	nearestGenMu = 
	  Common::nearestObject(removedMuonRefs[removedMuonRefs.size() - 1], muGenParticlePtrs, 
				nearestGenMuIndex);
      }

      //find nearest gen Z decay muon to soft muon
      int nearestGenZOrTTMuIndex = -1;
      const reco::Candidate* nearestGenZOrTTMu = NULL;
      if (MC_) {
	nearestGenZOrTTMu = 
	  Common::nearestObject(removedMuonRefs[removedMuonRefs.size() - 1], 
				ZOrTTMuGenParticlePtrs, nearestGenZOrTTMuIndex);
      }

      //find nearest gen Z decay muon that radiated to soft muon
      int nearestGenZOrTTMuFSRIndex = -1;
      const reco::Candidate* nearestGenZOrTTMuFSR = NULL;
      if (MC_) {
	nearestGenZOrTTMuFSR = Common::nearestObject(removedMuonRefs[removedMuonRefs.size() - 1], 
						     ZOrTTMuFSRGenParticlePtrs, 
						     nearestGenZOrTTMuFSRIndex);
      }

      //plot soft muon nearest muon properties
      double dRSoftMuNearestMu = -1.0;
      double dRSoftMuNearestGenMu = -1.0;
      int xBinVal = -1;
      int yBinVal = -1;
      if (nearestMuIndex != -1) {
	dRSoftMuNearestMu = reco::deltaR(*removedMuonRefs[removedMuonRefs.size() - 1], *nearestMu);
      }
      if (nearestGenMuIndex != -1) {
	dRSoftMuNearestGenMu = 
	  reco::deltaR(*removedMuonRefs[removedMuonRefs.size() - 1], *nearestGenMu);
      }
      if ((dRSoftMuNearestMu >= 0.0) && (dRSoftMuNearestMu < dR_)) {
	if (nearestMuIndex == (int)WMuonRefs[WMuonRefs.size() - 1].key()) {
	  xBinVal = 0;
	  if ((dRSoftMuNearestGenMu >= 0.0) && (dRSoftMuNearestGenMu < dR_)) yBinVal = 0;
	  else yBinVal = 1;
	}
	else {
	  xBinVal = 1;
	  if ((dRSoftMuNearestGenMu >= 0.0) && (dRSoftMuNearestGenMu < dR_)) yBinVal = 0;
	  else yBinVal = 1;
	}
      }
      else if (MC_) {
	xBinVal = 2;
	if ((dRSoftMuNearestGenMu >= 0.0) && (dRSoftMuNearestGenMu < dR_)) yBinVal = 0;
	else yBinVal = 1;
      }
      if (MC_) genMuExistsVsSoftMuNearestMuProperties_->Fill(xBinVal, yBinVal, PUWeight);

      /*plot charge product of W muon and third tight muon
	- one entry per additional muon past the W and tau muons)
	- 3rd muon not overlapping the tau muon or the W muon and passing tight ID
	- if no 3rd muon exists, fill one zero*/
      for (std::vector<double>::const_iterator iChargeProduct = chargeProduct.begin(); 
	   iChargeProduct != chargeProduct.end(); ++iChargeProduct) {
	WMu3rdTightMuChargeProduct_->Fill(*iChargeProduct, PUWeight);
      }
      if (chargeProduct.size() == 0) WMu3rdTightMuChargeProduct_->Fill(0.0, PUWeight);

      //plot hadronic tau eta
      tauHadEta_->Fill((*iTau)->eta(), PUWeight*tauHadPTWeight);

      //plot pT(soft muon)/(mu+had mass)
      softMuPTOverMuHadMass_->Fill(removedMuonRefs[removedMuonRefs.size() - 1]->pt()/muHadMass, 
				   PUWeight);

      //plot pT(mu+had)/(mu+had mass) for a given uncorrected jet pT cut
      double muHadPT = (removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).pt();
      const double muHadPTOverM = tauOldJetRef->pt()/tauOldJetRef->mass();
      if (tauOldJetRef->pt() > uncorrJetPTMin_) {
	muHadPTOverMuHadMass_->Fill(muHadPTOverM, PUWeight);
      }

      //plot dR(soft muon, nearest gen muon)
      if (MC_) dRSoftMuNearestGenMuHist_->Fill(dRSoftMuNearestGenMu, PUWeight);

      //plot mu+had pT
      muHadPT_->Fill(muHadPT, PUWeight);

      //plot invariant mass of W muon and tau muon
      mWMuTauMu_->Fill(mWMuTauMu, PUWeight);

      /*plot the PDG ID of the nearest status 1 particle to the soft muon and the PDG ID of the 
	source particle of the gen muon nearest to the soft muon*/
      int nearestStatus1GenParticleIndex = -1;
      const reco::GenParticle* nearestStatus1GenParticle = NULL;
      if (MC_) {
	nearestStatus1GenParticle = 
	  Common::nearestObject(removedMuonRefs[removedMuonRefs.size() - 1], 
				status1GenParticlePtrs, nearestStatus1GenParticleIndex);
      }
      if (nearestStatus1GenParticle != NULL) {
	if (reco::deltaR(*nearestStatus1GenParticle, 
			 *removedMuonRefs[removedMuonRefs.size() - 1]) < dR_) {
	  const unsigned int absNearestStatus1GenParticlePDGID = 
	    fabs(nearestStatus1GenParticle->pdgId());
	  unsigned int nearestGenParticleVal = 5;
	  unsigned int muSrcVal = 1;
	  const reco::Candidate* ancestor = nearestStatus1GenParticle->mother();
	  switch (absNearestStatus1GenParticlePDGID) {
	  case GenTauDecayID::EPDGID:
	    nearestGenParticleVal = 0;
	    break;
	  case GenTauDecayID::MUPDGID:
	    nearestGenParticleVal = 1;
	    while ((fabs(ancestor->pdgId()) != GenTauDecayID::B) && 
		   (fabs(ancestor->pdgId()) != GenTauDecayID::C) && 
		   (ancestor->numberOfMothers() > 0)) {
	      ancestor = ancestor->mother();
	    }
	    if ((fabs(ancestor->pdgId()) == GenTauDecayID::B) || 
		(fabs(ancestor->pdgId()) == GenTauDecayID::C)) muSrcVal = 0;
	    break;
	  case GenTauDecayID::ENEUTRINOPDGID:
	    nearestGenParticleVal = 2;
	    break;
	  case GenTauDecayID::MUNEUTRINOPDGID:
	    nearestGenParticleVal = 2;
	    break;
	  case GenTauDecayID::TAUNEUTRINOPDGID:
	    nearestGenParticleVal = 2;
	    break;
	  case GenTauDecayID::GAMMAPDGID:
	    nearestGenParticleVal = 3;
	    break;
	  default:
	    break;
	  }
	  if (absNearestStatus1GenParticlePDGID >= 111) nearestGenParticleVal = 4;
	  PDGIDNearestStatus1GenParticleToSoftMu_->Fill(nearestGenParticleVal/*, PUWeight*/);
	  PDGIDMuSrc_->Fill(muSrcVal/*, PUWeight*/);
	}
	else PDGIDNearestStatus1GenParticleToSoftMu_->Fill(6/*, PUWeight*/);
      }
      else if (MC_) PDGIDNearestStatus1GenParticleToSoftMu_->Fill(-1/*, PUWeight*/);

//       //debug
//       std::cerr << "All jets not overlapping the W muon\n";
//       for (std::vector<reco::PFJetRef>::const_iterator iJet = correctedOldJetRefs.begin(); 
//       	   iJet != correctedOldJetRefs.end(); ++iJet) {
//       	std::cerr << iJet->key() << " " << (*iJet)->pt() << " " << (*iJet)->eta() << " ";
//       	std::cerr << (*iJet)->phi() << std::endl;
//       }
//       std::cerr << "Tau jet\n";
//       std::cerr << tauOldJetRef.key() << " " << tauOldJetRef->pt() << " ";
//       std::cerr << tauOldJetRef->eta() << " " << tauOldJetRef->phi() << std::endl;

      //plot pT rank of mu+had parent uncleaned jet
      Common::sortByPT(correctedOldJetRefs);
      std::reverse(correctedOldJetRefs.begin(), correctedOldJetRefs.end());
      if (correctedOldJetRefs.size() > 0) {
	if (tauOldJetRef.key() == correctedOldJetRefs[0].key()) {
	  muHadUncleanedJetPTRank_->Fill(0.0, PUWeight); //compiler likes 0.0 but not 0
	}
	else if ((correctedOldJetRefs.size() > 1) && 
		 (tauOldJetRef.key() == correctedOldJetRefs[1].key())) {
	  muHadUncleanedJetPTRank_->Fill(0.0, PUWeight); //compiler likes 0.0 but not 0
	}
	else muHadUncleanedJetPTRank_->Fill(1, PUWeight);
      }
      else muHadUncleanedJetPTRank_->Fill(-1, PUWeight);

      /*plot the number of additional loose isolated tight muons with pT > 25 GeV in the event 
	(not counting the W muon)*/
      nAddlHardMuons_->Fill(WMuonRefs.size() - 1, PUWeight);

      /*fill vectors of the additional jets in the event
	vector #1
	- corrected pT > 30 GeV
	- |eta| < 2.4
	- not overlapping W muon and mu+had
	vector #2
	- corrected pT > 40 GeV
	- |eta| < 2.4
	- not overlapping W muon and mu+had*/
      std::vector<reco::PFJet*> addlJetsHBHEPTGeq30;
      std::vector<reco::PFJet*> addlJetsHBHEPTGeq40;
      for (reco::PFJetCollection::const_iterator iJet = pCorrJets->begin(); 
	   iJet != pCorrJets->end(); ++iJet) {
	if (iJet->pt() > 30.0/*GeV*/) {
	  addlJetsHBHEPTGeq30.push_back(const_cast<reco::PFJet*>(&*iJet));
	}
	if (iJet->pt() > 40.0/*GeV*/) {
	  addlJetsHBHEPTGeq40.push_back(const_cast<reco::PFJet*>(&*iJet));
	}
      }

      /*plot the number of additional jets in the event
	- corrected pT > 20 GeV
	- |eta| < 2.4
	- not overlapping W muon and mu+had*/
      nAddlJetsPTGeq20_->Fill(pCorrJets->size(), PUWeight);

      /*plot the number of additional jets in the event
	- corrected pT > 30 GeV
	- |eta| < 2.4
	- not overlapping W muon and mu+had*/
      nAddlJetsPTGeq30_->Fill(addlJetsHBHEPTGeq30.size(), PUWeight);

      /*plot the number of additional jets in the event
	- corrected pT > 40 GeV
	- |eta| < 2.4
	- not overlapping W muon and mu+had*/
      nAddlJetsPTGeq40_->Fill(addlJetsHBHEPTGeq40.size(), PUWeight);

      //plot hadronic tau decay mode
      int mode = -1;
      switch ((*iTau)->decayMode()) {
      case reco::PFTau::kOneProng0PiZero:
	mode = 0;
	break;
      case reco::PFTau::kOneProng1PiZero:
	mode = 1;
	break;
      case reco::PFTau::kOneProng2PiZero:
	mode = 2;
	break;
      case reco::PFTau::kThreeProng0PiZero:
	mode = 3;
	break;
      case reco::PFTau::kNull:
	mode = -1;
	break;
      default:
	mode = 4;
	break;
      }
      tauHadDecayMode_->Fill(mode, PUWeight);

      //plot dR(soft muon, nearest Z decay muon)
      double dRSoftMuNearestGenZOrTTMu = -1.0;
      if (nearestGenZOrTTMuIndex != -1) {
	dRSoftMuNearestGenZOrTTMu = reco::deltaR(*removedMuonRefs[removedMuonRefs.size() - 1], 
						 *nearestGenZOrTTMu);
      }
      dRSoftMuNearestGenZOrTTMu_->Fill(dRSoftMuNearestGenZOrTTMu, PUWeight);

      //plot dR(soft muon, nearest Z decay muon that radiated)
      double dRSoftMuNearestGenZOrTTMuFSR = -1.0;
      if (nearestGenZOrTTMuFSRIndex != -1) {
	dRSoftMuNearestGenZOrTTMuFSR = reco::deltaR(*removedMuonRefs[removedMuonRefs.size() - 1], 
						    *nearestGenZOrTTMuFSR);
      }
      dRSoftMuNearestGenZOrTTMuFSR_->Fill(dRSoftMuNearestGenZOrTTMuFSR, PUWeight);

      //calculate the hadronic tau vector excluding photons
      reco::LeafCandidate::LorentzVector tauVecExcludingPhotons;
      const reco::PFCandidateRefVector& tauHadSigChargedHadrons = 
	(*iTau)->signalPFChargedHadrCands();
      const reco::PFCandidateRefVector& tauHadSigNeutralHadrons = 
	(*iTau)->signalPFNeutrHadrCands();
      for (reco::PFCandidateRefVector::const_iterator 
	     iTauHadSigChargedHadron = tauHadSigChargedHadrons.begin(); 
	   iTauHadSigChargedHadron != tauHadSigChargedHadrons.end(); ++iTauHadSigChargedHadron) {
	tauVecExcludingPhotons+=(*iTauHadSigChargedHadron)->p4();
      }
      for (reco::PFCandidateRefVector::const_iterator 
	     iTauHadSigNeutralHadron = tauHadSigNeutralHadrons.begin(); 
	   iTauHadSigNeutralHadron != tauHadSigNeutralHadrons.end(); ++iTauHadSigNeutralHadron) {
	tauVecExcludingPhotons+=(*iTauHadSigNeutralHadron)->p4();
      }

      /*plot photon energy fraction of hadronic tau and the opening angle between the photon and 
	the other hadronic tau constituents*/
      const reco::PFCandidateRefVector& tauHadSigPhotons = (*iTau)->signalPFGammaCands();
      double tauHadPhotonEnergy = 0.0;
      for (reco::PFCandidateRefVector::const_iterator iTauHadSigPhoton = tauHadSigPhotons.begin(); 
	   iTauHadSigPhoton != tauHadSigPhotons.end(); ++iTauHadSigPhoton) {
	tauHadPhotonEnergy+=(*iTauHadSigPhoton)->et();
	dThetaPhotonOtherTauConstituents_->
	  Fill(angle((*iTauHadSigPhoton)->p4(), tauVecExcludingPhotons), PUWeight);
      }
      tauHadPhotonEnergyFraction_->Fill(tauHadPhotonEnergy/(*iTau)->et(), PUWeight);

      //plot dPhi(W muon, leading soft muon per jet)
//       double PhiW = (*WMuonRefs[WMuonRefs.size()-1]).phi();
//       double PhiTau = (*removedMuonRefs[removedMuonRefs.size() - 1]).phi();
      double dPhi = 
	/*PhiW - PhiTau*/fabs(reco::deltaPhi(WMuonRefs[WMuonRefs.size()-1]->phi(), 
					     removedMuonRefs[removedMuonRefs.size() - 1]->phi()));
//       const double Pi = 3.14159265359;
//       while (dPhi > Pi)
// 	dPhi -= 2.0*Pi;
//       while (dPhi < -1.0*Pi)
// 	dPhi += 2.0*Pi;
      dPhiWMuSoftMu_->Fill(dPhi, PUWeight);
//       double Mass_WMuTauMu = (WMuonRefs[WMuonRefs.size() - 1]->p4() + 
// 			      removedMuonRefs[removedMuonRefs.size() - 1]->p4()).M();
      if (/*Mass_WMuTauMu*/mWMuTauMu > 20.)
	dPhiWMuSoftMu_withCut_->Fill(dPhi, PUWeight);

      //plot cleaned jet pT vs. cleaned tau pT
      cleanedJetPTVsCleanedTauPT_->Fill((*iTau)->pt(), tauJetRef->pt(), PUWeight);

      //plot uncleaned jet pT vs. cleaned tau pT
      uncleanedJetPTVsCleanedTauPT_->Fill((*iTau)->pt(), tauOldJetRef->pt(), PUWeight);

      //calculate N-subjettiness of cleaned jet
            fastjet::Pruner pruner(fastjet::kt_algorithm, zCut_, RcutFactor_);
      vector<reco::PFCandidatePtr> pfCands = tauJetRef->getPFConstituents();
//       vector<reco::PFCandidatePtr> pfCands = tauOldJetRef->getPFConstituents();
      vector<const reco::PFCandidate*> all_particles;
      for (unsigned j = 0; j < pfCands.size(); j++){
	const reco::PFCandidate *thisPF = pfCands.at(j).get(); 
	all_particles.push_back( thisPF );	
      }
      vector<fastjet::PseudoJet> FJparticles;
      for (unsigned particle = 0; particle < all_particles.size(); particle++) {
	const reco::PFCandidate *thisParticle = all_particles.at(particle);
	FJparticles.push_back(fastjet::PseudoJet(thisParticle->px(),thisParticle->py(),
						 thisParticle->pz(),thisParticle->energy()));
      }
      fastjet::JetDefinition jet_def(fastjet::kt_algorithm, 0.5);
      fastjet::ClusterSequence thisClustering(FJparticles, jet_def);
      vector<fastjet::PseudoJet> FJjet = thisClustering.inclusive_jets();
      fastjet::PseudoJet thisMainJet = FJjet.at(0);
      fastjet::PseudoJet thisGroomedJet = pruner(thisMainJet); // jet grooming
      vector<fastjet::PseudoJet> FJparticles2 = thisGroomedJet.constituents();
//       vector<fastjet::PseudoJet> FJparticles2 = thisMainJet.constituents(); // if you don't want to prune
//       NsubParameters paraNsub = NsubParameters(1.0, 0.8);
      Nsubjettiness routine1(1, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
      Nsubjettiness routine2(2, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
      Nsubjettiness routine3(3, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
      Nsubjettiness routine4(4, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
      double tau1 = routine1.result(thisGroomedJet);
      double tau2 = routine2.result(thisGroomedJet);
      double tau3 = routine3.result(thisGroomedJet);
//       double tau4 = routine4.result(thisGroomedJet);
//       double tau1 = routine1.result(thisMainJet); // if you don't want to prune
//       double tau2 = routine2.result(thisMainJet); // if you don't want to prune
//       double tau3 = routine3.result(thisMainJet); // if you don't want to prune
//       double tau4 = routine4.result(thisMainJet); // if you don't want to prune
      if (tau1 != 0.)
	{
	  muHad_t3t1_->Fill(tau3/tau1, PUWeight);
	  muHad_t2t1_->Fill(tau2/tau1, PUWeight);
	  muHad_t3t1Vsptmj_->Fill(tauJetRef->pt()/tauJetRef->mass(), tau3/tau1, PUWeight);
	  //muHad_t3t1Vsptmj_->Fill(tauOldJetRef->pt()/tauOldJetRef->mass(), tau3/tau1, PUWeight);

	  if (tauJetRef->pt() >= 10. && tauJetRef->pt() < 20.)
	    muHad_t3t1_pT1020_->Fill(tau3/tau1, PUWeight);
	  else if (tauJetRef->pt() >= 20. && tauJetRef->pt() < 30.)
	    muHad_t3t1_pT2030_->Fill(tau3/tau1, PUWeight);
	  else if (tauJetRef->pt() >= 30. && tauJetRef->pt() < 40.)
	    muHad_t3t1_pT3040_->Fill(tau3/tau1, PUWeight);
	  else if (tauJetRef->pt() >= 40. && tauJetRef->pt() < 50.)
	    muHad_t3t1_pT4050_->Fill(tau3/tau1, PUWeight);
	  else if (tauJetRef->pt() > 50.)
	    muHad_t3t1_pT50Up_->Fill(tau3/tau1, PUWeight);

	  if (NJets == 0)
	    muHad_t3t1_0Jets_->Fill(tau3/tau1, PUWeight);
	  else if (NJets == 1)
	    muHad_t3t1_1Jets_->Fill(tau3/tau1, PUWeight);
	  else if (NJets == 2)
	    muHad_t3t1_2Jets_->Fill(tau3/tau1, PUWeight);
	  else if (NJets == 3)
	    muHad_t3t1_3Jets_->Fill(tau3/tau1, PUWeight);
	  else if (NJets == 4)
	    muHad_t3t1_4Jets_->Fill(tau3/tau1, PUWeight);
	  else if (NJets == 5)
	    muHad_t3t1_5Jets_->Fill(tau3/tau1, PUWeight);
	  else if (NJets > 5)
	    muHad_t3t1_MoreJets_->Fill(tau3/tau1, PUWeight);

	  if (tau3/tau1 < 0.02)
	    { // if tau3/tau1 < 0.02
// 	      cout << "FLAG: tau3/tau1 = " << tau3/tau1 << endl;
// 	      cout << "FLAG: tau3 = " << tau3 << endl;
// 	      cout << "FLAG: no. of tracks in jet = " << tauJetRef->getTrackRefs().size() << endl;
// 	      cout << "FLAG: no. of tracks in jet = " << tauOldJetRef->getTrackRefs().size();
// 	      cout << endl;
// 	      cout << "FLAG: no. of groomed jet constituents = " << FJparticles2.size() << endl;
	      std::vector<reco::PFCandidatePtr> tauJetCands = tauJetRef->getPFConstituents();
// 	      std::vector<reco::PFCandidatePtr> tauJetCands = tauOldJetRef->getPFConstituents();
	      for (std::vector<fastjet::PseudoJet>::const_iterator fjconst = FJparticles2.begin(); 
		   fjconst != FJparticles2.end(); fjconst++)
		{ // loop over FJparticles2
		  double minDelR = 10000.;
// 		  unsigned int tauJetCandID = 0;
		  for (std::vector<reco::PFCandidatePtr>::const_iterator jetcand = 
			 tauJetCands.begin(); jetcand != tauJetCands.end(); jetcand++)
		    {
		      double delR = deltaR(fjconst->eta(), fjconst->phi(), (**jetcand).eta(), 
					   (**jetcand).phi());
		      if (delR < minDelR)
			{
			  minDelR = delR;
// 			  tauJetCandID = (**jetcand).particleId();
			}
		    }
		  if (minDelR < 0.3)
		    {
// 		      cout << "delR between groomed jet constituent and jet PFCandidate = " << minDelR << endl;
// 		      cout << "Particle ID of jet PFCandidate = " << tauJetCandID << endl;
		    }
		} // loop over FJparticles2
	    } // if tau3/tau1 < 0.02
	  muHad_t3t1VsDecayMode_->Fill((*iTau)->decayMode(), tau3/tau1, PUWeight);
	}
      muHad_t3t1_->Fill(tau3/tau1, PUWeight);
      muHad_t2t1_->Fill(tau2/tau1, PUWeight);

      //plot number of charged tracks of cleaned jet
      reco::TrackRefVector jetTracks = tauJetRef->getTrackRefs();
      double Nchtrk_0 = 0.; // no cut on track pT
      double Nchtrk_1 = 0.; // track pT > 1 GeV
      double Nchtrk_10 = 0.;
      double Nchtrk_30 = 0.; // track pT > 30 GeV
      for (reco::TrackRefVector::const_iterator iTrack = jetTracks.begin(); 
	   iTrack != jetTracks.end(); ++iTrack)     
	{ // loop over tracks
	  if ((*iTrack)->charge() != 0.)
	    {
	      Nchtrk_0 += 1.;
	      if ((*iTrack)->pt() > 1.)
		{
		  Nchtrk_1 += 1.;
		}
	      if ((*iTrack)->pt() > 10.)
		{
		  Nchtrk_10 += 1.;
		}
	      if ((*iTrack)->pt() > 30.)
		{
		  Nchtrk_30 += 1.;
		}
	    }
	} // loop over tracks
      muHad_Nchtrk_0_->Fill(Nchtrk_0, PUWeight);
      muHad_Nchtrk_1_->Fill(Nchtrk_1, PUWeight);
      muHad_Nchtrk_10_->Fill(Nchtrk_10, PUWeight);
      muHad_Nchtrk_30_->Fill(Nchtrk_30, PUWeight);

      //plot mu+had mass vs. dR(tagged soft muon, tau axis)
      const double dRSoftMuTau = 
	reco::deltaR(*removedMuonRefs[removedMuonRefs.size() - 1], **iTau);
      muHadMassVsDRSoftMuTau_->Fill(dRSoftMuTau, muHadMass, PUWeight);

      //plot mu+had mass vs. tagged soft muon pT
      muHadMassVsSoftMuPT_->Fill(tauMuPT, muHadMass, PUWeight);

      //plot hadronic tau isolation vs. soft muon pT
      tauHadIsoVsSoftMuPT_->Fill(tauMuPT, (*pTauHadIso)[*iTau], PUWeight);

      //plot mu+had mass vs. hadronic tau eta
      muHadMassVsTauHadEta_->Fill((*iTau)->eta(), muHadMass, PUWeight);

      //plot mu+had mass vs. soft muon eta
      muHadMassVsSoftMuEta_->Fill(removedMuonRefs[removedMuonRefs.size() - 1]->eta(), muHadMass, 
				  PUWeight);

      //plot mu+had mass vs. hadronic tau isolation
      muHadMassVsTauHadIso_->Fill((*pTauHadIso)[*iTau], muHadMass, PUWeight);

      //plot mu+had mass vs. hadronic tau pT
      muHadMassVsTauHadPT_->Fill((*iTau)->pt(), muHadMass, PUWeight);

      //plot hadronic tau isolation vs. eta
      tauHadIsoVsEta_->Fill((*iTau)->eta(), (*pTauHadIso)[*iTau], PUWeight);

      //plot hadronic tau eta vs. soft muon eta (sanity check)
      tauHadEtaVsSoftMuEta_->Fill(removedMuonRefs[removedMuonRefs.size() - 1]->eta(), 
				  (*iTau)->eta(), PUWeight);

      //plot histogram of dEta(hadronic tau, soft muon) vs. dPhi(hadronic tau, soft muon)
      dEtaTauHadSoftMuVsDPhiTauHadSoftMu_->
	Fill(fabs(reco::deltaPhi(**iTau, *removedMuonRefs[removedMuonRefs.size() - 1])), 
	     fabs((*iTau)->eta() - removedMuonRefs[removedMuonRefs.size() - 1]->eta()), PUWeight);

      //plot histogram of pT(hadronic tau)/(mu+had mass) vs. hadronic tau isolation
      tauHadPTOverMuHadMassVsTauHadIso_->Fill((*pTauHadIso)[*iTau], (*iTau)->pt()/muHadMass, 
					      PUWeight);

      //plot histogram of pT(soft muon)/(mu+had mass) vs. hadronic tau isolation
      softMuPTOverMuHadMassVsTauHadIso_->
	Fill((*pTauHadIso)[*iTau], removedMuonRefs[removedMuonRefs.size() - 1]->pt()/muHadMass, 
	     PUWeight);

      /*plot histogram of (pT(soft muon) + pT(hadronic tau)/(2(mu+had mass)) vs. hadronic tau 
	isolation*/
      avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_->
	Fill((*pTauHadIso)[*iTau], 
	     (removedMuonRefs[removedMuonRefs.size() - 1]->pt() + (*iTau)->pt())/(2*muHadMass), 
	     PUWeight);

      /*plot histogram of pT(mu+had)/(mu+had mass) vs. hadronic tau isolation for pT(mu+had) > 30 
	GeV*/
      if (muHadPT > 30.0/*GeV*/) {
	muHadPTOverMuHadMassVsTauHadIso_->
	  Fill((*pTauHadIso)[*iTau], muHadPTOverM, PUWeight);
      }

      //plot soft muon pT vs. hadronic tau pT
      softMuPTVsTauHadPT_->Fill((*iTau)->pt(), removedMuonRefs[removedMuonRefs.size() - 1]->pt(), 
				PUWeight);

      //plot pT(mu+had parent jet)/m(mu+had parent jet) vs. m(W muon, soft muon)
      if (tauOldJetRef->pt() > uncorrJetPTMin_) {
	muHadPTOverMuHadMassVsMWMuSoftMu_->
	  Fill((WMuonRefs[WMuonRefs.size() - 1]->p4() + 
		removedMuonRefs[removedMuonRefs.size() - 1]->p4()).M(), 
	       muHadPTOverM, PUWeight);
      }

      //plot soft muon PF pT vs. RECO pT
      softMuPFPTVsRECOPT_->Fill(removedMuonRefs[removedMuonRefs.size() - 1]->pt(), 
				removedMuonRefs[removedMuonRefs.size() - 1]->pfP4().Pt(), 
				PUWeight);
      softMuPFEtaVsRECOEta_->Fill(removedMuonRefs[removedMuonRefs.size() - 1]->eta(), 
				  removedMuonRefs[removedMuonRefs.size() - 1]->pfP4().Eta(), 
				  PUWeight);

      //plot HT vs. MET
      tauMuTauHadJetWMuHTVsMET_->Fill(pMET->refAt(0)->et(), HT(tauMuTauHadJet), PUWeight);

      //plot dimuon invariant mass vs. mumutau invariant mass
      mWMuTauMuVsMWMuTauMuTauHad_->Fill(mWMuTauMuTauHad, mWMuTauMu, PUWeight);

      /*plot dimuon invariant mass vs. mumutau invariant mass for muons that radiated at the gen 
	level*/
      if (dRSoftMuNearestGenZOrTTMuFSR < dR_) {
	mWMuTauMuVsMWMuTauMuTauHadGenFSR_->Fill(mWMuTauMuTauHad, mWMuTauMu, PUWeight);
      }

      //plot mu+had mass vs. mumutau invariant mass
      muHadMassVsMWMuTauMuTauHad_->Fill(mWMuTauMuTauHad, muHadMass, PUWeight);

      //plot mu+had mass vs. Nj
      muHadMassVsNAddlJets_->Fill(pCorrJets->size(), muHadMass, PUWeight);

      //increment mu+had multiplicity counter
      ++nMuHad;

      //fill second jet plots
      if (oldJetRefsExcludingTau.size() > 0) {

	/*find the highest pT corrected AK5 PF jet in the original collection distinct from this 
	  tau and the W muon and passing |eta| and pT cuts*/
	const reco::PFJetRef& secondJetRef = 
	  oldJetRefsExcludingTau[oldJetRefsExcludingTau.size() - 1];

	// plot basic kinematic variables
	jet_eta->Fill(secondJetRef->eta(), PUWeight);
	jet_phi->Fill(secondJetRef->phi(), PUWeight);

	// fill etacut plots
	const double secondJetPTOverM = secondJetRef->pt()/secondJetRef->mass();
	jet_pt_etacut->Fill(secondJetRef->pt(), PUWeight);
	jet_mass_etacut->Fill(secondJetRef->mass(), PUWeight);
	jet_ptmj_etacut->Fill(secondJetPTOverM, PUWeight);

	//plot the invariant mass of the mu+had object and the hardest distinct jet
	const double mSecondJetMuHad = (removedMuonRefs[removedMuonRefs.size() - 1]->p4() + 
					(*iTau)->p4() + secondJetRef->p4()).M();
	mSecondJetMuHad_->Fill(mSecondJetMuHad, PUWeight);

	/*plot the invariant mass of the mu+had object and the hardest distinct jet vs. 
	  mu+had mass*/
	mSecondJetMuHadVsMuHadMass_->Fill(muHadMass, mSecondJetMuHad, PUWeight);

	//plot mu+had mass vs. second jet mass
	muHadMassVsSecondJetMass_->Fill(secondJetRef->mass(), muHadMass, PUWeight);

	//plot mu+had pT/m vs. second jet pT/m
	muHadPTOverMVsSecondJetPTOverM_->Fill(secondJetPTOverM, muHadPTOverM, PUWeight);

	//plot mu+had mass vs. second jet pT/m
	muHadMassVsSecondJetPTOverM_->Fill(secondJetPTOverM, muHadMass, PUWeight);

	/*plot Dphi(mu+had, second jet)
	  with fastjet stuff included for some reason deltaPhi with reco::PFJet arguments 
	  doesn't work*/
	dPhiMuHadSecondJet_->
	  Fill(fabs(reco::deltaPhi(tauOldJetRef->phi(), secondJetRef->phi())), PUWeight);

	// plot dPhi between W muon and second jet
// 	double Pi = 3.14159265359;
// 	double deltaphi = secondJetRef->phi() - WMuonRefs[WMuonRefs.size() - 1]->phi();
// 	while (deltaphi > Pi)
// 	  deltaphi -= 2.0*Pi;
// 	while (deltaphi < -1.0*Pi)
// 	  deltaphi += 2.0*Pi;
	dPhiWMuSecJet_->
	  Fill(/*deltaphi*/fabs(reco::deltaPhi(secondJetRef->phi(), 
					       WMuonRefs[WMuonRefs.size() - 1]->phi())), PUWeight);

	//       // plot N-subjettiness of second jet
	//       fastjet::Pruner pruner(fastjet::kt_algorithm, zCut_, RcutFactor_);
	//       vector<reco::PFCandidatePtr> pfCands = secondJetRef->getPFConstituents();
	//       vector<const reco::PFCandidate*> all_particles;
	//       for (unsigned j = 0; j < pfCands.size(); j++){
	// 	const reco::PFCandidate *thisPF = pfCands.at(j).get(); 
	// 	all_particles.push_back( thisPF );	
	//       }
	//       vector<fastjet::PseudoJet> FJparticles;
	//       for (unsigned particle = 0; particle < all_particles.size(); particle++) {
	// 	const reco::PFCandidate *thisParticle = all_particles.at(particle);
	// 	FJparticles.push_back(fastjet::PseudoJet(thisParticle->px(),thisParticle->py(),
	// 						 thisParticle->pz(),thisParticle->energy()));
	//       }
	//       fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.5);
	//       fastjet::ClusterSequence thisClustering(FJparticles, jet_def);
	//       vector<fastjet::PseudoJet> FJjet = thisClustering.inclusive_jets();
	//       fastjet::PseudoJet thisMainJet = FJjet.at(0);
	//       fastjet::PseudoJet thisGroomedJet = pruner(thisMainJet); // jet grooming
	//       vector<fastjet::PseudoJet> FJparticles2 = thisGroomedJet.constituents();
	//       NsubParameters paraNsub = NsubParameters(1.0, 0.8);
	//       Nsubjettiness routine1(1, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
	//       Nsubjettiness routine2(2, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
	//       Nsubjettiness routine3(3, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
	//       Nsubjettiness routine4(4, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
	//       double tau1 = routine1.result(thisGroomedJet);
	//       double tau2 = routine2.result(thisGroomedJet);
	//       double tau3 = routine3.result(thisGroomedJet);
	//       double tau4 = routine4.result(thisGroomedJet);
	      
	//       if (tau1 != 0.)
	// 	{
	// 	  muHad_t3t1_->Fill(tau3/tau1, PUWeight);
	// 	  muHad_t2t1_->Fill(tau2/tau1, PUWeight);
	// 	  muHad_t3t1Vsptmj_->Fill(tauJetRef->pt()/tauJetRef->mass(), tau3/tau1, PUWeight);
	// 	}

	// plot number of charged tracks in second jet
	reco::TrackRefVector secondjetTracks = secondJetRef->getTrackRefs();
	double second_Nchtrk_0 = 0.; // no cut on track pT
	double second_Nchtrk_1 = 0.; // track pT > 1 GeV
	double second_Nchtrk_10 = 0.;
	double second_Nchtrk_30 = 0.; // track pT > 30 GeV
	for (reco::TrackRefVector::const_iterator iTrack = secondjetTracks.begin(); 
	     iTrack != secondjetTracks.end(); ++iTrack)     
	  { // loop over tracks
	    if ((*iTrack)->charge() != 0.)
	      {
		second_Nchtrk_0 += 1.;
		if ((*iTrack)->pt() > 1.)
		  {
		    second_Nchtrk_1 += 1.;
		  }
		if ((*iTrack)->pt() > 10.)
		  {
		    second_Nchtrk_10 += 1.;
		  }
		if ((*iTrack)->pt() > 30.)
		  {
		    second_Nchtrk_30 += 1.;
		  }
	      }
	  } // loop over tracks
	second_Nchtrk_0_->Fill(second_Nchtrk_0, PUWeight);
	second_Nchtrk_1_->Fill(second_Nchtrk_1, PUWeight);
	second_Nchtrk_10_->Fill(second_Nchtrk_10, PUWeight);
	second_Nchtrk_30_->Fill(second_Nchtrk_30, PUWeight);
      }
    }

    //increment tau iterator
    ++iTau;
  }

  //fill an STL container of gen-matched muon refs
  std::vector<unsigned int> genMatchedMuonRefs;
  if (MC_) {
    for (reco::MuonRefVector::const_iterator iGenMatchedMu = pGenMatchedMuons->begin(); 
	 iGenMatchedMu != pGenMatchedMuons->end(); ++iGenMatchedMu) {
      genMatchedMuonRefs.push_back(iGenMatchedMu->key());
    }
  }

  //loop over cleaned jets
  for (reco::PFJetCollection::const_iterator iNewJet = pNewJets->begin(); 
       iNewJet != pNewJets->end(); ++iNewJet) {
    const reco::MuonRefVector& removedMuons = 
      (*pMuonJetMap)[reco::PFJetRef(pNewJets, iNewJet - pNewJets->begin())];

    //find the highest pT associated muon
    std::vector<reco::MuonRef> removedMuonRefs;
    for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
	 iMuon != removedMuons.end(); ++iMuon) { removedMuonRefs.push_back(*iMuon); }
    Common::sortByPT(removedMuonRefs);
    if (removedMuonRefs.size() > 0) {
      const double dR = reco::deltaR(*WMuonRefs[WMuonRefs.size() - 1], 
				     *removedMuonRefs[removedMuonRefs.size() - 1]);

      //plot dR(W muon, leading soft muon per jet)
      dRWMuSoftMu_->Fill(dR, PUWeight);

      //plot dR(W muon, leading soft gen-matched muon per jet)
      if (std::find(genMatchedMuonRefs.begin(), genMatchedMuonRefs.end(), 
		    removedMuonRefs[removedMuonRefs.size() - 1].key()) != 
	  genMatchedMuonRefs.end()) dRWMuSoftGenMatchedMu_->Fill(dR, PUWeight);
    }
  }

  //plot mu+had multiplicity
  muHadMultiplicity_->Fill(nMuHad, PUWeight);
} // analyzer


// ------------ method called once each job just before starting event loop  ------------
void TauAnalyzer::beginJob()
{
  //open output files
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book histograms
  MET_ = new TH1F("MET", ";#slash{E}_{T} (GeV);", 40, 0.0, 200.0);
  hadTauAssociatedMuMultiplicity_ = 
    new TH1F("hadTauAssociatedMuMultiplicity", ";N_{#mu}/#tau;", 2, 0.5, 2.5);
  muHadMass_ = 
    new TH1F("muHadMass", ";m_{#mu+had} (GeV);", muHadMassBins_.size() - 1, &muHadMassBins_[0]);
  muHadMassReweightErrSq_ = new TH1F("muHadMassReweightErrSq", ";m_{#mu+had} (GeV);", 
				     muHadMassBins_.size() - 1, &muHadMassBins_[0]);
  muHadCharge_ = new TH1F("muHadCharge", ";q_{#mu} + q_{had};", 5, -2.5, 2.5);
  bTagDiscrim_ = new TH1F("bTagDiscrim", ";b-tag discriminator;", 50, 0.0, 1.0);
  WMuIso_ = new TH1F("WMuIso", ";W muon PFRelIso;", 4000, 0.0, 40.0);
  WMuMT_ = new TH1F("WMuMT", ";W muon M_{T} (GeV);", 100, 0.0, 400.0);
  tauMuMT_ = new TH1F("tauMuMT", ";#tau muon M_{T} (GeV);", 50, 0.0, 200.0);
  tauHadMT_ = new TH1F("tauHadMT", ";#tau_{had} M_{T} (GeV);", 50, 0.0, 200.0);
  dPhiWMuMET_ = new TH1F("dPhiWMuMET", ";#Delta#phi(W muon, #slash{E}_{T});", 63, 0.0, 3.15);
  dPhiTauMuMET_ = 
    new TH1F("dPhiTauMuMET", ";#Delta#phi(#tau muon, #slash{E}_{T});", 63, 0.0, 3.15);
  tauMuTauHadJetHT_ = new TH1F("tauMuTauHadJetHT", ";H_{T} (GeV);", 40, 0.0, 1000.0);
  diJetHT_ = new TH1F("diJetHT", ";H_{T} (GeV);", 40, 0.0, 1000.0);
  jetTauJetHT_ = new TH1F("jetTauJetHT", ";H_{T} (GeV);", 40, 0.0, 1000.0);
  tauMuTauHadJetWMuHT_ = new TH1F("tauMuTauHadJetWMuHT", ";H_{T} (GeV);", 40, 0.0, 1000.0);
  tauMuTauHadJetWMuMETHT_ = new TH1F("tauMuTauHadJetWMuMETHT", ";H_{T} (GeV);", 40, 0.0, 1000.0);
  diJetWMuHT_ = new TH1F("diJetWMuHT", ";H_{T} (GeV);", 40, 0.0, 1000.0);
  jetTauJetWMuHT_ = new TH1F("jetTauJetWMuHT", ";H_{T} (GeV);", 40, 0.0, 1000.0);
  dRWMuSoftMu_ = new TH1F("dRWMuSoftMu", ";#DeltaR(W muon, soft muon);", 30, 0.0, 3.0);
  dPhiWMuSoftMu_ = new TH1F("dPhiWMuSoftMu", ";#Delta#phi(W muon, soft muon);", 400, -4.0, 4.0);
  dPhiWMuSoftMu_withCut_ = 
    new TH1F("dPhiWMuSoftMu_withCut", ";#Delta#phi(W muon, soft muon);", 400, -4.0, 4.0);
  dRWMuSoftGenMatchedMu_ = 
    new TH1F("dRWMuSoftGenMatchedMu", ";#DeltaR(W muon, soft muon);", 30, 0.0, 3.0);
  dRWMuSoftMuMuHadMassGe2_ = 
    new TH1F("dRWMuSoftMuMuHadMassGe2", ";#DeltaR(W muon, soft muon);", 30, 0.0, 3.0);
  tauMuPT_ = new TH1F("tauMuPT", ";p_{T} (GeV);", 20, 0.0, 100.0);
//   tauHadPT_ = new TH1F("tauHadPT", ";p_{T} (GeV);", 40, 0.0, 200.0);
  tauHadPT_ = new TH1F("tauHadPT", ";p_{T} (GeV);", tauHadPTBins_.size() - 1, &tauHadPTBins_[0]);
  tauHadIso_ = new TH1F("tauHadIso", ";Isolation energy (GeV);", 20, 0.0, 20.0);
  WMuLeptonRelIso_ = new TH1F("WMuLeptonRelIso", ";Relative isolation;", 100, 0.0, 2.0);
  tauHadEta_ = new TH1F("tauHadEta", ";#eta;", 46, -2.3, 2.3);
  softMuPTOverMuHadMass_ = 
    new TH1F("softMuPTOverMuHadMass", ";#frac{p_{T}^{#mu}}{m_{#mu+had}};", 80, 0.0, 80.0);
  muHadPTOverMuHadMass_ = 
    new TH1F("muHadPTOverMuHadMass", ";#frac{p_{T}}{m};", 80, 0.0, 80.0);
  dRSoftMuNearestGenMuHist_ = 
    new TH1F("dRSoftMuNearestGenMuHist", ";#DeltaR(soft #mu, gen #mu);", 40, 0.0, 2.0);
  muHadPT_ = new TH1F("muHadPT", ";p_{T}^{#mu+had} (GeV);", 50, 0.0, 100.0);
  muHadMultiplicity_ = new TH1F("muHadMultiplicity", ";N_{#mu+had}/event;", 4, 0.5, 4.5);
  nGoodVtx_ = new TH1F("nGoodVtx", ";No. good vertices;", 60, -0.5, 59.5);
  mWMuTauMu_ = new TH1F("mWMuTauMu", ";m_{#mu#mu} (GeV);", 36, 0.0, 180.0);
  PDGIDNearestStatus1GenParticleToSoftMu_ = 
    new TH1F("PDGIDNearestStatus1GenParticleToSoftMu", ";PDG ID;", 6, -0.5, 5.5);
  PDGIDMuSrc_ = new TH1F("PDGIDMuSrc", ";PDG ID;", 2, -0.5, 1.5);
  mSecondJetMuHad_ = new TH1F("mSecondJetMuHad", ";m_{#mu#tauj} (GeV);", 30, 0.0, 150.0);
  dPhiMuHadSecondJet_ = new TH1F("dPhiMuHadSecondJet", ";#Delta#phi(#mu+had,j);", 63, 0.0, 3.15);
  muHadUncleanedJetPTRank_ = new TH1F("muHadUncleanedJetPTRank", "", 2, -0.5, 1.5);
  nAddlHardMuons_ = new TH1F("nAddlHardMuons", ";N_{#mu};", 6, -0.5, 5.5);
  nAddlJetsPTGeq20_ = new TH1F("nAddlJetsPTGeq0", ";N_{j}(p_{T} > 20 GeV);", 76, -0.5, 75.5);
  nAddlJetsPTGeq30_ = new TH1F("nAddlJetsPTGeq20", ";N_{j}(p_{T} > 30 GeV);", 16, -0.5, 15.5);
  nAddlJetsPTGeq40_ = new TH1F("nAddlJetsPTGeq40", ";N_{j}(p_{T} > 40 GeV);", 16, -0.5, 15.5);
  tauHadDecayMode_ = new TH1F("tauHadDecayMode", ";Decay mode;", 5, -0.5, 4.5);
  dRSoftMuNearestGenZOrTTMu_ = 
    new TH1F("dRSoftMuNearestGenZOrTTMu", ";#DeltaR(soft #mu, nearest gen Z #mu);", 40, 0.0, 2.0);
  dRSoftMuNearestGenZOrTTMuFSR_ = new TH1F("dRSoftMuNearestGenZOrTTMuFSR", 
					   ";#DeltaR(soft #mu, nearest gen Z #mu from FSR);", 
					   40, 0.0, 2.0);
  genMuInfo_ = new TH1F("genMuInfo", ";;", 7, -0.5, 6.5);
  WMu3rdTightMuChargeProduct_ = 
    new TH1F("WMu3rdTightMuChargeProduct", ";q_{W#mu}#timesq_{#mu3};", 3, -1.5, 1.5);
  tauHadPhotonEnergyFraction_ = 
    new TH1F("tauHadPhotonEnergyFraction", ";#frac{E_{T#gamma}}{E_{T#tau}};", 20, 0.0, 1.0);
  dThetaPhotonOtherTauConstituents_ = 
    new TH1F("dThetaPhotonOtherTauConstituents", ";#Delta#theta;", 50, 0.0, 0.5);
  cleanedJetPTVsCleanedTauPT_ = 
    new TH2F("cleanedJetPTVsCleanedTauPT", ";#tau p_{T} (GeV);Jet p_{T} (GeV)", 
	     50, 0.0, 100.0, 50, 0.0, 100.0);
  uncleanedJetPTVsCleanedTauPT_ = 
    new TH2F("uncleanedJetPTVsCleanedTauPT", ";#tau p_{T} (GeV);Jet p_{T} (GeV)", 
	     50, 0.0, 100.0, 50, 0.0, 100.0);
  muHadMassVsDRSoftMuTau_ = 
    new TH2F("muHadMassVsDRSoftMuTau", ";#DeltaR(soft muon, #tau_{had});m_{#mu+had} (GeV)", 
	     30, 0.0, 3.0, 20, 0.0, 20.0);
  muHadMassVsSoftMuPT_ = 
    new TH2F("muHadMassVsSoftMuPT", ";p_{T} (GeV);m_{#mu+had} (GeV)", 
	     20, 0.0, 100.0, 20, 0.0, 20.0);
  tauHadIsoVsSoftMuPT_ = 
    new TH2F("tauHadIsoVsSoftMuPT", ";p_{T} (GeV);Isolation energy (GeV)", 
	     20, 0.0, 100.0, 20, 0.0, 20.0);
  genMuExistsVsSoftMuNearestMuProperties_ = 
    new TH2F("genMuExistsVsSoftMuNearestMuProperties", ";Nearest muon property;Gen muon property", 
	     3, -0.5, 2.5, 2, -0.5, 1.5);
  muHadMassVsTauHadEta_ = 
    new TH2F("muHadMassVsTauHadEta", ";#eta;m_{#mu+had} (GeV)", 46, -2.3, 2.3, 20, 0.0, 20.0);
  muHadMassVsSoftMuEta_ = 
    new TH2F("muHadMassVsSoftMuEta", ";#eta;m_{#mu+had} (GeV)", 46, -2.3, 2.3, 20, 0.0, 20.0);
  muHadMassVsTauHadIso_ = 
    new TH2F("muHadMassVsTauHadIso", ";Isolation energy (GeV);m_{#mu+had} (GeV)", 
	     20, 0.0, 20.0, 20, 0.0, 20.0);
  muHadMassVsTauHadPT_ = 
    new TH2F("muHadMassVsTauHadPT", ";p_{T} (GeV);m_{#mu+had} (GeV)", 
	     20, 0.0, 100.0, 20, 0.0, 20.0);
  tauHadIsoVsEta_ = 
    new TH2F("tauHadIsoVsEta", ";#eta;Isolation energy (GeV)", 46, -2.3, 2.3, 20, 0.0, 20.0);
  tauHadEtaVsSoftMuEta_ = 
    new TH2F("tauHadEtaVsSoftMuEta", ";#eta_{#mu};#eta_{#tau}", 46, -2.3, 2.3, 46, -2.3, 2.3);
  dEtaTauHadSoftMuVsDPhiTauHadSoftMu_ = 
    new TH2F("dEtaTauHadSoftMuVsDPhiTauHadSoftMu", ";#Delta#phi(#tau, #mu);#Delta#eta(#tau, #mu)", 
	     20, 0.0, 1.0, 20, 0.0, 1.0);
  tauHadPTOverMuHadMassVsTauHadIso_ = 
    new TH2F("tauHadPTOverMuHadMassVsTauHadIso", 
	     ";Isolation energy (GeV);#frac{p_{T}^{#tau}}{m_{#mu+had}}", 
	     20, 0.0, 20.0, 80, 0.0, 80.0);
  softMuPTOverMuHadMassVsTauHadIso_ = 
    new TH2F("softMuPTOverMuHadMassVsTauHadIso", 
	     ";Isolation energy (GeV);#frac{p_{T}^{#mu}}{m_{#mu+had}}", 
	     20, 0.0, 20.0, 80, 0.0, 80.0);
  avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_ = 
    new TH2F("avgTauHadSoftMuPTOverMuHadMassVsTauHadIso", 
	     ";Isolation energy (GeV);#frac{p_{T}^{#tau} + p_{T}^{#mu}}{2m_{#mu+had}}", 
	     20, 0.0, 20.0, 80, 0.0, 80.0);
  muHadPTOverMuHadMassVsTauHadIso_ = 
    new TH2F("muHadPTOverMuHadMassVsTauHadIso", 
	     ";Isolation energy (GeV);#frac{p_{T}^{#mu+had}}{m_{#mu+had}}", 
	     20, 0.0, 20.0, 80, 0.0, 80.0);
  WMuIsoVsTauHadIso_ = 
    new TH2F("WMuIsoVsTauHadIso", 
	     ";#tau_{h} Isolation energy (GeV);W_{#mu} PFRelIsolation", 
	     20, 0.0, 20.0, 400, 0.0, 40.0);
  softMuPTVsTauHadPT_ = 
    new TH2F("softMuPTVsTauHadPT", ";p_{T}^{#tau} (GeV);p_{T}^{#mu} (GeV)", 
	     20, 0.0, 100.0, 20, 0.0, 100.0);
  muHadPTOverMuHadMassVsMWMuSoftMu_ = 
    new TH2F("muHadPTOverMuHadMassVsMWMuSoftMu", ";m_{#mu#mu} (GeV);#frac{p_{T}}{m}", 
	     36, 0.0, 180.0, 80, 0.0, 80.0);
  softMuPFPTVsRECOPT_ = new TH2F("softMuPFPTVsRECOPT", ";RECO p_{T} (GeV);PF p_{T} (GeV)", 
				 20, 0.0, 100.0, 20, 0.0, 100.0);
  softMuPFEtaVsRECOEta_ = new TH2F("softMuPFEtaVsRECOEta", ";RECO #eta;PF #eta", 
				   46, -2.3, 2.3, 46, -2.3, 2.3);
  mSecondJetMuHadVsMuHadMass_ = 
    new TH2F("mSecondJetMuHadVsMuHadMass", ";m_{#mu+had} (GeV);m_{#mu#tauj} (GeV)", 
	     20, 0.0, 20.0, 60, 0.0, 300.0);
  muHadMassVsSecondJetMass_ = 
    new TH2F("muHadMassVsSecondJetMass", ";m_{j} (GeV);m_{#mu+had} (GeV)", 
	     20, 0.0, 200.0, 20, 0.0, 20.0);
  muHadPTOverMVsSecondJetPTOverM_ = 
    new TH2F("muHadPTOverMVsSecondJetPTOverM", ";p_{T}^{j}/m^{j};p_{T}^{#mu+had}/m^{#mu+had}", 
	     20, 0.0, 20.0, 20, 0.0, 20.0);
  muHadMassVsSecondJetPTOverM_ = 
    new TH2F("muHadMassVsSecondJetPTOverM", ";p_{T}^{j}/m^{j};m_{#mu+had} (GeV)", 
	     40, 0.0, 40.0, 20, 0.0, 20.0);
  tauMuTauHadJetWMuHTVsMET_ = 
    new TH2F("tauMuTauHadJetWMuHTVsMET", ";#slash{E}_{T} (GeV);H_{T} (GeV)", 
	     40, 0.0, 200.0, 40, 0.0, 1000.0);
  mWMuTauMuVsMWMuTauMuTauHad_ = 
    new TH2F("mWMuTauMuVsMWMuTauMuTauHad", ";m_{#mu#mu#tau} (GeV);m_{#mu#mu} (GeV)", 
	     36, 0.0, 180.0, 36, 0.0, 180.0);
  mWMuTauMuVsMWMuTauMuTauHadGenFSR_ = 
    new TH2F("mWMuTauMuVsMWMuTauMuTauHadGenFSR", ";m_{#mu#mu#tau} (GeV);m_{#mu#mu} (GeV)", 
	     36, 0.0, 180.0, 36, 0.0, 180.0);
  muHadMassVsMWMuTauMuTauHad_ = 
    new TH2F("muHadMassVsMWMuTauMuTauHad", ";m_{#mu#mu#tau} (GeV);m_{#mu+had} (GeV)", 
	     36, 0.0, 180.0, 20, 0.0, 20.0);
  jet_pt_etacut = new TH1F("jet_pt_etacut", ";p_{T} (GeV);", 20, 0., 200.);
  muHad_t3t1Vsptmj_ = new TH2F("muHad_t3t1Vsptmj", ";#frac{p_{T}}{m};#frac{#tau_{3}}{#tau_{1}}", 
			       80, 0.0, 80.0, 500, 0.0, 2.0);
  muHad_t3t1VsDecayMode_ = new TH2F("muHad_t3t1VsDecayMode", 
				    ";HPS #tau decay mode;#frac{#tau_{3}}{#tau_{1}}", 
				    16, -1.5, 14.5, 500, 0.0, 2.0);
  muHadMassVsNAddlJets_ = 
    new TH2F("muHadMassVsNAddlJets", ";N_{j};m_{#mu+had} (GeV)", 10, -0.5, 9.5, 20, 0.0, 20.0);
  jet_eta = new TH1F("jet_eta", "#eta", 70, -3.5, 3.5);
  jet_phi = new TH1F("jet_phi", "#phi", 14, -3.5, 3.5);
  jet_mass_etacut = new TH1F("jet_mass_etacut", "m (GeV)", 20, 0., 200.);
  jet_ptmj_etacut = new TH1F("jet_ptmj_etacut", "#frac{p_{T}}{m}", 80, 0., 80.);
  muHad_t3t1_ =
    new TH1F("muHad_t3t1", "{#tau}_{3}/{#tau}_{1} of cleaned parent jet", 50, 0.0, 2.0);
  muHad_t2t1_ =
    new TH1F("muHad_t2t1", "{#tau}_{2}/{#tau}_{1} of cleaned parent jet", 50, 0.0, 2.0);
  muHad_t3t1_pT1020_ = 
    new TH1F("muHad_t3t1_pT1020", 
	     "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with pT between 10 and 20", 
	     50, 0.0, 2.0);
  muHad_t3t1_pT2030_ =
    new TH1F("muHad_t3t1_pT2030", 
	     "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with pT between 20 and 30", 
	     50, 0.0, 2.0);
  muHad_t3t1_pT3040_ =
    new TH1F("muHad_t3t1_pT3040", 
	     "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with pT between 30 and 40", 
	     50, 0.0, 2.0);
  muHad_t3t1_pT4050_ =
    new TH1F("muHad_t3t1_pT4050", 
	     "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with pT between 40 and 50", 
	     50, 0.0, 2.0);
  muHad_t3t1_pT50Up_ =
    new TH1F("muHad_t3t1_pT50Up", "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with pT above 50", 
	     50, 0.0, 2.0);
  muHad_t3t1_0Jets_ =
    new TH1F("muHad_t3t1_0Jets", "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with NJets = 0", 
	     50, 0.0, 2.0);
  muHad_t3t1_1Jets_ =
    new TH1F("muHad_t3t1_1Jets", "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with NJets = 1", 
	     50, 0.0, 2.0);
  muHad_t3t1_2Jets_ =
    new TH1F("muHad_t3t1_2Jets", "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with NJets = 2", 
	     50, 0.0, 2.0);
  muHad_t3t1_3Jets_ =
    new TH1F("muHad_t3t1_3Jets", "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with NJets = 3", 
	     50, 0.0, 2.0);
  muHad_t3t1_4Jets_ =
    new TH1F("muHad_t3t1_4Jets", "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with NJets = 4", 
	     50, 0.0, 2.0);
  muHad_t3t1_5Jets_ =
    new TH1F("muHad_t3t1_5Jets", "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with NJets = 5", 
	     50, 0.0, 2.0);
  muHad_t3t1_MoreJets_ =
    new TH1F("muHad_t3t1_MoreJets", "{#tau}_{3}/{#tau}_{1} of cleaned parent jet with NJets > 5", 
	     50, 0.0, 2.0);
  muHad_Nchtrk_0_ = 
    new TH1F("muHad_Nchtrk_0", "charged track multiplicity of jet (pT > 0)", 40, -0.5, 39.5);
  muHad_Nchtrk_1_ = 
    new TH1F("muHad_Nchtrk_1", "charged track multiplicity of jet (pT > 1)", 40, -0.5, 39.5);
  muHad_Nchtrk_10_ = 
    new TH1F("muHad_Nchtrk_10", "charged track multiplicity of jet (pT > 10)", 40, -0.5, 39.5);
  muHad_Nchtrk_30_ = 
    new TH1F("muHad_Nchtrk_30", "charged track multiplicity of jet (pT > 30)", 40, -0.5, 39.5);
  second_Nchtrk_0_ = new TH1F("second_Nchtrk_0", 
			      "charged track multiplicity of second jet (pT > 0)", 40, -0.5, 39.5);
  second_Nchtrk_1_ = new TH1F("second_Nchtrk_1", 
			      "charged track multiplicity of second jet (pT > 1)", 40, -0.5, 39.5);
  second_Nchtrk_10_ = new TH1F("second_Nchtrk_10", 
			       "charged track multiplicity of second jet (pT > 10)", 
			       40, -0.5, 39.5);
  second_Nchtrk_30_ = new TH1F("second_Nchtrk_30", 
			       "charged track multiplicity of second jet (pT > 30)", 
			       40, -0.5, 39.5);
  dPhiWMuSecJet_ = 
    new TH1F("dPhiWMuSecJet", "#Delta#phi between W muon and second jet", 400, -4., 4.);

  //set bin labels
  std::stringstream dRStream;
  dRStream << dR_;
  genMuExistsVsSoftMuNearestMuProperties_->
    GetXaxis()->SetBinLabel(1, ("#DeltaR < " + dRStream.str() + " and W muon").c_str());
  genMuExistsVsSoftMuNearestMuProperties_->
    GetXaxis()->SetBinLabel(2, ("#DeltaR < " + dRStream.str() + " and not W muon").c_str());
  genMuExistsVsSoftMuNearestMuProperties_->
    GetXaxis()->SetBinLabel(3, ("#DeltaR #geq " + dRStream.str()).c_str());
  genMuExistsVsSoftMuNearestMuProperties_->
    GetYaxis()->SetBinLabel(1, ("#DeltaR < " + dRStream.str()).c_str());
  genMuExistsVsSoftMuNearestMuProperties_->
    GetYaxis()->SetBinLabel(2, ("#DeltaR #geq " + dRStream.str()).c_str());
  PDGIDNearestStatus1GenParticleToSoftMu_->GetXaxis()->SetBinLabel(1, "e");
  PDGIDNearestStatus1GenParticleToSoftMu_->GetXaxis()->SetBinLabel(2, "#mu");
  PDGIDNearestStatus1GenParticleToSoftMu_->GetXaxis()->SetBinLabel(3, "#nu");
  PDGIDNearestStatus1GenParticleToSoftMu_->GetXaxis()->SetBinLabel(4, "#gamma");
  PDGIDNearestStatus1GenParticleToSoftMu_->GetXaxis()->SetBinLabel(5, "h");
  PDGIDNearestStatus1GenParticleToSoftMu_->GetXaxis()->SetBinLabel(6, "Other");
  PDGIDMuSrc_->GetXaxis()->SetBinLabel(1, "b/c");
  PDGIDMuSrc_->GetXaxis()->SetBinLabel(2, "Other");
  muHadUncleanedJetPTRank_->GetXaxis()->SetBinLabel(1, "1 or 2");
  muHadUncleanedJetPTRank_->GetXaxis()->SetBinLabel(2, ">3");
  tauHadDecayMode_->GetXaxis()->SetBinLabel(1, "1-prong");
  tauHadDecayMode_->GetXaxis()->SetBinLabel(2, "1-prong + 1 #pi^{0}");
  tauHadDecayMode_->GetXaxis()->SetBinLabel(3, "1-prong + 2 #pi^{0}");
  tauHadDecayMode_->GetXaxis()->SetBinLabel(4, "3-prong");
  tauHadDecayMode_->GetXaxis()->SetBinLabel(5, "Other");
  genMuInfo_->GetXaxis()->SetBinLabel(1, "Correct");
  genMuInfo_->GetXaxis()->SetBinLabel(2, "W muon correct");
  genMuInfo_->GetXaxis()->SetBinLabel(3, "#tau muon correct");
  genMuInfo_->GetXaxis()->SetBinLabel(4, "Switched");
  genMuInfo_->GetXaxis()->SetBinLabel(5, "#tau muon #rightarrow W muon");
  genMuInfo_->GetXaxis()->SetBinLabel(6, "W muon #rightarrow #tau muon");
  genMuInfo_->GetXaxis()->SetBinLabel(7, "Neither");

  //set sumw2
  MET_->Sumw2();
  hadTauAssociatedMuMultiplicity_->Sumw2();
  muHadMass_->Sumw2();
  muHadMassReweightErrSq_->Sumw2();
  muHadCharge_->Sumw2();
  bTagDiscrim_->Sumw2();
  WMuIso_->Sumw2();
  WMuMT_->Sumw2();
  tauMuMT_->Sumw2();
  tauHadMT_->Sumw2();
  dPhiWMuMET_->Sumw2();
  dPhiTauMuMET_->Sumw2();
  tauMuTauHadJetHT_->Sumw2();
  diJetHT_->Sumw2();
  jetTauJetHT_->Sumw2();
  tauMuTauHadJetWMuHT_->Sumw2();
  tauMuTauHadJetWMuMETHT_->Sumw2();
  diJetWMuHT_->Sumw2();
  jetTauJetWMuHT_->Sumw2();
  dRWMuSoftMu_->Sumw2();
  dPhiWMuSoftMu_->Sumw2();
  dPhiWMuSoftMu_withCut_->Sumw2();
  dRWMuSoftGenMatchedMu_->Sumw2();
  dRWMuSoftMuMuHadMassGe2_->Sumw2();
  tauMuPT_->Sumw2();
  tauHadPT_->Sumw2();
  tauHadIso_->Sumw2();
  WMuLeptonRelIso_->Sumw2();
  tauHadEta_->Sumw2();
  softMuPTOverMuHadMass_->Sumw2();
  muHadPTOverMuHadMass_->Sumw2();
  dRSoftMuNearestGenMuHist_->Sumw2();
  muHadPT_->Sumw2();
  muHadMultiplicity_->Sumw2();
  nGoodVtx_->Sumw2();
  mWMuTauMu_->Sumw2();
  PDGIDNearestStatus1GenParticleToSoftMu_->Sumw2();
  PDGIDMuSrc_->Sumw2();
  mSecondJetMuHad_->Sumw2();
  dPhiMuHadSecondJet_->Sumw2();
  muHadUncleanedJetPTRank_->Sumw2();
  nAddlHardMuons_->Sumw2();
  nAddlJetsPTGeq20_->Sumw2();
  nAddlJetsPTGeq30_->Sumw2();
  nAddlJetsPTGeq40_->Sumw2();
  tauHadDecayMode_->Sumw2();
  dRSoftMuNearestGenZOrTTMu_->Sumw2();
  dRSoftMuNearestGenZOrTTMuFSR_->Sumw2();
  genMuInfo_->Sumw2();
  WMu3rdTightMuChargeProduct_->Sumw2();
  tauHadPhotonEnergyFraction_->Sumw2();
  dThetaPhotonOtherTauConstituents_->Sumw2();
  cleanedJetPTVsCleanedTauPT_->Sumw2();
  uncleanedJetPTVsCleanedTauPT_->Sumw2();
  muHadMassVsDRSoftMuTau_->Sumw2();
  muHadMassVsSoftMuPT_->Sumw2();
  tauHadIsoVsSoftMuPT_->Sumw2();
  genMuExistsVsSoftMuNearestMuProperties_->Sumw2();
  muHadMassVsTauHadEta_->Sumw2();
  muHadMassVsSoftMuEta_->Sumw2();
  muHadMassVsTauHadIso_->Sumw2();
  muHadMassVsTauHadPT_->Sumw2();
  tauHadIsoVsEta_->Sumw2();
  tauHadEtaVsSoftMuEta_->Sumw2();
  dEtaTauHadSoftMuVsDPhiTauHadSoftMu_->Sumw2();
  tauHadPTOverMuHadMassVsTauHadIso_->Sumw2();
  softMuPTOverMuHadMassVsTauHadIso_->Sumw2();
  avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_->Sumw2();
  muHadPTOverMuHadMassVsTauHadIso_->Sumw2();
  WMuIsoVsTauHadIso_->Sumw2();
  softMuPTVsTauHadPT_->Sumw2();
  muHad_t3t1Vsptmj_->Sumw2();
  muHad_t3t1VsDecayMode_->Sumw2();
  muHad_t3t1_->Sumw2();
  muHad_t2t1_->Sumw2();
  muHad_t3t1_pT1020_->Sumw2();
  muHad_t3t1_pT2030_->Sumw2();
  muHad_t3t1_pT3040_->Sumw2();
  muHad_t3t1_pT4050_->Sumw2();
  muHad_t3t1_pT50Up_->Sumw2();
  muHad_t3t1_0Jets_->Sumw2();
  muHad_t3t1_1Jets_->Sumw2();
  muHad_t3t1_2Jets_->Sumw2();
  muHad_t3t1_3Jets_->Sumw2();
  muHad_t3t1_4Jets_->Sumw2();
  muHad_t3t1_5Jets_->Sumw2();
  muHad_t3t1_MoreJets_->Sumw2();
  muHadPTOverMuHadMassVsMWMuSoftMu_->Sumw2();
  softMuPFPTVsRECOPT_->Sumw2();
  softMuPFEtaVsRECOEta_->Sumw2();
  mSecondJetMuHadVsMuHadMass_->Sumw2();
  muHadMassVsSecondJetMass_->Sumw2();
  muHadPTOverMVsSecondJetPTOverM_->Sumw2();
  muHadMassVsSecondJetPTOverM_->Sumw2();
  tauMuTauHadJetWMuHTVsMET_->Sumw2();
  mWMuTauMuVsMWMuTauMuTauHad_->Sumw2();
  mWMuTauMuVsMWMuTauMuTauHadGenFSR_->Sumw2();
  muHadMassVsMWMuTauMuTauHad_->Sumw2();
  muHadMassVsNAddlJets_->Sumw2();
  WMuMTVsMET_->Sumw2();
  jet_pt_etacut->Sumw2();
  jet_eta->Sumw2();
  jet_phi->Sumw2();
  jet_mass_etacut->Sumw2();
  jet_ptmj_etacut->Sumw2();
  muHad_Nchtrk_0_->Sumw2();
  muHad_Nchtrk_1_->Sumw2();
  muHad_Nchtrk_10_->Sumw2();
  muHad_Nchtrk_30_->Sumw2();
  second_Nchtrk_0_->Sumw2();
  second_Nchtrk_1_->Sumw2();
  second_Nchtrk_10_->Sumw2();
  second_Nchtrk_30_->Sumw2();
  dPhiWMuSecJet_->Sumw2();

  //maximum hadronic tau PT
  maxTauHadPT_ = 0.0;

  //maximum mu+had mass
  maxMuHadMass_ = 0.0;
}

// ------------ method called once each job just after ending the event loop  ------------
void TauAnalyzer::endJob() 
{
  //make canvases
  out_->cd();
  TCanvas METCanvas("METCanvas", "", 600, 600);
  TCanvas hadTauAssociatedMuMultiplicityCanvas("hadTauAssociatedMuMultiplicityCanvas", "", 
					       600, 600);
  TCanvas muHadMassCanvas("muHadMassCanvas", "", 600, 600);
  TCanvas muHadMassReweightErrSqCanvas("muHadMassReweightErrSqCanvas", "", 600, 600);
  TCanvas muHadChargeCanvas("muHadChargeCanvas", "", 600, 600);
  TCanvas bTagDiscrimCanvas("bTagDiscrimCanvas", "", 600, 600);
  TCanvas WMuIsoCanvas("WMuIsoCanvas", "", 600, 600);
  TCanvas WMuMTCanvas("WMuMTCanvas", "", 600, 600);
  TCanvas tauMuMTCanvas("tauMuMTCanvas", "", 600, 600);
  TCanvas tauHadMTCanvas("tauHadMTCanvas", "", 600, 600);
  TCanvas dPhiWMuMETCanvas("dPhiWMuMETCanvas", "", 600, 600);
  TCanvas dPhiTauMuMETCanvas("dPhiTauMuMETCanvas", "", 600, 600);
  TCanvas tauMuTauHadJetHTCanvas("tauMuTauHadJetHTCanvas", "", 600, 600);
  TCanvas diJetHTCanvas("diJetHTCanvas", "", 600, 600);
  TCanvas jetTauJetHTCanvas("jetTauJetHTCanvas", "", 600, 600);
  TCanvas tauMuTauHadJetWMuHTCanvas("tauMuTauHadJetWMuHTCanvas", "", 600, 600);
  TCanvas tauMuTauHadJetWMuMETHTCanvas("tauMuTauHadJetWMuMETHTCanvas", "", 600, 600);
  TCanvas diJetWMuHTCanvas("diJetWMuHTCanvas", "", 600, 600);
  TCanvas jetTauJetWMuHTCanvas("jetTauJetWMuHTCanvas", "", 600, 600);
  TCanvas dRWMuSoftMuCanvas("dRWMuSoftMuCanvas", "", 600, 600);
  TCanvas dPhiWMuSoftMuCanvas("dPhiWMuSoftMuCanvas", "", 600, 600);
  TCanvas dPhiWMuSoftMuWithCutCanvas("dPhiWMuSoftMuWithCutCanvas", "", 600, 600);
  TCanvas dRWMuSoftGenMatchedMuCanvas("dRWMuSoftGenMatchedMuCanvas", "", 600, 600);
  TCanvas dRWMuSoftMuMuHadMassGe2Canvas("dRWMuSoftMuMuHadMassGe2Canvas", "", 600, 600);
  TCanvas tauMuPTCanvas("tauMuPTCanvas", "", 600, 600);
  TCanvas tauHadPTCanvas("tauHadPTCanvas", "", 600, 600);
  TCanvas tauHadIsoCanvas("tauHadIsoCanvas", "", 600, 600);
  TCanvas WMuLeptonRelIsoCanvas("WMuLeptonRelIsoCanvas", "", 600, 600);
  TCanvas tauHadEtaCanvas("tauHadEtaCanvas", "", 600, 600);
  TCanvas softMuPTOverMuHadMassCanvas("softMuPTOverMuHadMassCanvas", "", 600, 600);
  TCanvas muHadPTOverMuHadMassCanvas("muHadPTOverMuHadMassCanvas", "", 600, 600);
  TCanvas dRSoftMuNearestGenMuHistCanvas("dRSoftMuNearestGenMuHistCanvas", "", 600, 600);
  TCanvas muHadPTCanvas("muHadPTCanvas", "", 600, 600);
  TCanvas muHadMultiplicityCanvas("muHadMultiplicityCanvas", "", 600, 600);
  TCanvas nGoodVtxCanvas("nGoodVtxCanvas", "", 600, 600);
  TCanvas mWMuTauMuCanvas("mWMuTauMuCanvas", "", 600, 600);
  TCanvas 
    PDGIDNearestStatus1GenParticleToSoftMuCanvas("PDGIDNearestStatus1GenParticleToSoftMuCanvas", 
						 "", 600, 600);
  TCanvas PDGIDMuSrcCanvas("PDGIDMuSrcCanvas", "", 600, 600);
  TCanvas mSecondJetMuHadCanvas("mSecondJetMuHadCanvas", "", 600, 600);
  TCanvas dPhiMuHadSecondJetCanvas("dPhiMuHadSecondJetCanvas", "", 600, 600);
  TCanvas muHadUncleanedJetPTRankCanvas("muHadUncleanedJetPTRankCanvas", "", 600, 600);
  TCanvas nAddlHardMuonsCanvas("nAddlHardMuonsCanvas", "", 600, 600);
  TCanvas nAddlJetsPTGeq20Canvas("nAddlJetsPTGeq0Canvas", "", 600, 600);
  TCanvas nAddlJetsPTGeq30Canvas("nAddlJetsPTGeq20Canvas", "", 600, 600);
  TCanvas nAddlJetsPTGeq40Canvas("nAddlJetsPTGeq40Canvas", "", 600, 600);
  TCanvas tauHadDecayModeCanvas("tauHadDecayModeCanvas", "", 600, 600);
  TCanvas dRSoftMuNearestGenZOrTTMuCanvas("dRSoftMuNearestGenZOrTTMuCanvas", "", 600, 600);
  TCanvas dRSoftMuNearestGenZOrTTMuFSRCanvas("dRSoftMuNearestGenZOrTTMuFSRCanvas", "", 600, 600);
  TCanvas genMuInfoCanvas("genMuInfoCanvas", "", 600, 600);
  TCanvas WMu3rdTightMuChargeProductCanvas("WMu3rdTightMuChargeProductCanvas", "", 600, 600);
  TCanvas tauHadPhotonEnergyFractionCanvas("tauHadPhotonEnergyFractionCanvas", "", 600, 600);
  TCanvas 
    dThetaPhotonOtherTauConstituentsCanvas("dThetaPhotonOtherTauConstituentsCanvas", "", 600, 600);
  TCanvas cleanedJetPTVsCleanedTauPTCanvas("cleanedJetPTVsCleanedTauPTCanvas", "", 600, 600);
  TCanvas uncleanedJetPTVsCleanedTauPTCanvas("uncleanedJetPTVsCleanedTauPTCanvas", "", 600, 600);
  TCanvas muHadMassVsDRSoftMuTauCanvas("muHadMassVsDRSoftMuTauCanvas", "", 600, 600);
  TCanvas muHadMassVsSoftMuPTCanvas("muHadMassVsSoftMuPTCanvas", "", 600, 600);
  TCanvas tauHadIsoVsSoftMuPTCanvas("tauHadIsoVsSoftMuPTCanvas", "", 600, 600);
  TCanvas 
    genMuExistsVsSoftMuNearestMuPropertiesCanvas("genMuExistsVsSoftMuNearestMuPropertiesCanvas", 
						 "", 600, 600);
  TCanvas muHadMassVsTauHadEtaCanvas("muHadMassVsTauHadEtaCanvas", "", 600, 600);
  TCanvas muHadMassVsSoftMuEtaCanvas("muHadMassVsSoftMuEtaCanvas", "", 600, 600);
  TCanvas muHadMassVsTauHadIsoCanvas("muHadMassVsTauHadIsoCanvas", "", 600, 600);
  TCanvas muHadMassVsTauHadPTCanvas("muHadMassVsTauHadPTCanvas", "", 600, 600);
  TCanvas tauHadIsoVsEtaCanvas("tauHadIsoVsEtaCanvas", "", 600, 600);
  TCanvas tauHadEtaVsSoftMuEtaCanvas("tauHadEtaVsSoftMuEtaCanvas", "", 600, 600);
  TCanvas 
    dEtaTauHadSoftMuVsDPhiTauHadSoftMuCanvas("dEtaTauHadSoftMuVsDPhiTauHadSoftMuCanvas", "", 
					     600, 600);
  TCanvas 
    tauHadPTOverMuHadMassVsTauHadIsoCanvas("tauHadPTOverMuHadMassVsTauHadIsoCanvas", "", 600, 600);
  TCanvas 
    softMuPTOverMuHadMassVsTauHadIsoCanvas("softMuPTOverMuHadMassVsTauHadIsoCanvas", "", 600, 600);
  TCanvas 
    avgTauHadSoftMuPTOverMuHadMassVsTauHadIsoCanvas("avgTauHadSoftMuPTOverMuHadMassVsTauHadIsoCanvas", "", 600, 600);
  TCanvas 
    muHadPTOverMuHadMassVsTauHadIsoCanvas("muHadPTOverMuHadMassVsTauHadIsoCanvas", "", 600, 600);
  TCanvas WMuIsoVsTauHadIsoCanvas("WMuIsoVsTauHadIsoCanvas", "", 600, 600);
  TCanvas softMuPTVsTauHadPTCanvas("softMuPTVsTauHadPTCanvas", "", 600, 600);
  TCanvas muHad_t3t1Canvas("muHad_t3t1Canvas", "", 600, 600);
  TCanvas muHad_t2t1Canvas("muHad_t2t1Canvas", "", 600, 600);
  TCanvas muHad_t3t1_pT1020Canvas("muHad_t3t1_pT1020Canvas", "", 600, 600);
  TCanvas muHad_t3t1_pT2030Canvas("muHad_t3t1_pT2030Canvas", "", 600, 600);
  TCanvas muHad_t3t1_pT3040Canvas("muHad_t3t1_pT3040Canvas", "", 600, 600);
  TCanvas muHad_t3t1_pT4050Canvas("muHad_t3t1_pT4050Canvas", "", 600, 600);
  TCanvas muHad_t3t1_pT50UpCanvas("muHad_t3t1_pT50UpCanvas", "", 600, 600);
  TCanvas muHad_t3t1_0JetsCanvas("muHad_t3t1_0JetsCanvas", "", 600, 600);
  TCanvas muHad_t3t1_1JetsCanvas("muHad_t3t1_1JetsCanvas", "", 600, 600);
  TCanvas muHad_t3t1_2JetsCanvas("muHad_t3t1_2JetsCanvas", "", 600, 600);
  TCanvas muHad_t3t1_3JetsCanvas("muHad_t3t1_3JetsCanvas", "", 600, 600);
  TCanvas muHad_t3t1_4JetsCanvas("muHad_t3t1_4JetsCanvas", "", 600, 600);
  TCanvas muHad_t3t1_5JetsCanvas("muHad_t3t1_5JetsCanvas", "", 600, 600);
  TCanvas muHad_t3t1_MoreJetsCanvas("muHad_t3t1_MoreJetsCanvas", "", 600, 600);
  TCanvas muHad_Nchtrk_0_Canvas("muHadNchtrk_0_Canvas", "", 600, 600);
  TCanvas muHad_Nchtrk_1_Canvas("muHadNchtrk_1_Canvas", "", 600, 600);
  TCanvas muHad_Nchtrk_10_Canvas("muHadNchtrk_10_Canvas", "", 600, 600);
  TCanvas muHad_Nchtrk_30_Canvas("muHadNchtrk_30_Canvas", "", 600, 600);
  TCanvas second_Nchtrk_0_Canvas("secondNchtrk_0_Canvas", "", 600, 600);
  TCanvas second_Nchtrk_1_Canvas("secondNchtrk_1_Canvas", "", 600, 600);
  TCanvas second_Nchtrk_10_Canvas("secondNchtrk_10_Canvas", "", 600, 600);
  TCanvas second_Nchtrk_30_Canvas("secondNchtrk_30_Canvas", "", 600, 600);
  TCanvas 
    muHadPTOverMuHadMassVsMWMuSoftMuCanvas("muHadPTOverMuHadMassVsMWMuSoftMuCanvas", "", 600, 600);
  TCanvas softMuPFPTVsRECOPTCanvas("softMuPFPTVsRECOPTCanvas", "", 600, 600);
  TCanvas softMuPFEtaVsRECOEtaCanvas("softMuPFEtaVsRECOEtaCanvas", "", 600, 600);
  TCanvas mSecondJetMuHadVsMuHadMassCanvas("mSecondJetMuHadVsMuHadMassCanvas", "", 600, 600);
  TCanvas muHadMassVsSecondJetMassCanvas("muHadMassVsSecondJetMassCanvas", "", 600, 600);
  TCanvas 
    muHadPTOverMVsSecondJetPTOverMCanvas("muHadPTOverMVsSecondJetPTOverMCanvas", "", 600, 600);
  TCanvas muHadMassVsSecondJetPTOverMCanvas("muHadMassVsSecondJetPTOverMCanvas", "", 600, 600);
  TCanvas tauMuTauHadJetWMuHTVsMETCanvas("tauMuTauHadJetWMuHTVsMETCanvas", "", 600, 600);
  TCanvas mWMuTauMuVsMWMuTauMuTauHadCanvas("mWMuTauMuVsMWMuTauMuTauHadCanvas", "", 600, 600);
  TCanvas 
    mWMuTauMuVsMWMuTauMuTauHadGenFSRCanvas("mWMuTauMuVsMWMuTauMuTauHadGenFSRCanvas", "", 600, 600);
  TCanvas muHadMassVsMWMuTauMuTauHadCanvas("muHadMassVsMWMuTauMuTauHadCanvas", "", 600, 600);
  TCanvas muHadMassVsNAddlJetsCanvas("muHadMassVsNAddlJetsCanvas", "", 600, 600);
  TCanvas WMuMTVsMETCanvas("WMuMTVsMETCanvas", "", 600, 600);
  TCanvas jet_pt_etacutCanvas("jet_pt_etacutCanvas", "", 600, 600);
  TCanvas jet_etaCanvas("jet_etaCanvas", "", 600, 600);
  TCanvas jet_phiCanvas("jet_phiCanvas", "", 600, 600);
  TCanvas jet_mass_etacutCanvas("jet_mass_etacutCanvas", "", 600, 600);
  TCanvas jet_ptmj_etacutCanvas("jet_ptmj_etacutCanvas", "", 600, 600);
  TCanvas dPhiWMuSecJetCanvas("dPhiWMuSecJetCanvas", "", 600, 600);
  TCanvas muHad_t3t1VsptmjCanvas("muHad_t3t1VsptmjCanvas", "", 600, 600);
  TCanvas muHad_t3t1VsDecayModeCanvas("muHad_t3t1VsDecayModeCanvas", "", 600, 600);

  //format and draw 1D plots
  Common::draw1DHistograms(METCanvas, MET_);
  Common::draw1DHistograms(hadTauAssociatedMuMultiplicityCanvas, hadTauAssociatedMuMultiplicity_);
  Common::draw1DHistograms(muHadMassCanvas, muHadMass_);
  Common::draw1DHistograms(muHadMassReweightErrSqCanvas, muHadMassReweightErrSq_);
  Common::draw1DHistograms(muHadChargeCanvas, muHadCharge_);
  Common::draw1DHistograms(bTagDiscrimCanvas, bTagDiscrim_);
  Common::draw1DHistograms(WMuIsoCanvas, WMuIso_);
  Common::draw1DHistograms(WMuMTCanvas, WMuMT_);
  Common::draw1DHistograms(tauMuMTCanvas, tauMuMT_);
  Common::draw1DHistograms(tauHadMTCanvas, tauHadMT_);
  Common::draw1DHistograms(dPhiWMuMETCanvas, dPhiWMuMET_);
  Common::draw1DHistograms(dPhiTauMuMETCanvas, dPhiTauMuMET_);
  Common::draw1DHistograms(tauMuTauHadJetHTCanvas, tauMuTauHadJetHT_);
  Common::draw1DHistograms(diJetHTCanvas, diJetHT_);
  Common::draw1DHistograms(jetTauJetHTCanvas, jetTauJetHT_);
  Common::draw1DHistograms(tauMuTauHadJetWMuHTCanvas, tauMuTauHadJetWMuHT_);
  Common::draw1DHistograms(tauMuTauHadJetWMuMETHTCanvas, tauMuTauHadJetWMuMETHT_);
  Common::draw1DHistograms(diJetWMuHTCanvas, diJetWMuHT_);
  Common::draw1DHistograms(jetTauJetWMuHTCanvas, jetTauJetWMuHT_);
  Common::draw1DHistograms(dRWMuSoftMuCanvas, dRWMuSoftMu_);
  Common::draw1DHistograms(dPhiWMuSoftMuCanvas, dPhiWMuSoftMu_);
  Common::draw1DHistograms(dPhiWMuSoftMuWithCutCanvas, dPhiWMuSoftMu_withCut_);
  Common::draw1DHistograms(dRWMuSoftGenMatchedMuCanvas, dRWMuSoftGenMatchedMu_);
  Common::draw1DHistograms(dRWMuSoftMuMuHadMassGe2Canvas, dRWMuSoftMuMuHadMassGe2_);
  Common::draw1DHistograms(tauMuPTCanvas, tauMuPT_);
  Common::draw1DHistograms(tauHadPTCanvas, tauHadPT_);
  Common::draw1DHistograms(tauHadIsoCanvas, tauHadIso_);
  Common::draw1DHistograms(WMuLeptonRelIsoCanvas, WMuLeptonRelIso_);
  Common::draw1DHistograms(tauHadEtaCanvas, tauHadEta_);
  Common::draw1DHistograms(softMuPTOverMuHadMassCanvas, softMuPTOverMuHadMass_);
  Common::draw1DHistograms(muHadPTOverMuHadMassCanvas, muHadPTOverMuHadMass_);
  Common::draw1DHistograms(dRSoftMuNearestGenMuHistCanvas, dRSoftMuNearestGenMuHist_);
  Common::draw1DHistograms(muHadPTCanvas, muHadPT_);
  Common::draw1DHistograms(muHadMultiplicityCanvas, muHadMultiplicity_);
  Common::draw1DHistograms(nGoodVtxCanvas, nGoodVtx_);
  Common::draw1DHistograms(mWMuTauMuCanvas, mWMuTauMu_);
  Common::draw1DHistograms(muHad_t3t1Canvas, muHad_t3t1_);
  Common::draw1DHistograms(muHad_t2t1Canvas, muHad_t2t1_);
  Common::draw1DHistograms(muHad_t3t1_pT1020Canvas, muHad_t3t1_pT1020_);
  Common::draw1DHistograms(muHad_t3t1_pT2030Canvas, muHad_t3t1_pT2030_);
  Common::draw1DHistograms(muHad_t3t1_pT3040Canvas, muHad_t3t1_pT3040_);
  Common::draw1DHistograms(muHad_t3t1_pT4050Canvas, muHad_t3t1_pT4050_);
  Common::draw1DHistograms(muHad_t3t1_pT50UpCanvas, muHad_t3t1_pT50Up_);
  Common::draw1DHistograms(muHad_t3t1_0JetsCanvas, muHad_t3t1_0Jets_);
  Common::draw1DHistograms(muHad_t3t1_1JetsCanvas, muHad_t3t1_1Jets_);
  Common::draw1DHistograms(muHad_t3t1_2JetsCanvas, muHad_t3t1_2Jets_);
  Common::draw1DHistograms(muHad_t3t1_3JetsCanvas, muHad_t3t1_3Jets_);
  Common::draw1DHistograms(muHad_t3t1_4JetsCanvas, muHad_t3t1_4Jets_);
  Common::draw1DHistograms(muHad_t3t1_5JetsCanvas, muHad_t3t1_5Jets_);
  Common::draw1DHistograms(muHad_t3t1_MoreJetsCanvas, muHad_t3t1_MoreJets_);
  Common::draw1DHistograms(muHad_Nchtrk_0_Canvas, muHad_Nchtrk_0_);
  Common::draw1DHistograms(muHad_Nchtrk_1_Canvas, muHad_Nchtrk_1_);
  Common::draw1DHistograms(muHad_Nchtrk_10_Canvas, muHad_Nchtrk_10_);
  Common::draw1DHistograms(muHad_Nchtrk_30_Canvas, muHad_Nchtrk_30_);
  Common::draw1DHistograms(second_Nchtrk_0_Canvas, second_Nchtrk_0_);
  Common::draw1DHistograms(second_Nchtrk_1_Canvas, second_Nchtrk_1_);
  Common::draw1DHistograms(second_Nchtrk_10_Canvas, second_Nchtrk_10_);
  Common::draw1DHistograms(second_Nchtrk_30_Canvas, second_Nchtrk_30_);
  Common::draw1DHistograms(PDGIDNearestStatus1GenParticleToSoftMuCanvas, 
			   PDGIDNearestStatus1GenParticleToSoftMu_);
  Common::draw1DHistograms(PDGIDMuSrcCanvas, PDGIDMuSrc_);
  Common::draw1DHistograms(jet_pt_etacutCanvas, jet_pt_etacut);
  Common::draw1DHistograms(jet_etaCanvas, jet_eta);
  Common::draw1DHistograms(jet_phiCanvas, jet_phi);
  Common::draw1DHistograms(jet_mass_etacutCanvas, jet_mass_etacut);
  Common::draw1DHistograms(jet_ptmj_etacutCanvas, jet_ptmj_etacut);
  Common::draw1DHistograms(mSecondJetMuHadCanvas, mSecondJetMuHad_);
  Common::draw1DHistograms(dPhiMuHadSecondJetCanvas, dPhiMuHadSecondJet_);
  Common::draw1DHistograms(muHadUncleanedJetPTRankCanvas, muHadUncleanedJetPTRank_);
  Common::draw1DHistograms(nAddlHardMuonsCanvas, nAddlHardMuons_);
  Common::draw1DHistograms(nAddlJetsPTGeq20Canvas, nAddlJetsPTGeq20_);
  Common::draw1DHistograms(nAddlJetsPTGeq30Canvas, nAddlJetsPTGeq30_);
  Common::draw1DHistograms(nAddlJetsPTGeq40Canvas, nAddlJetsPTGeq40_);
  Common::draw1DHistograms(tauHadDecayModeCanvas, tauHadDecayMode_);
  Common::draw1DHistograms(dRSoftMuNearestGenZOrTTMuCanvas, dRSoftMuNearestGenZOrTTMu_);
  Common::draw1DHistograms(dRSoftMuNearestGenZOrTTMuFSRCanvas, dRSoftMuNearestGenZOrTTMuFSR_);
  Common::draw1DHistograms(genMuInfoCanvas, genMuInfo_);
  Common::draw1DHistograms(WMu3rdTightMuChargeProductCanvas, WMu3rdTightMuChargeProduct_);
  Common::draw1DHistograms(dPhiWMuSecJetCanvas, dPhiWMuSecJet_);
  Common::draw1DHistograms(tauHadPhotonEnergyFractionCanvas, tauHadPhotonEnergyFraction_);
  Common::draw1DHistograms(dThetaPhotonOtherTauConstituentsCanvas, 
			   dThetaPhotonOtherTauConstituents_);

  //format and draw 2D plots
  Common::draw2DHistograms(cleanedJetPTVsCleanedTauPTCanvas, cleanedJetPTVsCleanedTauPT_);
  Common::draw2DHistograms(uncleanedJetPTVsCleanedTauPTCanvas, uncleanedJetPTVsCleanedTauPT_);
  Common::draw2DHistograms(muHadMassVsDRSoftMuTauCanvas, muHadMassVsDRSoftMuTau_);
  Common::draw2DHistograms(muHadMassVsSoftMuPTCanvas, muHadMassVsSoftMuPT_);
  Common::draw2DHistograms(tauHadIsoVsSoftMuPTCanvas, tauHadIsoVsSoftMuPT_);
  Common::draw2DHistograms(genMuExistsVsSoftMuNearestMuPropertiesCanvas, 
			   genMuExistsVsSoftMuNearestMuProperties_);
  Common::draw2DHistograms(muHadMassVsTauHadEtaCanvas, muHadMassVsTauHadEta_);
  Common::draw2DHistograms(muHadMassVsSoftMuEtaCanvas, muHadMassVsSoftMuEta_);
  Common::draw2DHistograms(muHadMassVsTauHadIsoCanvas, muHadMassVsTauHadIso_);
  Common::draw2DHistograms(muHadMassVsTauHadPTCanvas, muHadMassVsTauHadPT_);
  Common::draw2DHistograms(tauHadIsoVsEtaCanvas, tauHadIsoVsEta_);
  Common::draw2DHistograms(tauHadEtaVsSoftMuEtaCanvas, tauHadEtaVsSoftMuEta_);
  Common::draw2DHistograms(dEtaTauHadSoftMuVsDPhiTauHadSoftMuCanvas, 
			   dEtaTauHadSoftMuVsDPhiTauHadSoftMu_);
  Common::draw2DHistograms(tauHadPTOverMuHadMassVsTauHadIsoCanvas, 
			   tauHadPTOverMuHadMassVsTauHadIso_);
  Common::draw2DHistograms(softMuPTOverMuHadMassVsTauHadIsoCanvas, 
			   softMuPTOverMuHadMassVsTauHadIso_);
  Common::draw2DHistograms(avgTauHadSoftMuPTOverMuHadMassVsTauHadIsoCanvas, 
			   avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_);
  Common::draw2DHistograms(muHadPTOverMuHadMassVsTauHadIsoCanvas, 
			   muHadPTOverMuHadMassVsTauHadIso_);
  Common::draw2DHistograms(WMuIsoVsTauHadIsoCanvas, WMuIsoVsTauHadIso_);
  Common::draw2DHistograms(softMuPTVsTauHadPTCanvas, softMuPTVsTauHadPT_);
  Common::draw2DHistograms(muHadPTOverMuHadMassVsMWMuSoftMuCanvas, 
			   muHadPTOverMuHadMassVsMWMuSoftMu_);
  Common::draw2DHistograms(softMuPFPTVsRECOPTCanvas, softMuPFPTVsRECOPT_);
  Common::draw2DHistograms(softMuPFEtaVsRECOEtaCanvas, softMuPFEtaVsRECOEta_);
  Common::draw2DHistograms(mSecondJetMuHadVsMuHadMassCanvas, mSecondJetMuHadVsMuHadMass_);
  Common::draw2DHistograms(muHadMassVsSecondJetMassCanvas, muHadMassVsSecondJetMass_);
  Common::draw2DHistograms(muHadPTOverMVsSecondJetPTOverMCanvas, muHadPTOverMVsSecondJetPTOverM_);
  Common::draw2DHistograms(muHadMassVsSecondJetPTOverMCanvas, muHadMassVsSecondJetPTOverM_);
  Common::draw2DHistograms(tauMuTauHadJetWMuHTVsMETCanvas, tauMuTauHadJetWMuHTVsMET_);
  Common::draw2DHistograms(mWMuTauMuVsMWMuTauMuTauHadCanvas, mWMuTauMuVsMWMuTauMuTauHad_);
  Common::draw2DHistograms(mWMuTauMuVsMWMuTauMuTauHadGenFSRCanvas, 
			   mWMuTauMuVsMWMuTauMuTauHadGenFSR_);
  Common::draw2DHistograms(muHadMassVsMWMuTauMuTauHadCanvas, muHadMassVsMWMuTauMuTauHad_);
  Common::draw2DHistograms(WMuMTVsMETCanvas, WMuMTVsMET_);
  Common::draw2DHistograms(muHad_t3t1VsptmjCanvas, muHad_t3t1Vsptmj_);
  Common::draw2DHistograms(muHad_t3t1VsDecayModeCanvas, muHad_t3t1VsDecayMode_);
  Common::draw2DHistograms(muHadMassVsNAddlJetsCanvas, muHadMassVsNAddlJets_);

  //set custom options

  //write output files
  out_->cd();
  METCanvas.Write();
  hadTauAssociatedMuMultiplicityCanvas.Write();
  muHadMassCanvas.Write();
  muHadMassReweightErrSqCanvas.Write();
  muHadChargeCanvas.Write();
  bTagDiscrimCanvas.Write();
  WMuIsoCanvas.Write();
  WMuMTCanvas.Write();
  tauMuMTCanvas.Write();
  tauHadMTCanvas.Write();
  dPhiWMuMETCanvas.Write();
  dPhiTauMuMETCanvas.Write();
  tauMuTauHadJetHTCanvas.Write();
  diJetHTCanvas.Write();
  jetTauJetHTCanvas.Write();
  tauMuTauHadJetWMuHTCanvas.Write();
  tauMuTauHadJetWMuMETHTCanvas.Write();
  diJetWMuHTCanvas.Write();
  jetTauJetWMuHTCanvas.Write();
  dRWMuSoftMuCanvas.Write();
  dPhiWMuSoftMuCanvas.Write();
  dPhiWMuSoftMuWithCutCanvas.Write();
  dRWMuSoftGenMatchedMuCanvas.Write();
  dRWMuSoftMuMuHadMassGe2Canvas.Write();
  tauMuPTCanvas.Write();
  tauHadPTCanvas.Write();
  tauHadIsoCanvas.Write();
  WMuLeptonRelIsoCanvas.Write();
  tauHadEtaCanvas.Write();
  softMuPTOverMuHadMassCanvas.Write();
  muHadPTOverMuHadMassCanvas.Write();
  dRSoftMuNearestGenMuHistCanvas.Write();
  muHadPTCanvas.Write();
  muHadMultiplicityCanvas.Write();
  nGoodVtxCanvas.Write();
  mWMuTauMuCanvas.Write();
  PDGIDNearestStatus1GenParticleToSoftMuCanvas.Write();
  PDGIDMuSrcCanvas.Write();
  mSecondJetMuHadCanvas.Write();
  dPhiMuHadSecondJetCanvas.Write();
  muHadUncleanedJetPTRankCanvas.Write();
  nAddlHardMuonsCanvas.Write();
  nAddlJetsPTGeq20Canvas.Write();
  nAddlJetsPTGeq30Canvas.Write();
  nAddlJetsPTGeq40Canvas.Write();
  tauHadDecayModeCanvas.Write();
  dRSoftMuNearestGenZOrTTMuCanvas.Write();
  dRSoftMuNearestGenZOrTTMuFSRCanvas.Write();
  genMuInfoCanvas.Write();
  WMu3rdTightMuChargeProductCanvas.Write();
  tauHadPhotonEnergyFractionCanvas.Write();
  dThetaPhotonOtherTauConstituentsCanvas.Write();
  cleanedJetPTVsCleanedTauPTCanvas.Write();
  uncleanedJetPTVsCleanedTauPTCanvas.Write();
  muHadMassVsDRSoftMuTauCanvas.Write();
  muHadMassVsSoftMuPTCanvas.Write();
  tauHadIsoVsSoftMuPTCanvas.Write();
  genMuExistsVsSoftMuNearestMuPropertiesCanvas.Write();
  muHadMassVsTauHadEtaCanvas.Write();
  muHadMassVsSoftMuEtaCanvas.Write();
  muHadMassVsTauHadIsoCanvas.Write();
  muHadMassVsTauHadPTCanvas.Write();
  tauHadIsoVsEtaCanvas.Write();
  tauHadEtaVsSoftMuEtaCanvas.Write();
  dEtaTauHadSoftMuVsDPhiTauHadSoftMuCanvas.Write();
  tauHadPTOverMuHadMassVsTauHadIsoCanvas.Write();
  softMuPTOverMuHadMassVsTauHadIsoCanvas.Write();
  avgTauHadSoftMuPTOverMuHadMassVsTauHadIsoCanvas.Write();
  muHadPTOverMuHadMassVsTauHadIsoCanvas.Write();
  WMuIsoVsTauHadIsoCanvas.Write();
  softMuPTVsTauHadPTCanvas.Write();
  muHad_t3t1Canvas.Write();
  muHad_t2t1Canvas.Write();
  muHad_t3t1_pT1020Canvas.Write();
  muHad_t3t1_pT2030Canvas.Write();
  muHad_t3t1_pT3040Canvas.Write();
  muHad_t3t1_pT4050Canvas.Write();
  muHad_t3t1_pT50UpCanvas.Write();
  muHad_t3t1_0JetsCanvas.Write();
  muHad_t3t1_1JetsCanvas.Write();
  muHad_t3t1_2JetsCanvas.Write();
  muHad_t3t1_3JetsCanvas.Write();
  muHad_t3t1_4JetsCanvas.Write();
  muHad_t3t1_5JetsCanvas.Write();
  muHad_t3t1_MoreJetsCanvas.Write();
  muHad_Nchtrk_0_Canvas.Write();
  muHad_Nchtrk_1_Canvas.Write();
  muHad_Nchtrk_10_Canvas.Write();
  muHad_Nchtrk_30_Canvas.Write();
  second_Nchtrk_0_Canvas.Write();
  second_Nchtrk_1_Canvas.Write();
  second_Nchtrk_10_Canvas.Write();
  second_Nchtrk_30_Canvas.Write();
  muHadPTOverMuHadMassVsMWMuSoftMuCanvas.Write();
  softMuPFPTVsRECOPTCanvas.Write();
  softMuPFEtaVsRECOEtaCanvas.Write();
  mSecondJetMuHadVsMuHadMassCanvas.Write();
  muHadMassVsSecondJetMassCanvas.Write();
  muHadPTOverMVsSecondJetPTOverMCanvas.Write();
  muHadMassVsSecondJetPTOverMCanvas.Write();
  tauMuTauHadJetWMuHTVsMETCanvas.Write();
  mWMuTauMuVsMWMuTauMuTauHadCanvas.Write();
  mWMuTauMuVsMWMuTauMuTauHadGenFSRCanvas.Write();
  muHadMassVsMWMuTauMuTauHadCanvas.Write();
  muHadMassVsNAddlJetsCanvas.Write();
  WMuMTVsMETCanvas.Write();
  jet_pt_etacutCanvas.Write();
  jet_etaCanvas.Write();
  jet_phiCanvas.Write();
  jet_mass_etacutCanvas.Write();
  jet_ptmj_etacutCanvas.Write();
  dPhiWMuSecJetCanvas.Write();
  muHad_t3t1VsptmjCanvas.Write();
  muHad_t3t1VsDecayModeCanvas.Write();
  out_->Write();
  out_->Close();

  //print maximum hadronic tau PT and mu+had mass
  cout << "Maximum HPS tau pT: " << maxTauHadPT_ << " GeV\n";
  cout << "Maximum mu+had mass: " << maxMuHadMass_ << " GeV\n\n";
}

// ------------ method called when starting to processes a run  ------------
void TauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void TauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void TauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
				       edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void TauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
				     edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

double TauAnalyzer::HT(const std::vector<reco::Candidate*>& cands) const
{
  double theHT = 0.0;
  for (std::vector<reco::Candidate*>::const_iterator iCand = cands.begin(); iCand != cands.end(); 
       ++iCand) {
    theHT+=(*iCand)->pt();
  }
  return theHT;
}

void TauAnalyzer::makePTRankCanvas(TCanvas& canvas, TLegend& legend, 
				   const std::string& header, 
				   std::vector<TH1F*>& hists)
{
  drawMultiplePTHistograms(canvas, hists, pTRankColors_, pTRankStyles_, legend, 
			   pTRankEntries_, header);
}

void TauAnalyzer::drawMultiplePTHistograms(TCanvas& canvas, 
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

void TauAnalyzer::plotDPhiMuMet(const reco::MuonRef& muonRef, 
				const edm::Handle<edm::View<reco::PFMET> >& pMET, TH1F* hist, 
				const double weight)
{
  hist->Fill(fabs(reco::deltaPhi(muonRef->phi(), pMET->refAt(0)->phi())), weight);
}

void TauAnalyzer::plotMT(const reco::MuonRef& muonRef, 
			 const edm::Handle<edm::View<reco::PFMET> >& pMET, TH1F* hist, 
			 const double weight)
{
  edm::RefToBase<reco::PFMET> METRefToBase = pMET->refAt(0);
  hist->Fill(sqrt(2*muonRef->pt()*METRefToBase->et()*
		  (1.0 - cos(reco::deltaPhi(muonRef->phi(), METRefToBase->phi())))), weight);
}

void TauAnalyzer::plotHT(const std::vector<reco::Candidate*>& cands, TH1F* hist, 
			 const double weight) { hist->Fill(HT(cands), weight); }

void TauAnalyzer::reset(const bool doDelete)
{
  if (doDelete && (out_ != NULL)) delete out_;
  out_ = NULL;
  if (doDelete && (MET_ != NULL)) delete MET_;
  MET_ = NULL;
  if (doDelete && (hadTauAssociatedMuMultiplicity_ != NULL)) {
    delete hadTauAssociatedMuMultiplicity_;
  }
  hadTauAssociatedMuMultiplicity_ = NULL;
  if (doDelete && (muHadMass_ != NULL)) delete muHadMass_;
  muHadMass_ = NULL;
  if (doDelete && (muHadMassReweightErrSq_ != NULL)) delete muHadMassReweightErrSq_;
  muHadMassReweightErrSq_ = NULL;
  if (doDelete && (muHadCharge_ != NULL)) delete muHadCharge_;
  muHadCharge_ = NULL;
  if (doDelete && (WMuMT_ != NULL)) delete WMuMT_;
  WMuMT_ = NULL;
  if (doDelete && (tauMuMT_ != NULL)) delete tauMuMT_;
  tauMuMT_ = NULL;
  if (doDelete && (dPhiWMuMET_ != NULL)) delete dPhiWMuMET_;
  dPhiWMuMET_ = NULL;
  if (doDelete && (dPhiTauMuMET_ != NULL)) delete dPhiTauMuMET_;
  dPhiTauMuMET_ = NULL;
  if (doDelete && (tauMuTauHadJetHT_ != NULL)) delete tauMuTauHadJetHT_;
  tauMuTauHadJetHT_ = NULL;
  if (doDelete && (diJetHT_ != NULL)) delete diJetHT_;
  diJetHT_ = NULL;
  if (doDelete && (jetTauJetHT_ != NULL)) delete jetTauJetHT_;
  jetTauJetHT_ = NULL;
  if (doDelete && (tauMuTauHadJetWMuHT_ != NULL)) delete tauMuTauHadJetWMuHT_;
  tauMuTauHadJetWMuHT_ = NULL;
  if (doDelete && (diJetWMuHT_ != NULL)) delete diJetWMuHT_;
  diJetWMuHT_ = NULL;
  if (doDelete && (jetTauJetWMuHT_ != NULL)) delete jetTauJetWMuHT_;
  jetTauJetWMuHT_ = NULL;
  if (doDelete && (dRWMuSoftMu_ != NULL)) delete dRWMuSoftMu_;
  dRWMuSoftMu_ = NULL;
  if (doDelete && (dRWMuSoftGenMatchedMu_ != NULL)) delete dRWMuSoftGenMatchedMu_;
  dRWMuSoftGenMatchedMu_ = NULL;
  if (doDelete && (dRWMuSoftMuMuHadMassGe2_ != NULL)) delete dRWMuSoftMuMuHadMassGe2_;
  dRWMuSoftMuMuHadMassGe2_ = NULL;
  if (doDelete && (tauMuPT_ != NULL)) delete tauMuPT_;
  tauMuPT_ = NULL;
  if (doDelete && (tauHadPT_ != NULL)) delete tauHadPT_;
  tauHadPT_ = NULL;
  if (doDelete && (tauHadIso_ != NULL)) delete tauHadIso_;
  tauHadIso_ = NULL;
  if (doDelete && (WMuLeptonRelIso_ != NULL)) delete WMuLeptonRelIso_;
  WMuLeptonRelIso_ = NULL;
  if (doDelete && (tauHadEta_ != NULL)) delete tauHadEta_;
  tauHadEta_ = NULL;
  if (doDelete && (softMuPTOverMuHadMass_ != NULL)) delete softMuPTOverMuHadMass_;
  softMuPTOverMuHadMass_ = NULL;
  if (doDelete && (muHadPTOverMuHadMass_ != NULL)) delete muHadPTOverMuHadMass_;
  muHadPTOverMuHadMass_ = NULL;
  if (doDelete && (dRSoftMuNearestGenMuHist_ != NULL)) delete dRSoftMuNearestGenMuHist_;
  dRSoftMuNearestGenMuHist_ = NULL;
  if (doDelete && (muHadPT_ != NULL)) delete muHadPT_;
  muHadPT_ = NULL;
  if (doDelete && (muHadMultiplicity_ != NULL)) delete muHadMultiplicity_;
  muHadMultiplicity_ = NULL;
  if (doDelete && (nGoodVtx_ != NULL)) delete nGoodVtx_;
  nGoodVtx_ = NULL;
  if (doDelete && (mWMuTauMu_ != NULL)) delete mWMuTauMu_;
  mWMuTauMu_ = NULL;
  if (doDelete && (PDGIDNearestStatus1GenParticleToSoftMu_ != NULL)) {
    delete PDGIDNearestStatus1GenParticleToSoftMu_;
  }
  PDGIDNearestStatus1GenParticleToSoftMu_ = NULL;
  if (doDelete && (PDGIDMuSrc_ != NULL)) delete PDGIDMuSrc_;
  PDGIDMuSrc_ = NULL;
  if (doDelete && (mSecondJetMuHad_ != NULL)) delete mSecondJetMuHad_;
  mSecondJetMuHad_ = NULL;
  if (doDelete && (dPhiMuHadSecondJet_ != NULL)) delete dPhiMuHadSecondJet_;
  dPhiMuHadSecondJet_ = NULL;
  if (doDelete && (muHadUncleanedJetPTRank_ != NULL)) delete muHadUncleanedJetPTRank_;
  muHadUncleanedJetPTRank_ = NULL;
  if (doDelete && (nAddlHardMuons_ != NULL)) delete nAddlHardMuons_;
  nAddlHardMuons_ = NULL;
  if (doDelete && (nAddlJetsPTGeq20_ != NULL)) delete nAddlJetsPTGeq20_;
  nAddlJetsPTGeq20_ = NULL;
  if (doDelete && (nAddlJetsPTGeq30_ != NULL)) delete nAddlJetsPTGeq30_;
  nAddlJetsPTGeq30_ = NULL;
  if (doDelete && (nAddlJetsPTGeq40_ != NULL)) delete nAddlJetsPTGeq40_;
  nAddlJetsPTGeq40_ = NULL;
  if (doDelete && (tauHadDecayMode_ != NULL)) delete tauHadDecayMode_;
  tauHadDecayMode_ = NULL;
  if (doDelete && (dRSoftMuNearestGenZOrTTMu_ != NULL)) delete dRSoftMuNearestGenZOrTTMu_;
  dRSoftMuNearestGenZOrTTMu_ = NULL;
  if (doDelete && (dRSoftMuNearestGenZOrTTMuFSR_ != NULL)) delete dRSoftMuNearestGenZOrTTMuFSR_;
  dRSoftMuNearestGenZOrTTMuFSR_ = NULL;
  if (doDelete && (genMuInfo_ != NULL)) delete genMuInfo_;
  genMuInfo_ = NULL;
  if (doDelete && (WMu3rdTightMuChargeProduct_ != NULL)) delete WMu3rdTightMuChargeProduct_;
  WMu3rdTightMuChargeProduct_ = NULL;
  if (doDelete && (tauHadPhotonEnergyFraction_ != NULL)) delete tauHadPhotonEnergyFraction_;
  tauHadPhotonEnergyFraction_ = NULL;
  if (doDelete && (dThetaPhotonOtherTauConstituents_ != NULL)) {
    delete dThetaPhotonOtherTauConstituents_;
  }
  dThetaPhotonOtherTauConstituents_ = NULL;
  if (doDelete && (cleanedJetPTVsCleanedTauPT_ != NULL)) delete cleanedJetPTVsCleanedTauPT_;
  cleanedJetPTVsCleanedTauPT_ = NULL;
  if (doDelete && (uncleanedJetPTVsCleanedTauPT_ != NULL)) delete uncleanedJetPTVsCleanedTauPT_;
  uncleanedJetPTVsCleanedTauPT_ = NULL;
  if (doDelete && (muHadMassVsDRSoftMuTau_ != NULL)) delete muHadMassVsDRSoftMuTau_;
  muHadMassVsDRSoftMuTau_ = NULL;
  if (doDelete && (muHadMassVsSoftMuPT_ != NULL)) delete muHadMassVsSoftMuPT_;
  muHadMassVsSoftMuPT_ = NULL;
  if (doDelete && (tauHadIsoVsSoftMuPT_ != NULL)) delete tauHadIsoVsSoftMuPT_;
  tauHadIsoVsSoftMuPT_ = NULL;
  if (doDelete && (genMuExistsVsSoftMuNearestMuProperties_ != NULL)) {
    delete genMuExistsVsSoftMuNearestMuProperties_;
  }
  genMuExistsVsSoftMuNearestMuProperties_ = NULL;
  if (doDelete && (muHadMassVsTauHadEta_ != NULL)) delete muHadMassVsTauHadEta_;
  muHadMassVsTauHadEta_ = NULL;
  if (doDelete && (muHadMassVsSoftMuEta_ != NULL)) delete muHadMassVsSoftMuEta_;
  muHadMassVsSoftMuEta_ = NULL;
  if (doDelete && (muHadMassVsTauHadIso_ != NULL)) delete muHadMassVsTauHadIso_;
  muHadMassVsTauHadIso_ = NULL;
  if (doDelete && (muHadMassVsTauHadPT_ != NULL)) delete muHadMassVsTauHadPT_;
  muHadMassVsTauHadPT_ = NULL;
  if (doDelete && (tauHadIsoVsEta_ != NULL)) delete tauHadIsoVsEta_;
  tauHadIsoVsEta_ = NULL;
  if (doDelete && (tauHadEtaVsSoftMuEta_ != NULL)) delete tauHadEtaVsSoftMuEta_;
  tauHadEtaVsSoftMuEta_ = NULL;
  if (doDelete && (dEtaTauHadSoftMuVsDPhiTauHadSoftMu_ != NULL)) {
    delete dEtaTauHadSoftMuVsDPhiTauHadSoftMu_;
  }
  dEtaTauHadSoftMuVsDPhiTauHadSoftMu_ = NULL;
  if (doDelete && (tauHadPTOverMuHadMassVsTauHadIso_ != NULL)) {
    delete tauHadPTOverMuHadMassVsTauHadIso_;
  }
  tauHadPTOverMuHadMassVsTauHadIso_ = NULL;
  if (doDelete && (softMuPTOverMuHadMassVsTauHadIso_ != NULL)) {
    delete softMuPTOverMuHadMassVsTauHadIso_;
  }
  softMuPTOverMuHadMassVsTauHadIso_ = NULL;
  if (doDelete && (avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_ != NULL)) {
    delete avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_;
  }
  avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_ = NULL;
  if (doDelete && (muHadPTOverMuHadMassVsTauHadIso_ != NULL)) {
    delete muHadPTOverMuHadMassVsTauHadIso_;
  }
  muHadPTOverMuHadMassVsTauHadIso_ = NULL;
  if (doDelete && (softMuPTVsTauHadPT_ != NULL)) delete softMuPTVsTauHadPT_;
  softMuPTVsTauHadPT_ = NULL;
  if (doDelete && (muHadPTOverMuHadMassVsMWMuSoftMu_ != NULL)) {
    delete muHadPTOverMuHadMassVsMWMuSoftMu_;
  }
  muHadPTOverMuHadMassVsMWMuSoftMu_ = NULL;
  if (doDelete && (softMuPFPTVsRECOPT_ != NULL)) delete softMuPFPTVsRECOPT_;
  softMuPFPTVsRECOPT_ = NULL;
  if (doDelete && (softMuPFEtaVsRECOEta_ != NULL)) delete softMuPFEtaVsRECOEta_;
  softMuPFEtaVsRECOEta_ = NULL;
  if (doDelete && (mSecondJetMuHadVsMuHadMass_ != NULL)) delete mSecondJetMuHadVsMuHadMass_;
  mSecondJetMuHadVsMuHadMass_ = NULL;
  if (doDelete && (muHadMassVsSecondJetMass_ != NULL)) delete muHadMassVsSecondJetMass_;
  muHadMassVsSecondJetMass_ = NULL;
  if (doDelete && (muHadPTOverMVsSecondJetPTOverM_ != NULL)) {
    delete muHadPTOverMVsSecondJetPTOverM_;
  }
  muHadPTOverMVsSecondJetPTOverM_ = NULL;
  if (doDelete && (muHadMassVsSecondJetPTOverM_ != NULL)) delete muHadMassVsSecondJetPTOverM_;
  muHadMassVsSecondJetPTOverM_ = NULL;
  if (doDelete && (tauMuTauHadJetWMuHTVsMET_ != NULL)) delete tauMuTauHadJetWMuHTVsMET_;
  tauMuTauHadJetWMuHTVsMET_ = NULL;
  if (doDelete && (mWMuTauMuVsMWMuTauMuTauHad_ != NULL)) delete mWMuTauMuVsMWMuTauMuTauHad_;
  mWMuTauMuVsMWMuTauMuTauHad_ = NULL;
  if (doDelete && (mWMuTauMuVsMWMuTauMuTauHadGenFSR_ != NULL)) {
    delete mWMuTauMuVsMWMuTauMuTauHadGenFSR_;
  }
  mWMuTauMuVsMWMuTauMuTauHadGenFSR_ = NULL;
  if (doDelete && (muHadMassVsMWMuTauMuTauHad_ != NULL)) delete muHadMassVsMWMuTauMuTauHad_;
  muHadMassVsMWMuTauMuTauHad_ = NULL;
  if (doDelete && (muHadMassVsNAddlJets_ != NULL)) delete muHadMassVsNAddlJets_;
  muHadMassVsNAddlJets_ = NULL;
  if (doDelete && (jet_pt_etacut != NULL)) delete jet_pt_etacut;
  jet_pt_etacut = NULL;
  if (doDelete && (jet_eta != NULL)) delete jet_eta;
  jet_eta = NULL;
  if (doDelete && (jet_phi != NULL)) delete jet_phi;
  jet_phi = NULL;
  if (doDelete && (jet_mass_etacut != NULL)) delete jet_mass_etacut;
  jet_mass_etacut = NULL;
  if (doDelete && (jet_ptmj_etacut != NULL)) delete jet_ptmj_etacut;
  jet_ptmj_etacut = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
