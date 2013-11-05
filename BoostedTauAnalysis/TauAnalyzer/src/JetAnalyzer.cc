// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      JetAnalyzer
// 
/**\class JetAnalyzer JetAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/JetAnalyzer.cc

   Description: analyze jet activity related to taus

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: JetAnalyzer.cc,v 1.7 2012/11/08 16:44:24 yohay Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <sstream>
#include <fstream>

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
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

//
// class declaration
//

class JetAnalyzer : public edm::EDAnalyzer {
public:
  explicit JetAnalyzer(const edm::ParameterSet&);
  ~JetAnalyzer();

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
  void fillETHistogram(edm::Handle<edm::View<U> >& pView, TH1F* hist)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->et());
  }

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

  //make visible gen tau pT canvas for 1 decay mode, multiple pT ranks
  void makePTRankCanvas(TCanvas&, TLegend&, const std::string&, std::vector<TH1F*>&);

  //format and draw multiple pT histograms on one canvas
  void drawMultiplePTHistograms(TCanvas&, std::vector<TH1F*>&, const std::vector<unsigned int>&, 
				const std::vector<unsigned int>&, TLegend&, 
				const std::vector<std::string>&, const std::string&);

  // ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //text output stream
  ofstream textOut_;

  //name of output file
  std::string outFileName_;

  //name of text output file
  std::string textOutFileName_;

  //set of parameters for GenTauDecayID class
  edm::ParameterSet genTauDecayIDPSet_;

  //flag indicating whether pT cuts should be applied in determining valid gen objects
  bool applyPTCuts_;

  //flag indicating whether KShorts should be counted as neutral hadrons
  bool countKShort_;

  //sister fabs(PDG ID) to match
  unsigned int sisterAbsMatchPDGID_;

  //vector of jet tags
  std::vector<edm::InputTag> jetTags_;

  //MET tag
  edm::InputTag METTag_;

  //input tag for base gen particle collection
  edm::InputTag genParticleTag_;

  //selected gen object tag
  edm::InputTag selectedGenObjTag_;

  //muon tag
  edm::InputTag muonTag_;

  //dR matching cut
  double dR_;

  /*flag indicating whether gen tau matching should be performed (it's time intensive if you don't 
    need it)*/
  bool doGenTauFinding_;

  //marker colors for histograms with different pT rank
  std::vector<unsigned int> pTRankColors_;

  //marker styles for histograms with different pT rank
  std::vector<unsigned int> pTRankStyles_;

  //legend entries for histograms with different pT rank
  std::vector<std::string> pTRankEntries_;

  //bins of (DR(gen mu, jet), gen mu pT)
  std::vector<std::pair<std::string, std::string> > muDRJetPTBins_;

  //histogram of eta for the single collection
  TH1F* eta_;

  //vector of jet pT histograms, highest pT rank first
  std::vector<TH1F*> jetPTHists_;

  //histogram of MET
  TH1F* MET_;

  //histogram of HT
  TH1F* HT_;

  //histogram of muon energy fraction
  TH1F* muEnergyFraction_;

  //histogram of muon energy fraction, zoomed in on [0, 0.02)
  TH1F* muEnergyFractionZoom_;

  //histogram of gen muon pT in mu+had objects matched to jets with 0 muons
  TH1F* muHadGenMuPT_;

  //histogram of DR(gen muon, jet) in mu+had objects matched to jets with 0 muons
  TH1F* muHadDRMuJet_;

  //histogram of muon multiplicity in events with 1 gen-matched jet
  TH1F* muMultiplicityJet1_;

  //histogram of muon multiplicity vs. energy fraction
  TH2F* muMultiplicityVsEnergyFraction_;

  //histogram of gen muon pT vs. DR(gen muon, jet) in mu+had objects matched to jets with 0 muons
  TH2F* muHadGenMuPTVsDRMuJet_;

  //histogram of muon multiplicity in events with 2 gen-matched jets
  TH2F* muMultiplicityJet2VsJet1_;

  /*histogram of reco muon match vs. gen muon in detector acceptance for events with 1 gen-matched 
    jet*/
  std::vector<TH2F*> recoMuExistsVsGenMuInAcceptanceJet1_;

  //histogram of the parent parton of the jet
  TH1F* jetParentParton_;

  //selected gen object counter
  unsigned int nSelectedGenObjs_;

  //vector of reco jet counters
  std::vector<unsigned int> nRecoJets_;

  //counter of events with all jet collections filled once
  unsigned int nEvts1RecoJetPerType_;
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
JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  textOutFileName_(iConfig.getParameter<std::string>("textOutFileName")),
  genTauDecayIDPSet_(iConfig.getParameter<edm::ParameterSet>("genTauDecayIDPSet")),
  applyPTCuts_(iConfig.getParameter<bool>("applyPTCuts")),
  countKShort_(iConfig.getParameter<bool>("countKShort")),
  sisterAbsMatchPDGID_(iConfig.getParameter<unsigned int>("sisterAbsMatchPDGID")),
  jetTags_(iConfig.getParameter<std::vector<edm::InputTag> >("jetTags")),
  METTag_(iConfig.getParameter<edm::InputTag>("METTag")),
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  selectedGenObjTag_(iConfig.getParameter<edm::InputTag>("selectedGenObjTag")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  dR_(iConfig.getParameter<double>("dR")),
  doGenTauFinding_(iConfig.getParameter<bool>("doGenTauFinding")),
  pTRankColors_(iConfig.getParameter<std::vector<unsigned int> >("pTRankColors")),
  pTRankStyles_(iConfig.getParameter<std::vector<unsigned int> >("pTRankStyles")),
  pTRankEntries_(iConfig.getParameter<std::vector<std::string> >("pTRankEntries")),
  nSelectedGenObjs_(0),
  nRecoJets_(std::vector<unsigned int>(jetTags_.size(), 0)),
  nEvts1RecoJetPerType_(0)
{
  //now do what ever initialization is needed
  muDRJetPTBins_.push_back(std::make_pair("LessThan0p5", "LessThan10"));
  muDRJetPTBins_.push_back(std::make_pair("LessThan0p5", "GreaterThanOrEqualTo10"));
  muDRJetPTBins_.push_back(std::make_pair("GreaterThanOrEqualTo0p5", "LessThan10"));
  muDRJetPTBins_.push_back(std::make_pair("GreaterThanOrEqualTo0p5", "GreaterThanOrEqualTo10"));
  reset(false);
}

JetAnalyzer::~JetAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void JetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get jet tags
  std::map<TH1F*, edm::Handle<edm::View<reco::PFJet> > > jetCollMap;
  for (std::vector<edm::InputTag>::const_iterator iTag = jetTags_.begin(); 
       iTag != jetTags_.end(); ++iTag) {
    edm::Handle<edm::View<reco::PFJet> > pView;
    iEvent.getByLabel(*iTag, pView);
    jetCollMap[jetPTHists_[iTag - jetTags_.begin()]] = pView;
  }
  
  //get MET tag
  edm::Handle<edm::View<reco::PFMET> > pMET;
  iEvent.getByLabel(METTag_, pMET);

  //get base gen particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //get selected gen object tag
  edm::Handle<reco::GenParticleRefVector> pSelectedGenObjs;
  iEvent.getByLabel(selectedGenObjTag_, pSelectedGenObjs);

  //get muons
  edm::Handle<reco::MuonRefVector> pMuons;
  iEvent.getByLabel(muonTag_, pMuons);
  
  //fill pT histograms for jets, 1 per pT rank
  for (std::map<TH1F*, edm::Handle<edm::View<reco::PFJet> > >::iterator i = jetCollMap.begin(); 
       i != jetCollMap.end(); ++i) { fillPTHistogramArbitrary(i->second, i->first); }
  
  //plot MET distribution
  fillETHistogram(pMET, MET_);
  
  /*plot HT distribution (sum ET of jets matched to di-tau objects in the case where there are 2 
    di-tau objects with jet matches)*/
  fillHTHistogram(jetCollMap, 2);

  //size of selected gen object collection
  const unsigned int nSelectedGenObjs = pSelectedGenObjs->size();

  //STL container of muon pointers
  std::vector<reco::Muon*> recoMuPtrVec;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) { recoMuPtrVec.push_back(const_cast<reco::Muon*>((*iMuon).get())); }
  Common::sortByPT(recoMuPtrVec);

  //fill an STL container of pointers to the gen particles
  std::vector<reco::GenParticle*> genParticlePtrs;
  for (reco::GenParticleRefVector::const_iterator iGenParticle = pSelectedGenObjs->begin(); 
       iGenParticle != pSelectedGenObjs->end(); ++iGenParticle) {
    genParticlePtrs.push_back(const_cast<reco::GenParticle*>((*iGenParticle).get()));
  }

  //quantities to store for each gen-matched jet
  std::vector<bool> muMultiplicity0(2, true);
  std::vector<bool> genMuInsideMuonAcceptance(2, true);
  std::vector<bool> genMuRecoMatchExists(2, true);
  std::vector<bool> muHadDRMuJetLessThan0p5(2, true);
  std::vector<bool> genMuPTLessThan10(2, true);
  
  //plot muon energy fraction of jets
  const unsigned int nGenMatchedRecoJets = jetCollMap.begin()->second->size();
  for (unsigned int iJet = 0; iJet < nGenMatchedRecoJets; ++iJet) {
    edm::RefToBase<reco::PFJet> jetRefToBase(jetCollMap.begin()->second->refAt(iJet));
    float muEnergyFraction = jetRefToBase->muonEnergyFraction();
    muEnergyFraction_->Fill(muEnergyFraction);

    //zoom in on jets with muon energy fraction < 0.02
    if (muEnergyFraction < 0.02) muEnergyFractionZoom_->Fill(muEnergyFraction);

    //plot muon multiplicity in jets with muon energy fraction =0.0
    int muMultiplicity = jetRefToBase->muonMultiplicity();
    if (muEnergyFraction == 0.0) {
      muMultiplicityVsEnergyFraction_->Fill(0.0, muMultiplicity);
      muMultiplicity0[iJet] = true;

      //find the matching gen tau (1st within dR_ of jet)
      GenTauDecayID genHadTau;
      std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> genHadTauDecayType;
      bool foundMatch = false;
      unsigned int iSelectedGenObj = 0;
      try {
	while (doGenTauFinding_ && (iSelectedGenObj < nSelectedGenObjs) && !foundMatch) {
	  GenTauDecayID tauDecay(genTauDecayIDPSet_, pGenParticles, 
				 Common::getStatus3Key(pSelectedGenObjs, pGenParticles, 
						       iSelectedGenObj));
	  std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> decayType = 
	    tauDecay.tauDecayType(applyPTCuts_, countKShort_);
	  reco::LeafCandidate::LorentzVector visibleGenP4 = tauDecay.getVisibleTauP4();
	  std::vector<reco::LeafCandidate> 
	    visibleGenParticle(1, reco::LeafCandidate(0.0, visibleGenP4));
	  edm::Ref<std::vector<reco::LeafCandidate> > 
	    visibleGenParticleRef(&visibleGenParticle, 0);
	  if (reco::deltaR(*jetRefToBase, *visibleGenParticleRef) < dR_) {
	    foundMatch = true;
	    genHadTau = tauDecay;
	    genHadTauDecayType = decayType;
	  }
	  ++iSelectedGenObj;
	}
      }
      catch (std::string& ex) { throw cms::Exception("JetAnalyzer") << ex; }

      //if the matched tau decayed hadronically...
      if ((foundMatch == true) && (genHadTauDecayType.second == GenTauDecayID::HAD)) {

	//...and its sister decayed muonically...
	try {
	  genHadTau.findSister();
	  if (genHadTau.sisterDecayType(applyPTCuts_, countKShort_).second == GenTauDecayID::MU) {
	    reco::LeafCandidate::LorentzVector genMuP4 = genHadTau.getVisibleTauSisterP4();
	    const double genMuPT = genMuP4.Pt();
	    const double genMuEta = genMuP4.Eta();
	    const double dRMuJet = 
	      reco::deltaR(jetRefToBase->eta(), jetRefToBase->phi(), genMuP4.Eta(), genMuP4.Phi());
	    std::vector<reco::LeafCandidate> genMuVec(1, reco::LeafCandidate(0.0, genMuP4));
	    edm::Ref<std::vector<reco::LeafCandidate> > genMuRef(&genMuVec, 0);

	    //plot the matched mu+had gen muon pT and DR(muon, jet)
	    muHadGenMuPT_->Fill(genMuPT);
	    muHadDRMuJet_->Fill(dRMuJet);
	    muHadGenMuPTVsDRMuJet_->Fill(dRMuJet, genMuPT);

	    //is the gen muon inside the muon detector acceptance?
	    if (fabs(genMuEta) < 2.1) genMuInsideMuonAcceptance[iJet] = true;
	    else genMuInsideMuonAcceptance[iJet] = false;

	    //can a reco muon match the the gen muon be found?
	    int iNearestRecoMuon = -1;
	    const reco::Muon* pNearestRecoMuon = 
	      Common::nearestObject(genMuRef, recoMuPtrVec, iNearestRecoMuon);
	    if ((pNearestRecoMuon != NULL) && 
		(reco::deltaR(genMuEta, genMuP4.Phi(), 
			      pNearestRecoMuon->eta(), pNearestRecoMuon->phi()) < dR_)) {
	      genMuRecoMatchExists[iJet] = true;
	    }
	    else genMuRecoMatchExists[iJet] = false;

	    //is the gen muon ostensibly inside the jet cone?
	    if (dRMuJet < 0.5) muHadDRMuJetLessThan0p5[iJet] = true;
	    else muHadDRMuJetLessThan0p5[iJet] = false;

	    //is the gen muon within the muon pT acceptance?
	    if (genMuPT < 10.0) genMuPTLessThan10[iJet] = true;
	    else genMuPTLessThan10[iJet] = false;

// 	    //print event and jet numbers when gen muon appears to be inside jet cone
// 	    if ((dRMuJet >= 0.2) && (dRMuJet < 0.4) && (genMuPT > 10.0)) {
// 	      textOut_ << "inside high " << iEvent.run() << " " << iEvent.id().event() << " ";
// 	      textOut_ << iEvent.id().luminosityBlock() << " " << iJet << std::endl;
// 	    }

// 	    //print event and jet numbers when gen muon appears to be outside jet cone
// 	    if ((dRMuJet >= 0.5) && (dRMuJet < 0.7)) {
// 	      textOut_ << "outside ";
// 	      if (genMuPT > 10.0) textOut_ << "high ";
// 	      else textOut_ << "low ";
// 	      textOut_ << iEvent.run() << " " << iEvent.id().event() << " ";
// 	      textOut_ << iEvent.id().luminosityBlock() << " " << iJet << std::endl;
// 	    }
	  }
	}
	catch (std::string& ex) {
	  if (sisterAbsMatchPDGID_ != GenTauDecayID::ANY_PDGID) { 
	    throw cms::Exception("JetAnalyzer") << ex; /*assuming you wanted a sister, so if one 
							 isn't found an exception should be 
							 thrown*/
	  }
	  //else assume particle really didn't have sister and you knew that, so just move on
	}
      }
    }

    //plot muon multiplicity in jets with muon energy fraction >0.0
    else if (muEnergyFraction > 0.0) {
      muMultiplicityVsEnergyFraction_->Fill(1.0, muMultiplicity);
      muMultiplicity0[iJet] = false;
    }

    //plot the parent parton of the jet, excluding selected muons from the jet collection
    int nearestPartonIndex = -1;
    const reco::GenParticle* nearestParton = 
      Common::nearestObject(jetRefToBase, genParticlePtrs, nearestPartonIndex);
    if (nearestParton != NULL) {
      if ((reco::deltaR(*recoMuPtrVec[recoMuPtrVec.size() - 1], *jetRefToBase) >= dR_) && 
	  (reco::deltaR(*nearestParton, *jetRefToBase) < dR_)) {
	const unsigned int absNearestPartonPDGID = fabs(nearestParton->pdgId());
	unsigned int bin = absNearestPartonPDGID == GenTauDecayID::G ? 0 : absNearestPartonPDGID;
	jetParentParton_->Fill(bin);
      }
      else jetParentParton_->Fill(-1);
    }
    else jetParentParton_->Fill(7);
  }

  //sanity check
  if ((genMuInsideMuonAcceptance.size() != genMuRecoMatchExists.size()) || 
      (genMuRecoMatchExists.size() != muHadDRMuJetLessThan0p5.size()) || 
      (muHadDRMuJetLessThan0p5.size() != genMuPTLessThan10.size())) {
    throw cms::Exception("JetAnalyzer") << "Unequal vector sizes.\n";
  }

  //this event has exactly 1 gen-matched jet...
  if (nGenMatchedRecoJets == 1) {

    //fill muon multiplicity plot
    muMultiplicityJet1_->Fill(!muMultiplicity0[0]);

    //fill reco mu exists vs. gen mu in acceptance plot for jets with 0 muons
    if (muMultiplicity0[0]) {
      unsigned int iHist = 0;
      if (muHadDRMuJetLessThan0p5[0]) {
	if (genMuPTLessThan10[0]) iHist = 0;
	else iHist = 1;
      }
      else {
	if (genMuPTLessThan10[0]) iHist = 2;
	else iHist = 3;
      }
      recoMuExistsVsGenMuInAcceptanceJet1_[iHist]->Fill(!genMuInsideMuonAcceptance[0], 
							!genMuRecoMatchExists[0]);
      textOut_ << iHist << " " << !genMuInsideMuonAcceptance[0] << " " << !genMuRecoMatchExists[0];
      textOut_ << " " << iEvent.run() << " " << iEvent.id().event() << std::endl;
    }
  }

  //this event has exactly 2 gen-matched jets...
  else if (nGenMatchedRecoJets == 2) {

    //fill muon multiplicity plot
    muMultiplicityJet2VsJet1_->Fill(!muMultiplicity0[0], !muMultiplicity0[1]);
  }

  //count selected gen objects
  nSelectedGenObjs_+=nSelectedGenObjs;

  //count reco jets
  unsigned int count = 0;
  unsigned int nRecoJets = 0;
  for (std::map<TH1F*, edm::Handle<edm::View<reco::PFJet> > >::iterator i = jetCollMap.begin(); 
       i != jetCollMap.end(); ++i) {
    nRecoJets_[count]+=i->second->size();
    ++count;
    nRecoJets+=i->second->size();
  }

  //count events with 1 reco jet per type
  if (nRecoJets == nSelectedGenObjs) ++nEvts1RecoJetPerType_;
}


// ------------ method called once each job just before starting event loop  ------------
void JetAnalyzer::beginJob()
{
  //open output files
  out_ = new TFile(outFileName_.c_str(), "RECREATE");
  textOut_.open(textOutFileName_.c_str());

  //book jet histograms
  for (std::vector<edm::InputTag>::const_iterator iTag = jetTags_.begin(); 
       iTag != jetTags_.end(); ++iTag) {
    jetPTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
				   .c_str(), "", 20, 0.0, 100.0));
  }
  MET_ = new TH1F("MET", "", 20, 0.0, 100.0);
  HT_ = new TH1F("HT", "", 100, 0.0, 500.0);
  muEnergyFraction_ = new TH1F("muEnergyFraction", ";Muon energy fraction;", 50, 0.0, 1.0);
  muEnergyFractionZoom_ = 
    new TH1F("muEnergyFractionZoom", ";Muon energy fraction;", 50, 0.0, 0.02);
  muHadGenMuPT_ = new TH1F("muHadGenMuPT", ";p_{T} (GeV);", 20, 0.0, 100.0);
  muHadDRMuJet_ = new TH1F("muHadDRMuJet", ";#DeltaR_{#muj};", 50, 0.0, 5.0);
  muMultiplicityJet1_ = new TH1F("muMultiplicityJet1", ";Muon multiplicity;", 2, -0.5, 1.5);
  muMultiplicityVsEnergyFraction_ = 
    new TH2F("muMultiplicityVsEnergyFraction", ";Muon energy fraction;Muon multiplicity", 
	     2, -0.5, 1.5, 4, -0.5, 3.5);
  muHadGenMuPTVsDRMuJet_ = 
    new TH2F("muHadGenMuPTVsDRMuJet", ";#DeltaR_{#muj};p_{T} (GeV)", 50, 0.0, 5.0, 20, 0.0, 100.0);
  muMultiplicityJet2VsJet1_ = 
    new TH2F("muMultiplicityJet2VsJet1", ";Jet 1 muon multiplicity;Jet 2 muon multiplicity", 
	     2, -0.5, 1.5, 2, -0.5, 1.5);
  for (std::vector<std::pair<std::string, std::string> >::const_iterator iMuDRJetPTBin = 
	 muDRJetPTBins_.begin(); iMuDRJetPTBin != muDRJetPTBins_.end(); ++iMuDRJetPTBin) {
    std::stringstream title;
    title << ";Gen muon |#eta| < 2.1;#DeltaR_{reco #mu, gen #mu} < " << dR_;
    recoMuExistsVsGenMuInAcceptanceJet1_.
      push_back(new TH2F(("recoMuExistsVsGenMuInAcceptanceJet1_dRGenMuJet" + 
			  iMuDRJetPTBin->first + "_genMuPT" + iMuDRJetPTBin->second).c_str(), 
			 title.str().c_str(), 2, -0.5, 1.5, 2, -0.5, 1.5));
  }
  jetParentParton_ = new TH1F("jetParentParton", ";Jet parent parton;", 7, -0.5, 6.5);

  //set bin labels
  muMultiplicityJet1_->GetXaxis()->SetBinLabel(1, "0");
  muMultiplicityJet1_->GetXaxis()->SetBinLabel(2, ">0");
  muMultiplicityVsEnergyFraction_->GetXaxis()->SetBinLabel(1, "0.0");
  muMultiplicityVsEnergyFraction_->GetXaxis()->SetBinLabel(2, ">0.0");
  muMultiplicityJet2VsJet1_->GetXaxis()->SetBinLabel(1, "0");
  muMultiplicityJet2VsJet1_->GetXaxis()->SetBinLabel(2, ">0");
  muMultiplicityJet2VsJet1_->GetYaxis()->SetBinLabel(1, "0");
  muMultiplicityJet2VsJet1_->GetYaxis()->SetBinLabel(2, ">0");
  for (unsigned int iMuDRJetPTBin = 0; iMuDRJetPTBin < muDRJetPTBins_.size(); ++iMuDRJetPTBin) {
    recoMuExistsVsGenMuInAcceptanceJet1_[iMuDRJetPTBin]->GetXaxis()->SetBinLabel(1, "Yes");
    recoMuExistsVsGenMuInAcceptanceJet1_[iMuDRJetPTBin]->GetXaxis()->SetBinLabel(2, "No");
    recoMuExistsVsGenMuInAcceptanceJet1_[iMuDRJetPTBin]->GetYaxis()->SetBinLabel(1, "Yes");
    recoMuExistsVsGenMuInAcceptanceJet1_[iMuDRJetPTBin]->GetYaxis()->SetBinLabel(2, "No");
  }
  jetParentParton_->GetXaxis()->SetBinLabel(1, "g");
  jetParentParton_->GetXaxis()->SetBinLabel(2, "d");
  jetParentParton_->GetXaxis()->SetBinLabel(3, "u");
  jetParentParton_->GetXaxis()->SetBinLabel(4, "s");
  jetParentParton_->GetXaxis()->SetBinLabel(5, "c");
  jetParentParton_->GetXaxis()->SetBinLabel(6, "b");
  jetParentParton_->GetXaxis()->SetBinLabel(7, "t");
}

// ------------ method called once each job just after ending the event loop  ------------
void JetAnalyzer::endJob() 
{
  //make the jet canvases
  out_->cd();
  TCanvas jetPTRankCanvas("jetPTRankCanvas", "", 600, 600);
  TCanvas muEnergyFractionCanvas("muEnergyFractionCanvas", "", 600, 600);
  TCanvas muEnergyFractionZoomCanvas("muEnergyFractionZoomCanvas", "", 600, 600);
  TCanvas muHadGenMuPTCanvas("muHadGenMuPTCanvas", "", 600, 600);
  TCanvas muHadDRMuJetCanvas("muHadDRMuJetCanvas", "", 600, 600);
  TCanvas muMultiplicityJet1Canvas("muMultiplicityJet1Canvas", "", 600, 600);
  TCanvas muMultiplicityVsEnergyFractionCanvas("muMultiplicityVsEnergyFractionCanvas", "", 600, 
					       600);
  TCanvas muHadGenMuPTVsDRMuJetCanvas("muHadGenMuPTVsDRMuJetCanvas", "", 600, 600);
  TCanvas muMultiplicityJet2VsJet1Canvas("muMultiplicityJet2VsJet1Canvas", "", 600, 600);
  std::vector<TCanvas*> recoMuExistsVsGenMuInAcceptanceJet1Canvases;
  for (std::vector<TH2F*>::const_iterator iHist = recoMuExistsVsGenMuInAcceptanceJet1_.begin(); 
       iHist != recoMuExistsVsGenMuInAcceptanceJet1_.end(); ++iHist) {
    recoMuExistsVsGenMuInAcceptanceJet1Canvases.
      push_back(new TCanvas((std::string((*iHist)->GetName()) + "_canvas").c_str(), "", 600, 600));
  }
  TCanvas jetParentPartonCanvas("jetParentPartonCanvas", "", 600, 600);

  //make legends
  TLegend jetPTRankLegend(0.4, 0.6, 0.8, 0.8);

  //make jet pT rank canvases
  makePTRankCanvas(jetPTRankCanvas, jetPTRankLegend, 
		   "gg fusion NMSSM Higgs-matched AK5 jets", jetPTHists_);

  //format and draw 1D plots
  Common::draw1DHistograms(muEnergyFractionCanvas, muEnergyFraction_);
  Common::draw1DHistograms(muEnergyFractionZoomCanvas, muEnergyFractionZoom_);
  Common::draw1DHistograms(muHadGenMuPTCanvas, muHadGenMuPT_);
  Common::draw1DHistograms(muHadDRMuJetCanvas, muHadDRMuJet_);
  Common::draw1DHistograms(muMultiplicityJet1Canvas, muMultiplicityJet1_);
  Common::draw1DHistograms(jetParentPartonCanvas, jetParentParton_);

  //format and draw 2D plots
  Common::draw2DHistograms(muMultiplicityVsEnergyFractionCanvas, muMultiplicityVsEnergyFraction_);
  Common::draw2DHistograms(muHadGenMuPTVsDRMuJetCanvas, muHadGenMuPTVsDRMuJet_);
  Common::draw2DHistograms(muMultiplicityJet2VsJet1Canvas, muMultiplicityJet2VsJet1_);
  for (std::vector<TCanvas*>::iterator iCanvas = 
	 recoMuExistsVsGenMuInAcceptanceJet1Canvases.begin(); 
       iCanvas != recoMuExistsVsGenMuInAcceptanceJet1Canvases.end(); ++iCanvas) {
    const unsigned int i = iCanvas - recoMuExistsVsGenMuInAcceptanceJet1Canvases.begin();
    Common::draw2DHistograms(**iCanvas, recoMuExistsVsGenMuInAcceptanceJet1_[i]);
  }

  //write output files
  out_->cd();
  for (std::vector<TH1F*>::iterator iHist = jetPTHists_.begin(); 
       iHist != jetPTHists_.end(); ++iHist) { (*iHist)->Write(); }
  muEnergyFractionCanvas.Write();
  muEnergyFractionZoomCanvas.Write();
  muHadGenMuPTCanvas.Write();
  muHadDRMuJetCanvas.Write();
  muMultiplicityJet1Canvas.Write();
  muMultiplicityVsEnergyFractionCanvas.Write();
  muHadGenMuPTVsDRMuJetCanvas.Write();
  muMultiplicityJet2VsJet1Canvas.Write();
  for (std::vector<TCanvas*>::iterator iCanvas = 
	 recoMuExistsVsGenMuInAcceptanceJet1Canvases.begin(); 
       iCanvas != recoMuExistsVsGenMuInAcceptanceJet1Canvases.end(); ++iCanvas) {
    (*iCanvas)->Write();
  }
  jetParentPartonCanvas.Write();
  MET_->Write();
  HT_->Write();
  out_->Write();
  out_->Close();
  for (std::vector<TH2F*>::iterator iHist = recoMuExistsVsGenMuInAcceptanceJet1_.begin(); 
       iHist != recoMuExistsVsGenMuInAcceptanceJet1_.end(); ++iHist) { delete *iHist; }
  for (std::vector<TCanvas*>::iterator iCanvas = 
	 recoMuExistsVsGenMuInAcceptanceJet1Canvases.begin(); 
       iCanvas != recoMuExistsVsGenMuInAcceptanceJet1Canvases.end(); ++iCanvas) {
    delete *iCanvas;
  }
  textOut_.close();

  //print counters
  std::cout << "No. selected gen objects: " << nSelectedGenObjs_ << std::endl;
  std::cout << "No. reco jets:\n";
  for (std::vector<unsigned int>::const_iterator iCounter = nRecoJets_.begin(); 
       iCounter != nRecoJets_.end(); ++iCounter) {
    std::cout << " In collection " << iCounter - nRecoJets_.begin() << ": " << *iCounter;
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "No. events with 1 reco jet per type: " << nEvts1RecoJetPerType_ << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void JetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void JetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void JetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
				       edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void JetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
				     edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void JetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void JetAnalyzer::makePTRankCanvas(TCanvas& canvas, TLegend& legend, 
				   const std::string& header, 
				   std::vector<TH1F*>& hists)
{
  drawMultiplePTHistograms(canvas, hists, pTRankColors_, pTRankStyles_, legend, 
			   pTRankEntries_, header);
}

void JetAnalyzer::drawMultiplePTHistograms(TCanvas& canvas, 
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

void JetAnalyzer::reset(const bool doDelete)
{
  if (doDelete && (out_ != NULL)) delete out_;
  out_ = NULL;
  for (std::vector<TH1F*>::iterator iHist = jetPTHists_.begin(); 
       iHist != jetPTHists_.end(); ++iHist) {
    if (doDelete && (*iHist != NULL)) delete *iHist;
    *iHist = NULL;
  }
  if (doDelete && (MET_ != NULL)) delete MET_;
  MET_ = NULL;
  if (doDelete && (HT_ != NULL)) delete HT_;
  HT_ = NULL;
  if (doDelete && (muEnergyFraction_ != NULL)) delete muEnergyFraction_;
  muEnergyFraction_ = NULL;
  if (doDelete && (muEnergyFractionZoom_ != NULL)) delete muEnergyFractionZoom_;
  muEnergyFractionZoom_ = NULL;
  if (doDelete && (muHadGenMuPT_ != NULL)) delete muHadGenMuPT_;
  muHadGenMuPT_ = NULL;
  if (doDelete && (muHadDRMuJet_ != NULL)) delete muHadDRMuJet_;
  muHadDRMuJet_ = NULL;
  if (doDelete && (muMultiplicityJet1_ != NULL)) delete muMultiplicityJet1_;
  muMultiplicityJet1_ = NULL;
  if (doDelete && (muMultiplicityVsEnergyFraction_ != NULL)) {
    delete muMultiplicityVsEnergyFraction_;
  }
  muMultiplicityVsEnergyFraction_ = NULL;
  if (doDelete && (muHadGenMuPTVsDRMuJet_ != NULL)) delete muHadGenMuPTVsDRMuJet_;
  muHadGenMuPTVsDRMuJet_ = NULL;
  if (doDelete && (muMultiplicityJet2VsJet1_ != NULL)) delete muMultiplicityJet2VsJet1_;
  muMultiplicityJet2VsJet1_ = NULL;
  if (doDelete && (jetParentParton_ != NULL)) delete jetParentParton_;
  jetParentParton_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetAnalyzer);
