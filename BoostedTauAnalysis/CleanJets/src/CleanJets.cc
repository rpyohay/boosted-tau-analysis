// -*- C++ -*-
//
// Package:    CleanJets
// Class:      CleanJets
// 
/**\class CleanJets CleanJets.cc BoostedTauAnalysis/CleanJets/src/CleanJets.cc

 Description: Removes PF muons from PFJet candidates and reconstructs the jets
              Matches those PF muons to muons from a --> tau --> muon decays
              Studies the kinematic properties of the discarded muons
	      Associates those muons to the jets from which they were removed

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//         Created:  Fri Aug 31 13:01:48 CEST 2012
// $Id: CleanJets.cc,v 1.6 2012/12/06 17:44:58 yohay Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <sstream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTauTagInfo.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/TauReco/interface/PFTauDecayModeAssociation.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

using namespace edm;
using namespace reco;
using namespace std;


//
// class declaration
//

class CleanJets : public edm::EDProducer {
   public:
      explicit CleanJets(const edm::ParameterSet&);
      ~CleanJets();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      // source of the jets to be cleaned of muons
      edm::InputTag jetSrc_;

      // source of muons that, if found within jet, should be removed
      edm::InputTag muonSrc_;

      // source of PF candidates
      edm::InputTag PFCandSrc_;

      //input, output
      TFile* out_;
      std::string outFileName_;
      bool cutOnGenMatches_;
      edm::InputTag thisTag_;
      unsigned int momPDGID_;
      double genMuTauPTMin_;
      double genMuPTMin_;
      edm::ParameterSet* cfg_;

      //histograms
      TH1D *fracElePizeroDecays;
      TH1D *fracPFMuMatched;
      TH1D *fracGenMuMatched;
      TH1D *genMuonatm_pT;
      TH1D *genMuonZtm_pT;
      TH1D *genMuonZm_pT;
      TH1D *genMuonOther_pT;
      TH1D *genMuonatm_matched_pT;
      TH1D *genMuonatm_unmatched_pT;
      TH1D *genMuonZtm_matched_pT;
      TH1D *genMuonZm_matched_pT;
      TH1D *pfMuonatm_pT;
      TH1D *pfMuonZtm_pT;
      TH1D *pfMuonZm_pT;
      TH1D *pfMuonOther_pT;
      TH1D *genEleZ_pT;
      TH1D *genElePiZero_pT;
      TH1D *genElePhoton_pT;
      TH1D *genEleOther_pT;
      // histograms for candidates of PFJets matched to tmth
      TH1D *jetMuonatm_ET;
      TH1D *jeth_ET;
      TH1D *jete_ET;
      TH1D *jetgamma_ET;
      TH1D *jeth0_ET;
      TH1D *jetMuonnotatm_ET;
      TH1D *jetOther_ET;
      ofstream EleStuff;
      ofstream MuStuff;
      ofstream PhotonStuff;
  
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
CleanJets::CleanJets(const edm::ParameterSet& iConfig)
{
  jetSrc_ = iConfig.getParameter<edm::InputTag>("jetSrc");
  muonSrc_ = iConfig.existsAs<edm::InputTag>("muonSrc") ? 
    iConfig.getParameter<edm::InputTag>("muonSrc") : edm::InputTag();
  PFCandSrc_ = iConfig.getParameter<edm::InputTag>("PFCandSrc");
  outFileName_ = iConfig.getParameter<std::string>("outFileName");
  cutOnGenMatches_ = iConfig.getParameter<bool>("cutOnGenMatches");
  thisTag_ = iConfig.getParameter<edm::InputTag>("thisTag");
  momPDGID_ = iConfig.getParameter<unsigned int>("momPDGID");
  genMuTauPTMin_ = iConfig.getParameter<double>("genMuTauPTMin");
  genMuPTMin_ = iConfig.getParameter<double>("genMuPTMin");
  cfg_ = const_cast<edm::ParameterSet*>(&iConfig);

  //warn if user doesn't want to cut on gen matches but did not supply a muon tag
  if (!cutOnGenMatches_ && (muonSrc_ == InputTag())) {
    std::cerr << "Warning: you don't want to cut on gen matches but you did not supply a muon ";
    std::cerr << "tag.  No muons will be removed.\n";
  }

  //register your products
  
  produces<PFJetCollection>( "ak5PFJetsNoMu" );
  produces<edm::ValueMap<bool> >();
  produces<edm::ValueMap<MuonRefVector> >();
  produces<edm::ValueMap<PFJetRef> >();
  produces<PFCandidateCollection>();
  
}


CleanJets::~CleanJets()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CleanJets::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   Handle<GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles", genParticles);

   Handle<PFJetCollection> PFJets;
   iEvent.getByLabel(jetSrc_, PFJets);
   auto_ptr<PFJetCollection> SetOfJets( new PFJetCollection );

   Handle<MuonRefVector> muons;
   if (muonSrc_ == InputTag()) {}
   else iEvent.getByLabel(muonSrc_, muons);

   Handle<PFCandidateCollection> PFCands;
   iEvent.getByLabel(PFCandSrc_, PFCands);
   auto_ptr<PFCandidateCollection> PFCandsExcludingSoftMuons(new PFCandidateCollection);

   //fill an STL container with muon ref keys
   std::vector<unsigned int> muonRefKeys;
   if (muons.isValid()) {
     for (MuonRefVector::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon) {
       muonRefKeys.push_back(iMuon->key());
     }
   }

   //vector of bools holding the signal muon tag decision for each jet
   std::vector<bool> muonTagDecisions;

   //map between new jet and refs to muons in the original collection that were removed
   std::vector<MuonRefVector> removedMuonMap;

   //map between new jet and ref to original jet
   std::vector<PFJetRef> oldJets;

   std::vector<reco::GenParticle*> genMusatm; //gen muons from a -> tau -> mu (atm)
   std::vector<reco::GenParticle*> genMusatm2; // for use in PFJet analysis
   std::vector<reco::GenParticle*> genMusZtm; //gen muons from Z -> tau -> mu (Ztm)
   std::vector<reco::GenParticle*> genMusZm; //gen muons from Z -> mu (Zm)
   std::vector<reco::GenParticle*> genMusOther; //gen muons from other sources
   std::vector<GenTauDecayID> tauDecays; // vector of taus from a decays
   std::vector<reco::PFJet> pfjetVector;
   for (PFJetCollection::const_iterator j = PFJets->begin(); j != PFJets->end(); ++ j)
     {
       pfjetVector.push_back(*j);
     }

   double total_number_of_pizeros = 0.;
   double pizero_electron_decays = 0.;
   bool pizeroelecdecay = false;
   /////// Getting vector of gen mus that come from tau decays ////////
   if (cutOnGenMatches_) {
     for (GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen)
       { // loop over genparticles to get mus
	 if ((fabs(iGen->pdgId()) == 13) && (fabs(iGen->mother()->pdgId()) != 13))
	   { // if it's a mu
	     //if(iGen->status() == 1 && fabs(iGen->mother()->pdgId() == 15))
	     if(fabs(iGen->mother()->pdgId()) == 15)
	       { // if its mother was a tau
		 if ((fabs(iGen->mother()->mother()->mother()->pdgId() == 36))&&(iGen->mother()->mother()->mother()->status() == 3))
		   { // if its grandmother was an a
		     genMusatm.push_back(const_cast<reco::GenParticle*>(&*iGen));
		     genMusatm2.push_back(const_cast<reco::GenParticle*>(&*iGen));
		     genMuonatm_pT->Fill((&*iGen)->p4().pt());
		   } // if its grandmother was an a
		 else if (iGen->status() == 1 && (fabs(iGen->mother()->mother()->mother()->pdgId() == 24)))
		   { // if its grandmother was a Z/W
		     genMusZtm.push_back(const_cast<reco::GenParticle*>(&*iGen));
		     genMuonZtm_pT->Fill((&*iGen)->p4().pt());
		   } // if its grandmother was a Z/W
	       } // if its mother was a tau
	     else if(iGen->mother()->pdgId() == 24 && iGen->status()==3)
	       { // if its mother was a Z/W
		 genMusZm.push_back(const_cast<reco::GenParticle*>(&*iGen));
		 genMuonZm_pT->Fill((&*iGen)->p4().pt());
	       } // if its mother was a Z/W
	     else
	       { // if its mother was not a tau or a Z/W
		 genMusOther.push_back(const_cast<reco::GenParticle*>(&*iGen));
		 genMuonOther_pT->Fill((&*iGen)->p4().pt());
	       } // if its mother was not a tau or a Z/W
	   } // if it's a mu
	 if (fabs(iGen->pdgId()) == 11) 
	   { // if it's an electron not descended from an electron
	     if(fabs(iGen->mother()->pdgId())!= 11)
	       { // if not from an electron
		 if(fabs(iGen->mother()->pdgId()) == 24)
		   { // e from Z/W decay
		     genEleZ_pT->Fill(iGen->pt());
		   } // e from Z/W decay
		 if(fabs(iGen->mother()->pdgId()) == 15)
		   { // e from tau decay
		     genEleZ_pT->Fill(iGen->pt());
		   } // e from tau decay
		 else if(fabs(iGen->mother()->pdgId()) == 111)
		   { // e from pizero decay
		     genElePiZero_pT->Fill(iGen->pt());
		   } // e from pizero decay
		 else if(fabs(iGen->mother()->pdgId()) == 22)
		   { // e from photon
		     genElePhoton_pT->Fill(iGen->pt());
		   } // e from photon
		 else
		   {
		     genEleOther_pT->Fill(iGen->pt());
		   }
	       } // if not from an electron
	   }// if it's an electron
	 else if((fabs(iGen->pdgId()) == 111) && (fabs(iGen->mother()->pdgId()) != 111))
	   { // if it's a pizero
	     total_number_of_pizeros += 1.;
	     if(iGen->numberOfDaughters() != 0)
	       { // if it has daughters
		 for(unsigned int i = 0; i < iGen->numberOfDaughters(); ++i)
		   { // loop over daughters
		     if(fabs(iGen->daughter(i)->pdgId()) == 11)
		       pizeroelecdecay = true;
		   } // loop over daughters
	       } // if it has daughters
	     if(pizeroelecdecay)
	       {
		 pizero_electron_decays += 1.;
		 pizeroelecdecay = false;
	       }
	   } // if it's a pizero
       } // loop over genparticles to get mus

     if(total_number_of_pizeros != 0.)
       {
	 double pizerofrac = pizero_electron_decays/total_number_of_pizeros;
	 //cout << "Fraction of pizero electron decays = " << pizerofrac << endl;
	 fracElePizeroDecays->Fill(pizerofrac);
       }
     /////////// Getting gen tau decays ///////////////////////////

     for (GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen)
       { // loop over genparticles to get tau decays
       
	 if((fabs(iGen->pdgId()) == 15) && (iGen->status() == 3))
	   { // if it is a status-3 tau
	     if((fabs(iGen->mother()->pdgId()) == 36))
	       { // if its mother was an a
		 //cout << "A status 3 tau from an a has been found!" << endl;
		 //cout << "status of a = " << iGen->mother()->status() << endl;
		 try
		   { // try loop
		     GenTauDecayID tauDecay(*cfg_, genParticles, iGen - genParticles->begin());
		     //cout << "Is this tau a status 3 decay product? Boolean answer: " << tauDecay.tauIsStatus3DecayProduct() << endl;
		     tauDecays.push_back(tauDecay);
		     //cout << "We have collected this tauDecay for a status-3 tau descended from an a" << endl;
		   } // try loop
		 catch (std::string& ex) { throw cms::Exception("CheckTauSigCands") << ex; }
	       } // if its mother was an a
	   } // if it is a status-3 tau
       } // loop over genparticles to get tau decays


     ////// Look at a -> tmth decays and match PFJets


     GenParticleRefVector tausFromTMTH; // hadronically decaying taus whose sisters decayed to muons

     //loop over taus from boson decay
     std::vector<unsigned int> keysToIgnore;
     for (std::vector<GenTauDecayID>::iterator iTau = tauDecays.begin();iTau != tauDecays.end();++iTau) {
       //cout << "We have entered the tau decay loop" << endl;

       //try for exceptions
       try {
       
	 //look for a hadronic tau decay from the signal sample
	 if ((iTau->tauDecayType() == GenTauDecayID::HAD)) {
	   //cout << "Hadronically decaying tau!" << endl;
	   //look for the other tau decay product of the a: is it a tau-->mu decay?
	   iTau->findSister();
	   const unsigned int iSister = iTau->getSisterIndex();
	   if ((std::find(keysToIgnore.begin(), keysToIgnore.end(), iSister) == 
		keysToIgnore.end()) && (iTau->sisterDecayType() == GenTauDecayID::MU)) {
	     //cout << "Its sister decayed to a muon!" << endl;
	     //does the sister decay mu pass the pT threshold?
	     reco::GenParticleRef sisterRef(genParticles, iSister);
	     if (sisterRef->pt() > genMuTauPTMin_) {
	       //cout << "sister tau_mu passes threshold" << endl;
	       //get the visible had 4-vector
	       reco::GenParticleRef tauRef(genParticles, iTau->getTauIndex());
	       reco::GenParticleRef status2HadTauRef = tauRef->daughterRef(0);
	       reco::LeafCandidate::LorentzVector visibleHadTauP4;
	       for (unsigned int iDaughter = 0; iDaughter < status2HadTauRef->numberOfDaughters(); 
		    ++iDaughter) {
		 reco::GenParticleRef hadTauDaughterRef = status2HadTauRef->daughterRef(iDaughter);
		 if ((fabs(hadTauDaughterRef->pdgId()) != GenTauDecayID::ENEUTRINOPDGID) && 
		     (fabs(hadTauDaughterRef->pdgId()) != GenTauDecayID::MUNEUTRINOPDGID) && 
		     (fabs(hadTauDaughterRef->pdgId()) != GenTauDecayID::TAUNEUTRINOPDGID)) {
		   visibleHadTauP4+=hadTauDaughterRef->p4();
		 }
	       }
	       const double visibleHadTauPT = visibleHadTauP4.Pt();
	       const double visibleHadTauEta = visibleHadTauP4.Eta();
	       double muTauEta;
	       double muTauPhi;
	     
	       //get the muon 4-vector
	       int iGenMu = -1;
	       unsigned int iDaughter = 0;
	       reco::GenParticleRef status2MuTauRef = sisterRef->daughterRef(0);
	       while ((iDaughter < status2MuTauRef->numberOfDaughters()) && (iGenMu == -1)) {
		 reco::GenParticleRef muTauDaughter = status2MuTauRef->daughterRef(iDaughter);
		 if (fabs(muTauDaughter->pdgId()) == GenTauDecayID::MUPDGID) {
		   iGenMu = muTauDaughter.key();
		   muTauEta = muTauDaughter->p4().eta();
		   muTauPhi = muTauDaughter->p4().phi();
		 }
		 ++iDaughter;
	       }
	       if (iGenMu == -1) throw cms::Exception("TauAnalyzer") << "Muon not found.\n";
	       reco::GenParticleRef muRef(genParticles, iGenMu);
	       const double genMuPT = muRef->pt();
	       const double genMuEta = muRef->eta();
	    
	       //get dR between visible parts of hadronic tau and muonic tau
	       const double muHadDR = reco::deltaR(visibleHadTauEta, visibleHadTauP4.Phi(), genMuEta, muRef->phi());
	    
	       // collect this gen tau_had if it satisfies mu and tau pT cutoffs
	       // and match it to an HPS tau
	       if (genMuPT > genMuPTMin_)
		 { // if it satisfies the cutoffs
		   tausFromTMTH.push_back(tauRef);
		   //cout << "We have found a tau_had whose sister was a tau_mu" << endl;
		   double mindelR_tm = 100.;
		   double delR_tm = 100.;
		   double mindelR_th = 100.;
		   double delR_th = 100.;
		   unsigned int jetIterator = 0;
		   unsigned int jetPointer = 0;
		   bool tmMatched = false;
		   bool thMatched = false;
		
		   // loop over PFJets
		   for(std::vector<reco::PFJet>::const_iterator iJet = pfjetVector.begin(); iJet != pfjetVector.end(); ++iJet)
		     { // loop over jets
		       delR_th = deltaR((*iJet).p4().eta(), (*iJet).p4().phi(), visibleHadTauEta, visibleHadTauP4.Phi());
		       if (delR_th < mindelR_th)
			 {
			   mindelR_th = delR_th;
			   jetPointer = jetIterator;
			 }
		       jetIterator += 1;
		     } // loop over jets
		   if (mindelR_th < 0.3)
		     { // if mindelR_th < 0.3
		       //cout << "Jet matched to this tau_had! Match to tau_mu now." << endl;
		       thMatched = true;
		       delR_tm = deltaR(pfjetVector[jetPointer].p4().eta(), pfjetVector[jetPointer].p4().phi(),muTauEta, muTauPhi);
		       if (delR_tm < 0.3)
			 { // if delR_tm < 0.3
			   //cout << "Jet matched to sister tau_mu too!" << endl;
			   tmMatched = true;
			   //look at pfjet candidates before erasing the jet
			   std::vector<reco::PFCandidatePtr> jetcands = pfjetVector[jetPointer].getPFConstituents();
			   for(std::vector<reco::PFCandidatePtr>::const_iterator cand = jetcands.begin(); cand != jetcands.end(); ++cand)
			     { // loop over jet pfcands
			       if((*cand)->particleId() == PFCandidate::h)
				 { // if it's a charged hadron
				   jeth_ET->Fill((**cand).et());
				 } // if it's a charged hadron
			       else if((*cand)->particleId() == PFCandidate::e)
				 { // if it's an electron
				   jete_ET->Fill((**cand).et());
				 } // if it's an electron
			       else if((*cand)->particleId() == PFCandidate::mu)
				 { // if it's a muon
				   if(genMusatm2.size() != 0)
				     { // if there are gen atm muons
				       double mindelR = 100.;
				       double delR = 100.;
				       unsigned int muonIterator = 0;
				       unsigned int iPointer = 0;
				    
				       for(std::vector<GenParticle*>::iterator iMuon = genMusatm2.begin(); iMuon != genMusatm2.end(); ++iMuon)
					 { // loop over gen muons (a to tau to mu) to find match
					
					   /* Calculate delR for PF muon and gen muon */
					   delR = deltaR((*cand)->eta(), (*cand)->phi(), (*iMuon)->eta(), (*iMuon)->phi());
					
					   /* See if it is less than then minimum*/
					   if (delR < mindelR)
					     {
					       mindelR = delR;
					       iPointer = muonIterator;
					     }
					   muonIterator += 1;
					 } // loop over gen muons to find match
				    
				       if(mindelR < 0.3)
					 {
					   jetMuonatm_ET->Fill((**cand).et());
					   genMusatm2.erase(genMusatm2.begin()+iPointer);
					 }
				       else
					 { jetMuonnotatm_ET->Fill((**cand).et()); }
				    
				     } // if there are gen atm muons
				   else
				     { // if there are no gen atm muons
				       jetMuonnotatm_ET->Fill((**cand).et());
				     } // if there are no gen atm muons
				 } // if it's a muon
			       else if((*cand)->particleId() == PFCandidate::gamma)
				 { // if it's a gamma
				   jetgamma_ET->Fill((**cand).et());
				 } // if it's a gamma
			       else if((*cand)->particleId() == PFCandidate::h0)
				 { // if it's a neutral hadron
				   jeth0_ET->Fill((**cand).et());
				 } // if it's a neutral hadron
			       else
				 { // if it's something else
				   jetOther_ET->Fill((**cand).et());
				 } // if it's something else
			    
			     } // loop over jet pfcands
			   pfjetVector.erase(pfjetVector.begin()+jetPointer);
			 } // if delR_tm < 0.3
		     } // if mindelR_th < 0.3
		   else
		     {
		       //cout << "Jet not matched to this tau_had" << endl;
		     }
		 } // if it satisfies the cutoffs
	    
	     } // if pT of sisterRef satisfies minimum tau_mu pt cutoff
	   
	     //add this tau's key to the list of keys to ignore when we get to its sister so the 
	     //histograms don't get filled twice
	     keysToIgnore.push_back(iTau->getTauIndex());
	   }
	   //else { cout << "Its sister was not a tau_mu" << endl; }
	 } // if hadronic tau decay
       } // matches "try"
       catch (std::string& ex) { throw cms::Exception("TauAnalyzer") << ex; }
     } // end loop over tauDecays
   } // matches "ifcutongenmatches"


   ///////////// NOW DO THE JET-CLEANING //////////////////////////////////////

   unsigned int numgenmuonsatm = genMusatm.size();
   //cout << "Number of atm gen muons in this event = " << numgenmuonsatm << endl;
   int PFCandMuons = 0;
   int PFCandMuonsMatched = 0;
   
   
   for(reco::PFJetCollection::const_iterator iJet = PFJets->begin(); iJet != PFJets->end(); ++iJet)
     { // loop over jets
       
       std::vector<reco::PFCandidatePtr> JetPFCands = iJet->getPFConstituents();
       reco::PFJet::Specific specs = iJet->getSpecific();
       /*       cout << "JET INFO" << endl;
       cout << "jet mMuonEnergy = " << iJet->getSpecific().mMuonEnergy << endl;
       cout << "jet mMuonMultiplicity = " << iJet->getSpecific().mMuonMultiplicity << endl;
       cout << "jet mChargedMuEnergy = " << iJet->getSpecific().mChargedMuEnergy << endl;
       cout << "jet mChargedMultiplicity = " << iJet->getSpecific().mChargedMultiplicity << endl;
       cout << "END JET INFO" << endl;*/
       
       math::XYZTLorentzVector pfmomentum;
       std::vector<edm::Ptr<Candidate> > jetConstituents;
       jetConstituents.clear();

       //flag indicating whether >=0 muons were tagged for removal
       bool taggedMuonForRemoval = false;

       //vector of removed muons for this jet
       MuonRefVector removedMuons;

       for(std::vector<edm::Ptr<reco::PFCandidate> >::iterator i = JetPFCands.begin(); i != JetPFCands.end(); ++i)
	 { // loop over PF candidates
	   reco::PFCandidate pfCand = *i;
	 
	   /* Is the PF Candidate a muon? */
	   if (pfCand.particleId() == 3) 
	     { // if it's a muon
	       
	       // get the ref to the corresponding muon
	       // and count one more PF muon

	       reco::MuonRef theRecoMuon = pfCand.muonRef();
	       PFCandMuons += 1;
	       bool matchedToA = false;
	       // loop over atm gen muons and
	       // look for the smallest delR
	       
	       unsigned int genMuonAtmVectorSize = genMusatm.size();
	       if(genMuonAtmVectorSize!=0)
		 { // if there are atm gen mus in event

// 		   cout << "There are gen mus in this event. Look for a match!" << endl;
		   
		   double mindelR = 100.;
		   double delR = 100.;
		   unsigned int muonIterator = 0;
		   unsigned int iPointer = 0;

		   for(std::vector<GenParticle*>::iterator iMuon = genMusatm.begin(); iMuon != genMusatm.end(); ++iMuon)
		     { // loop over gen muons (a to tau to mu) to find match
		       
		       /* Calculate delR for PF muon and gen muon */
		       delR = deltaR(theRecoMuon->eta(), theRecoMuon->phi(), (*iMuon)->eta(), (*iMuon)->phi());
		      
		       /* See if it is less than then minimum*/
		       if (delR < mindelR)
			 {
			   mindelR = delR;
			   iPointer = muonIterator;
			 }
		       muonIterator += 1;
		     } // loop over gen muons to find match

		   //cout << "The winning delR = " << mindelR << endl;

		   if(mindelR < 0.3)
		     {
		       //cout << "We have found a match!" << endl;
		       matchedToA = true;
		       PFCandMuonsMatched += 1;
		       genMuonatm_matched_pT->Fill((genMusatm[iPointer])->p4().pt());
		       pfMuonatm_pT->Fill(theRecoMuon->p4().pt());
		       genMusatm.erase(genMusatm.begin()+iPointer);
		       genMuonAtmVectorSize -= 1;
		       
		       ////// subtract muon specs ///////
		       
		       specs.mMuonEnergy -= pfCand.p4().e();
		       specs.mMuonMultiplicity -= 1;
		       specs.mChargedMuEnergy -= pfCand.p4().e();
		       specs.mChargedMultiplicity -= 1;
		       
		       /////// and don't keep it in the jet //////
		       
		       //save tag decision for this muon
		       taggedMuonForRemoval = true;
		       
		       //does this muon pass the desired muon ID?
		       std::vector<unsigned int>::const_iterator iSoftMuon = 
			 std::find(muonRefKeys.begin(), muonRefKeys.end(), theRecoMuon.key());
		       
		       
		       /*add this muon ref to the vector of removed muons for this jet
			 iSoftMuon - muonRefKeys.begin() is the index into muonRefKeys of the soft 
			 muon
			 since muonRefKeys was filled in order of muons, it is also the index into 
			 muons of the soft muon*/
		       if (iSoftMuon != muonRefKeys.end()) {
			 removedMuons.push_back(muons->at(iSoftMuon - muonRefKeys.begin()));
		       }
		     }
		   else
		     {
		       //cout << "No match found within delR = 0.3" << endl;
		       //cout << "Since we can't match it to a gen muon from a tau decay from an a decay, we'll keep it in the jet" << endl;
		       pfmomentum += pfCand.p4(); // total p4()
		       jetConstituents.push_back((*i));
		     }
		   
		 } // if there are atm gen mus in event
	       else
		 { // if there are no atm gen mus in the event

		   //cout << "There are no atm gen mus in this event" << endl;
		   // it is a PF muon, but...
		   // since it's not matched to an atm gen mu, keep it

		   //does this muon pass the desired muon ID?
		   std::vector<unsigned int>::const_iterator iSoftMuon = 
		     std::find(muonRefKeys.begin(), muonRefKeys.end(), theRecoMuon.key());
		   
		   if (cutOnGenMatches_)
		     {
		       pfmomentum += pfCand.p4(); // total p4()
		       jetConstituents.push_back((*i));
		     }

		   /*if we're not requiring gen matching but instead looking for muons in jets 
		     that pass the desired muon ID...*/
		   else if (iSoftMuon != muonRefKeys.end())
		     {
		       
		       //remove this muon
		       specs.mMuonEnergy -= pfCand.p4().e();
		       specs.mMuonMultiplicity -= 1;
		       specs.mChargedMuEnergy -= pfCand.p4().e();
		       specs.mChargedMultiplicity -= 1;
		       
		       //save tag decision for this muon
		       taggedMuonForRemoval = true;
		       
		       /*add this muon ref to the vector of removed muons for this jet
			 iSoftMuon - muonRefKeys.begin() is the index into muonRefKeys of the soft 
			 muon
			 since muonRefKeys was filled in order of muons, it is also the index into 
			 muons of the soft muon*/
		       removedMuons.push_back(muons->at(iSoftMuon - muonRefKeys.begin()));
		     }

		   //this muon doesn't pass the soft ID, so keep it in the jet
		   else {
		     pfmomentum += pfCand.p4(); // total p4()
		     jetConstituents.push_back((*i));
		   }
		   
		 } // if there are no atm gen mus in the event

	       if(matchedToA == false)
		 { // if no match to an atm gen muon was found

		   // plot any unmatched atm gen muons left
		   if (genMusatm.size()!=0)
		     { // if there are still atm gen muons
		       for(unsigned int i = 0; i < genMusatm.size(); ++i)
			 {
			   genMuonatm_unmatched_pT->Fill((genMusatm[i])->p4().pt());

			 }
		     } // if there are still atm gen muons

		   // Check if the pf muon can be matched to other muons of interest

		   unsigned int genMuonZtmVectorSize = genMusZtm.size();
		   unsigned int genMuonZmVectorSize = genMusZm.size();
		   bool matchedToZtm = false;
		   bool matchedToZm = false;

		   if(genMuonZtmVectorSize!=0 && matchedToZtm == false)
		     { // if there are Ztm gen muons

		       //cout << "There are Ztm gen mus in this event. Look for a match!" << endl;
		       
		       double mindelR = 100.;
		       double delR = 100.;
		       unsigned int muonIterator = 0;
		       unsigned int iPointer = 0;
		       
		       for(std::vector<GenParticle*>::iterator iMuon = genMusZtm.begin(); iMuon != genMusZtm.end(); ++iMuon)
			 { // loop over gen muons (a to tau to mu) to find match
			   
			   /* Calculate delR for PF muon and gen muon */
			   delR = deltaR(theRecoMuon->eta(), theRecoMuon->phi(), (*iMuon)->eta(), (*iMuon)->phi());
			   
			   /* See if it is less than then minimum*/
			   if (delR < mindelR)
			     {
			       mindelR = delR;
			       iPointer = muonIterator;
			     }
			   muonIterator += 1;
			 } // loop over gen muons to find match
		       		       
		       if(mindelR < 0.3)
			 {
			   //cout << "We have found a match!" << endl;
			   matchedToZtm = true;
			   genMuonZtm_matched_pT->Fill((genMusZtm[iPointer])->p4().pt());
			   genMusZtm.erase(genMusZtm.begin()+iPointer);
			   genMuonZtmVectorSize -= 1;			   
			   pfMuonZtm_pT->Fill(theRecoMuon->p4().pt());
			 }
		     } // if there are Ztm gen muons
		   if(genMuonZmVectorSize!=0 && matchedToZtm == false && matchedToZm == false)
		     { // if there are Zm gen muons
		       
		       //cout << "There are Zm gen mus in this event. Look for a match!" << endl;
		       
		       double mindelR = 100.;
		       double delR = 100.;
		       unsigned int muonIterator = 0;
		       unsigned int iPointer = 0;
		       
		       for(std::vector<GenParticle*>::iterator iMuon = genMusZm.begin(); iMuon != genMusZm.end(); ++iMuon)
			 { // loop over gen muons (a to tau to mu) to find match
			   
			   /* Calculate delR for PF muon and gen muon */
			   delR = deltaR(theRecoMuon->eta(), theRecoMuon->phi(), (*iMuon)->eta(), (*iMuon)->phi());
			   
			   /* See if it is less than then minimum*/
			   if (delR < mindelR)
			     {
			       mindelR = delR;
			       iPointer = muonIterator;
			     }
			   muonIterator += 1;
			 } // loop over gen muons to find match
		       
		       //cout << "The winning delR = " << mindelR << endl;
		       
		       if(mindelR < 0.3)
			 {
			   //cout << "We have found a match!" << endl;
			   matchedToZm = true;
			   genMuonZm_matched_pT->Fill((genMusZm[iPointer])->p4().pt());
			   genMusZm.erase(genMusZm.begin()+iPointer);
			   genMuonZmVectorSize -= 1;
			   pfMuonZm_pT->Fill(theRecoMuon->p4().pt());
			 }
		       
		     } // if there are Zm gen muons
		   if (matchedToZtm == false && matchedToZm == false)
		     { // if there are gen muons, they are Other
		       pfMuonOther_pT->Fill(theRecoMuon->p4().pt());
		     } // if there are gen muons, they are Other

		 } // if no match to an atm gen muon was found
	       
	     } // if it's a muon
	   else // if it's not a muon
	     { // get p4 and constituents
	       
	       pfmomentum += pfCand.p4(); // total p4()
	       jetConstituents.push_back((*i));

	     } //get p4 and constituents

	   
	 } // loop over PF candidates
       
       ////// Build a new jet without the muon /////////////

       PFJet muonfreePFJet(pfmomentum, specs, jetConstituents);
       SetOfJets->push_back( muonfreePFJet );

       //if at least 1 muon was tagged for removal, save a positive muon tag decision for this jet
       muonTagDecisions.push_back(taggedMuonForRemoval);

       //save the ref vector of removed muons
       removedMuonMap.push_back(removedMuons);

       //ref to this (old) jet
       oldJets.push_back(PFJetRef(PFJets, iJet - PFJets->begin()));
      
     } // loop over jets
   
   if(numgenmuonsatm != 0)
     {
       //cout << "Printing the matching fractions:" << endl;
       double ratio1 = (PFCandMuonsMatched*1.)/(1.*numgenmuonsatm);
       //cout << "fracGenMuMatched = " << ratio1 << endl;
       fracGenMuMatched->Fill(ratio1);
       double ratio2 = (PFCandMuonsMatched*1.)/(1.*PFCandMuons);
       fracPFMuMatched->Fill(ratio2);
       if(PFCandMuonsMatched != 0)
	 {
	   //cout << "fracPFMuMatched = " << ratio2 << endl;
	 }
     }

   //fill an STL container of keys of removed muons
   std::vector<unsigned int> removedMuRefKeys;
   for (std::vector<MuonRefVector>::const_iterator iJet = removedMuonMap.begin(); 
	iJet != removedMuonMap.end(); ++iJet) {
     for (MuonRefVector::const_iterator iRemovedMu = iJet->begin(); iRemovedMu != iJet->end(); 
	  ++iRemovedMu) { removedMuRefKeys.push_back(iRemovedMu->key()); }
   }

   /*build a collection of PF candidates excluding soft muons
     we will still tag the jet as signal-like by the presence of a soft muon IN the jet, but this 
     ensures that such jets also cannot have the soft muon enter the isolation candidate 
     collection
     right now only remove muons that are inside jets; later expand to muons within dR = X.X of 
     the jets*/
   for (PFCandidateCollection::const_iterator iPFCand = PFCands->begin(); 
	iPFCand != PFCands->end(); ++iPFCand) {
     MuonRef removedMuRef = iPFCand->muonRef();
     if ((removedMuRef.isNonnull() && 
	  (std::find(removedMuRefKeys.begin(), removedMuRefKeys.end(), removedMuRef.key()) == 
	   removedMuRefKeys.end())) || 
	 removedMuRef.isNull()) PFCandsExcludingSoftMuons->push_back(*iPFCand);
   }

//    //debug
//    for (PFJetCollection::const_iterator iJet = SetOfJets->begin(); iJet != SetOfJets->end(); 
// 	++iJet) {
//      if (muonTagDecisions.at(iJet - SetOfJets->begin()) == 1) {
//        std::cerr << "Selected jet pT: " << iJet->pt() << " GeV\n";
//        std::cerr << "Selected jet eta: " << iJet->eta() << std::endl;
//        std::cerr << "Selected jet phi: " << iJet->phi() << std::endl;
//        std::cerr << "Selected jet ref key: " << iJet - SetOfJets->begin() << std::endl;
//      }
//    }

   const OrphanHandle<PFJetCollection> cleanedJetsRefProd = 
     iEvent.put( SetOfJets, "ak5PFJetsNoMu"  );

   //fill the value map of muon tag decision for each cleaned jet
   std::auto_ptr<edm::ValueMap<bool> > valMap(new edm::ValueMap<bool>());
   edm::ValueMap<bool>::Filler filler(*valMap);
   filler.insert(cleanedJetsRefProd, muonTagDecisions.begin(), muonTagDecisions.end());
   filler.fill();
   iEvent.put(valMap);

   //fill the value map of removed muon refs for each cleaned jet
   std::auto_ptr<edm::ValueMap<MuonRefVector> > muonValMap(new edm::ValueMap<MuonRefVector>());
   edm::ValueMap<MuonRefVector>::Filler muonFiller(*muonValMap);
   muonFiller.insert(cleanedJetsRefProd, removedMuonMap.begin(), removedMuonMap.end());
   muonFiller.fill();
   iEvent.put(muonValMap);

   //fill the value map of old jet refs for each cleaned jet
   std::auto_ptr<edm::ValueMap<PFJetRef> > jetValMap(new edm::ValueMap<PFJetRef>());
   edm::ValueMap<PFJetRef>::Filler jetFiller(*jetValMap);
   jetFiller.insert(cleanedJetsRefProd, oldJets.begin(), oldJets.end());
   jetFiller.fill();
   iEvent.put(jetValMap);

   //put the soft-muon-free PF cands into the event
   iEvent.put(PFCandsExcludingSoftMuons);

}

// ------------ method called once each job just before starting event loop  ------------
void 
CleanJets::beginJob()
{
  out_ = new TFile(outFileName_.c_str(), "RECREATE");
  //  EleStuff.open ("ElectronAncestry.txt");
  //  MuStuff.open("MuAncestry.txt");
  //  PhotonStuff.open("PhotonAncestry.txt");

  fracElePizeroDecays = new TH1D("fracElePizeroDecays", "Fraction of electrons that came from pizero", 11., 0., 1.1);
  fracPFMuMatched = new TH1D("fracPFMuMatched", "Fraction of PFCand muons from a jet that were matched to an atm gen muon", 11.,0., 1.1);
  fracGenMuMatched = new TH1D("fracGenMuMatched", "Fraction of atm gen muons that were matched to a PFCand muon from a jet", 11.,0., 1.1);
  genMuonatm_pT = new TH1D("genMuonatm_pT", "pT of gen muons from a --> tau decays", 100., 0., 200.);
  genMuonZtm_pT = new TH1D("genMuonZtm_pT", "pT of gen muons from Wtm decays", 100., 0., 200.);
  genMuonZm_pT = new TH1D("genMuonZm_pT", "pT of gen muons from Wm decays", 100., 0., 200.);
  genMuonOther_pT = new TH1D("genMuonOther_pT", "pT of gen muons from other sources", 100., 0., 200.);
  genMuonatm_matched_pT = new TH1D("genMuonatm_matched_pT", "pT of gen muons from a --> tau decays matched to PF muons", 100., 0., 200.);
  genMuonatm_unmatched_pT = new TH1D("genMuonatm_unmatched_pT", "pT of gen muons from a --> tau decays NOT matched to PF muons", 100., 0., 200.);
  genMuonZtm_matched_pT = new TH1D("genMuonZtm_matched_pT", "pT of gen muons from Wtm decays matched to PF muons", 100., 0., 200.);
  genMuonZm_matched_pT = new TH1D("genMuonZm_matched_pT", "pT of gen muons from Wm decays matched to PF muons", 100., 0., 200.);
  pfMuonatm_pT = new TH1D("pfMuonatm_pT", "pT of PF muons matched to atm gen muons", 100., 0., 200.);
  pfMuonZtm_pT = new TH1D("pfMuonZtm_pT", "pT of PF muons matched to Wtm gen muons", 100., 0., 200.);
  pfMuonZm_pT = new TH1D("pfMuonZm_pT", "pT of PF muons matched to Wm gen muons", 100., 0., 200.);
  pfMuonOther_pT = new TH1D("pfMuonOther_pT", "pT of PF muons not matched to atm, Ztm, or Zm gen muons", 100., 0., 200.);
  genEleZ_pT = new TH1D("genEleZ_pT", "pT of gen electrons from W decays", 100., 0., 200.);
  genElePiZero_pT = new TH1D("genElePiZero_pT", "pT of gen electrons from neutral pion decays", 100., 0., 200.);
  genElePhoton_pT = new TH1D("genElePhoton_pT", "pT of gen electrons from photons", 100., 0., 200.);
  genEleOther_pT = new TH1D("genEleOther_pT", "pT of gen electrons from other sources", 100., 0., 200.); // e.g.: D*0 -> D0 -> e, rho -> pizero -> e
  jetMuonatm_ET = new TH1D("jetMuonatm_ET", "ET of PFJet muon candidates matched to gen atm muons", 100., 0., 200.);
  jetMuonnotatm_ET = new TH1D("jetMuonnotatm_ET", "ET of PFJet muon candidates not matched to gen atm muons", 100., 0., 200.);
  jeth_ET = new TH1D("jeth_ET", "ET of PFJet charged hadron candidates", 100., 0., 200.);
  jete_ET = new TH1D("jete_ET", "ET of PFJet electron candidates", 100., 0., 200.);
  jetgamma_ET = new TH1D("jetgamma_ET", "ET of PFJet gamma candidates", 100., 0., 200.);
  jeth0_ET = new TH1D("jeth0_ET", "ET of PFJet neutral hadron candidates", 100., 0., 200.);
  jetOther_ET = new TH1D("jetOther_ET", "ET of PFJet Other candidates", 100., 0., 200.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
CleanJets::endJob()
{

  out_->cd();
  fracElePizeroDecays->Write();
  fracPFMuMatched->Write();
  fracGenMuMatched->Write();
  genMuonatm_pT->Write();
  genMuonZtm_pT->Write();
  genMuonZm_pT->Write();
  genMuonOther_pT->Write();
  genMuonatm_matched_pT->Write();
  genMuonatm_unmatched_pT->Write();
  genMuonZtm_matched_pT->Write();
  genMuonZm_matched_pT->Write();
  pfMuonatm_pT->Write();
  pfMuonZtm_pT->Write();
  pfMuonZm_pT->Write();
  pfMuonOther_pT->Write();
  genEleZ_pT->Write();
  genElePiZero_pT->Write();
  genElePhoton_pT->Write();
  genEleOther_pT->Write();
  jeth_ET->Write();
  jete_ET->Write();
  jetMuonatm_ET->Write();
  jetMuonnotatm_ET->Write();
  jetgamma_ET->Write();
  jeth0_ET->Write();
  jetOther_ET->Write();
  out_->Write();
  out_->Close();
  //  EleStuff.close();
  //  MuStuff.close();
  //  PhotonStuff.close();
}

// ------------ method called when starting to processes a run  ------------
void 
CleanJets::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CleanJets::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CleanJets::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CleanJets::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CleanJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CleanJets);
