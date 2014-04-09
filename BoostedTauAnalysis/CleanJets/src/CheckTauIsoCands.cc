// -*- C++ -*-
//
// Package:    CheckTauIsoCands
// Class:      CheckTauIsoCands
// 
/**\class CheckTauIsoCands CheckTauIsoCands.cc BoostedTauAnalysis/CheckTauIsoCands/src/CheckTauIsoCands.cc

 Description: study the properties of signal and isolation
              candidates of HPS PFTau objects

 Implementation:
     Remember to specify whether you are looking
     at the signal candidates or the isolation
     candidates of the HPS PFTaus

*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//         Created:  Mon Sep 10 10:28:20 CEST 2012
// $Id: CheckTauIsoCands.cc,v 1.3 2012/11/08 16:51:51 friccita Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
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
#include "TLegend.h"

using namespace edm;
using namespace reco;
using namespace std;



//
// class declaration
//

class CheckTauIsoCands : public edm::EDAnalyzer {
   public:
      explicit CheckTauIsoCands(const edm::ParameterSet&);
      ~CheckTauIsoCands();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  
  // source of taus to be examined
  edm::InputTag tauSrc_;

  //input,output
  TFile* out_;
  unsigned int momPDGID_;
  double genMuTauPTMin_;
  double genMuPTMin_;
  edm::ParameterSet* cfg_;
  std::string outFileName_;

  //histograms
  TH1D *genEleZ_pT;
  TH1D *genElePizero_pT;
  TH1D *genEleTau_pT;
  TH1D *genEleOther_pT;
  TH1D* pfElesZee_pT;
  TH1D* pfElesPizero_pT;
  TH1D* pfElesTau_pT;
  TH1D* pfElesOther_pT;
  TH1D *tmthDR;
  TH1D *IsoCandFrequency;
  //TH1D *X_pT;
  TH1D *h_pT;
  TH1D *e_pT;
  TH1D *mu_pT;
  TH1D *gamma_pT;
  TH1D *h0_pT;
  //TH1D *hHF_pT;
  //TH1D *egammaHF_pT;
  //TH1D *X_ET;
  TH1D *h_ET;
  TH1D *e_ET;
  TH1D *mu_ET;
  TH1D *gamma_ET;
  TH1D *h0_ET;
  //TH1D *hHF_ET;
  //TH1D *egammaHF_ET;
  //TH1D *X_dxy;
  TH1D *h_dxy;
  TH1D *e_dxy;
  TH1D *mu_dxy;
  //TH1D *gamma_dxy;
  //TH1D *h0_dxy;
  //TH1D *hHF_dxy;
  //TH1D *egammaHF_dxy;
  //TH1D *X_dz;
  TH1D *h_dz;
  TH1D *e_dz;
  TH1D *mu_dz;
  TH1D *gamma_dz;
  //TH1D *h0_dz;
  //TH1D *hHF_dz;
  //TH1D *egammaHF_dz;
  TH1D *pfMu_atmMatch_pT;
  TH1D *pfMu_ZtmMatch_pT;
  TH1D *pfMu_ZmMatch_pT;
  TH1D *pfMu_OtherSource_pT;
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
CheckTauIsoCands::CheckTauIsoCands(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  tauSrc_ = iConfig.getParameter<edm::InputTag>("tauSrc");
  outFileName_ = iConfig.getParameter<std::string>("outFileName");
  momPDGID_ = iConfig.getParameter<unsigned int>("momPDGID");
  genMuTauPTMin_ = iConfig.getParameter<double>("genMuTauPTMin");
  genMuPTMin_ = iConfig.getParameter<double>("genMuPTMin");
  cfg_ = const_cast<edm::ParameterSet*>(&iConfig);
}


CheckTauIsoCands::~CheckTauIsoCands()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CheckTauIsoCands::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  /* Get gen particle collection */
   Handle<GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles", genParticles);

  /* Get tau collection */
  cout << "Getting taus" << endl;
  Handle<PFTauCollection> HPSPFTaus;
  iEvent.getByLabel(tauSrc_, HPSPFTaus); // tauSrc_ = hpsPFTauProducerNoMuons
  cout << "Got taus" << endl;

  // Loop over genparticles to get mus

   std::vector<reco::GenParticle*> genMusatm; //gen muons from a -> tau -> mu (atm)
   std::vector<reco::GenParticle*> genMusZtm; //gen muons from Z -> tau -> mu (Ztm)
   std::vector<reco::GenParticle*> genMusZm; //gen muons from Z -> mu (Zm)
   std::vector<reco::GenParticle*> genMusOther; //gen muons from other sources
   std::vector<reco::GenParticle*> genElesZee; // gen e's from Z -> ee
   std::vector<reco::GenParticle*> genElesPizero; // gen e's from pizero -> e's + gammas
   std::vector<reco::GenParticle*> genElesTau; // gen e's from Z -> ee
   std::vector<GenTauDecayID> tauDecays; // vector of taus from a decays

   // make vector of gen mus that come from tau decays
   for (GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen)
     { // loop over genparticles to get mus
       if ((fabs(iGen->pdgId()) == 13) && (fabs(iGen->mother()->pdgId()) != 13))
	 { // if it's a mu
	   cout << "Status of mu = " << iGen->status() << endl;
	   cout << "Mother pdgid = " << iGen->mother()->pdgId() << endl;
	   if(fabs(iGen->mother()->pdgId()) == 15)
	     { // if its mother was a tau
	       if ((fabs(iGen->mother()->mother()->mother()->pdgId() == 36))&&(iGen->mother()->mother()->mother()->status() == 3))
		 { // if its grandmother was a status-3 a
		   genMusatm.push_back(const_cast<reco::GenParticle*>(&*iGen));
		  
		 } // if its grandmother was a status-3 a
	       else if (iGen->status() == 1 && (fabs(iGen->mother()->mother()->mother()->pdgId()) == 24))
		 { // if its grandmother was a Z/W
		   genMusZtm.push_back(const_cast<reco::GenParticle*>(&*iGen));
		 } // if its grandmother was a Z/W
	     } // if its mother was a tau
	   else if(fabs(iGen->mother()->pdgId()) == 24 && iGen->status()==3)
	     { // if its mother was a Z/W
	       genMusZm.push_back(const_cast<reco::GenParticle*>(&*iGen));
	     } // if its mother was a Z/W
	   else
	     { // if its mother was not a tau or a Z/W
		 genMusOther.push_back(const_cast<reco::GenParticle*>(&*iGen));
	     } // if its mother was not a tau or a Z/W
	 } // if it's a mu
       if ((fabs(iGen->pdgId()) == 11)) 
	 { // if it's an electron
	   if(fabs(iGen->mother()->pdgId())!= 11)
	     { // if not from an electron
	       if(fabs(iGen->mother()->pdgId()) == 24)
		 { // e from Z/W decay
		   genElesZee.push_back(const_cast<reco::GenParticle*>(&*iGen));
		   genEleZ_pT->Fill(iGen->pt());
		 } // e from Z/W decay
	       else if(fabs(iGen->mother()->pdgId()) == 111)
		 { // e from pizero decay
		   genElesPizero.push_back(const_cast<reco::GenParticle*>(&*iGen));
		   genElePizero_pT->Fill(iGen->pt());
		 } // e from pizero decay
	       else if(fabs(iGen->mother()->pdgId()) == 15)
		 { // e from tau decay
		   genElesTau.push_back(const_cast<reco::GenParticle*>(&*iGen));
		   genEleTau_pT->Fill(iGen->pt());
		 } // e from tau decay
	       else
		 {
		   genEleOther_pT->Fill(iGen->pt());
		 }
	     } // if not from an electron
	 }// if it's an electron
       
     } // loop over genparticles to get mus

   ///////////////////////////////////////////////////////////////////////////
   
   // Collect hadronically decaying taus

  for (GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen)
     { // loop over genparticles to get tau decays
       
       if((fabs(iGen->pdgId()) == 15) && (iGen->status() == 3))
	 { // if it is a status-3 tau
	   if((fabs(iGen->mother()->pdgId()) == 36))
	     { // if its mother was an a
	       cout << "A status 3 tau from an a has been found" << endl;
	       cout << "status of a = " << iGen->mother()->status() << endl;
	       try
		     { // try loop
		       GenTauDecayID tauDecay(*cfg_, genParticles, iGen - genParticles->begin());
		       cout << "Is this tau a status 3 decay product? Boolean answer: " << tauDecay.tauIsStatus3DecayProduct() << endl;
		       tauDecays.push_back(tauDecay);
		       cout << "We have collected this tauDecay for a status-3 tau descended from an a" << endl;
		     } // try loop
	       catch (std::string& ex) { throw cms::Exception("CheckTauIsoCands") << ex; }
	     } // if its mother was an a
	 } // if it is a status-3 tau
       
     } // loop over genparticles to get tau decays



   GenParticleRefVector tausFromTMTH; // hadronically decaying taus whose sisters decayed to muons
   PFTauRefVector matchedHPSPFTaus; // HPS PFTaus matched to the above

     //loop over taus from boson decay
   std::vector<unsigned int> keysToIgnore;
   for (std::vector<GenTauDecayID>::iterator iTau = tauDecays.begin();iTau != tauDecays.end();++iTau) {
     cout << "We have entered the tau decay loop" << endl;

     //try for exceptions
     try {
       
       //look for a hadronic tau decay from the signal sample
       if ((iTau->tauDecayType() == GenTauDecayID::HAD)) {
	 cout << "Hadronically decaying tau!" << endl;
	 //look for the other tau decay product of the a: is it a tau-->mu decay?
	 iTau->findSister();
	 const unsigned int iSister = iTau->getSisterIndex();
	 if ((std::find(keysToIgnore.begin(), keysToIgnore.end(), iSister) == 
	      keysToIgnore.end()) && (iTau->sisterDecayType() == GenTauDecayID::MU)) {
	   cout << "Its sister decayed to a muon!" << endl;
	   //does the sister decay mu pass the pT threshold?
	   reco::GenParticleRef sisterRef(genParticles, iSister);
	   if (sisterRef->pt() > genMuTauPTMin_) {
	     cout << "sister tau_mu passes threshold" << endl;
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
// 	     const double visibleHadTauPT = visibleHadTauP4.Pt();
	     const double visibleHadTauEta = visibleHadTauP4.Eta();
	     
	     //get the muon 4-vector
	     int iGenMu = -1;
	     unsigned int iDaughter = 0;
	    reco::GenParticleRef status2MuTauRef = sisterRef->daughterRef(0);
	    while ((iDaughter < status2MuTauRef->numberOfDaughters()) && (iGenMu == -1)) {
	      reco::GenParticleRef muTauDaughter = status2MuTauRef->daughterRef(iDaughter);
	      if (fabs(muTauDaughter->pdgId()) == GenTauDecayID::MUPDGID) {
		iGenMu = muTauDaughter.key();
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
		tmthDR->Fill(muHadDR);
		cout << "We have found a tau_had whose sister was a tau_mu" << endl;
	      } // if it satisfies the cutoffs
	    
	   }
	   
	   //add this tau's key to the list of keys to ignore when we get to its sister so the 
	   //histograms don't get filled twice
	   keysToIgnore.push_back(iTau->getTauIndex());
	 }
	 else { cout << "Its sister was not a tau_mu" << endl; }
       }
     }
     catch (std::string& ex) { throw cms::Exception("TauAnalyzer") << ex; }
   }

   for(PFTauCollection::const_iterator pfTau = HPSPFTaus->begin(); pfTau != HPSPFTaus->end(); ++pfTau)
     { // loop over HPS PFTaus
       
       double mindelR = 100.;
       double delR = 100.;
       unsigned int tauIterator = 0;
       unsigned int iPointer = 0;
       
       for (GenParticleRefVector::const_iterator gent = tausFromTMTH.begin(); gent != tausFromTMTH.end(); ++gent)
	 { // loop over gen taus
	   reco::GenParticleRef status2HadTauRef = (*gent)->daughterRef(0);
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
	   delR = deltaR(pfTau->p4().eta(), pfTau->p4().phi(), visibleHadTauP4.eta(), visibleHadTauP4.phi());
	   //See if it is less than then minimum
	   if (delR < mindelR)
	     {
	       mindelR = delR;
	       iPointer = tauIterator;
	     }
	   tauIterator += 1;
	 } // loop over gen taus
       if(mindelR < 0.3)
	 {
	   cout << "We have matched a pftau to an atm gen tau" << endl;
	   matchedHPSPFTaus.push_back(PFTauRef(HPSPFTaus, pfTau - HPSPFTaus->begin()));
	   tausFromTMTH.erase(tausFromTMTH.begin()+iPointer);
	 }
       
     } // loop over HPS PFTaus
   
   
/////////////////////////////////////////////////////////////////////////


  /* Loop through matched taus to study sig/iso candidate properties */

   for(PFTauCollection::const_iterator iTau = HPSPFTaus->begin(); iTau != HPSPFTaus->end(); ++iTau)
     //for(PFTauRefVector::const_iterator iTau = matchedHPSPFTaus.begin(); iTau != matchedHPSPFTaus.end(); ++iTau) // after you've matched to TMTH
    { // loop over matched HPS PFTaus
      cout << "Looping over taus" << endl;

      /* Get isolation candidates */
      PFCandidateRefVector IsoCands = iTau->isolationPFCands();
      //PFCandidateRefVector IsoCands = (*iTau)->isolationPFCands(); // after you've matched to TMTH

      for(PFCandidateRefVector::const_iterator cand = IsoCands.begin(); cand != IsoCands.end(); ++cand)
	{
	  unsigned int Id = (**cand).particleId();
	  IsoCandFrequency->Fill(Id);

	  // Plot the kinematic parameters for various particle types
	  if(Id == PFCandidate::h)
	    { // if it's an h
	      cout << "Charged hadron" << endl;
	      h_pT->Fill((*((**cand).trackRef())).pt());
	      h_ET->Fill((**cand).et());
	      h_dxy->Fill((*((**cand).trackRef())).dxy());
	      h_dz->Fill((*((**cand).trackRef())).dz());
	    } // if it's an h
	  else if(Id == PFCandidate::e)
	    { // if it's an e
	      cout << "Electron" << endl;
	      if((**cand).gsfTrackRef().isNonnull())
		{
		  e_pT->Fill((*((**cand).gsfTrackRef())).pt());
		  e_ET->Fill((**cand).et());
		  e_dxy->Fill((*((**cand).gsfTrackRef())).dxy());
		  e_dz->Fill((*((**cand).gsfTrackRef())).dz());
		}
	      else
		{ cout << "there are no gsftracks, strangely" << endl; }
	      /// PERFORM ELECTRON MATCHING ///

	      reco::GsfTrackRef theRecoEle = (**cand).gsfTrackRef();
	      bool matchedToZ = false;
	      bool matchedToPizero = false;
	      bool matchedToTau = false;
	      // loop over atm gen electrons and
	      // look for the smallest delR
	      
	      if(genElesZee.size()!=0)
		{ // if there are Z->ee electrons, loop over them
		  double mindelR = 100.;
		  double delR = 100.;
		  unsigned int eleIterator = 0;
		  unsigned int iPointer = 0;
		  
		  for(std::vector<GenParticle*>::iterator iEle = genElesZee.begin(); iEle != genElesZee.end(); ++iEle)
		    { // loop over Zees to find match
		      delR = deltaR(theRecoEle->eta(), theRecoEle->phi(), (*iEle)->eta(), (*iEle)->phi());
		      if (delR < mindelR)
			{
			  mindelR = delR;
			  iPointer = eleIterator;
			}
		      eleIterator += 1;
		    } // loop over Zees to find match
		  cout << "The winning delR = " << mindelR << endl;
		  if (mindelR < 0.3)
		    {
		      cout << "Zee match found!" << endl;
		      matchedToZ = true;
		      pfElesZee_pT->Fill(theRecoEle->pt());
		      genElesZee.erase(genElesZee.begin()+iPointer);
		    }
		  else
		    { cout << "No Zee match found" << endl; }
		} // if there are Z->ee electrons, loop over them

	      if(genElesPizero.size()!=0 && matchedToZ == false)
		{ // if there are electrons from pizero decays, loop over them
		  double mindelR = 100.;
		  double delR = 100.;
		  unsigned int eleIterator = 0;
		  unsigned int iPointer = 0;
		  
		  for(std::vector<GenParticle*>::iterator iEle = genElesPizero.begin(); iEle != genElesPizero.end(); ++iEle)
		    { // loop over pizero e's to find match
		      delR = deltaR(theRecoEle->eta(), theRecoEle->phi(), (*iEle)->eta(), (*iEle)->phi());
		      if (delR < mindelR)
			{
			  mindelR = delR;
			  iPointer = eleIterator;
			}
		      eleIterator += 1;
		    } // loop over pizero e's to find match
		  cout << "The winning delR = " << mindelR << endl;
		  if (mindelR < 0.3)
		    {
		      cout << "pizero match found!" << endl;
		      matchedToPizero = true;
		      pfElesPizero_pT->Fill(theRecoEle->pt());
		      genElesPizero.erase(genElesPizero.begin()+iPointer);
		    }
		  else
		    { cout << "No pizero match found" << endl; }
		} // if there are electrons from pizero decays, loop over them
	      
	      if(genElesTau.size()!=0 && matchedToZ == false && matchedToPizero == false)
		{ // if there are electrons from tau decays, loop over them
		  double mindelR = 100.;
		  double delR = 100.;
		  unsigned int eleIterator = 0;
		  unsigned int iPointer = 0;
		  unsigned int pos = 0;
		  unsigned int counter = 0;
		  cout << "number of gen tau_e's = " << genElesTau.size() << endl;
		  for(std::vector<GenParticle*>::iterator iEle = genElesTau.begin(); iEle != genElesTau.end(); ++iEle)
		    { // loop over tau e's to find match
		      delR = deltaR(theRecoEle->eta(), theRecoEle->phi(), (*iEle)->eta(), (*iEle)->phi());
		      if (delR < mindelR)
			{
			  mindelR = delR;
			  iPointer = eleIterator;
			  pos = counter;
			}
		      eleIterator += 1;
		      counter += 1;
		    } // loop over tau e's to find match
		  cout << "The winning delR = " << mindelR << endl;
		  if (mindelR < 0.3)
		    {
		      cout << "tau match found!" << endl;
		      matchedToTau = true;
		      cout << "IMPORTANT NOTE: the mother of the tau is " << genElesTau[pos]->mother()->mother()->mother()->pdgId() << endl;
		      pfElesTau_pT->Fill(theRecoEle->pt());
		      genElesTau.erase(genElesTau.begin()+iPointer);
		    }
		  else
		    { cout << "No tau match found" << endl; }
		} // if there are electrons from tau decays, loop over them

	      if((matchedToZ == false) && (matchedToPizero == false) && (matchedToTau == false))
		{
		  pfElesOther_pT->Fill(theRecoEle->pt());
		}
	    } // if it's an e
	  else if(Id == PFCandidate::mu)
	    { // if it's a mu
	      cout << "Muon" << endl;
	      mu_pT->Fill((*((**cand).trackRef())).pt());
	      mu_ET->Fill((**cand).et());
	      mu_dxy->Fill((*((**cand).trackRef())).dxy());
	      mu_dz->Fill((*((**cand).trackRef())).dz());
	      /* Matching to gen mu */
	      
	      // get the ref to the corresponding muon
	      // and count one more PF muon
	      
	      reco::MuonRef theRecoMuon = (**cand).muonRef();
	      bool matchedToA = false;
	      // loop over atm gen muons and
	      // look for the smallest delR
	      
	      unsigned int genMuonAtmVectorSize = genMusatm.size();
	      if(genMuonAtmVectorSize!=0)
		{ // if there are atm gen mus in event
		  
		  cout << "There are gen mus in this event. Look for a match!" << endl;
		  
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
		  
		  cout << "The winning delR = " << mindelR << endl;
		  
		  if(mindelR < 0.3)
		    {
		      cout << "We have found a match!" << endl;
		      matchedToA = true;
		      genMusatm.erase(genMusatm.begin()+iPointer);
		      genMuonAtmVectorSize -= 1;
		      pfMu_atmMatch_pT->Fill(theRecoMuon->p4().pt());
		    }
		  else
		    {
		      cout << "No match found within delR = 0.3" << endl;
		    }
		} // if there are atm gen mus in event
	      else
		{ // if there are no atm gen mus in the event
		  cout << "There are no atm gen mus in this event" << endl;
		} // if there are no atm gen mus in the event
	      
	      if(matchedToA == false)
		{ // if no match to an atm gen muon was found
		  
		  cout << "We couldn't find a match to an atm gen mu," << endl;
		  cout << "so let's see if we can match this PF mu to something else." << endl;
		  
		  // Check if it can be matched to other muons of interest
		  
		  unsigned int genMuonZtmVectorSize = genMusZtm.size();
		  unsigned int genMuonZmVectorSize = genMusZm.size();
		  bool matchedToZtm = false;
		  bool matchedToZm = false;
		  
		  if(genMuonZtmVectorSize!=0 && matchedToZtm == false)
		    { // if there are Ztm gen muons
		      
		      cout << "There are Ztm gen mus in this event. Look for a match!" << endl;
		      
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
			  cout << "We have found a match!" << endl;
			  matchedToZtm = true;
			  genMusZtm.erase(genMusZtm.begin()+iPointer);
			  genMuonZtmVectorSize -= 1;
			  pfMu_ZtmMatch_pT->Fill(theRecoMuon->p4().pt());			   
			}
		    } // if there are Ztm gen muons
		  if(genMuonZmVectorSize!=0 && matchedToZtm == false && matchedToZm == false)
		    { // if there are Zm gen muons
		      
		      cout << "There are Zm gen mus in this event. Look for a match!" << endl;
		      
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
		      
		      cout << "The winning delR = " << mindelR << endl;
		      
		      if(mindelR < 0.3)
			{
			  cout << "We have found a match!" << endl;
			  matchedToZm = true;
			  genMusZm.erase(genMusZm.begin()+iPointer);
			  genMuonZmVectorSize -= 1;
			  pfMu_ZmMatch_pT->Fill(theRecoMuon->p4().pt());
			}
		      
		    } // if there are Zm gen muons
		  if (matchedToZtm == false && matchedToZm == false)
		    { // if there are gen muons, they are Other
		      pfMu_OtherSource_pT->Fill(theRecoMuon->p4().pt());
		    } // if there are gen muons, they are Other
		  
		} // if no match to an atm gen muon was found
	      
	      /* Matching to gen mu */
	    } // if it's a mu
	  else if(Id == PFCandidate::gamma)
	    { // if it's a gamma
	      cout << "Photon" << endl;
	      gamma_ET->Fill((**cand).et());
	    } // if it's a gamma
	  else if(Id == PFCandidate::h0)
	    { // if it's an h0
	      cout << "Neutral hadron" << endl;
	      h0_ET->Fill((**cand).et());
	    } // if it's an h0
	  else
	    { // if it's none of the above
	      cout << "Something completely different!" << endl;
	    } // if it's none of the above
	    
	    /*case PFCandidate::X:
	    cout << "Unknown" << endl;
	    // X_pT->Fill((*((**cand).trackRef())).pt());
	    X_ET->Fill((**cand).et());
	    // X_dxy->Fill((*((**cand).trackRef())).dxy());
	    //X_dz->Fill((*((**cand).trackRef())).dz());
	    break;*/
	 
	}

    } // loop over matched HPS PFTaus
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
CheckTauIsoCands::beginJob()
{
  cout << "Beginning job" << endl;
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  IsoCandFrequency = new TH1D("IsoCandFrequency", "Frequency of particle IDs among the isolation PF Candidates of the taus", 10., 0., 10.);
  tmthDR = new TH1D("tmthDR", "deltaR of gen tau_mu and gen tau_had from a decay", 50., 0., 5.);
  genEleZ_pT = new TH1D("genEleZ_pT", "pT of gen e from W->e", 100., 0., 200.);
  genElePizero_pT = new TH1D("genElePizero_pT", "pT of gen e directly from pizero decay", 100., 0., 200.);
  genEleTau_pT = new TH1D("genEleTau_pT", "pT of gen e from tau decay", 100., 0., 200.);
  genEleOther_pT = new TH1D("genEleOther_pT", "pT of gen e from other sources", 100., 0., 200.);
  pfElesZee_pT = new TH1D("pfElesZee_pT", "pT of PF e matched to gen e from W->e", 100., 0., 200.);
  pfElesPizero_pT = new TH1D("pfElesPizero_pT", "pT of PF e matched to gen e from a pizero decay", 100., 0., 200.);
  pfElesTau_pT = new TH1D("pfElesTau_pT", "pT of PF e matched to gen e from tau decay", 100., 0., 200.);  
  pfElesOther_pT = new TH1D("pfElesOther_pT", "pT of PF e from other sources", 100., 0., 200.);
  //X_pT = new TH1D("X_pT", "pT of X PFCands", 100., 0., 200.);
  h_pT = new TH1D("h_pT", "pT of h PFCands", 100., 0., 200.);
  e_pT = new TH1D("e_pT", "pT of e PFCands", 100., 0., 200.);
  mu_pT = new TH1D("mu_pT", "pT of mu PFCands", 100., 0., 200.);
  gamma_pT = new TH1D("gamma_pT", "pT of gamma PFCands", 100., 0., 200.);
  h0_pT = new TH1D("h0_pT", "pT of h0 PFCands", 100., 0., 200.);
  //hHF_pT = new TH1D("hHF_pT", "pT of hHF PFCands", 100., 0., 200.);
  //egammaHF_pT = new TH1D("egammaHF_pT", "pT of egammaHF PFCands", 100., 0., 200.);
  //X_ET = new TH1D("X_ET", "ET of X PFCands", 100., 0., 200.);
  h_ET = new TH1D("h_ET", "ET of h PFCands", 100., 0., 200.);
  e_ET = new TH1D("e_ET", "ET of e PFCands", 100., 0., 200.);
  mu_ET = new TH1D("mu_ET", "ET of mu PFCands", 100., 0., 200.);
  gamma_ET = new TH1D("gamma_ET", "ET of gamma PFCands", 100., 0., 200.);
  h0_ET = new TH1D("h0_ET", "ET of h0 PFCands", 100., 0., 200.);
  //hHF_ET = new TH1D("hHF_ET", "ET of hHF PFCands", 100., 0., 200.);
  //egammaHF_ET = new TH1D("egammaHF_ET", "ET of egammaHF PFCands", 100., 0., 200.);
  //X_dxy = new TH1D("X_dxy", "dxy of X PFCands", 10., 0., .05);
  h_dxy = new TH1D("h_dxy", "dxy of h PFCands", 10., 0., .05);
  e_dxy = new TH1D("e_dxy", "dxy of e PFCands", 10., 0., .05);
  mu_dxy = new TH1D("mu_dxy", "dxy of mu PFCands", 10., 0., .05);
  //gamma_dxy = new TH1D("gamma_dxy", "dxy of gamma PFCands", 10., 0., .05);
  //h0_dxy = new TH1D("h0_dxy", "dxy of h0 PFCands", 10., 0., .05);
  //hHF_dxy = new TH1D("hHF_dxy", "dxy of hHF PFCands", 10., 0., .05);
  //egammaHF_dxy = new TH1D("egammaHF_dxy", "dxy of egammaHF PFCands", 10., 0., .05);
  //X_dz = new TH1D("X_dz", "dz of X PFCands", 10., 0., .5);
  h_dz = new TH1D("h_dz", "dz of h PFCands", 10., 0., .5);
  e_dz = new TH1D("e_dz", "dz of e PFCands", 10., 0., .5);
  mu_dz = new TH1D("mu_dz", "dz of mu PFCands", 10., 0., .5);
  //gamma_dz = new TH1D("gamma_dz", "dz of gamma PFCands", 10., 0., .5);
  //h0_dz = new TH1D("h0_dz", "dz of h0 PFCands", 10., 0., .5);
  //hHF_dz = new TH1D("hHF_dz", "dz of hHF PFCands", 10., 0., .5);
  //egammaHF_dz = new TH1D("egammaHF_dz", "dz of egammaHF PFCands", 10., 0., .5);
  pfMu_atmMatch_pT = new TH1D("pfMu_atmMatch_pT", "pT of PF mus matched to atm gen mus", 100., 0., 200.);
  pfMu_ZtmMatch_pT = new TH1D("pfMu_ZtmMatch_pT", "pT of PF mus matched to Ztm gen mus", 100., 0., 200.);
  pfMu_ZmMatch_pT = new TH1D("pfMu_ZmMatch_pT", "pT of PF mus matched to Zm gen mus", 100., 0., 200.);
  pfMu_OtherSource_pT = new TH1D("pfMu_OtherSource_pT", "pT of PF mus from other sources", 100., 0., 200.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
CheckTauIsoCands::endJob() 
{
  out_->cd();
  tmthDR->Write();
  IsoCandFrequency->Write();
  genEleZ_pT->Write();
  genElePizero_pT->Write();
  genEleTau_pT->Write();
  genEleOther_pT->Write();
  pfElesZee_pT->Write();
  pfElesPizero_pT->Write();
  pfElesTau_pT->Write();
  pfElesOther_pT->Write();
  TCanvas *candpT = new TCanvas("candpT", "pT's of different kinds of isolation PF candidates", 750, 750);
  candpT->cd();
  //X_pT->Write();
  if(h_pT->Integral()!=0.)
    h_pT->Scale(1./h_pT->Integral());
  h_pT->Draw();
  if(e_pT->Integral()!=0.)
     e_pT->Scale(1./e_pT->Integral());
  e_pT->SetLineColor(4);
   e_pT->Draw("same");
  if(mu_pT->Integral()!=0.)
    mu_pT->Scale(1./mu_pT->Integral());
  mu_pT->SetLineColor(2);
  mu_pT->Draw("same");
  if(gamma_pT->Integral()!=0.)
    gamma_pT->Scale(1./gamma_pT->Integral());
  gamma_pT->SetLineColor(3);
  gamma_pT->Draw("same");
  if(h0_pT->Integral()!=0.)
    h0_pT->Scale(1./h0_pT->Integral());
  h0_pT->SetLineColor(6);
  h0_pT->Draw("same");
  TLegend *pTLegend = new TLegend(0.4,0.6,0.89,0.89);
  pTLegend->AddEntry(h_pT, "h");
  pTLegend->AddEntry(e_pT, "e");
  pTLegend->AddEntry(mu_pT, "mu");
  pTLegend->AddEntry(gamma_pT, "gamma");
  pTLegend->AddEntry(h0_pT, "h0");
  pTLegend->Draw("same");
  //hHF_pT->Write();
  //egammaHF_pT->Write();
  //X_ET->Write();
  out_->cd();
  candpT->Write();
  TCanvas *candET = new TCanvas("candET", "ET's of different kinds of isolation PF candidates", 750, 750);
  candET->cd();
  if(h_ET->Integral()!=0.)
    h_ET->Scale(1./h_ET->Integral());
  h_ET->Draw();
  if(e_ET->Integral()!=0.)
    e_ET->Scale(1./e_ET->Integral());
  e_ET->SetLineColor(4);
  e_ET->Draw("same");
  if(mu_ET->Integral()!=0.)
    mu_ET->Scale(1./mu_ET->Integral());
  mu_ET->SetLineColor(2);
  mu_ET->Draw("same");
  if(gamma_ET->Integral()!=0.)
    gamma_ET->Scale(1./gamma_ET->Integral());
  e_ET->SetLineColor(4);
  gamma_ET->SetLineColor(3);
  gamma_ET->Draw("same");
  if(h0_ET->Integral()!=0.)
    h0_ET->Scale(1./h0_ET->Integral());
  h0_ET->SetLineColor(6);
  h0_ET->Draw("same");
  TLegend *ETLegend = new TLegend(0.4,0.6,0.89,0.89);
  ETLegend->AddEntry(h_pT, "h");
  ETLegend->AddEntry(e_pT, "e");
  ETLegend->AddEntry(mu_pT, "mu");
  ETLegend->AddEntry(gamma_pT, "gamma");
  ETLegend->AddEntry(h0_pT, "h0");
  ETLegend->Draw("same"); //hHF_ET->Write();
  //egammaHF_ET->Write();
  //X_dxy->Write();
  out_->cd();
  candET->Write();
  h_dxy->Write();
  e_dxy->Write();
  mu_dxy->Write();
  //gamma_dxy->Write();
  //h0_dxy->Write();
  //hHF_dxy->Write();
  //egammaHF_dxy->Write();
  //X_dz->Write();
  h_dz->Write();
  e_dz->Write();
  mu_dz->Write();
  //gamma_dz->Write();
  //h0_dz->Write();
  //hHF_dz->Write();
  //egammaHF_dz->Write();
  TCanvas *pfMu_pT = new TCanvas("pfMu_pT", "pT distributions of PF muon isolation candidates matched to gen muons", 750, 750);
  pfMu_pT->cd();
if(pfMu_atmMatch_pT->Integral()!=0.)
    pfMu_atmMatch_pT->Scale(1./pfMu_atmMatch_pT->Integral());
  pfMu_atmMatch_pT->Draw();
  if(pfMu_ZtmMatch_pT->Integral()!=0.)
    pfMu_ZtmMatch_pT->Scale(1./pfMu_ZtmMatch_pT->Integral());
  pfMu_ZtmMatch_pT->SetLineColor(2);
  pfMu_ZtmMatch_pT->Draw("same");
  if(pfMu_ZmMatch_pT->Integral()!=0.)
    pfMu_ZmMatch_pT->Scale(1./pfMu_ZmMatch_pT->Integral());
  pfMu_ZmMatch_pT->SetLineColor(3);
  pfMu_ZmMatch_pT->Draw("same");
  if(pfMu_OtherSource_pT->Integral()!=0.)
    pfMu_OtherSource_pT->Scale(1./pfMu_OtherSource_pT->Integral());
  pfMu_OtherSource_pT->SetLineColor(4);
  pfMu_OtherSource_pT->Draw("same");
  TLegend *pfpTLegend = new TLegend(0.4,0.6,0.89,0.89);
  pfpTLegend->AddEntry(pfMu_atmMatch_pT, "pfMu_atmMatch");
  pfpTLegend->AddEntry(pfMu_ZtmMatch_pT, "pfMu_ZtmMatch");
  pfpTLegend->AddEntry(pfMu_ZmMatch_pT, "pfMu_ZmMatch");
  pfpTLegend->AddEntry(pfMu_OtherSource_pT, "pfMu_OtherSource");
  pfpTLegend->Draw("same");
  out_->cd();
  pfMu_pT->Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
CheckTauIsoCands::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CheckTauIsoCands::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CheckTauIsoCands::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CheckTauIsoCands::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CheckTauIsoCands::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CheckTauIsoCands);
