#ifndef BoostedTauAnalysis_Common_interface_GenTauDecayID_h
#define BoostedTauAnalysis_Common_interface_GenTauDecayID_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TauReco/interface/PFTau.h"

class GenTauDecayID {

 public:

  //PDG IDs
  enum PDGID {        EPDGID = 11,         MUPDGID = 13,         TAUPDGID = 15, 
	      ENEUTRINOPDGID = 12, MUNEUTRINOPDGID = 14, TAUNEUTRINOPDGID = 16, 
		      ZPDGID = 23,          APDGID = 36,       GAMMAPDGID = 22, 
		      TPDGID = 6,        ANY_PDGID = 0};

  //decay types
  enum DecayType {HAD = 0, MU, E, UNKNOWN, NOT_TAU};

  //pT rank
  enum PTRank {ANY_PT_RANK = -1};

  //parton PDG ID codes
  enum PartonPDGID {D = 1, U = 2, S = 3, C = 4, B = 5, T = 6, G = 21};

  //default constructor
  GenTauDecayID();

  //constructor that takes a parameter set of arguments
  GenTauDecayID(edm::ParameterSet& parSet, 
		const edm::Handle<reco::GenParticleCollection>& pGenParticles, 
		const unsigned int iTau, const bool warn = true);

  //copy constructor
  GenTauDecayID(const GenTauDecayID&);

  //destructor
  ~GenTauDecayID();

  //assignment operator
  GenTauDecayID& operator=(const GenTauDecayID&);

  //getter for parameter set
  edm::ParameterSet* getParSet() const;

  //getter for gen particle handle
  edm::Handle<reco::GenParticleCollection> getGenParticleHandle() const;

  //getter for tau index
  unsigned int getTauIndex() const;

  //getter for tau sister index; warning suppressed for copy constructor and assignment operator
  unsigned int getSisterIndex(const bool warn = true) const;

  //getter for visible tau 4-vector
  reco::LeafCandidate::LorentzVector getVisibleTauP4(const bool warn = true) const;

  //getter for visible tau sister 4-vector
  reco::LeafCandidate::LorentzVector getVisibleTauSisterP4(const bool warn = true) const;

  //getter for overall tau decay type
  DecayType getTauDecayType(const bool warn = true) const;

  //getter for overall tau sister decay type
  DecayType getTauSisterDecayType(const bool warn = true) const;

  //getter for hadronic tau decay type
  reco::PFTau::hadronicDecayMode getTauHadronicDecayType(const bool warn = true) const;

  //getter for hadronic tau sister decay type
  reco::PFTau::hadronicDecayMode getTauSisterHadronicDecayType(const bool warn = true) const;

  //getter for charged hadron min pT
  double getChargedHadronPTMin() const;

  //getter for neutral hadron min pT
  double getNeutralHadronPTMin() const;

  //getter for charged lepton min pT
  double getChargedLeptonPTMin() const;

  //getter for total min pT of visible decay products
  double getTotalPTMin() const;

  //true if parameter set has been unpacked
  bool unpacked() const;

  //true if sister has been found
  bool foundSister() const;

  //true if visible decay products have been found
  bool foundVisibleDecayProducts() const;

  //true if sister visible decay products have been found
  bool foundSisterVisibleDecayProducts() const;

  //true if pT rank has been set
  bool pTRankSet() const;

  //find sister
  void findSister();

  //is tau a status 3 decay product?
  bool tauIsStatus3DecayProduct() const;

  //is particle a status 3 decay product, optionally only consider particles with given PDG ID
  bool isStatus3DecayProduct(const int PDGID = 0) const;

  //get tau decay type
  DecayType tauDecayType() const;

  //get tau sister decay type
  DecayType sisterDecayType() const;

  //get tau hadronic decay type
  reco::PFTau::hadronicDecayMode tauHadronicDecayType(const bool countKShort = false) const;

  //get tau sister hadronic decay type
  reco::PFTau::hadronicDecayMode sisterHadronicDecayType(const bool countKShort = false) const;

  /*get fully specified tau decay type, with option to only consider a decay type as valid if pT 
    cuts on the decay products and on the total visible pT are passed*/
  std::pair<reco::PFTau::hadronicDecayMode, DecayType> tauDecayType(const bool, const bool);

  /*get fully specified tau sister decay type, with option to only consider a decay type as valid 
    if pT cuts on the decay products and on the total visible pT are passed*/
  std::pair<reco::PFTau::hadronicDecayMode, DecayType> sisterDecayType(const bool, const bool);

  /*get pT rank (i.e. 0 means highest pT object with specified mother, 1 means 2nd highest pT 
    object with specified mother, etc.)*/
  unsigned int getPTRank(const bool warn = true) const;

  /*set pT rank (i.e. 0 means highest pT object with specified mother, 1 means 2nd highest pT 
    object with specified mother, etc.)*/
  void setPTRank(const unsigned int);

 private:

  //parameters that might be needed in the decay ID (e.g. eta)
  edm::ParameterSet* pParSet_;

  //gen particle collection handle
  edm::Handle<reco::GenParticleCollection> pGenParticles_;

  //index of the tau
  unsigned int iTau_;

  //index of the tau sister
  unsigned int iSister_;

  //vector of PDG IDs of tau mothers
  std::vector<int> momPDGID_;

  //visible 4-vector of tau
  reco::LeafCandidate::LorentzVector visibleTauP4_;

  //visible 4-vector of tau sister
  reco::LeafCandidate::LorentzVector visibleTauSisterP4_;

  //overall decay type of tau
  DecayType tauDecayType_;

  //overall decay type of tau sister
  DecayType tauSisterDecayType_;

  //hadronic decay type of tau
  reco::PFTau::hadronicDecayMode tauHadronicDecayType_;

  //hadronic decay type of tau sister
  reco::PFTau::hadronicDecayMode tauSisterHadronicDecayType_;

  //minimum pT to be counted as a charged hadron in hadronic tau decay
  double chargedHadronPTMin_;

  //minimum pT to be counted as a neutral hadron in hadronic tau decay
  double neutralHadronPTMin_;

  //minimum pT to be counted as a charged lepton in leptonic tau decay or mother decay
  double chargedLeptonPTMin_;

  //minimum pT of the visible tau decay products or gen lepton from mother decay
  double totalPTMin_;

  //pT rank of this object in the event
  unsigned int pTRank_;

  //true if parameter set has been unpacked
  bool unpacked_;

  //true if pGenParticles_ points to a valid handle
  bool validHandle_;

  //true if tau sister has been found
  bool foundSister_;

  //true if visible decay products have been found
  bool foundVisibleDecayProducts_;

  //true if sister visible decay products have been found
  bool foundSisterVisibleDecayProducts_;

  //true if pT rank was set
  bool pTRankSet_;

  //unpack parameter set
  void unpackParSet(const bool warn = true);

  //find sister decay product
  int sister() const;

  //classify decay for given particle index
  DecayType decayType(const unsigned int) const;

  /*recursively find the total number of charged and neutral hadrons in a tau decay, optionally 
    counting Kshorts*/
  void numChargedAndNeutralHadronsInTauDecay(const reco::GenParticleRef&, unsigned int&, 
					     unsigned int&, const bool) const;

  //classify hadronic decay
  reco::PFTau::hadronicDecayMode hadronicDecayType(const unsigned int, const bool) const;

  /*recursively find the total number of charged and neutral hadrons and charged leptons in a tau 
    decay, optionally counting Kshorts as neutral hadrons and applying pT cuts, and return a 
    4-vector sum of the visible decay products*/
  reco::LeafCandidate::LorentzVector numAndP4VisibleDecayProducts(const reco::GenParticleRef&, 
								  unsigned int&, unsigned int&, 
								  unsigned int&, unsigned int&, 
								  const bool, const bool) const;

  //get decay type and p4 of visible decay products
  void decayTypeAndVisibleP4(const unsigned int, reco::LeafCandidate::LorentzVector&, 
			     std::pair<reco::PFTau::hadronicDecayMode, DecayType>&, const bool, 
			     const bool) const;

  //classify decay with pT cuts on visible decay products
  //need to add a flag to indicate when iParticle does not refer to a tau
  std::pair<reco::PFTau::hadronicDecayMode, DecayType>
    decayTypeWithPTCuts(const unsigned int, const bool, const bool, 
			reco::LeafCandidate::LorentzVector&);

  //uniquely classify tau decay
  std::pair<reco::PFTau::hadronicDecayMode, DecayType> decayType(const unsigned int, 
							    const unsigned int, 
							    const unsigned int, 
							    const unsigned int) const;

  //if particle is not a tau, set visible 4-vector to its 4-vector and set the necessary flags true
  void checkIfTau(const unsigned int, reco::LeafCandidate::LorentzVector&, DecayType&, 
		  reco::PFTau::hadronicDecayMode&, bool&);

  //warning that sister was never found
  std::string warnSisterNeverFound(const std::string&) const;

  //warning that parameter could not be found
  std::string warnParameterNotFound(const std::string&, const std::string&) const;

  //warning that tau visible decay products were never found
  std::string warnVisibleDecayProductsNeverFound(const std::string&) const;

  //warning that pT rank was not set
  std::string warnPTRankNotSet(const std::string&) const;

  //error that sister could not be found
  std::string errorSisterNotFound(const std::string&, const unsigned int) const;

  //error that gen particle handle is invalid
  std::string errorInvalidGenParticleHandle(const std::string&) const;

  //error that parameter set is not unpacked
  std::string errorParameterSetNotUnpacked(const std::string&) const;

  //error that supplied parameter set pointer is null
  std::string errorNullParameterSetPointer(const std::string&) const;

  //error that number of daughters is unexpected
  std::string errorUnexpectedNumDaughters(const std::string&, const unsigned int) const;

  //error that numbers of types of decay particles are inconsistent with a tau decay
  std::string errorInvalidTauDecay(const std::string&, const unsigned int, const unsigned int, 
				   const unsigned int, const unsigned int) const;

  //error that 2 4-vectors that are supposed to be equal aren't
  std::string errorUnequalP4(const std::string&) const;
};

#endif
