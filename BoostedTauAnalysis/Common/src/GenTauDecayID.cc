#include <string>
#include <sstream>
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"

//default constructor
GenTauDecayID::GenTauDecayID() :
  pParSet_(NULL),
  iTau_(0),
  iSister_(0),
  tauDecayType_(UNKNOWN),
  tauSisterDecayType_(UNKNOWN),
  tauHadronicDecayType_(reco::PFTau::kNull),
  tauSisterHadronicDecayType_(reco::PFTau::kNull),
  chargedHadronPTMin_(0.0),
  neutralHadronPTMin_(0.0),
  chargedLeptonPTMin_(0.0),
  totalPTMin_(0.0),
  pTRank_(0),
  unpacked_(false),
  validHandle_(false),
  foundSister_(false),
  foundVisibleDecayProducts_(false),
  foundSisterVisibleDecayProducts_(false),
  pTRankSet_(false) {}

//constructor that takes a parameter set of arguments
GenTauDecayID::GenTauDecayID(edm::ParameterSet& parSet, 
			     const edm::Handle<reco::GenParticleCollection>& pGenParticles, 
			     const unsigned int iTau, const bool warn) :
  pParSet_(&parSet),
  pGenParticles_(pGenParticles),
  iTau_(iTau),
  iSister_(0),
  tauDecayType_(UNKNOWN),
  tauSisterDecayType_(UNKNOWN),
  tauHadronicDecayType_(reco::PFTau::kNull),
  tauSisterHadronicDecayType_(reco::PFTau::kNull),
  chargedHadronPTMin_(0.0),
  neutralHadronPTMin_(0.0),
  chargedLeptonPTMin_(0.0),
  totalPTMin_(0.0),
  pTRank_(0),
  unpacked_(false),
  validHandle_(pGenParticles_.isValid()),
  foundSister_(false),
  foundVisibleDecayProducts_(false),
  foundSisterVisibleDecayProducts_(false),
  pTRankSet_(false)
{
  try { unpackParSet(warn); }
  catch (std::string& ex) { throw; }
}

//copy constructor
GenTauDecayID::GenTauDecayID(const GenTauDecayID& other) :
  pParSet_(other.getParSet()),
  pGenParticles_(other.getGenParticleHandle()),
  iTau_(other.getTauIndex()),
  iSister_(other.getSisterIndex(false)),
  visibleTauP4_(other.getVisibleTauP4(false)),
  visibleTauSisterP4_(other.getVisibleTauSisterP4(false)),
  tauDecayType_(other.getTauDecayType(false)),
  tauSisterDecayType_(other.getTauSisterDecayType(false)),
  tauHadronicDecayType_(other.getTauHadronicDecayType(false)),
  tauSisterHadronicDecayType_(other.getTauSisterHadronicDecayType(false)),
  chargedHadronPTMin_(0.0),
  neutralHadronPTMin_(0.0),
  chargedLeptonPTMin_(0.0),
  totalPTMin_(0.0),
  pTRank_(other.getPTRank(false)),
  unpacked_(other.unpacked()),
  validHandle_(pGenParticles_.isValid()),
  foundSister_(other.foundSister()),
  foundVisibleDecayProducts_(other.foundVisibleDecayProducts()),
  foundSisterVisibleDecayProducts_(other.foundSisterVisibleDecayProducts()),
  pTRankSet_(other.pTRankSet())
{
  try { unpackParSet(); }
  catch (std::string& ex) { throw; }
}

//destructor
GenTauDecayID::~GenTauDecayID()
{
  pParSet_ = NULL;
  pGenParticles_.clear();
  iTau_ = 0;
  iSister_ = 0;
  momPDGID_.clear();
  visibleTauP4_ = reco::LeafCandidate::LorentzVector();
  visibleTauSisterP4_ = reco::LeafCandidate::LorentzVector();
  tauDecayType_ = UNKNOWN;
  tauSisterDecayType_ = UNKNOWN;
  tauHadronicDecayType_ = reco::PFTau::kNull;
  tauSisterHadronicDecayType_ = reco::PFTau::kNull;
  chargedHadronPTMin_ = 0.0;
  neutralHadronPTMin_ = 0.0;
  chargedLeptonPTMin_ = 0.0;
  totalPTMin_ = 0.0;
  pTRank_ = 0;
  unpacked_ = false;
  validHandle_ = false;
  foundSister_ = false;
  foundVisibleDecayProducts_ = false;
  foundSisterVisibleDecayProducts_ = false;
  pTRankSet_ = false;
}

//assignment operator
GenTauDecayID& GenTauDecayID::operator=(const GenTauDecayID& rhs)
{
  if (this != &rhs) {
    pParSet_ = rhs.getParSet();
    pGenParticles_ = rhs.getGenParticleHandle();
    iTau_ = rhs.getTauIndex();
    iSister_ = rhs.getSisterIndex(false);
    visibleTauP4_ = rhs.getVisibleTauP4(false);
    visibleTauSisterP4_ = rhs.getVisibleTauSisterP4(false);
    tauDecayType_ = rhs.getTauDecayType(false);
    tauSisterDecayType_ = rhs.getTauSisterDecayType(false);
    tauHadronicDecayType_ = rhs.getTauHadronicDecayType(false);
    tauSisterHadronicDecayType_ = rhs.getTauSisterHadronicDecayType(false);
    chargedHadronPTMin_ = 0.0;
    neutralHadronPTMin_ = 0.0;
    chargedLeptonPTMin_ = 0.0;
    totalPTMin_ = 0.0;
    pTRank_ = rhs.getPTRank(false);
    unpacked_ = rhs.unpacked();
    validHandle_ = pGenParticles_.isValid();
    foundSister_ = rhs.foundSister();
    foundVisibleDecayProducts_ = rhs.foundVisibleDecayProducts();
    foundSisterVisibleDecayProducts_ = rhs.foundSisterVisibleDecayProducts();
    pTRankSet_ = rhs.pTRankSet();
    try { unpackParSet(); }
    catch (std::string& ex) { throw; }
  }
  return *this;
}

//getter for parameter set
edm::ParameterSet* GenTauDecayID::getParSet() const { return pParSet_; }

//getter for gen particle handle
edm::Handle<reco::GenParticleCollection> GenTauDecayID::getGenParticleHandle() const
{
  return pGenParticles_;
}

//getter for tau index
unsigned int GenTauDecayID::getTauIndex() const { return iTau_; }

//getter for tau sister index
unsigned int GenTauDecayID::getSisterIndex(const bool warn) const
{
  std::string fnName("const unsigned int GenTauDecayID::getSisterIndex(const bool warn) const");
  if (warn && !foundSister_) std::cerr << warnSisterNeverFound(fnName);
  return iSister_;
}

//getter for visible tau 4-vector
reco::LeafCandidate::LorentzVector GenTauDecayID::getVisibleTauP4(const bool warn) const
{
  std::string fnName("reco::LeafCandidate::LorentzVector ");
  fnName+="GenTauDecayID::getVisibleTauP4(const bool warn) const";
  if (warn && !foundVisibleDecayProducts_) std::cerr << warnVisibleDecayProductsNeverFound(fnName);
  return visibleTauP4_;
}

//getter for visible tau sister 4-vector
reco::LeafCandidate::LorentzVector GenTauDecayID::getVisibleTauSisterP4(const bool warn) const
{
  std::string fnName("reco::LeafCandidate::LorentzVector ");
  fnName+="GenTauDecayID::getVisibleTauSisterP4(const bool warn) const";
  if (warn && !foundSister_) std::cerr << warnSisterNeverFound(fnName);
  if (warn && !foundVisibleDecayProducts_) std::cerr << warnVisibleDecayProductsNeverFound(fnName);
  return visibleTauSisterP4_;
}

//getter for overall tau decay type
GenTauDecayID::DecayType GenTauDecayID::getTauDecayType(const bool warn) const
{
  std::string fnName("GenTauDecayID::DecayType ");
  fnName+="GenTauDecayID::getTauDecayType(const bool warn) const";
  if (warn && !foundVisibleDecayProducts_) std::cerr << warnVisibleDecayProductsNeverFound(fnName);
  return tauDecayType_;
}

//getter for overall tau sister decay type
GenTauDecayID::DecayType GenTauDecayID::getTauSisterDecayType(const bool warn) const
{
  std::string fnName("GenTauDecayID::DecayType ");
  fnName+="GenTauDecayID::getTauSisterDecayType(const bool warn) const";
  if (warn && !foundSister_) std::cerr << warnSisterNeverFound(fnName);
  if (warn && !foundVisibleDecayProducts_) std::cerr << warnVisibleDecayProductsNeverFound(fnName);
  return tauSisterDecayType_;
}

//getter for hadronic tau decay type
reco::PFTau::hadronicDecayMode GenTauDecayID::getTauHadronicDecayType(const bool warn) const
{
  std::string fnName("reco::PFTau::hadronicDecayMode ");
  fnName+="GenTauDecayID::getTauHadronicDecayType(const bool warn) const";
  if (warn && !foundVisibleDecayProducts_) std::cerr << warnVisibleDecayProductsNeverFound(fnName);
  return tauHadronicDecayType_;
}

//getter for hadronic tau sister decay type
reco::PFTau::hadronicDecayMode GenTauDecayID::getTauSisterHadronicDecayType(const bool warn) const
{
  std::string fnName("reco::PFTau::hadronicDecayMode ");
  fnName+="GenTauDecayID::getTauSisterHadronicDecayType(const bool warn) const";
  if (warn && !foundSister_) std::cerr << warnSisterNeverFound(fnName);
  if (warn && !foundVisibleDecayProducts_) std::cerr << warnVisibleDecayProductsNeverFound(fnName);
  return tauSisterHadronicDecayType_;
}

//getter for charged hadron min pT
double GenTauDecayID::getChargedHadronPTMin() const { return chargedHadronPTMin_; }

//getter for neutral hadron min pT
double GenTauDecayID::getNeutralHadronPTMin() const { return neutralHadronPTMin_; }

//getter for charged lepton min pT
double GenTauDecayID::getChargedLeptonPTMin() const { return chargedLeptonPTMin_; }

//getter for total min pT of visible decay products
double GenTauDecayID::getTotalPTMin() const { return totalPTMin_; }

//true if parameter set has been unpacked
bool GenTauDecayID::unpacked() const { return unpacked_; }

//true if sister has been found
bool GenTauDecayID::foundSister() const { return foundSister_; }

//true if visible decay products have been found
bool GenTauDecayID::foundVisibleDecayProducts() const { return foundVisibleDecayProducts_; }

//true if sister visible decay products have been found
bool GenTauDecayID::foundSisterVisibleDecayProducts() const
{
  return foundSisterVisibleDecayProducts_;
}

//true if pT rank has been set
bool GenTauDecayID::pTRankSet() const { return pTRank_; }

//find sister
void GenTauDecayID::findSister()
{
  if (!foundSister_) {
    int iSister = -1;
    try { iSister = sister(); }
    catch (std::string& ex) { throw; }
    if (iSister != -1) {
      iSister_ = iSister;
      foundSister_ = true;
    }
    else throw errorSisterNotFound("void GenTauDecayID::findSister()", iTau_);
  }
}

//is tau a status 3 decay product?
bool GenTauDecayID::tauIsStatus3DecayProduct() const
{
  std::string fnName("bool GenTauDecayID::tauIsStatus3DecayProduct() const");
  bool ret = false;
  try { ret = isStatus3DecayProduct(TAUPDGID); }
  catch (std::string& ex) { throw; }
  return ret;
}

//is particle a status 3 decay product, optionally only consider particles with given PDG ID
bool GenTauDecayID::isStatus3DecayProduct(const int PDGID) const
{
  std::string fnName("bool GenTauDecayID::isStatus3DecayProduct() const");
  bool ret = false;
  if (unpacked_) {
    if (validHandle_) {
      reco::GenParticleRef particleRef(pGenParticles_, iTau_);
      ret = (((PDGID == ANY_PDGID) || (fabs(particleRef->pdgId()) == PDGID)) && 
	     (particleRef->status() == 3)); 
      bool rightMom = false;
      std::vector<int>::const_iterator iPDGID = momPDGID_.begin();
      while ((iPDGID != momPDGID_.end()) && !rightMom) {
	rightMom = (*iPDGID == ANY_PDGID) || //don't consider mother PDG ID in selection
	  ((particleRef->numberOfMothers() > 0) && /*insure that the mothers of particles with 0 
						     mothers (i.e. the protons) are not accessed*/
	   (fabs(particleRef->mother()->pdgId()) == fabs(*iPDGID))); /*require that the 0th mother 
								     match 1 of these PDG IDs
								     in the case of >1 mother, the 
								     sister should share both 
								     mothers anyway so it's enough 
								     to just check one of them*/
	++iPDGID;
      }
      ret = ret && rightMom;
    }
    else throw errorInvalidGenParticleHandle(fnName);
  }
  else throw errorParameterSetNotUnpacked(fnName);
  return ret;
}

//get tau decay type
GenTauDecayID::DecayType GenTauDecayID::tauDecayType() const
{
  DecayType decayTypeCode = UNKNOWN;
  try { decayTypeCode = decayType(iTau_); }
  catch (std::string& ex) { throw; }
  return decayTypeCode;
}

//get tau sister decay type
GenTauDecayID::DecayType GenTauDecayID::sisterDecayType() const
{
  std::string fnName("GenTauDecayID::DecayType GenTauDecayID::sisterDecayType() const");
  if (!foundSister_) std::cerr << warnSisterNeverFound(fnName);
  DecayType decayTypeCode = UNKNOWN;
  try { decayTypeCode = decayType(iSister_); }
  catch (std::string& ex) { throw; }
  return decayTypeCode;
}

//get tau hadronic decay type
reco::PFTau::hadronicDecayMode GenTauDecayID::tauHadronicDecayType(const bool countKShort) const
{
  reco::PFTau::hadronicDecayMode hadronicDecayTypeCode = reco::PFTau::kNull;
  try { hadronicDecayTypeCode = hadronicDecayType(iTau_, countKShort); }
  catch (std::string& ex) { throw; }
  return hadronicDecayTypeCode;
}

//get tau hadronic sister decay type
reco::PFTau::hadronicDecayMode GenTauDecayID::sisterHadronicDecayType(const bool countKShort) const
{
  std::string fnName("reco::PFTau::hadronicDecayMode GenTauDecayID::sisterHadronicDecayType(");
  fnName+="const bool countKShort) const";
  if (!foundSister_) std::cerr << warnSisterNeverFound(fnName);
  reco::PFTau::hadronicDecayMode hadronicDecayTypeCode = reco::PFTau::kNull;
  try { hadronicDecayTypeCode = hadronicDecayType(iSister_, countKShort); }
  catch (std::string& ex) { throw; }
  return hadronicDecayTypeCode;
}

/*get pT rank (i.e. 0 means highest pT object with specified mother, 1 means 2nd highest pT object 
  with specified mother, etc.)*/
unsigned int GenTauDecayID::getPTRank(const bool warn) const {
  std::string fnName("unsigned int GenTauDecayID::getPTRank(const bool warn) const");
  if (warn && !pTRankSet_) std::cerr << warnPTRankNotSet(fnName);
  return pTRank_;
}

/*set pT rank (i.e. 0 means highest pT object with specified mother, 1 means 2nd highest pT object 
  with specified mother, etc.)*/
void GenTauDecayID::setPTRank(const unsigned int pTRank)
{
  pTRank_ = pTRank;
  pTRankSet_ = true;
}

/*get fully specified tau decay type, with option to only consider a decay type as valid if pT 
  cuts on the decay products and on the total visible pT are passed*/
std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType>
GenTauDecayID::tauDecayType(const bool applyPTCuts, const bool countKShort)
{
  checkIfTau(iTau_, visibleTauP4_, tauDecayType_, tauHadronicDecayType_, 
	     foundVisibleDecayProducts_);
  std::pair<reco::PFTau::hadronicDecayMode, DecayType> decayTypeCode(reco::PFTau::kNull, UNKNOWN);
  if (foundVisibleDecayProducts_) {
    decayTypeCode = 
      std::pair<reco::PFTau::hadronicDecayMode, DecayType>(tauHadronicDecayType_, tauDecayType_);
  }
  else {
    try { decayTypeCode = decayTypeWithPTCuts(iTau_, countKShort, applyPTCuts, visibleTauP4_); }
    catch (std::string& ex) { throw; }
    tauDecayType_ = decayTypeCode.second;
    tauHadronicDecayType_ = decayTypeCode.first;
    foundVisibleDecayProducts_ = true;
  }
  return decayTypeCode;
}

/*get fully specified tau sister decay type, with option to only consider a decay type as valid if 
  pT cuts on the decay products and on the total visible pT are passed*/
std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType>
GenTauDecayID::sisterDecayType(const bool applyPTCuts, const bool countKShort)
{
  std::string fnName("std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType>");
  fnName+="GenTauDecayID::sisterDecayType(const bool applyPTCuts, ";
  fnName+="const bool countKShort) const";
  if (!foundSister_) std::cerr << warnSisterNeverFound(fnName);
  checkIfTau(iSister_, visibleTauSisterP4_, tauSisterDecayType_, tauSisterHadronicDecayType_, 
	     foundSisterVisibleDecayProducts_);
  std::pair<reco::PFTau::hadronicDecayMode, DecayType> decayTypeCode(reco::PFTau::kNull, UNKNOWN);
  if (foundSisterVisibleDecayProducts_) {
    decayTypeCode = 
      std::pair<reco::PFTau::hadronicDecayMode, DecayType>(tauSisterHadronicDecayType_, 
							   tauSisterDecayType_);
  }
  else {
    try {
      decayTypeCode = decayTypeWithPTCuts(iSister_, countKShort, applyPTCuts, visibleTauSisterP4_);
    }
    catch (std::string& ex) { throw; }
    tauSisterDecayType_ = decayTypeCode.second;
    tauSisterHadronicDecayType_ = decayTypeCode.first;
    foundSisterVisibleDecayProducts_ = true;
  }
  return decayTypeCode;
}

//unpack parameter set
void GenTauDecayID::unpackParSet(const bool warn)
{
  std::string fnName("void GenTauDecayID::unpackParSet()");
  if (pParSet_ == NULL) {
    throw errorNullParameterSetPointer(fnName);
  }
  if (pParSet_->existsAs<std::vector<int> >("momPDGID")) {
    momPDGID_ = pParSet_->getParameter<std::vector<int> >("momPDGID");
  }
  else {
    std::cerr << warnParameterNotFound(fnName, "momPDGID");
    momPDGID_.clear();
  }
  if (pParSet_->existsAs<double>("chargedHadronPTMin")) {
    chargedHadronPTMin_ = pParSet_->getParameter<double>("chargedHadronPTMin");
  }
  else {
    if (warn) std::cerr << warnParameterNotFound(fnName, "chargedHadronPTMin");
    chargedHadronPTMin_ = 0.0;
  }
  if (pParSet_->existsAs<double>("neutralHadronPTMin")) {
    neutralHadronPTMin_ = pParSet_->getParameter<double>("neutralHadronPTMin");
  }
  else {
    if (warn) std::cerr << warnParameterNotFound(fnName, "neutralHadronPTMin");
    neutralHadronPTMin_ = 0.0;
  }
  if (pParSet_->existsAs<double>("chargedLeptonPTMin")) {
    chargedLeptonPTMin_ = pParSet_->getParameter<double>("chargedLeptonPTMin");
  }
  else {
    if (warn) std::cerr << warnParameterNotFound(fnName, "chargedLeptonPTMin");
    chargedLeptonPTMin_ = 0.0;
  }
  if (pParSet_->existsAs<double>("totalPTMin")) {
    totalPTMin_ = pParSet_->getParameter<double>("totalPTMin");
  }
  else {
    if (warn) std::cerr << warnParameterNotFound(fnName, "totalPTMin");
    totalPTMin_ = 0.0;
  }
  unpacked_ = true;
}

//find first sister decay product
int GenTauDecayID::sister() const
{
  int iSister = -1;
  if (validHandle_) {
    reco::GenParticleRef tauRef(pGenParticles_, iTau_);
    reco::GenParticleRef momRef = tauRef->motherRef();
    unsigned int iDaughter = 0;
    while ((iDaughter < momRef->numberOfDaughters()) && (iSister == -1)) {
      reco::GenParticleRef childRef = momRef->daughterRef(iDaughter);
      unsigned int childRefKey = childRef.key();
      if ((childRefKey != iTau_) && (childRef->status() == tauRef->status())) {
	iSister = childRefKey;
      }
      ++iDaughter;
    }
  }
  else throw errorInvalidGenParticleHandle("int GenTauDecayID::sister() const");
  return iSister;
}

//classify decay
GenTauDecayID::DecayType GenTauDecayID::decayType(const unsigned int iParticle) const
{
  std::string fnName("unsigned int GenTauDecayID::decayType() const");
  DecayType ret = UNKNOWN;
  unsigned int leptonPDGID = 0;
  unsigned int neutrinoPDGID = 0;
  if (validHandle_) {
    reco::GenParticleRef tauRef(pGenParticles_, iParticle);
    const size_t numDaughters = tauRef->numberOfDaughters();
    if (numDaughters == 1) {
      const reco::Candidate* daughter = tauRef->daughter(0);
      reco::Candidate::const_iterator iDaughter = daughter->begin();
      while ((iDaughter != daughter->end()) && 
	     ((leptonPDGID == 0) || (neutrinoPDGID == 0))) {
	if (leptonPDGID == 0) {
	  if (((neutrinoPDGID == 0) || (neutrinoPDGID == ENEUTRINOPDGID)) && 
	      (fabs(iDaughter->pdgId()) == EPDGID)) leptonPDGID = EPDGID;
	  if (((neutrinoPDGID == 0) || (neutrinoPDGID == MUNEUTRINOPDGID)) && 
	      (fabs(iDaughter->pdgId()) == MUPDGID)) leptonPDGID = MUPDGID;
	}
	if (neutrinoPDGID == 0) {
	  if (((leptonPDGID == 0) || (leptonPDGID == EPDGID)) && 
	      (fabs(iDaughter->pdgId()) == ENEUTRINOPDGID)) neutrinoPDGID = ENEUTRINOPDGID;
	  if (((leptonPDGID == 0) || (leptonPDGID == MUPDGID)) && 
	      (fabs(iDaughter->pdgId()) == MUNEUTRINOPDGID)) {
	    neutrinoPDGID = MUNEUTRINOPDGID;
	  }
	}
	++iDaughter;
      }
    }
    else throw errorUnexpectedNumDaughters(fnName, numDaughters);
  }
  else throw errorInvalidGenParticleHandle(fnName);
  if ((leptonPDGID == 0) || (neutrinoPDGID == 0)) ret = HAD;
  if ((leptonPDGID == EPDGID) || (neutrinoPDGID == ENEUTRINOPDGID)) ret = E;
  if ((leptonPDGID == MUPDGID) || (neutrinoPDGID == MUNEUTRINOPDGID)) ret = MU;
  return ret;
}

/*recursively find the total number of charged and neutral hadrons in a tau decay, optionally 
  counting Kshorts*/
void GenTauDecayID::numChargedAndNeutralHadronsInTauDecay(const reco::GenParticleRef& momRef, 
							  unsigned int& nChargedHadrons, 
							  unsigned int& nNeutralHadrons, 
							  const bool countKShort) const
{
  for (unsigned int iDaughter = 0; iDaughter < momRef->numberOfDaughters(); 
       ++iDaughter) {
    reco::GenParticleRef kidRef = momRef->daughterRef(iDaughter);
    const unsigned int absDaughterPDGID = fabs(kidRef->pdgId());
    if ((absDaughterPDGID == 111) || (countKShort && (absDaughterPDGID == 310)) || 
	(absDaughterPDGID == 221)) ++nNeutralHadrons;
    else if ((absDaughterPDGID == 211) || (absDaughterPDGID == 321)) ++nChargedHadrons;
    else if (kidRef->status() != 1) {
      numChargedAndNeutralHadronsInTauDecay(kidRef, nChargedHadrons, nNeutralHadrons, countKShort);
    }
  }
}

//classify hadronic decay
reco::PFTau::hadronicDecayMode GenTauDecayID::hadronicDecayType(const unsigned int iParticle, 
								const bool countKShort) const
{
  std::string fnName("reco::PFTau::hadronicDecayMode GenTauDecayID::hadronicDecayType(const ");
  fnName+="unsigned int iParticle, const bool countKShort) const";
  reco::PFTau::hadronicDecayMode type = reco::PFTau::kNull;
  unsigned int nChargedHadrons = 0;
  unsigned int nPi0s = 0;
  if (validHandle_) {
    numChargedAndNeutralHadronsInTauDecay(reco::GenParticleRef(pGenParticles_, 
							       iParticle)->daughterRef(0), 
					  nChargedHadrons, nPi0s, countKShort);
  }
  else throw errorInvalidGenParticleHandle(fnName);
  if (nChargedHadrons > 0) {
    const unsigned int maxPiZeros = reco::PFTau::kOneProngNPiZero;
    unsigned int trackIndex = (nChargedHadrons - 1)*(maxPiZeros + 1);
    if (trackIndex >= reco::PFTau::kRareDecayMode) type = reco::PFTau::kRareDecayMode;
    else {
      nPi0s = (nPi0s <= maxPiZeros) ? nPi0s : maxPiZeros;
      type = static_cast<reco::PFTau::hadronicDecayMode>(trackIndex + nPi0s);
    }
  }
  return type;
}

/*recursively find the total number of charged and neutral hadrons and charged leptons in a tau 
  decay, optionally counting Kshorts as neutral hadrons and applying pT cuts, and return a 
  4-vector sum of the visible decay products*/
reco::LeafCandidate::LorentzVector
GenTauDecayID::numAndP4VisibleDecayProducts(const reco::GenParticleRef& momRef, 
					    unsigned int& nChargedHadrons, 
					    unsigned int& nNeutralHadrons, 
					    unsigned int& nElectrons, unsigned int& nMuons, 
					    const bool countKShort, const bool applyPTCuts) const
{
  std::string fnName("reco::LeafCandidate::LorentzVector ");
  fnName+="GenTauDecayID::numAndP4VisibleDecayProducts(const reco::GenParticleRef& momRef, ";
  fnName+="unsigned int& nChargedHadrons, unsigned int& nNeutralHadrons, ";
  fnName+="unsigned int& nElectrons, unsigned int& nMuons, const bool countKShort, ";
  fnName+="const bool applyPTCuts) const)";
  reco::LeafCandidate::LorentzVector p4;
  for (unsigned int iDaughter = 0; iDaughter < momRef->numberOfDaughters(); 
       ++iDaughter) {
    const reco::GenParticleRef kidRef = momRef->daughterRef(iDaughter);
    const unsigned int absDaughterPDGID = fabs(kidRef->pdgId());
    const reco::LeafCandidate::LorentzVector daughterP4 = kidRef->p4();
    const double daughterPT = kidRef->pt();
    if (unpacked_) {
      if (((absDaughterPDGID == 111) || (countKShort && (absDaughterPDGID == 310)) || 
	   (absDaughterPDGID == 221)) && (!applyPTCuts || (daughterPT > neutralHadronPTMin_))) {
	++nNeutralHadrons;
	p4+=daughterP4;
      }
      else if (((absDaughterPDGID == 211) || (absDaughterPDGID == 321)) && 
	       (!applyPTCuts || (daughterPT > chargedHadronPTMin_))) {
	++nChargedHadrons;
	p4+=daughterP4;
      }
      else if ((absDaughterPDGID == EPDGID) && 
	       (!applyPTCuts || (daughterPT > chargedLeptonPTMin_))) {
	++nElectrons;
	p4+=daughterP4;
      }
      else if ((absDaughterPDGID == MUPDGID) && 
	       (!applyPTCuts || (daughterPT > chargedLeptonPTMin_))) {
	++nMuons;
	p4+=daughterP4;
      }
      else if (kidRef->status() != 1) {
	p4+=numAndP4VisibleDecayProducts(kidRef, nChargedHadrons, nNeutralHadrons, nElectrons, 
					 nMuons, countKShort, applyPTCuts);
      }
    }
    else throw errorParameterSetNotUnpacked(fnName);
  }
  return p4;
}

//get decay type and p4 of visible decay products
void GenTauDecayID::decayTypeAndVisibleP4(const unsigned int iParticle, 
					  reco::LeafCandidate::LorentzVector& visibleP4, 
					  std::pair<reco::PFTau::hadronicDecayMode, 
					  DecayType>& decayTypePair, const bool countKShort, 
					  const bool applyPTCuts) const
{
  std::string fnName("void ");
  fnName+="GenTauDecayID::decayTypeAndVisibleP4(const unsigned int iParticle, ";
  fnName+="reco::LeafCandidate::LorentzVector& visibleP4, ";
  fnName+="std::pair<reco::PFTau::hadronicDecayMode, DecayType>& decayTypePair, ";
  fnName+="const bool countKShort, const bool applyPTCuts) const";
  unsigned int nChargedHadrons = 0;
  unsigned int nNeutralHadrons = 0;
  unsigned int nElectrons = 0;
  unsigned int nMuons = 0;
  if (validHandle_) {
    try {
      //does the daughterRef(0) cut apply only to taus?
      visibleP4 = numAndP4VisibleDecayProducts(reco::GenParticleRef(pGenParticles_, 
								    iParticle)->daughterRef(0), 
					       nChargedHadrons, nNeutralHadrons, nElectrons, 
					       nMuons, countKShort, applyPTCuts);
      decayTypePair = decayType(nChargedHadrons, nNeutralHadrons, nElectrons, nMuons);
    }
    catch (std::string& ex) { throw; }
  }
  else throw errorInvalidGenParticleHandle(fnName);
}

//classify decay with pT cuts on visible decay products
//need to add a flag to indicate when iParticle does not refer to a tau
std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType>
GenTauDecayID::decayTypeWithPTCuts(const unsigned int iParticle, const bool countKShort, 
				   const bool applyPTCuts, 
				   reco::LeafCandidate::LorentzVector& visibleP4)
{
  std::string fnName("std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> ");
  fnName+="GenTauDecayID::decayTypeWithPTCuts(const unsigned int iParticle, ";
  fnName+="const bool countKShort, const bool applyPTCuts) const";
  reco::LeafCandidate::LorentzVector visibleP4NoPTCuts;
  std::pair<reco::PFTau::hadronicDecayMode, DecayType> 
    decayTypePairNoPTCuts(reco::PFTau::kNull, UNKNOWN);
  reco::LeafCandidate::LorentzVector visibleP4PTCuts;
  std::pair<reco::PFTau::hadronicDecayMode, DecayType> 
    decayTypePairPTCuts(reco::PFTau::kNull, UNKNOWN);
  try {
    decayTypeAndVisibleP4(iParticle, visibleP4NoPTCuts, decayTypePairNoPTCuts, countKShort, false);
    if (applyPTCuts) {
      decayTypeAndVisibleP4(iParticle, visibleP4PTCuts, decayTypePairPTCuts, countKShort, true);
    }
  }
  catch (std::string& ex) { throw; }
  const double visiblePTNoPTCuts = visibleP4NoPTCuts.Pt();
  const double visiblePTPTCuts = visibleP4PTCuts.Pt();
  if ((applyPTCuts && ((decayTypePairNoPTCuts != decayTypePairPTCuts) || 
		       (visiblePTPTCuts < totalPTMin_))) || (visiblePTNoPTCuts < totalPTMin_)) {
    decayTypePairNoPTCuts = 
      std::pair<reco::PFTau::hadronicDecayMode, DecayType>(reco::PFTau::kNull, UNKNOWN);
  }
  else if (applyPTCuts && (visibleP4NoPTCuts != visibleP4PTCuts)) throw errorUnequalP4(fnName);
  visibleP4 = visibleP4NoPTCuts;
  return decayTypePairNoPTCuts;
}

//uniquely classify tau decay
std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType>
GenTauDecayID::decayType(const unsigned int nChargedHadrons, const unsigned int nNeutralHadrons, 
			 const unsigned int nElectrons, const unsigned int nMuons) const
{
  std::string fnName("std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> ");
  fnName+="GenTauDecayID::decayType(const unsigned int nChargedHadrons, ";
  fnName+="const unsigned int nNeutralHadrons, const unsigned int nElectrons, ";
  fnName+="const unsigned int nMuons) const";
  reco::PFTau::hadronicDecayMode hadronicDecayType = reco::PFTau::kNull;
  DecayType overallDecayType = UNKNOWN;
  if (((nChargedHadrons > 0) && ((nElectrons + nMuons) > 0)) || 
      ((nElectrons > 1) || (nMuons > 1)) || 
      ((nChargedHadrons == 0) && (nElectrons == 1) && (nMuons == 1))) {
    throw errorInvalidTauDecay(fnName, nChargedHadrons, nNeutralHadrons, nElectrons, nMuons);
  }
  if (nChargedHadrons > 0) {
    overallDecayType = HAD;
    const unsigned int maxPiZeros = reco::PFTau::kOneProngNPiZero;
    unsigned int trackIndex = (nChargedHadrons - 1)*(maxPiZeros + 1);
    if (trackIndex >= reco::PFTau::kRareDecayMode) hadronicDecayType = reco::PFTau::kRareDecayMode;
    else {
      unsigned int validNNeutralHadrons = 
	(nNeutralHadrons <= maxPiZeros) ? nNeutralHadrons : maxPiZeros;
      hadronicDecayType = 
	static_cast<reco::PFTau::hadronicDecayMode>(trackIndex + validNNeutralHadrons);
    }
  }
  else if (nElectrons == 1) overallDecayType = E;
  else if (nMuons == 1) overallDecayType = MU;
  return std::pair<reco::PFTau::hadronicDecayMode, DecayType>
    (hadronicDecayType, overallDecayType);
}

//if particle is not a tau, set visible 4-vector to its 4-vector and set the necessary flags true
void GenTauDecayID::checkIfTau(const unsigned int iParticle, 
			       reco::LeafCandidate::LorentzVector& p4, 
			       DecayType& decayType, 
			       reco::PFTau::hadronicDecayMode& hadronicDecayType, 
			       bool& foundVisibleDecayProducts)
{
  std::string fnName("void GenTauDecayID::checkIfTau(const unsigned int iParticle, ");
  fnName+="reco::LeafCandidate::LorentzVector& p4, DecayType& decayType, ";
  fnName+="reco::PFTau::hadronicDecayMode& hadronicDecayType, bool& foundVisibleDecayProducts)";
  if (validHandle_) {
    reco::GenParticleRef particleRef(pGenParticles_, iParticle);
    if (fabs(particleRef->pdgId()) != TAUPDGID) {
      p4 = particleRef->p4();
      decayType = NOT_TAU;
      hadronicDecayType = reco::PFTau::kNull;
      foundVisibleDecayProducts = true;
    }
  }
  else throw errorInvalidGenParticleHandle(fnName);
}

//warning that sister was never found
std::string GenTauDecayID::warnSisterNeverFound(const std::string& fnName) const
{
  return ("Warning in " + fnName + 
	  ":\nSister index should not be trusted because sister was never found.\n");
}

//warning that parameter could not be found
std::string GenTauDecayID::warnParameterNotFound(const std::string& fnName, 
						 const std::string& par) const
{
  return ("Warning in " + fnName + ":\n" + par + " not found.\n");
}

//warning that tau visible decay products were never found
std::string GenTauDecayID::warnVisibleDecayProductsNeverFound(const std::string& fnName) const
{
  return ("Warning in " + fnName + ":\nTau visible decay products never found.\n");
}

//warning that pT rank was not set
std::string GenTauDecayID::warnPTRankNotSet(const std::string& fnName) const
{
  return ("Error in " + fnName + ":\npT rank not set.\n");
}

//error that sister could not be found
std::string GenTauDecayID::errorSisterNotFound(const std::string& fnName, 
					       const unsigned int iTau) const
{
  std::stringstream err;
  err << "Error in " << fnName << ":\nCould not find sister of tau with index " << iTau_ << ".\n";
  return (err.str());
}

//error that gen particle handle is invalid
std::string GenTauDecayID::errorInvalidGenParticleHandle(const std::string& fnName) const
{
  return ("Error in " + fnName + ":\nInvalid handle.\n");
}

//error that parameter set is not unpacked
std::string GenTauDecayID::errorParameterSetNotUnpacked(const std::string& fnName) const
{
  return ("Error in " + fnName + ":\nParameter set not unpacked.\n");
}

//error that supplied parameter set pointer is null
std::string GenTauDecayID::errorNullParameterSetPointer(const std::string& fnName) const
{
  return ("Error in " + fnName + ":\nedm::ParameterSet* parSet_ is NULL.\n");
}

//error that number of daughters is unexpected
std::string GenTauDecayID::errorUnexpectedNumDaughters(const std::string& fnName, 
						       const unsigned int numDaughters) const
{
  std::stringstream err;
  err << "Error in " << fnName << ":\npTau_ has " << numDaughters << " daughters.\n";
  return err.str();
}

//error that numbers of types of decay particles are inconsistent with a tau decay
std::string GenTauDecayID::errorInvalidTauDecay(const std::string& fnName, 
						const unsigned int nChargedHadrons, 
						const unsigned int nNeutralHadrons, 
						const unsigned int nElectrons, 
						const unsigned int nMuons) const
{
  std::stringstream err;
  err << "Error in " << fnName << ":\n" << nChargedHadrons << " charged hadrons, ";
  err << nNeutralHadrons << " neutral hadrons, " << nElectrons << " electrons, and " << nMuons;
  err << " muons do not constitute a valid tau decay.\n";
  return err.str();
}

//error that 2 4-vectors that are supposed to be equal aren't
std::string GenTauDecayID::errorUnequalP4(const std::string& fnName) const
{
  return ("Error in " + fnName + ":\nUnequal 4-vectors.\n");
}
