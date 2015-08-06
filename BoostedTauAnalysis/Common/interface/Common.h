#ifndef BoostedTauAnalysis_Common_interface_Common_h
#define BoostedTauAnalysis_Common_interface_Common_h

#include <vector>
#include <algorithm>
#include <string>
#include <boost/assign/list_of.hpp>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TH2F.h"

class minDR {

 public:

  void setCandidate(const reco::Candidate*);

  void deleteCandidate();

  reco::Candidate* getCandidate() const;

  template<typename T>
    bool operator()(const T* cand1, const T* cand2)
    {
      return (deltaR(*cand_, *cand1) < deltaR(*cand_, *cand2));
    }

 private:

  reco::Candidate* cand_;

};

class maxM {

 public:

  void setMuonJetMap(const edm::Handle<edm::ValueMap<reco::MuonRefVector> >*);

  void deleteMuonJetMap();

  edm::Handle<edm::ValueMap<reco::MuonRefVector> >* getMuonJetMap() const;

  reco::MuonRef highestPTMu(const reco::PFTauRef& tau) const;

  bool operator()(const reco::PFTauRef&, const reco::PFTauRef&) const;

 private:

  edm::Handle<edm::ValueMap<reco::MuonRefVector> >* pMuonJetMap_;

};

class Common {

 public:

  template<typename T, typename U>
    static std::vector<double> sortByProximity(const std::vector<T*>& objectsToSort, 
					       const U& referenceObject)
    {
      std::vector<double> sortedDR2;
      for (typename std::vector<T*>::const_iterator iObject = objectsToSort.begin(); 
	   iObject != objectsToSort.end(); ++iObject) {
	sortedDR2.push_back(reco::deltaR2(referenceObject, **iObject));
      }
      std::sort(sortedDR2.begin(), sortedDR2.end());
      return sortedDR2;
    }

  static void sortByPT(std::vector<reco::Candidate*>&);

  static void sortByPT(std::vector<GenTauDecayID>&);

  static void sortByPT(std::vector<reco::PFJet*>&);

  static void sortByPT(std::vector<reco::PFTau*>&);

  static void sortByPT(std::vector<reco::Muon*>&);

  static void sortByPT(std::vector<reco::MuonRef>&);

  static void sortByPT(std::vector<reco::PFJetRef>&);

  static void sortByPT(std::vector<reco::PFTauRef>&);

  static void sortByPT(std::vector<reco::PhotonRef>&);

  template<typename T>
    static void sortByPT(std::vector<edm::Ref<T> >& objects)
    {
      std::sort(objects.begin(), objects.end(), comparePT);
    }

  template<typename T, typename U>
    static const T* nearestObject(const U& obj, const std::vector<T*>& objs, int& index)
    {
      minDR comp;
      comp.setCandidate(dynamic_cast<const reco::Candidate*>(obj.get()));
      typename std::vector<T*>::const_iterator iMinElement = 
	min_element(objs.begin(), objs.end(), comp);
      const T* nearestObj = NULL;
      index = -1;
      if (iMinElement != objs.end()) {
	nearestObj = *iMinElement;
	index = iMinElement - objs.begin();
      }
      comp.deleteCandidate();
      return nearestObj;
    }

  static void sortByMass(const edm::Handle<edm::ValueMap<reco::MuonRefVector> >&, 
			 std::vector<reco::PFTauRef>&);

  //return true if the vertex passes the standard selection
  static bool isGoodVertex(const reco::Vertex*);

  //identify the first good vertex (the "primary" (?))
  static reco::Vertex* getPrimaryVertex(edm::Handle<reco::VertexCollection>&);

  //get the number of vertices passing the standard selection
  static unsigned int numGoodVertices(edm::Handle<reco::VertexCollection>&);

  //get muon combined particle isolation with adjustable PU subtraction
  static float getMuonCombPFIso(const reco::Muon&, const double);

  //get muon combined particle isolation with nearby tau candidates' pT subtracted
  static float getMuonCombPFIsoModified(const reco::Muon&, const reco::PFTauRef&, const double);

  //get muon isolation computed from leptons only
  static float getMuonLeptonPFIso(const reco::Muon&);

  //get muon isolation with tau pt subtracted
  static float getMuonCombPFIsoMinusTau(const reco::Muon&, const reco::LeafCandidate::LorentzVector, const double);

  /*fill STL container with muons passing the 2012 tight selection, PF isolation, and |eta|
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId and 
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1)*/
  static std::vector<reco::MuonRef> 
    getTightPFIsolatedRecoMuons(const edm::Handle<reco::MuonCollection>&, const reco::Vertex*, 
				const double, const double, const double, const bool);

  /*fill STL container with muons passing the 2012 tight selection, PF isolation, and |eta|
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId and 
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1)*/
  static std::vector<reco::MuonRef> 
    getTightPFIsolatedRecoMuons(const edm::Handle<reco::MuonRefVector>&, 
				const edm::Handle<reco::MuonCollection>&, const reco::Vertex*, 
				const double, const double, const double, const bool);

  /*fill STL container with muons passing the 2012 tight selection, detector isolation, and |eta|
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId and 
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1)*/
  static std::vector<reco::MuonRef>
    getTightDetectorIsolatedRecoMuons(const edm::Handle<reco::MuonCollection>&, 
				      const reco::Vertex*, const double, const double, const bool);

  /*fill STL container with muons passing the 2012 tight selection, detector isolation, and |eta|
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId and 
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1)*/
  static std::vector<reco::MuonRef>
    getTightDetectorIsolatedRecoMuons(const edm::Handle<reco::MuonRefVector>&, 
				      const edm::Handle<reco::MuonCollection>&, 
				      const reco::Vertex*, const double, const double, const bool);

  /*fill STL container with muons passing the 2012 soft selection and |eta|
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId)*/
  static std::vector<reco::MuonRef> getSoftRecoMuons(const edm::Handle<reco::MuonCollection>&, 
						     const reco::Vertex*, const double);

  /*fill STL container with muons passing the 2012 soft selection and |eta|
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId)*/
  static std::vector<reco::MuonRef> getSoftRecoMuons(const edm::Handle<reco::MuonRefVector>&, 
						     const edm::Handle<reco::MuonCollection>&, 
						     const reco::Vertex*, const double);

  /*fill STL container with muons passing the 2012 soft selection and |eta| along with dz < 0.5 cm
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId)*/
  static std::vector<reco::MuonRef> 
    getSoftRecoMuonsFromPV(const edm::Handle<reco::MuonCollection>&, const reco::Vertex*, 
			   const double);

  //fill STL container with taus passing specified discriminators in specified eta range
  static std::vector<reco::PFTauRef>
    getRecoTaus(const edm::Handle<reco::PFTauCollection>&,
		const std::vector<edm::Handle<reco::PFTauDiscriminator> >&, 
		const edm::Handle<reco::PFTauDiscriminator>&, const double, const double, 
		const bool, const double);

  //fill STL container with taus passing specified discriminators in specified eta range
  static std::vector<reco::PFTauRef>
    getRecoTaus(const edm::Handle<reco::PFTauRefVector>&, 
		const edm::Handle<reco::PFTauCollection>&,
		const std::vector<edm::Handle<reco::PFTauDiscriminator> >&, 
		const edm::Handle<reco::PFTauDiscriminator>&, const double, const double, 
		const bool, const double);
  
  //fill STL container with photons passing specified cuts in specified eta range
  static std::vector<reco::PhotonRef> 
    getRecoPhotons(const edm::Handle<reco::PhotonCollection>&, 
		   const edm::Handle<reco::BeamSpot>&, 
		   const edm::Handle<reco::ConversionCollection>&, 
		   const edm::Handle<reco::GsfElectronCollection>&, 
		   const edm::Handle<edm::ValueMap<double> >&, 
		   const edm::Handle<edm::ValueMap<double> >&, 
		   const edm::Handle<edm::ValueMap<double> >&, 
		   const edm::Handle<double>&, 
		   const double, const double, const double, 
		   const bool, 
		   const double, const bool, 
		   const double, const bool, 
		   const double, const bool, 
		   const double, const double, const bool, 
		   const double, const double, const bool);

  //fill STL container with photons passing specified cuts in specified eta range
  static std::vector<reco::PhotonRef> 
    getRecoPhotons(const edm::Handle<reco::PhotonRefVector>&, 
		   const edm::Handle<reco::PhotonCollection>&, 
		   const edm::Handle<reco::BeamSpot>&, 
		   const edm::Handle<reco::ConversionCollection>&, 
		   const edm::Handle<reco::GsfElectronCollection>&, 
		   const edm::Handle<edm::ValueMap<double> >&, 
		   const edm::Handle<edm::ValueMap<double> >&, 
		   const edm::Handle<edm::ValueMap<double> >&, 
		   const edm::Handle<double>&, 
		   const double, const double, const double, 
		   const bool, 
		   const double, const bool, 
		   const double, const bool, 
		   const double, const bool, 
		   const double, const double, const bool, 
		   const double, const double, const bool);

  //set canvas drawing options
  static void setCanvasOptions(TCanvas&, const Int_t, const Int_t, const Int_t);

  //set canvas margins
  static void setCanvasMargins(TCanvas&, const float, const float, const float, const float);

  //draw 1D histograms
  static void draw1DHistograms(TCanvas&, TH1F*);

  //draw 2D histograms
  static void draw2DHistograms(TCanvas&, TH2F*);

  //set axis drawing options
  static void setAxisOptions(TAxis*, const Float_t, const Float_t, const Float_t, const char*);

  //set axis drawing options excluding title
  static void setAxisOptions(TAxis*, const Float_t, const Float_t, const Float_t);

  //set axis labels
  static void setAxisLabels(TAxis*, const std::vector<std::string>&);

  //set graph drawing options
  static void setGraphOptions(TGraphAsymmErrors&, const Color_t, const Size_t, const Style_t, 
			      const char*, const char*);

  //set 1D histogram drawing options
  static void setHistogramOptions(TH1*, const Color_t, const Size_t, const Style_t, 
				  const Double_t, const char*, const char*, const Double_t);

  //set 1D histogram drawing options, excluding axis titles
  static void setHistogramOptions(TH1F*, const Color_t, const Size_t, const Style_t, 
				  const Double_t, const Double_t);

  //set 2D histogram drawing options
  static void setHistogramOptions(TH2F*, const Color_t, const Size_t, const Style_t, 
				  const Float_t, const Double_t, const char*, const char*);

  //set 2D histogram drawing options, excluding axis titles
  static void setHistogramOptions(TH2F*, const Color_t, const Size_t, const Style_t, 
				  const Float_t, const Double_t);

  //set legend options
  static void setLegendOptions(TLegend&, const char*);

  /*tauDecay function expects a status 3 tau, so if this is a status 1 light lepton from tau 
    decay, pass in the key of the status 3 tau ref*/
  static unsigned int getStatus3Key(const edm::Handle<reco::GenParticleRefVector>&, 
				    const edm::Handle<reco::GenParticleCollection>&, 
				    const unsigned int);

  //determine if a particle with a given PDG ID is an ancestor of the given particle
  static bool hasAncestor(const reco::GenParticleRef&, const std::vector<int>&);

  //get rho-corrected charged hadron isolation
  static double rhoCorrChargedHadronIso(const edm::Handle<edm::ValueMap<double> >&, 
					const edm::Handle<double>&, const reco::PhotonRef&);

  //get rho-corrected neutral hadron isolation
  static double rhoCorrNeutralHadronIso(const edm::Handle<edm::ValueMap<double> >&, 
					const edm::Handle<double>&, const reco::PhotonRef&);

  //get rho-corrected photon isolation
  static double rhoCorrPhotonIso(const edm::Handle<edm::ValueMap<double> >&, 
				 const edm::Handle<double>&, const reco::PhotonRef&);

  //decide if muon passes tight ID or not
  static bool isTightIsolatedRecoMuon(const reco::Muon* iMuon, 
				      const reco::Vertex* pPV, const bool usePFIso, 
				      const double PUSubtractionCoeff, const double isoMax, 
				      const double etaMax, const bool passIso);

  //decide if muon passes tight ID or not
  static bool isTightIsolatedRecoMuon(const edm::RefToBase<reco::Muon> iMuon, 
				      const reco::Vertex* pPV, const bool usePFIso, 
				      const double PUSubtractionCoeff, const double isoMax, 
				      const double etaMax, const bool passIso);
    
 private:

  static bool compareCandidatePT(reco::Candidate*, reco::Candidate*);

  static bool compareGenTauDecayIDPT(GenTauDecayID, GenTauDecayID);

  static bool comparePFJetPT(reco::PFJet*, reco::PFJet*);

  static bool comparePFTauPT(reco::PFTau*, reco::PFTau*);

  static bool compareMuonPT(reco::Muon*, reco::Muon*);

  static bool compareMuonRefPT(reco::MuonRef, reco::MuonRef);

  static bool comparePFJetRefPT(reco::PFJetRef, reco::PFJetRef);

  static bool comparePFTauRefPT(reco::PFTauRef, reco::PFTauRef);

  static bool comparePhotonRefPT(reco::PhotonRef, reco::PhotonRef);

  template<typename T>
    bool comparePT(edm::Ref<T> object1, edm::Ref<T> object2)
    {
      return (object1->pt() < object2->pt());
    }

  static std::vector<reco::MuonRef> 
    getTightIsolatedRecoMuons(const edm::Handle<reco::MuonCollection>&, const reco::Vertex*, 
			      const bool, const double, const double, const double, const bool);

  static std::vector<reco::MuonRef> 
    getTightIsolatedRecoMuons(const edm::Handle<reco::MuonRefVector>&, 
			      const edm::Handle<reco::MuonCollection>&, const reco::Vertex*, 
			      const bool, const double, const double, const double, const bool);

  static double rhoCorrIso(const edm::Handle<edm::ValueMap<double> >&, 
			   const edm::Handle<double>&, const reco::PhotonRef&, 
			   const std::map<double, double>&);

  static const std::vector<double> absEtaBinEdges;

  static const std::map<double, double> chargedHadronIsoEAs;

  static const std::map<double, double> neutralHadronIsoEAs;

  static const std::map<double, double> photonIsoEAs;

};

#endif
