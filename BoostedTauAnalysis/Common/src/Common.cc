#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TH1F.h"
#include "TLegend.h"

void minDR::setCandidate(const reco::Candidate* cand)
{
  cand_ = const_cast<reco::Candidate*>(cand);
}

void minDR::deleteCandidate() { cand_ = NULL; }

reco::Candidate* minDR::getCandidate() const { return cand_; }

void maxM::setMuonJetMap(const edm::Handle<edm::ValueMap<reco::MuonRefVector> >* pMuonJetMap)
{
  pMuonJetMap_ = const_cast<edm::Handle<edm::ValueMap<reco::MuonRefVector> >*>(pMuonJetMap);
}

void maxM::deleteMuonJetMap() { pMuonJetMap_ = NULL; }

edm::Handle<edm::ValueMap<reco::MuonRefVector> >* maxM::getMuonJetMap() const
{
  return pMuonJetMap_;
}

reco::MuonRef maxM::highestPTMu(const reco::PFTauRef& tau) const
{
  const reco::PFJetRef& tauJetRef = tau->jetRef();
  const reco::MuonRefVector& removedMuons = (**pMuonJetMap_)[tauJetRef];
  std::vector<reco::MuonRef> removedMuonRefs;
  for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
       iMuon != removedMuons.end(); ++iMuon) { removedMuonRefs.push_back(*iMuon); }
  Common::sortByPT(removedMuonRefs);
  std::reverse(removedMuonRefs.begin(), removedMuonRefs.end());
  if (removedMuonRefs.size() > 0) return removedMuonRefs[0];
  else return reco::MuonRef();
}

bool maxM::operator()(const reco::PFTauRef& tau1, const reco::PFTauRef& tau2) const
{
  reco::MuonRef highestPTMu1 = highestPTMu(tau1);
  reco::MuonRef highestPTMu2 = highestPTMu(tau2);
  bool retVal = false;
  if (highestPTMu1.isNull() && highestPTMu2.isNonnull()) retVal = false;
  if (highestPTMu2.isNull() && highestPTMu1.isNonnull()) retVal = true;
  if (highestPTMu1.isNull() && highestPTMu2.isNull()) retVal = true;
  if (highestPTMu1.isNonnull() && highestPTMu2.isNonnull()) {
    retVal = ((highestPTMu1->p4() + tau1->p4()).M() > (highestPTMu2->p4() + tau2->p4()).M());
  }
  return retVal;
}

void Common::sortByPT(std::vector<reco::Candidate*>& objects)
{
  std::sort(objects.begin(), objects.end(), compareCandidatePT);
}

void Common::sortByPT(std::vector<GenTauDecayID>& objects)
{
  std::sort(objects.begin(), objects.end(), compareGenTauDecayIDPT);
}

void Common::sortByPT(std::vector<reco::PFJet*>& objects)
{
  std::sort(objects.begin(), objects.end(), comparePFJetPT);
}

void Common::sortByPT(std::vector<reco::PFTau*>& objects)
{
  std::sort(objects.begin(), objects.end(), comparePFTauPT);
}

void Common::sortByPT(std::vector<reco::Muon*>& objects)
{
  std::sort(objects.begin(), objects.end(), compareMuonPT);
}

void Common::sortByPT(std::vector<reco::MuonRef>& objects)
{
  std::sort(objects.begin(), objects.end(), compareMuonRefPT);
}

void Common::sortByPT(std::vector<reco::PFJetRef>& objects)
{
  std::sort(objects.begin(), objects.end(), comparePFJetRefPT);
}

void Common::sortByPT(std::vector<reco::PFTauRef>& objects)
{
  std::sort(objects.begin(), objects.end(), comparePFTauRefPT);
}

void Common::sortByMass(const edm::Handle<edm::ValueMap<reco::MuonRefVector> >& muonJetMap, 
			std::vector<reco::PFTauRef>& taus)
{
  maxM comp;
  comp.setMuonJetMap(&muonJetMap);
  sort(taus.begin(), taus.end(), comp);
  comp.deleteMuonJetMap();
}

bool Common::compareCandidatePT(reco::Candidate* object1, reco::Candidate* object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::compareGenTauDecayIDPT(GenTauDecayID object1, GenTauDecayID object2)
{
  return (object1.getVisibleTauP4().Pt() < object2.getVisibleTauP4().Pt());
}

bool Common::comparePFJetPT(reco::PFJet* object1, reco::PFJet* object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::comparePFTauPT(reco::PFTau* object1, reco::PFTau* object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::compareMuonPT(reco::Muon* object1, reco::Muon* object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::compareMuonRefPT(reco::MuonRef object1, reco::MuonRef object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::comparePFJetRefPT(reco::PFJetRef object1, reco::PFJetRef object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::comparePFTauRefPT(reco::PFTauRef object1, reco::PFTauRef object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::isGoodVertex(const reco::Vertex* pVtx)
{
  bool retVal = false;
  if (!pVtx->isFake() && 
      (pVtx->ndof() > 4) && 
      (fabs(pVtx->x()) <= 24.0/*cm*/) && 
      (fabs(pVtx->position().Rho()) <= 2.0/*cm*/)) retVal = true;
  return retVal;
}

reco::Vertex* Common::getPrimaryVertex(edm::Handle<reco::VertexCollection>& pVertices)
{
  reco::VertexCollection::const_iterator iVtx = pVertices->begin();
  reco::Vertex* pPV = NULL;
  while ((iVtx != pVertices->end()) && (pPV == NULL)) {
    if (isGoodVertex(&*iVtx)) pPV = const_cast<reco::Vertex*>(&*iVtx);
    ++iVtx;
  }
  return pPV;
}

unsigned int Common::numGoodVertices(edm::Handle<reco::VertexCollection>& pVertices)
{
  unsigned int nGoodVtx = 0;
  for (reco::VertexCollection::const_iterator iVtx = pVertices->begin(); 
       iVtx != pVertices->end(); ++iVtx) {
    if (isGoodVertex(&*iVtx)) ++nGoodVtx;
  }
  return nGoodVtx;
}

float Common::getMuonCombPFIso(const reco::Muon& muon, const double PUSubtractionCoeff)
{
  const reco::MuonPFIsolation isoBlock = muon.pfIsolationR04();
  return (isoBlock.sumChargedHadronPt + 
	  std::max(0.0, (double)(isoBlock.sumNeutralHadronEt + isoBlock.sumPhotonEt - 
				 PUSubtractionCoeff*isoBlock.sumPUPt)));
}

std::vector<reco::MuonRef>
Common::getTightPFIsolatedRecoMuons(const edm::Handle<reco::MuonCollection>& pMuons, 
				    const reco::Vertex* pPV, const double PUSubtractionCoeff, 
				    const double PFIsoMax, const double etaMax, const bool passIso)
{
  return getTightIsolatedRecoMuons(pMuons, pPV, true, PUSubtractionCoeff, PFIsoMax, etaMax, 
				   passIso);
}

std::vector<reco::MuonRef>
Common::getTightPFIsolatedRecoMuons(const edm::Handle<reco::MuonRefVector>& pMuons, 
				    const edm::Handle<reco::MuonCollection>& pBaseMuons, 
				    const reco::Vertex* pPV, const double PUSubtractionCoeff, 
				    const double PFIsoMax, const double etaMax, const bool passIso)
{
  return getTightIsolatedRecoMuons(pMuons, pBaseMuons, pPV, true, PUSubtractionCoeff, PFIsoMax, 
				   etaMax, passIso);
}

std::vector<reco::MuonRef>
Common::getTightDetectorIsolatedRecoMuons(const edm::Handle<reco::MuonCollection>& pMuons, 
					  const reco::Vertex* pPV, const double detectorIsoMax, 
					  const double etaMax, const bool passIso)
{
  return getTightIsolatedRecoMuons(pMuons, pPV, false, 0.0, detectorIsoMax, etaMax, passIso);
}

std::vector<reco::MuonRef>
Common::getTightDetectorIsolatedRecoMuons(const edm::Handle<reco::MuonRefVector>& pMuons, 
					  const edm::Handle<reco::MuonCollection>& pBaseMuons, 
					  const reco::Vertex* pPV, const double detectorIsoMax, 
					  const double etaMax, const bool passIso)
{
  return getTightIsolatedRecoMuons(pMuons, pBaseMuons, pPV, false, 0.0, detectorIsoMax, etaMax, 
				   passIso);
}

std::vector<reco::MuonRef>
Common::getSoftRecoMuons(const edm::Handle<reco::MuonCollection>& pMuons, const reco::Vertex* pPV, 
			 const double etaMax)
{
  std::vector<reco::MuonRef> softMuons;
  for (reco::MuonCollection::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) {
    const reco::TrackRef innerTrack = iMuon->innerTrack();
    const reco::Vertex::Point PVPos = pPV->position();
    if (muon::isGoodMuon(*iMuon, muon::TMOneStationTight) && 
	(iMuon->track()->hitPattern().trackerLayersWithMeasurement() > 5) && 
	(innerTrack->hitPattern().pixelLayersWithMeasurement() > 1) && 
	(innerTrack->normalizedChi2() < 1.8) && 
	(fabs(innerTrack->dxy(PVPos)) < 3.0) && 
	(fabs(innerTrack->dz(PVPos)) < 30.0) && 
	((etaMax == -1.0) || (fabs(iMuon->eta()) < etaMax))) {
      softMuons.push_back(reco::MuonRef(pMuons, iMuon - pMuons->begin()));
    }
  }
  return softMuons;
}

std::vector<reco::MuonRef>
Common::getSoftRecoMuons(const edm::Handle<reco::MuonRefVector>& pMuons, 
			 const edm::Handle<reco::MuonCollection>& pBaseMuons, 
			 const reco::Vertex* pPV, const double etaMax)
{
  std::vector<reco::MuonRef> softMuons;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) {
    const reco::TrackRef innerTrack = (*iMuon)->innerTrack();
    const reco::Vertex::Point PVPos = pPV->position();
    if (muon::isGoodMuon(**iMuon, muon::TMOneStationTight) && 
	((*iMuon)->track()->hitPattern().trackerLayersWithMeasurement() > 5) && 
	(innerTrack->hitPattern().pixelLayersWithMeasurement() > 1) && 
	(innerTrack->normalizedChi2() < 1.8) && 
	(fabs(innerTrack->dxy(PVPos)) < 3.0) && 
	(fabs(innerTrack->dz(PVPos)) < 30.0) && 
	((etaMax == -1.0) || (fabs((*iMuon)->eta()) < etaMax))) {
      softMuons.push_back(reco::MuonRef(pBaseMuons, iMuon->key()));
    }
  }
  return softMuons;
}

std::vector<reco::PFTauRef> 
Common::getRecoTaus(const edm::Handle<reco::PFTauCollection>& pTaus, 
		    const std::vector<edm::Handle<reco::PFTauDiscriminator> >& pTauDiscriminators, 
		    const double pTMin, const double etaMax, const bool passIso)
{
  std::vector<reco::PFTauRef> taus;
  for (reco::PFTauCollection::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau) {
    reco::PFTauRef tauRef(pTaus, iTau - pTaus->begin());
    bool passTauDiscriminators = true;
    std::vector<edm::Handle<reco::PFTauDiscriminator> >::const_iterator iDiscriminator = 
      pTauDiscriminators.begin();
    while ((iDiscriminator != pTauDiscriminators.end()) && passTauDiscriminators) {
      if ((**iDiscriminator)[tauRef] != 1.0) passTauDiscriminators = false;
      ++iDiscriminator;
    }
    if (((passIso && passTauDiscriminators) || (!passIso && !passTauDiscriminators)) && 
	((etaMax == -1.0) || (fabs(iTau->eta()) < etaMax)) && 
	((pTMin == -1.0) || (iTau->pt() > pTMin))) {
      taus.push_back(tauRef);
    }
  }
  return taus;
}

std::vector<reco::PFTauRef> 
Common::getRecoTaus(const edm::Handle<reco::PFTauRefVector>& pTaus, 
		    const edm::Handle<reco::PFTauCollection>& pBaseTaus, 
		    const std::vector<edm::Handle<reco::PFTauDiscriminator> >& pTauDiscriminators, 
		    const double pTMin, const double etaMax, const bool passIso)
{
  std::vector<reco::PFTauRef> taus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau) {
    reco::PFTauRef tauRef(pBaseTaus, iTau->key());
    bool passTauDiscriminators = true;
    std::vector<edm::Handle<reco::PFTauDiscriminator> >::const_iterator iDiscriminator = 
      pTauDiscriminators.begin();
    while ((iDiscriminator != pTauDiscriminators.end()) && passTauDiscriminators) {
      if ((**iDiscriminator)[tauRef] != 1.0) passTauDiscriminators = false;
      ++iDiscriminator;
    }
    if (((passIso && passTauDiscriminators) || (!passIso && !passTauDiscriminators)) && 
	((etaMax == -1.0) || (fabs((*iTau)->eta()) < etaMax)) && 
	((pTMin == -1.0) || ((*iTau)->pt() > pTMin))) {
      taus.push_back(tauRef);
    }
  }
  return taus;
}

void Common::setCanvasOptions(TCanvas& canvas, const Int_t grid, const Int_t logY, 
			      const Int_t logZ)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
}

void Common::setCanvasMargins(TCanvas& canvas, const float left, const float top, 
			      const float right, const float bottom)
{
  canvas.cd()->SetLeftMargin(left);
  canvas.cd()->SetTopMargin(top);
  canvas.cd()->SetRightMargin(right);
  canvas.cd()->SetBottomMargin(bottom);
}

void Common::draw1DHistograms(TCanvas& canvas, TH1F* hist)
{
  setCanvasOptions(canvas, 1, 0, 0);
  setHistogramOptions(hist, kBlack, 0.7, 20, 1.0, 0.05);
  hist->SetLineWidth(2);
  canvas.cd();
  hist->Draw();
}

void Common::draw2DHistograms(TCanvas& canvas, TH2F* hist)
{
  setCanvasOptions(canvas, 0, 0, 0);
  setCanvasMargins(canvas, 0.2, 0.2, 0.2, 0.2);
  setHistogramOptions(hist, kBlack, 0.7/*4.2*/, 20, 1.6, 1.0);
  canvas.cd();
  hist->Draw("COLZ");
}

void Common::setAxisOptions(TAxis* axis, const Float_t labelSize, const Float_t titleSize, 
			    const Float_t titleOffset, const char* title)
{
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetTitleSize(titleSize);
  axis->SetTitleOffset(titleOffset);
  axis->SetTitle(title);
}

void Common::setAxisOptions(TAxis* axis, const Float_t labelSize, const Float_t titleSize, 
			    const Float_t titleOffset)
{
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetTitleSize(titleSize);
  axis->SetTitleOffset(titleOffset);
}

void Common::setAxisLabels(TAxis* axis, const std::vector<std::string>& binLabels)
{
  for (Int_t iBin = 1; iBin <= axis->GetNbins(); ++iBin) {
    if (iBin <= (int)binLabels.size()) axis->SetBinLabel(iBin, binLabels[iBin - 1].c_str());
  }
}

void Common::setGraphOptions(TGraphAsymmErrors& graph, const Color_t color, const Size_t size, 
			     const Style_t style, const char* xAxisTitle, const char* yAxisTitle)
{
  graph.SetMarkerColor(color);
  graph.SetMarkerSize(size);
  graph.SetMarkerStyle(style);
  graph.SetLineColor(color);
  graph.SetTitle("");
  setAxisOptions(graph.GetXaxis(), 0.05, 0.05, 0.9, xAxisTitle);
  setAxisOptions(graph.GetYaxis(), 0.05, 0.05, 1.05, yAxisTitle);
}

void Common::setHistogramOptions(TH1F* histogram, const Color_t color, const Size_t size, 
				 const Style_t style, const Double_t scale, 
				 const char* xAxisTitle, const char* yAxisTitle, 
				 const Double_t xAxisLabelSize)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), xAxisLabelSize, 0.05, 0.9, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), 0.05, 0.05, 1.05, yAxisTitle);
  histogram->Scale(scale);
}

void Common::setHistogramOptions(TH1F* histogram, const Color_t color, const Size_t size, 
				 const Style_t style, const Double_t scale, 
				 const Double_t xAxisLabelSize)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), xAxisLabelSize, 0.05, 0.9);
  setAxisOptions(histogram->GetYaxis(), 0.05, 0.05, 1.05);
  histogram->Scale(scale);
}

void Common::setHistogramOptions(TH2F* histogram, const Color_t color, const Size_t size, 
				 const Style_t style, const Float_t yAxisTitleOffset, 
				 const Double_t scale, const char* xAxisTitle, 
				 const char* yAxisTitle)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), 0.04, 0.04, 1.1, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), 0.04, 0.04, yAxisTitleOffset, yAxisTitle);
  setAxisOptions(histogram->GetZaxis(), 0.04, 0.04, histogram->GetZaxis()->GetTitleOffset(), "");
  histogram->Scale(scale);
}

void Common::setHistogramOptions(TH2F* histogram, const Color_t color, const Size_t size, 
				 const Style_t style, const Float_t yAxisTitleOffset, 
				 const Double_t scale)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), 0.04, 0.04, 1.1);
  setAxisOptions(histogram->GetYaxis(), 0.04, 0.04, yAxisTitleOffset);
  setAxisOptions(histogram->GetZaxis(), 0.04, 0.04, histogram->GetZaxis()->GetTitleOffset(), "");
  histogram->Scale(scale);
}

void Common::setLegendOptions(TLegend& legend, const char* header)
{
  legend.SetFillColor(0);
  legend.SetTextFont(42);
  legend.SetHeader(header);
}

unsigned int
Common::getStatus3Key(const edm::Handle<reco::GenParticleRefVector>& pSelectedGenParticles, 
		      const edm::Handle<reco::GenParticleCollection>& pGenParticles, 
		      const unsigned int i)
{
  reco::GenParticleRef genParticleRef = pSelectedGenParticles->at(i);
  unsigned int key = genParticleRef.key();
  if (genParticleRef->status() == 1) {
    key = reco::GenParticleRef(pGenParticles, key)->motherRef()->motherRef().key();
  }
  return key;
}

std::vector<reco::MuonRef>
Common::getTightIsolatedRecoMuons(const edm::Handle<reco::MuonCollection>& pMuons, 
				  const reco::Vertex* pPV, const bool usePFIso, 
				  const double PUSubtractionCoeff, const double isoMax, 
				  const double etaMax, const bool passIso)
{
  std::vector<reco::MuonRef> tightMuons;
  for (reco::MuonCollection::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) {
    if ((pPV != NULL) && 
	muon::isTightMuon(*iMuon, *pPV) && 
	iMuon->isPFMuon() && 
	(fabs(iMuon->innerTrack()->dz(pPV->position())) < 0.5) && 
	(iMuon->track()->hitPattern().trackerLayersWithMeasurement() > 5) && 
	((etaMax == -1.0) || (fabs(iMuon->eta()) < etaMax))) {
      float iso = 0.0;
      if (usePFIso) {
	iso = getMuonCombPFIso(*iMuon, PUSubtractionCoeff)/iMuon->pt();
      }
      else {
	const reco::MuonIsolation isoBlock = iMuon->isolationR03();
	iso = (isoBlock.sumPt + isoBlock.emEt + isoBlock.hadEt)/iMuon->pt();
      }
      if ((isoMax == -1.0) || (passIso && (iso < isoMax)) || (!passIso && (iso >= isoMax))) {
	tightMuons.push_back(reco::MuonRef(pMuons, iMuon - pMuons->begin()));
      }
    }
  }
  return tightMuons;
}

std::vector<reco::MuonRef>
Common::getTightIsolatedRecoMuons(const edm::Handle<reco::MuonRefVector>& pMuons, 
				  const edm::Handle<reco::MuonCollection>& pBaseMuons, 
				  const reco::Vertex* pPV, const bool usePFIso, 
				  const double PUSubtractionCoeff, const double isoMax, 
				  const double etaMax, const bool passIso)
{
  std::vector<reco::MuonRef> tightMuons;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) {
    if ((pPV != NULL) && 
	muon::isTightMuon(**iMuon, *pPV) && 
	(*iMuon)->isPFMuon() && 
	(fabs((*iMuon)->innerTrack()->dz(pPV->position())) < 0.5) && 
	((*iMuon)->track()->hitPattern().trackerLayersWithMeasurement() > 5) && 
	((etaMax == -1.0) || (fabs((*iMuon)->eta()) < etaMax))) {
      float iso = 0.0;
      if (usePFIso) {
	iso = getMuonCombPFIso(**iMuon, PUSubtractionCoeff)/(*iMuon)->pt();
      }
      else {
	const reco::MuonIsolation isoBlock = (*iMuon)->isolationR03();
	iso = (isoBlock.sumPt + isoBlock.emEt + isoBlock.hadEt)/(*iMuon)->pt();
      }
      if ((isoMax == -1.0) || (passIso && (iso < isoMax)) || (!passIso && (iso >= isoMax))) {
	tightMuons.push_back(reco::MuonRef(pBaseMuons, iMuon->key()));
      }
    }
  }
  return tightMuons;
}
