#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

typedef SingleObjectSelector<
  reco::CandidateView, 
  StringCutObjectSelector<reco::Candidate>,
  reco::CandidateCollection
  > CandidateViewCloneSelector;

DEFINE_FWK_MODULE(CandidateViewCloneSelector);
