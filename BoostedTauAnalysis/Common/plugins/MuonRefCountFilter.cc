#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/ObjectCountFilter.h"
#include "CommonTools/UtilAlgos/interface/MinNumberSelector.h"
#include "CommonTools/UtilAlgos/interface/MaxNumberSelector.h"
#include "CommonTools/UtilAlgos/interface/AndSelector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

typedef ObjectCountFilter<
  reco::MuonRefVector, 
  AnySelector, 
  AndSelector<MinNumberSelector, MaxNumberSelector> 
  >::type MuonRefCountFilter;

DEFINE_FWK_MODULE(MuonRefCountFilter);
