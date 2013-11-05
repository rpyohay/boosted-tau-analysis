#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

typedef SingleObjectSelector<
  reco::PFJetCollection, 
  StringCutObjectSelector<reco::PFJet>,
  reco::PFJetRefVector
  > PFJetRefSelector;

DEFINE_FWK_MODULE(PFJetRefSelector);
