#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

typedef SingleObjectSelector<
  edm::View<reco::PFJet>, 
  StringCutObjectSelector<reco::PFJet>,
  reco::PFJetCollection
  > PFJetViewCloneSelector;

DEFINE_FWK_MODULE(PFJetViewCloneSelector);
