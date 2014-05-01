import FWCore.ParameterSet.Config as cms

looseEBPhotonSelector = cms.EDFilter(
    'CustomPhotonSelector',
    basePhotonTag = cms.InputTag('photons'),
    beamSpotTag = cms.InputTag('offlineBeamSpot'),
    conversionTag = cms.InputTag('allConversions'),
    electronTag = cms.InputTag('gsfElectrons'),
    chargedHadronIsoTag = cms.InputTag('phPFIsoValueCharged03PFIdPFIso'),
    neutralHadronIsoTag = cms.InputTag('phPFIsoValueNeutral03PFIdPFIso'),
    photonIsoTag = cms.InputTag('phPFIsoValueGamma03PFIdPFIso'),
    rhoTag = cms.InputTag('kt6PFJets', 'rho', 'RECO'),
    pTMin = cms.double(31.0), #GeV
    etaMin = cms.double(0.0), #EB only
    etaMax = cms.double(1.4442), #EB only
    passCSEV = cms.bool(True),
    HOverEMax = cms.double(0.05),
    passHOverE = cms.bool(True),
    sigmaIetaIetaMax = cms.double(0.012),
    passSigmaIetaIeta = cms.bool(True),
    chargedHadronIsoMax = cms.double(2.6), #GeV
    passChargedHadronIso = cms.bool(True),
    neutralHadronIsoConst = cms.double(3.5), #GeV
    neutralHadronIsoPTMultiplier = cms.double(0.04),
    passNeutralHadronIso = cms.bool(True),
    photonIsoConst = cms.double(1.3), #GeV
    photonIsoPTMultiplier = cms.double(0.005),
    passPhotonIso = cms.bool(True),
    minNumObjsToPassFilter = cms.uint32(1)
    )

mediumEBPhotonSelector = looseEBPhotonSelector.clone()
mediumEBPhotonSelector.sigmaIetaIetaMax = cms.double(0.011)
mediumEBPhotonSelector.chargedHadronIsoMax = cms.double(1.5) #GeV
mediumEBPhotonSelector.neutralHadronIsoConst = cms.double(1.0) #GeV
mediumEBPhotonSelector.photonIsoConst = cms.double(0.7) #GeV

tightEBPhotonSelector = mediumEBPhotonSelector.clone()
tightEBPhotonSelector.chargedHadronIsoMax = cms.double(0.7) #GeV
tightEBPhotonSelector.neutralHadronIsoConst = cms.double(0.4) #GeV
tightEBPhotonSelector.photonIsoConst = cms.double(0.5) #GeV

looseEEPhotonSelector = looseEBPhotonSelector.clone()
looseEEPhotonSelector.etaMin = cms.double(1.566) #EE only
looseEEPhotonSelector.etaMax = cms.double(2.5) #EE only
looseEEPhotonSelector.sigmaIetaIetaMax = cms.double(0.034)
looseEEPhotonSelector.chargedHadronIsoMax = cms.double(2.3) #GeV
looseEEPhotonSelector.neutralHadronIsoConst = cms.double(2.9) #GeV
looseEEPhotonSelector.photonIsoConst = cms.double(-1.0)

mediumEEPhotonSelector = looseEEPhotonSelector.clone()
mediumEEPhotonSelector.sigmaIetaIetaMax = cms.double(0.033)
mediumEEPhotonSelector.chargedHadronIsoMax = cms.double(1.2) #GeV
mediumEEPhotonSelector.neutralHadronIsoConst = cms.double(1.5) #GeV
mediumEEPhotonSelector.photonIsoConst = cms.double(1.0)

tightEEPhotonSelector = mediumEEPhotonSelector.clone()
tightEEPhotonSelector.sigmaIetaIetaMax = cms.double(0.031)
tightEEPhotonSelector.chargedHadronIsoMax = cms.double(0.5) #GeV
