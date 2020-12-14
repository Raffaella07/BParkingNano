import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

electronPair = cms.EDProducer(
    'DiElectronBuilder',
    src = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    transientTracksSrc = cms.InputTag('electronsForAnalysis', 'SelectedTransientElectrons'),
    lep1Selection = cms.string('pt > 1.3'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string(
        
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1.  '
        '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03 && userInt("nlowpt")<2'
        
    ),
    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'),
)

electronPairsForKee = cms.EDProducer(
    'DiElectronBuilder',
    src = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    transientTracksSrc = cms.InputTag('electronsForAnalysis', 'SelectedTransientElectrons'),
    lep1Selection = cms.string('pt > 1.3'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string(
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1. &&  mass() < 5  '
        '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03 && userInt("nlowpt")<2'
        
    ),
    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'),
)
BToKee = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('electronPairsForKee'),
    leptonTransientTracks = electronPairsForKee.transientTracksSrc,
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    kaonSelection = cms.string(''),
    isoTracksSelection = cms.string('pt > 0.5 && abs(eta)<2.5'),
    preVtxSelection = cms.string(
        'pt > 1.75 && userFloat("min_dr") > 0.03 '
        '&& mass < 5.4 && mass > 4.8'
        ),
    postVtxSelection = cms.string(
         ''
    )
)

muonPairsForKmumu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 1.5'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
                                 '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'),
    postVtxSelection = electronPairsForKee.postVtxSelection,
)

BToKmumu = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('muonPairsForKmumu'),
    leptonTransientTracks = muonPairsForKmumu.transientTracksSrc,
    kaons = BToKee.kaons,
    kaonsTransientTracks = BToKee.kaonsTransientTracks,
    beamSpot = cms.InputTag("offlineBeamSpot"),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    kaonSelection = cms.string(''),
    isoTracksSelection = BToKee.isoTracksSelection,
    # This in principle can be different between electrons and muons
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& mass < 7. && mass > 4.'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.'
    )
)

electronPairTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("electronPair"),
    cut = cms.string(""),
    name = cms.string("electronPair"),
    doc = cms.string("electronPair Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),
       # vtx_ez = ufloat('vtx_ez'),
        l1_SC_rawEnergy = ufloat('l1_SC_Eraw'),
        l2_SC_rawEnergy = ufloat('l2_SC_Eraw'),
        l1_SC_Energy = ufloat('l1_SC_E'),
        l2_SC_Energy = ufloat('l2_SC_E'),
        l1_SC_RegEnergy = ufloat('l1_SC_Ereg'),
        l2_SC_RegEnergy = ufloat('l2_SC_Ereg'),
        l1_SCeta = ufloat('l1_SC_eta'),
        l2_SCeta= ufloat('l2_SC_eta'),
        l1_SCphi = ufloat('l1_SC_phi'),
        l2_SCphi = ufloat('l2_SC_phi'),
        #l1_SC_EOverP = ufloat("l1_Esc/p_track")',float),
        #l1_SC_EOverP = ufloat('l1_Esc/p_track'),
       # l2_SC_EOverP = ufloat("l2_Esc/p_track")',float),
        # fit and vtx info
        #chi2 = ufloat('sv_chi2'),
        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = ufloat('fitted_mass'), # this might not work
      #  mllErr_llfit = Var('userCand("dilepton").userFloat("fitted_massErr")', float), # this might not work
     #   mll_fullfit = ufloat('fitted_mll'),
        # Cos(theta)
        #cos2D = ufloat('cos_theta_2D'),
       # fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
      # fit_mass = ufloat('fitted_mass'),
      # fit_massErr = ufloat('fitted_massErr'),
      # fit_pt = ufloat('fitted_pt'),
      # fit_eta = ufloat('fitted_eta'),
      # fit_phi = ufloat('fitted_phi'),
      # fit_l1_pt = ufloat('fitted_l1_pt'),
      # fit_l1_eta = ufloat('fitted_l1_eta'),
      # fit_l1_phi = ufloat('fitted_l1_phi'),
      # fit_l2_pt = ufloat('fitted_l2_pt'),
      # fit_l2_eta = ufloat('fitted_l2_eta'),
      # fit_l2_phi = ufloat('fitted_l2_phi'),
      # fit_k_pt = ufloat('fitted_k_pt'),
      # fit_k_eta = ufloat('fitted_k_eta'),
      # fit_k_phi = ufloat('fitted_k_phi'),
   #    l1_iso03 = ufloat('l1_iso03'),
   #    l1_iso04 = ufloat('l1_iso04'),
   #    l2_iso03 = ufloat('l2_iso03'),
   #    l2_iso04 = ufloat('l2_iso04'),
   #    n_l1_used = uint('n_l1_used'),
   #    n_l2_used = uint('n_l2_used'),
    )
)

BToKeeTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToKee"),
    cut = cms.string(""),
    name = cms.string("electronPair"),
    doc = cms.string("electronPair from B Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
#	Bmass = ufloat('mass'),
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),
       # vtx_ez = ufloat('vtx_ez'),
        l1_SC_rawEnergy = Var('userCand("dilepton").userFloat("l1_SC_Eraw")',float),
        l2_SC_rawEnergy = Var('userCand("dilepton").userFloat("l2_SC_Eraw")',float),
        l1_SC_Energy = Var('userCand("dilepton").userFloat("l1_SC_E")',float),
        l2_SC_Energy = Var('userCand("dilepton").userFloat("l2_SC_E")',float),
        l1_SC_RegEnergy = Var('userCand("dilepton").userFloat("l1_SC_Ereg")',float),
        l2_SC_RegEnergy = Var('userCand("dilepton").userFloat("l2_SC_Ereg")',float),
        l1_SCeta = Var('userCand("dilepton").userFloat("l1_SC_eta")',float),
        l2_SCeta= Var('userCand("dilepton").userFloat("l2_SC_eta")',float),
        l1_SCphi = Var('userCand("dilepton").userFloat("l1_SC_phi")',float),
        l2_SCphi = Var('userCand("dilepton").userFloat("l2_SC_phi")',float),
        l1_SmearE = Var('userCand("dilepton").userFloat("l1_SmearedE")',float),
        l2_SmearE = Var('userCand("dilepton").userFloat("l2_SmearedE")',float),
        l1_Smearing = Var('userCand("dilepton").userFloat("l1_Smearing")',float),
        l2_Smearing = Var('userCand("dilepton").userFloat("l2_Smearing")',float),
        l1_Scale = Var('userCand("dilepton").userFloat("l1_Scale")',float),
        l2_Scale = Var('userCand("dilepton").userFloat("l2_Scale")',float),
        #l1_SC_EOverP = ufloat("l1_Esc/p_track")',float),
        #l1_SC_EOverP = ufloat('l1_Esc/p_track'),
       # l2_SC_EOverP = ufloat("l2_Esc/p_track")',float),
        # fit and vtx info
        #chi2 = ufloat('sv_chi2'),
        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")',float) # this might not work
    )
)

BToKmumuTable = BToKeeTable.clone(
    src = cms.InputTag("BToKmumu"),
    name = cms.string("BToKMuMu"),
    doc = cms.string("BToKMuMu Variable")
)


CountBToKee = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToKee")
)    
CountBToKmumu = CountBToKee.clone(
    minNumber = cms.uint32(1),
    src = cms.InputTag("BToKmumu")
)

CountElePair = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("electronPair")
)
BToKMuMuSequence = cms.Sequence(
    (muonPairsForKmumu * BToKmumu)
)
BToKEESequence = cms.Sequence(
    (electronPairsForKee * BToKee)
)

BToKLLSequence = cms.Sequence(
    (electronPairsForKee * BToKee) 
#    (muonPairsForKmumu * BToKmumu)
)
BToKLLTables = cms.Sequence(BToKeeTable)# + BToKmumuTable)

electronPairSequence = cms.Sequence(
    electronPair
)

electronPairTables = cms.Sequence(electronPairTable)
