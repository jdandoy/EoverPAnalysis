from calculation import Calculation, CalculationDataMC
from variables import trkEtaPhiECAL
import numpy as np

def NoSelection(trk):
    return np.ones(len(trk)) > 0.0
branches = []
sel_NoSelection = Calculation(NoSelection, branches)

def hasHADExtrapolation(trk):
    return (trk["trk_phiHEC1"] > -100) | (trk["trk_phiTileBar2"] > -100) | (trk["trk_phiTileExt1"] > -100)
branches =["trk_phiHEC1", "trk_phiTileBar2", "trk_phiTileExt1"]
sel_hasHADExtrapolation = Calculation(hasHADExtrapolation, branches)

#A selection for plotting event-based variables
def Event(trk):
    '''get only the variables associated with the first track in the event.
    This is so that we can reweight variables with an event, and not with individual tracks'''
    return (trk["trkIndex"] < 0.5) & (trk["trkIndex"] > -0.5)
branches = ["trkIndex"]
sel_Event = Calculation(Event, branches)

def SubleadingTrack(trk):
    '''get only the variables associated with the first track in the event.
    This is so that we can reweight variables with an event, and not with individual tracks'''
    return (trk["trkIndex"] < 1.5) & (trk["trkIndex"] > 0.5)
branches = ["trkIndex"]
sel_SubleadingTrack = Calculation(SubleadingTrack, branches)

#These are the definitions of the different types of track classifications according to tracking CP
def TruthLink(trk):
    return trk["trk_hasTruthParticle"] == 1

def TruthMatched(trk):
    return (trk["trk_truthProb"] > 0.5) 

def Fake(trk):
    return np.logical_not(TruthMatched(trk))
branches = ["trk_truthProb"]
sel_Fake = Calculation(Fake, branches)

def Pileup(trk):
    return TruthMatched(trk) & np.logical_not(TruthLink(trk))
branches = ["trk_truthProb", "trk_hasTruthParticle"]
sel_Pileup = Calculation(Pileup, branches)

def HardScatter(trk):
    return TruthMatched(trk) & TruthLink(trk)
branches = ["trk_truthProb", "trk_hasTruthParticle"]
sel_HardScatter = Calculation(HardScatter, branches)

def ParticlePDGID_ABS(trk, pdgID):
    return ((np.abs(trk["trk_truthPdgId"]) - abs(pdgID)) < 0.1)

def ParticlePDGID(trk, pdgID):
    return (np.abs(trk["trk_truthPdgId"] - pdgID)) < 0.1

#Define a significant hadron energy deposit with the amount of energy in the HAD calorimeter
def EHadBetween30And90OfMomentum(trk):
    '''At least 30 % of the track momentum and no more than 90 % of the track momentum was found the HAD calorimeter'''
    E_HAD_frac = trk["trk_ClusterEnergy_HAD_200"]/trk["trk_p"]
    return (E_HAD_frac > 0.3) & (E_HAD_frac < 0.9) #This selection only works for EM-scale
branches = ["trk_ClusterEnergy_HAD_200", "trk_p"]
sel_EHadBetween30And90OfMomentum = Calculation(EHadBetween30And90OfMomentum, branches)

def Lar1_1GeV(trk):
    return trk["trk_ClusterEnergy_EM_100"] < 1.1
branches = ["trk_ClusterEnergy_EM_100"]
sel_Lar1_1GeV = Calculation(Lar1_1GeV, branches)

def Z0SinThetaLess1_5(trk):
    return trk["trk_z0sintheta"] < 1.5
branches = ["trk_z0sintheta"]
sel_Z0SinThetaLess1_5 = Calculation(Z0SinThetaLess1_5, branches)

def d0Less1_5(trk):
    return trk["trk_d0"] < 1.5
branches = ["trk_d0"]
sel_d0Less1_5 = Calculation(d0Less1_5, branches)

def ExtrapolAcceptanceCalculator(trk, min_cut, max_cut):
    trk_eta = trkEtaPhiECAL(trk)[0]

    #check that both tracks are in the acceptance
    in_acceptance = (trk_eta < max_cut) & (trk_eta > min_cut)

    return in_acceptance

def NonZeroEnergy(trk):
    return (trk["trk_nclusters_EM_200"] + trk["trk_nclusters_HAD_200"]) > 0.5 #there was at least one cluster assocated with the track
branches = ["trk_nclusters_EM_200", "trk_nclusters_HAD_200"]
sel_NonZeroEnergy = Calculation(NonZeroEnergy, branches)

def ELessEqual0(trk):
    return (trk["trk_ClusterEnergy_EM_200"] + trk["trk_ClusterEnergy_HAD_200"]) <= 1e-6
branches = ["trk_ClusterEnergy_EM_200", "trk_ClusterEnergy_HAD_200"]
sel_ELessEqual0 = Calculation(ELessEqual0, branches)

def EHadFracAbove70(trk):
    return ((trk["trk_ClusterEnergy_HAD_200"]) / (trk["trk_ClusterEnergy_EM_200"] + trk["trk_ClusterEnergy_HAD_200"]) >= 0.7) & ((trk["trk_ClusterEnergy_EM_200"] + trk["trk_ClusterEnergy_HAD_200"]) > 0.0)
branches = ["trk_ClusterEnergy_EM_200", "trk_ClusterEnergy_HAD_200"]
sel_EHadFracAbove70 = Calculation(EHadFracAbove70, branches)

#A general function to pick different regions of the atlas detector based on track eta in the ID
def IDAcceptanceCalculator(trk, min_cut, max_cut):
    trk_etaID = np.abs(trk["trk_etaID"])
    upper_in_acceptance = trk_etaID < max_cut
    lower_in_acceptance = trk_etaID > min_cut
    return (upper_in_acceptance & lower_in_acceptance)

#You can do fancy things with lambda functions here
def EtaBin(trk, min_cut, max_cut):
    return ExtrapolAcceptanceCalculator(trk, min_cut, max_cut)

def PBin(trk, min_cut, max_cut):
    return (trk["trk_p"] > min_cut) & (trk["trk_p"] <= max_cut)

def NTRTX(trk, nTRT_low, nTRT_high):
    return ((trk["trk_nTRT"] >= nTRT_low) & (trk["trk_nTRT"] < nTRT_high)) | (np.abs(trk["trk_etaID"]) > 2.0)

def NPVBin(trk, npv_low, npv_high):
    return (trk["trk_NPV2"] >= npv_low) & (trk["trk_NPV2"] < npv_high)
