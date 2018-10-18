from calculation.calculation import calculation, calculationDataMC
import numpy as np

def NoSelection(trk):
    return np.ones(len(trk)) > 0.0
branches = []
sel_NoSelection = calculation(NoSelection, branches)

def Lar1GeV(trk):
    return trk["trk_sumEPos_Lar_100"] < 1.1
branches = ["trk_sumEPos_Lar_100"]
sel_Lar1GeV = calculation(Lar1GeV, branches)

def Z0SinThetaLess1_5(trk):
    return trk["trk_z0sintheta"] < 1.5
branches = ["trk_z0sintheta"]
sel_Z0SinThetaLess1_5 = calculation(Z0SinThetaLess1_5, branches)

def d0Less1_5(trk):
    return trk["trk_d0"] < 1.5
branches = ["trk_d0"]
sel_d0Less1_5 = calculation(d0Less1_5, branches)

def EM2AcceptanceCalculator(trk, min_cut, max_cut):
    trk_etaEMB = np.abs(trk["trk_etaEMB2"])
    trk_etaEME = np.abs(trk["trk_etaEME2"])

    upper_in_acceptance_EMB = trk_etaEMB < max_cut
    upper_in_acceptance_EME = trk_etaEME < max_cut

    lower_in_acceptance_EMB = trk_etaEMB > min_cut
    lower_in_acceptance_EME = trk_etaEME > min_cut

    return (upper_in_acceptance_EMB | upper_in_acceptance_EME) & (lower_in_acceptance_EMB | lower_in_acceptance_EME)

def NonZeroEnergy(trk):
    return np.logical_not(np.abs(trk["trk_E_Total_200"]) < 1e-10)
branches = ["trk_E_Total_200"]
sel_NonZeroEnergy = calculation(NonZeroEnergy, branches)

def ELessEqual0(trk):
    return trk["trk_E_Total_200"] <= 1e-10
branches = ["trk_E_Total_200"]
sel_ELessEqual0 = calculation(ELessEqual0, branches)

def IDAcceptanceCalculator(trk, min_cut, max_cut):
    trk_etaID = np.abs(trk["trk_etaID"])
    upper_in_acceptance = trk_etaID < max_cut
    lower_in_acceptance = trk_etaID > min_cut
    return (upper_in_acceptance & lower_in_acceptance)

def ECALEta0_6(trk):
    return EM2AcceptanceCalculator(trk, 0.0, 0.6)
branches = ["trk_etaEMB2", "trk_etaEME2"]
sel_ECALEta0_6 = calculation(ECALEta0_6, branches)

#### Here are the eta track selections ###
def IDEta0_6(trk):
    return IDAcceptanceCalculator(trk, 0.0, 0.6)
branches = ["trk_etaID"]
sel_IDEta0_6 = calculation(IDEta0_6, branches)

def IDEta06_11(trk):
    return IDAcceptanceCalculator(trk, 0.6, 1.1)
branches = ["trk_etaID"]
sel_IDEta06_11 = calculation(IDEta06_11, branches)

def IDEta11_14(trk):
    return IDAcceptanceCalculator(trk, 1.1, 1.4)
branches = ["trk_etaID"]
sel_IDEta11_14 = calculation(IDEta11_14, branches)

def IDEta14_15(trk):
    return IDAcceptanceCalculator(trk, 1.4, 1.5)
branches = ["trk_etaID"]
sel_IDEta14_15 = calculation(IDEta14_15, branches)

def IDEta15_18(trk):
    return IDAcceptanceCalculator(trk, 1.5, 1.8)
branches = ["trk_etaID"]
sel_IDEta15_18 = calculation(IDEta15_18, branches)

def IDEta18_23(trk):
    return IDAcceptanceCalculator(trk, 1.8, 2.3)
branches = ["trk_etaID"]
sel_IDEta18_23 = calculation(IDEta18_23, branches)

def IDEta19_23(trk):
    return IDAcceptanceCalculator(trk, 1.9, 2.3)
branches = ["trk_etaID"]
sel_IDEta19_23 = calculation(IDEta19_23, branches)

def PBetween12_18(trk):
    return (1.2 < trk["trk_p"]) & (trk["trk_p"] < 1.8)
branches = ["trk_p"]
sel_PBetween12_18 = calculation(PBetween12_18, branches)

def PBetween22_28(trk):
    return (2.2 < trk["trk_p"]) & (trk["trk_p"] < 2.8)
branches = ["trk_p"]
sel_PBetween22_28 = calculation(PBetween22_28, branches)

def PBetween28_36(trk):
    return (2.8 < trk["trk_p"]) & (trk["trk_p"] < 3.6)
branches = ["trk_p"]
sel_PBetween28_36 = calculation(PBetween28_36, branches)

def NTRT15(trk):
    return trk["trk_nTRT"] >= 15
branches = ["trk_nTRT"]
sel_NTRT15 = calculation(NTRT15, branches)

def NTRT20(trk):
    return trk["trk_nTRT"] >= 20
branches = ["trk_nTRT"]
sel_NTRT20 = calculation(NTRT20, branches)

def NTRT25(trk):
    return trk["trk_nTRT"] >= 25
branches = ["trk_nTRT"]
sel_NTRT25 = calculation(NTRT25, branches)

def NTRT30(trk):
    return trk["trk_nTRT"] >= 30
branches = ["trk_nTRT"]
sel_NTRT30 = calculation(NTRT30, branches)

def NTRT35(trk):
    return trk["trk_nTRT"] >= 35
branches = ["trk_nTRT"]
sel_NTRT35 = calculation(NTRT35, branches)

def PGreater1(trk):
    return trk["trk_p"] > 1.0
branches =["trk_p"]
sel_PGreater1 = calculation(PGreater1, branches)

def PGreater1_5(trk):
    return trk["trk_p"] > 1.5
branches =["trk_p"]
sel_PGreater1_5 = calculation(PGreater1_5, branches)

def PGreater2(trk):
    return trk["trk_p"] > 2.0
branches =["trk_p"]
sel_PGreater2 = calculation(PGreater2, branches)

def PGreater2_5(trk):
    return trk["trk_p"] > 2.5
branches =["trk_p"]
sel_PGreater2_5 = calculation(PGreater2_5, branches)
