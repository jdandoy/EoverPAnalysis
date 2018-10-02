from calculation.calculation import calculation, calculationDataMC
import numpy as np

def NoSelection(trk):
    return np.ones(len(trk)) > 0.0
branches = []
sel_NoSelection = calculation(NoSelection, branches)

def Lar1GeV(trk):
    return trk["trk_sumEPos_Lar_200"] > 1.0
branches = ["trk_sumEPos_Lar_200"]
sel_Lar1GeV = calculation(Lar1GeV, branches)

def acceptanceCalculator(trk, min_cut, max_cut):
    trk_etaEMB = np.abs(trk["trk_etaEMB2"])
    trk_etaEME = np.abs(trk["trk_etaEME2"])

    upper_in_acceptance_EMB = trk_etaEMB < max_cut
    upper_in_acceptance_EME = trk_etaEME < max_cut

    lower_in_acceptance_EMB = trk_etaEMB > min_cut
    lower_in_acceptance_EME = trk_etaEME > min_cut

    return (upper_in_acceptance_EMB | upper_in_acceptance_EME) & (lower_in_acceptance_EMB | lower_in_acceptance_EME)

def Eta0_6(trk):
    return acceptanceCalculator(trk, 0.0, 0.6)
branches = ["trk_etaEMB2", "trk_etaEME2"]
sel_Eta0_6 = calculation(Eta0_6, branches)

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
