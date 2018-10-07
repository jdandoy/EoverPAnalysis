from calculation.calculation import calculation, calculationDataMC
import numpy as np

def trkPt(trk):
    return trk["trk_pt"]
branches = ["trk_pt"]
calc_trkPt = calculation(trkPt, branches)

def trkEtaID(trk):
    return trk["trk_etaID"]
branches = ["trk_etaID"]
calc_trkEta = calculation(trkEtaID, branches)

def trkEtaID_ABS(trk):
    return np.abs(trk["trk_etaID"])
branches = ["trk_etaID"]
calc_trkEta_ABS = calculation(trkEtaID_ABS, branches)

def trkEtaECAL(trk):
    trk_etaEMB = np.abs(trk["trk_etaEMB2"])
    trk_etaEME = np.abs(trk["trk_etaEME2"])

    trkEtaECAL = np.ones(len(trk)) * 1000.0

    hasEMB = trk_etaEMB < 100.0
    hasEME = np.logical_not(hasEMB) & (trk_etaEME < 100.0)

    print "This many tracks had extensions to both the barrel and the emec " + str(np.sum(1*(hasEMB & hasEME)))

    trkEtaECAL[hasEMB] = trk["trk_etaEMB2"][hasEMB]
    trkEtaECAL[hasEME] = trk["trk_etaEME2"][hasEME]

    print(np.logical_not(hasEMB | hasEME))

    if np.any(np.logical_not(hasEMB | hasEME)):
        print(np.sum( 1 * np.logical_not(hasEMB | hasEME)))
        raise ValueError("one of the trakc did not have an extrapolation to the electromagnetic calorimeter")
    return trkEtaECAL

branches = ["trk_etaEMB2","trk_etaEME2"]
calc_trkEtaECAL = calculation(trkEtaECAL, branches)

def trkNearestNeighbourEM2(trk):
    return trk["trk_nearest_dR"]
branches = ["trk_nearest_dR"]
calc_trkNearestNeighbourEM2 =  calculation(trkNearestNeighbourEM2, branches)

def trkNPV2(trk):
    return trk["trk_NPV_2"]
branches = ["trk_NPV_2"]
calc_trkNPV2 = calculation(trkNPV2, branches)

def trkNPV4(trk):
    return trk["trk_NPV_4"]
branches = ["trk_NPV_4"]
calc_trkNPV4 = calculation(trkNPV4, branches)

def trkAverageMu(trk):
    return trk["trk_averagemu"]
branches = ["trk_averagemu"]
calc_trkAverageMu = calculation(trkAverageMu, branches)

def trkP(trk):
    return trk["trk_p"]
branches = ["trk_p"]
calc_trkP = calculation(trkP, branches)

def trkEtaID(trk):
    return trk["trk_etaID"]
branches = ["trk_etaID"]
calc_trkEtaID = calculation(trkEtaID, branches)

def trkEtaEME2(trk):
    return trk["trk_etaEME2"]
branches = ["trk_etaEME2"]
calc_trkEtaEME2 = calculation(trkEtaEME2, branches)

def trkEtaEMB2(trk):
    return trk["trk_etaEMB2"]
branches = ["trk_etaEMB2"]
calc_trkEtaEMB2 = calculation(trkEtaEMB2, branches)

def EOP(trk):
    return trk["trk_sumE_Total_200"]/trk["trk_p"]
branches = ["trk_sumE_Total_200", "trk_p"]
calc_EOP = calculation(EOP, branches)

def weight(trk, isData):
    if not isData:
        return trk["trkWeight"]
    return np.ones(len(trk))
branches = ["trkWeight"]
calc_weight = calculationDataMC(weight, branches)
