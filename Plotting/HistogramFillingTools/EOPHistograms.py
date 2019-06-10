
from HistogramFillingTools.HistogramFiller import HistogramFiller, get_p, get_bins, get_log_bins, create_selection_function
import ROOT

import os
from math import pi
import pickle
import numpy as np
from selections.selections import EtaBin, PBin, sel_SubleadingTrack
from variables.variables import calc_trkP


def PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description):
    #create std::vectors for the eta ranges of std::vectors
    VectorForEtaLow = ROOT.std.vector("float")()
    VectorForEtaHigh = ROOT.std.vector("float")()

    VectorForPBinsLowList = []
    VectorForPBinsHighList = []

    for i in range(0, len(eta_ranges)):
        VectorForPBinsLowList.append(0)
        VectorForPBinsHighList.append(0)

    if not outFile.GetListOfKeys().Contains(description + "BinningTree"):
        outFile.cd()
        binningTree = ROOT.TTree(description + "BinningTree", description + "BinningTree")

        binningTree.Branch(description + "EtaBinsLow", VectorForEtaLow)
        binningTree.Branch(description + "EtaBinsHigh", VectorForEtaHigh)

        eta_count = -1

        for eta_range in eta_ranges:
            eta_count += 1

            vlow = ROOT.std.vector("float")()
            vhigh = ROOT.std.vector("float")()

            VectorForPBinsLowList[eta_count] = vlow
            VectorForPBinsHighList[eta_count] = vhigh

            binningTree.Branch(description + "PBinsLow_Eta" + str(eta_count), vlow)
            binningTree.Branch(description + "PBinsHigh_Eta" + str(eta_count), vhigh)

        for eta_range, p_bins, VectorForPBinsLow, VectorForPBinsHigh in zip(eta_ranges, p_bins_for_eta_range, VectorForPBinsLowList, VectorForPBinsHighList):
            p_ranges = [ (p_bins[i], p_bins[i+1])  for i in range(0, len(p_bins)-1) ]
            for p_range in p_ranges:
                VectorForPBinsLow.push_back(p_range[0])
                VectorForPBinsHigh.push_back(p_range[1])

            VectorForEtaLow.push_back(eta_range[0])
            VectorForEtaHigh.push_back(eta_range[1])

        binningTree.Fill()
        binningTree.Write()
    

def CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges,p_bins_for_eta_range, description, doClusterPlots=False, doCalibHitPlots=False, optimalBinningFile=None):

    #define a set of eta bins
    eta_count = -1
    for eta_range, p_bins in zip(eta_ranges, p_bins_for_eta_range):
        eta_count += 1
        #get the function that selects tracks in that bin
        eta_bin_selection = create_selection_function(EtaBin, ["trk_etaEMB2","trk_etaEME2","trk_phiEME2", "trk_phiEMB2"], eta_range[0], eta_range[1])

        selections = base_selection + [eta_bin_selection]

        NPtBins = len(p_bins)
        Pt_low = 0.5
        Pt_high = max(p_bins)
        ptbins = get_log_bins(Pt_low, Pt_high, NPtBins)
        eop_bins = get_bins(-1.0, +3.0, 800) # use a super fine binning

        histogramName = "EOPProfileVsMomentum"
        histogramName = histogramName + "_" + "_" + description + "_Eta_" + str(eta_count)

        from variables.variables import calc_EOP
        AverageEOP  =  hist_filler.BookTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )

        histogramName = "2DHist_EOPVsMomentum"
        histogramName = histogramName + "_" + "_" + description + "_Eta_" + str(eta_count)
        AverageEOP  =  hist_filler.Book2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p",\
                                                  )

        histogramName = "EnergyAnulusProfileVsMomentum"
        histogramName = histogramName + "_" + "_" + description + "_Eta_" + str(eta_count)
        from variables.variables import calc_EnergyAnulus
        AverageAnulus =  hist_filler.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )

        histogramName = "2DHist_EnergyAnulusVsMomentum"
        histogramName = histogramName + "_" + "_" + description + "_Eta_" + str(eta_count)
        AverageAnulus =  hist_filler.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )

        histogramName = "EnergyBkgProfileVsMomentum"
        histogramName = histogramName + "_" + "_" + description + "_Eta_" + str(eta_count)
        from variables.variables import calc_EOPBkg
        AverageAnulus =  hist_filler.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )


        histogramName = "2DHist_EnergyBkgVsMomentum"
        histogramName = histogramName + "_" + "_" + description + "_Eta_" + str(eta_count)
        AverageAnulus =  hist_filler.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p BKG",\
                                                  )

        p_count = -1
        p_ranges = [(p_bins[i],p_bins[i+1]) for i in range(0, len(p_bins)-1)]
        for p_range in p_ranges:
            p_count += 1
            print("The prange is " + str(p_range))
            p_bin_selection = create_selection_function(PBin, ["trk_p"], p_range[0], p_range[1])
            selections = base_selection + [eta_bin_selection] + [p_bin_selection]

            histogramName = "EOPDistribution" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            EOPDist  =  hist_filler.BookHistograms(histogramName,
                                                      calc_EOP,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )

            if doCalibHitPlots:
              from variables.variables import calc_CalibHitFrac, calc_PhotonCalibHitFrac, calc_HadronCalibHitFrac, sel_HasCalibHit
              from variables.variables import calc_EMCalibHitFrac, calc_PhotonEMCalibHitFrac, calc_HadronEMCalibHitFrac, sel_HasEMCalibHit
              from variables.variables import calc_HADCalibHitFrac, calc_PhotonHADCalibHitFrac, calc_HadronHADCalibHitFrac, sel_HasHADCalibHit

              histogramName = "CalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.Book2DHistograms(histogramName,
                                                    calc_EOP,\
                                                    calc_CalibHitFrac,\
                                                    list_selections = selections + [sel_HasCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    range_low_y = 0.0,\
                                                    range_high_y = 1.0,\
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{Calib}_{Track}/E^{Calib}_{Total}",\
                                                    )

              histogramName = "PhotonCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.Book2DHistograms(histogramName,
                                                    calc_EOP,\
                                                    calc_PhotonCalibHitFrac,\
                                                    list_selections = selections + [sel_HasCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{Calib}_{Photons}/E^{Calib}_{Total}",\
                                                    )

              histogramName = "HadronicCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.Book2DHistograms(histogramName,
                                                    calc_EOP,\
                                                    calc_HadronCalibHitFrac,\
                                                    list_selections = selections + [sel_HasCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{Calib}_{Neutral Hadrons}/E^{Calib}_{Total}",\
                                                    )

              histogramName = "EMCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.Book2DHistograms(histogramName,
                                                    calc_EOP,\
                                                    calc_EMCalibHitFrac,\
                                                    list_selections = selections + [sel_HasEMCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{EM Calib}_{Track}/E^{EM Calib}_{Total}",\
                                                    )

              histogramName = "PhotonEMCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.Book2DHistograms(histogramName,
                                                    calc_EOP,\
                                                    calc_PhotonEMCalibHitFrac,\
                                                    list_selections = selections + [sel_HasEMCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{EM Calib}_{Photons}/E^{EM Calib}_{Total}",\
                                                    )

              histogramName = "HadronicEMCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.Book2DHistograms(histogramName,
                                                    calc_EOP,\
                                                    calc_HadronEMCalibHitFrac,\
                                                    list_selections = selections + [sel_HasEMCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{HAD Calib}_{Neutral Hadrons}/E^{HAD Calib}_{Total}",\
                                                    )

              histogramName = "HADCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.Book2DHistograms(histogramName,
                                                    calc_EOP,\
                                                    calc_HADCalibHitFrac,\
                                                    list_selections = selections + [sel_HasHADCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{HAD Calib}_{Track}/E^{HAD Calib}_{Total}",\
                                                    )

              histogramName = "PhotonHADCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.Book2DHistograms(histogramName,
                                                    calc_EOP,\
                                                    calc_PhotonHADCalibHitFrac,\
                                                    list_selections = selections + [sel_HasHADCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{HAD Calib}_{Photons}/E^{HAD Calib}_{Total}",\
                                                    )

              histogramName = "HadronicHADCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.Book2DHistograms(histogramName,
                                                    calc_EOP,\
                                                    calc_HadronHADCalibHitFrac,\
                                                    list_selections = selections + [sel_HasHADCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{HAD Calib}_{Neutral Hadrons}/E^{HAD Calib}_{Neutral Hadrons}",\
                                                    )

            histogramName = "EOPBkgDistribution" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            EOPBkgDist  =  hist_filler.BookHistograms(histogramName,
                                                      calc_EOPBkg,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p Bkg",\
                                                      )
            histogram_name = "trkTRTHits" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            from variables.variables import calc_nTRT
            hist_filler.BookHistograms(histogram_name,\
                                   calc_nTRT,\
                                   list_selections = selections,\
                                   range_low = -0.5,\
                                   range_high = 59.5,\
                                   bins = 60,\
                                   xlabel = "Number of TRT Hits",\
                                   ylabel = "Number of Tracks")

            histogram_name = "trkEMDR100" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            from variables.variables import calc_EnergyEMDR100
            hist_filler.BookHistograms(histogram_name,\
                                   calc_EnergyEMDR100,\
                                   list_selections = selections,\
                                   range_low = -2.0,\
                                   range_high = + 10.0,\
                                   bins = 48,\
                                   xlabel = "E_{EM}^{#DeltaR<0.1}[GeV]",\
                                   ylabel = "Number of Tracks")

            histogram_name = "MomentumHadFrac" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            from variables.variables import calc_MomentumHadFrac
            hist_filler.BookHistograms(histogram_name,\
                                   calc_MomentumHadFrac,\
                                   list_selections = selections,\
                                   range_low = -1.0,\
                                   range_high = + 5.0,\
                                   bins = 48,\
                                   xlabel = "E^{HAD}/P",\
                                   ylabel = "Number of Tracks")

            histogram_name = "HadFrac" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            from variables.variables import calc_HadFrac
            hist_filler.BookHistograms(histogram_name,\
                                   calc_HadFrac,\
                                   list_selections = selections,\
                                   range_low = -1.0,\
                                   range_high = + 2.0,\
                                   bins = 48,\
                                   xlabel = "E^{HAD}/E^{Total}",\
                                   ylabel = "Number of Tracks")

            if doClusterPlots:

               from variables.variables import calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD,  calc_trkNClusters_emlike, calc_trkNClusters_hadlike
               histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
               xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
               variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]

               for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
                   histogram_name = histogram_name + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
                   hist_filler.BookHistograms(histogram_name,\
                                          variable,\
                                          list_selections = selections,\
                                          bins = 10,\
                                          range_low = -0.5,\
                                          range_high = 9.5,\
                                          xlabel=xlabel,\
                                          ylabel="Number of Tracks")
