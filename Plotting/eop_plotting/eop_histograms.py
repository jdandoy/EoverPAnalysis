from eop_plotting.histogram_filling import HistogramFiller, get_p, get_bins, get_log_bins, create_selection_function
import ROOT
import os
from math import pi
import pickle
import numpy as np
from eop_plotting.selections import EtaBin, PBin, sel_SubleadingTrack
from eop_plotting.variables import calc_trkP


def put_binning_vectors_in_file(outFile, eta_ranges, p_bins_for_eta_range, description):
    #create std::vectors for the eta ranges of std::vectors
    eta_low_vector = ROOT.std.vector("float")()
    eta_high_vector = ROOT.std.vector("float")()

    p_bins_low_list = []
    p_bins_high_list = []

    for i in range(0, len(eta_ranges)):
        p_bins_low_list.append(0)
        p_bins_high_list.append(0)

    if not outFile.GetListOfKeys().Contains(description + "BinningTree"):
        outFile.cd()
        binningTree = ROOT.TTree(description + "BinningTree", description + "BinningTree")

        binningTree.Branch(description + "EtaBinsLow", eta_low_vector)
        binningTree.Branch(description + "EtaBinsHigh", eta_high_vector)

        eta_count = -1

        for eta_range in eta_ranges:
            eta_count += 1

            vlow = ROOT.std.vector("float")()
            vhigh = ROOT.std.vector("float")()

            p_bins_low_list[eta_count] = vlow
            p_bins_high_list[eta_count] = vhigh

            binningTree.Branch(description + "PBinsLow_Eta" + str(eta_count), vlow)
            binningTree.Branch(description + "PBinsHigh_Eta" + str(eta_count), vhigh)

        for eta_range, p_bins, VectorForPBinsLow, VectorForPBinsHigh in zip(eta_ranges, p_bins_for_eta_range, p_bins_low_list, p_bins_high_list):
            p_ranges = [ (p_bins[i], p_bins[i+1])  for i in range(0, len(p_bins)-1) ]
            for p_range in p_ranges:
                VectorForPBinsLow.push_back(p_range[0])
                VectorForPBinsHigh.push_back(p_range[1])

            eta_low_vector.push_back(eta_range[0])
            eta_high_vector.push_back(eta_range[1])

        binningTree.Fill()
        binningTree.Write()

def create_eop_histograms(hist_filler, base_selection, eta_ranges,p_bins_for_eta_range, description, do_cluster_plots=False, do_calib_hit_plots=False):
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

        histogram_name = "EOPProfileVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)

        from eop_plotting.variables import calc_EOP
        hist_filler.book_tprofile_fill(histogram_name,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )

        histogram_name = "2DHist_EOPVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        hist_filler.book_2dhistogram_fill(histogram_name,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p",\
                                                  )
        from variables import cone_strings, total_energy_annulus_template
        from math import pi
        def annulus_area(low,high):
            return pi * ((high ** 2) - (low ** 2))
        for (low, high) in zip(cone_strings[:-1], cone_strings[1:]):
            #this function calculates the energy in a cone from low to high
            func_e = lambda x, y = low, z = high : total_energy_annulus_template(x,y,z)
            func_e.__name__ = "energy_in_cone_{}_{}".format(low, high)
            #this function calculates the energy in a cone from low to high
            func_e_area = lambda x, y = low, z = high, l=low, h=high : total_energy_annulus_template(x,y,z)/(annulus_area(float(l)/1000.0, float(h)/1000.0))
            func_e_area.__name__ = "energy_in_cone_area_{}_{}".format(low, high)
            #this function calculates the eop in a cone from low to high
            func_eop = lambda x, y = low, z = high : total_energy_annulus_template(x,y,z)/x["trk_p"]
            func_eop.__name__ = "eop_in_cone_{}_{}".format(low, high)
            #this function calculates the eop in a cone from low to high
            func_eop_area = lambda x, y = low, z = high, l=low, h=high : (total_energy_annulus_template(x,y,z)/x["trk_p"])/(annulus_area(float(l)/1000.0, float(h)/1000.0))
            func_eop_area.__name__ = "eop_in_cone_area_{}_{}".format(low, high)

            #define the calculation and the branches that we need for this
            branches = ["trk_ClusterEnergy_EM_{}".format(low), "trk_ClusterEnergy_HAD_{}".format(low)]
            branches += ["trk_ClusterEnergy_EM_{}".format(high), "trk_ClusterEnergy_HAD_{}".format(high)]
            clean_branches = []
            for b in branches:
                if "000" not in b:
                    clean_branches.append(b)
            branches = clean_branches
            from calculation import Calculation
            calc_cone_e = Calculation(func_e, branches)
            calc_cone_e_area = Calculation(func_e_area, branches)
            branches = branches + ["trk_p"]
            calc_cone_eop = Calculation(func_eop, branches)
            calc_cone_eop_area = Calculation(func_eop_area, branches)

            #book the histograms
            low_descr = float(low)/10.0
            high_descr = float(high)/10.0
            histogram_name = "EnergyAnnulusProfileVsMomentum_{}_{}".format(low, high)
            histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
            hist_filler.book_tprofile_fill(histogram_name,\
                                                      calc_trkP,\
                                                      calc_cone_e,\
                                                      selections = selections,\
                                                      bins = p_bins,\
                                                      xlabel ="P[GeV]",\
                                                      ylabel = "<E_{r#in[" + "{},{}".format(low_descr, high_descr) + "]}>[GeV]",\
                                                      )

            histogram_name = "EOPProfileVsMomentum_{}_{}".format(low, high)
            histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
            hist_filler.book_tprofile_fill(histogram_name,\
                                                      calc_trkP,\
                                                      calc_cone_eop,\
                                                      selections = selections,\
                                                      bins = p_bins,\
                                                      xlabel ="P[GeV]",\
                                                      ylabel = "<E_{r#in[" + "{},{}".format(low_descr, high_descr) + "]}/p>",\
                                                      )

            #book the histograms
            low_descr = float(low)/10.0
            high_descr = float(high)/10.0
            histogram_name = "EnergyAnnulusProfileVsMomentum_{}_{}_Area".format(low, high)
            histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
            hist_filler.book_tprofile_fill(histogram_name,\
                                                      calc_trkP,\
                                                      calc_cone_e_area,\
                                                      selections = selections,\
                                                      bins = p_bins,\
                                                      xlabel ="P[GeV]",\
                                                      ylabel = "<E_{r#in[" + "{},{}".format(low_descr, high_descr) + "]}>/Area [GeV]",\
                                                      )

            histogram_name = "EOPProfileVsMomentum_{}_{}_Area".format(low, high)
            histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
            hist_filler.book_tprofile_fill(histogram_name,\
                                                      calc_trkP,\
                                                      calc_cone_eop_area,\
                                                      selections = selections,\
                                                      bins = p_bins,\
                                                      xlabel ="P[GeV]",\
                                                      ylabel = "<E_{r#in[" + "{},{}".format(low_descr, high_descr) + "]}/p>/ Area",\
                                                      )

        histogram_name = "EnergyAnulusProfileVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        from eop_plotting.variables import calc_EnergyAnulus
        hist_filler.book_tprofile_fill(histogram_name,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )

#       from variables import calc_EOTotalEMCalibHitEnergy
#       histogram_name = "EOTotalEMCalibHitEnergyProfileVsMomentum"  + "_" + "_" + description + "_Eta_" + str(eta_count)
#       hist_filler.book_tprofile_fill(histogram_name,\
#                                      calc_trkP,\
#                                      calc_EOTotalEMCalibHitEnergy,\
#                                      selections = selections,\
#                                      bins = p_bins,\
#                                      xlabel = "P[GeV]",\
#                                      ylabel = "<E_{Total Calibration}/P>",\
#                                      )


        histogram_name = "2DHist_EnergyAnulusVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        hist_filler.book_2dhistogram_fill(histogram_name,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )

        histogram_name = "EnergyBkgProfileVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        from eop_plotting.variables import calc_EOPBkg
        hist_filler.book_tprofile_fill(histogram_name,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogram_name = "EnergyBkgUpProfileVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        from eop_plotting.variables import calc_EOPBkgUp
        hist_filler.book_tprofile_fill(histogram_name,\
                                                  calc_trkP,\
                                                  calc_EOPBkgUp,\
                                                  selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogram_name = "EnergyBkgDownProfileVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        from eop_plotting.variables import calc_EOPBkgDown
        hist_filler.book_tprofile_fill(histogram_name,\
                                                  calc_trkP,\
                                                  calc_EOPBkgDown,\
                                                  selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogram_name = "EnergyBigBkgProfileVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        from eop_plotting.variables import calc_EOPBigBkg
        hist_filler.book_tprofile_fill(histogram_name,\
                                                  calc_trkP,\
                                                  calc_EOPBigBkg,\
                                                  selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogram_name = "EnergyBigBkgUpProfileVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        from eop_plotting.variables import calc_EOPBigBkgUp
        hist_filler.book_tprofile_fill(histogram_name,\
                                                  calc_trkP,\
                                                  calc_EOPBigBkgUp,\
                                                  selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogram_name = "EnergyBigBkgDownProfileVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        from eop_plotting.variables import calc_EOPBigBkgDown
        hist_filler.book_tprofile_fill(histogram_name,\
                                                  calc_trkP,\
                                                  calc_EOPBigBkgDown,\
                                                  selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )


        histogram_name = "2DHist_EnergyBkgVsMomentum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        hist_filler.book_2dhistogram_fill(histogram_name,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  selections = selections,\
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

            histogram_name = "EOPDistribution" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            EOPDist  =  hist_filler.book_histogram_fill(histogram_name,
                                                      calc_EOP,\
                                                      selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )

            if do_calib_hit_plots:
              from eop_plotting.variables import calc_CalibHitFrac, calc_PhotonCalibHitFrac, calc_HadronCalibHitFrac, sel_HasCalibHit
              from eop_plotting.variables import calc_EMCalibHitFrac, calc_PhotonEMCalibHitFrac, calc_HadronEMCalibHitFrac, sel_HasEMCalibHit
              from eop_plotting.variables import calc_HADCalibHitFrac, calc_PhotonHADCalibHitFrac, calc_HadronHADCalibHitFrac, sel_HasHADCalibHit

              histogram_name = "CalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.book_2dhistogram_fill(histogram_name,
                                                    calc_EOP,\
                                                    calc_CalibHitFrac,\
                                                    selections = selections + [sel_HasCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    range_low_y = 0.0,\
                                                    range_high_y = 1.0,\
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{Calib}_{Track}/E^{Calib}_{Total}",\
                                                    )

              histogram_name = "PhotonCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.book_2dhistogram_fill(histogram_name,
                                                    calc_EOP,\
                                                    calc_PhotonCalibHitFrac,\
                                                    selections = selections + [sel_HasCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{Calib}_{Photons}/E^{Calib}_{Total}",\
                                                    )

              histogram_name = "HadronicCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.book_2dhistogram_fill(histogram_name,
                                                    calc_EOP,\
                                                    calc_HadronCalibHitFrac,\
                                                    selections = selections + [sel_HasCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{Calib}_{Neutral Hadrons}/E^{Calib}_{Total}",\
                                                    )

              histogram_name = "EMCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.book_2dhistogram_fill(histogram_name,
                                                    calc_EOP,\
                                                    calc_EMCalibHitFrac,\
                                                    selections = selections + [sel_HasEMCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{EM Calib}_{Track}/E^{EM Calib}_{Total}",\
                                                    )

              histogram_name = "PhotonEMCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.book_2dhistogram_fill(histogram_name,
                                                    calc_EOP,\
                                                    calc_PhotonEMCalibHitFrac,\
                                                    selections = selections + [sel_HasEMCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{EM Calib}_{Photons}/E^{EM Calib}_{Total}",\
                                                    )

              histogram_name = "HadronicEMCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.book_2dhistogram_fill(histogram_name,
                                                    calc_EOP,\
                                                    calc_HadronEMCalibHitFrac,\
                                                    selections = selections + [sel_HasEMCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{HAD Calib}_{Neutral Hadrons}/E^{HAD Calib}_{Total}",\
                                                    )

              histogram_name = "HADCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.book_2dhistogram_fill(histogram_name,
                                                    calc_EOP,\
                                                    calc_HADCalibHitFrac,\
                                                    selections = selections + [sel_HasHADCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{HAD Calib}_{Track}/E^{HAD Calib}_{Total}",\
                                                    )

              histogram_name = "PhotonHADCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.book_2dhistogram_fill(histogram_name,
                                                    calc_EOP,\
                                                    calc_PhotonHADCalibHitFrac,\
                                                    selections = selections + [sel_HasHADCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{HAD Calib}_{Photons}/E^{HAD Calib}_{Total}",\
                                                    )

              histogram_name = "HadronicHADCalibrationHitTwoDHist_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
              hist_filler.book_2dhistogram_fill(histogram_name,
                                                    calc_EOP,\
                                                    calc_HadronHADCalibHitFrac,\
                                                    selections = selections + [sel_HasHADCalibHit],\
                                                    bins_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],\
                                                    bins_y =list( np.linspace(0.0,1.0,20)),
                                                    xlabel ="E/p",\
                                                    ylabel = "E^{HAD Calib}_{Neutral Hadrons}/E^{HAD Calib}_{Neutral Hadrons}",\
                                                    )

            histogram_name = "EOPBkgDistribution" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            EOPBkgDist  =  hist_filler.book_histogram_fill(histogram_name,
                                                      calc_EOPBkg,\
                                                      selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p Bkg",\
                                                      )
            histogram_name = "trkTRTHits" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            from eop_plotting.variables import calc_nTRT
            hist_filler.book_histogram_fill(histogram_name,\
                                   calc_nTRT,\
                                   selections = selections,\
                                   range_low = -0.5,\
                                   range_high = 59.5,\
                                   bins = 60,\
                                   xlabel = "Number of TRT Hits",\
                                   ylabel = "Number of Tracks")

            histogram_name = "trkEMDR100" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            from eop_plotting.variables import calc_EnergyEMDR100
            hist_filler.book_histogram_fill(histogram_name,\
                                   calc_EnergyEMDR100,\
                                   selections = selections,\
                                   range_low = -2.0,\
                                   range_high = + 10.0,\
                                   bins = 48,\
                                   xlabel = "E_{EM}^{#DeltaR<0.1}[GeV]",\
                                   ylabel = "Number of Tracks")

            histogram_name = "MomentumHadFrac" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            from eop_plotting.variables import calc_MomentumHadFrac
            hist_filler.book_histogram_fill(histogram_name,\
                                   calc_MomentumHadFrac,\
                                   selections = selections,\
                                   range_low = -1.0,\
                                   range_high = + 5.0,\
                                   bins = 48,\
                                   xlabel = "E^{HAD}/P",\
                                   ylabel = "Number of Tracks")

            histogram_name = "HadFrac" + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            from eop_plotting.variables import calc_HadFrac
            hist_filler.book_histogram_fill(histogram_name,\
                                   calc_HadFrac,\
                                   selections = selections,\
                                   range_low = -1.0,\
                                   range_high = + 2.0,\
                                   bins = 48,\
                                   xlabel = "E^{HAD}/E^{Total}",\
                                   ylabel = "Number of Tracks")

            if do_cluster_plots:

               from eop_plotting.variables import calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD,  calc_trkNClusters_emlike, calc_trkNClusters_hadlike
               histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
               xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
               variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]

               for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
                   histogram_name = histogram_name + "_" + description + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
                   hist_filler.book_histogram_fill(histogram_name,\
                                          variable,\
                                          selections = selections,\
                                          bins = 10,\
                                          range_low = -0.5,\
                                          range_high = 9.5,\
                                          xlabel=xlabel,\
                                          ylabel="Number of Tracks")
