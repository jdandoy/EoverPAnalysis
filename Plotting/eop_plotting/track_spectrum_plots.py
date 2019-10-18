from histogram_filling import HistogramFiller, get_p, get_bins, get_log_bins, create_selection_function
import ROOT
import os
from math import pi
import pickle
import numpy as np
from selections import EtaBin, PBin, sel_SubleadingTrack
from variables import calc_trkP, calc_trkPt, calc_TruthMomentum, calc_trkEta

def create_spectrum_plots(hist_filler, base_selection, eta_ranges,  p_bins_for_eta_range, description):
    eta_count = -1
    for eta_range, p_bins in zip(eta_ranges, p_bins_for_eta_range):
        eta_count += 1
        #get the function that selects tracks in that bin
        eta_bin_selection = create_selection_function(EtaBin, ["trk_etaEMB2","trk_etaEME2","trk_phiEME2", "trk_phiEMB2"], eta_range[0], eta_range[1])

        selections = base_selection + [eta_bin_selection]

        NPtBins = len(p_bins) - 1
        Pt_low = 0.5
        Pt_high = max(p_bins)
        ptbins = get_log_bins(Pt_low, Pt_high, NPtBins)

        ################################################################################
        histogram_name = "TrackPtSpectrum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        hist_filler.book_histogram_fill(histogram_name,\
                                           calc_trkPt,\
                                           selections = selections,\
                                           bins = ptbins,\
                                           xlabel ="Track P_{T} [GeV]",\
                                           ylabel = "Number of Tracks")

        ################################################################################
        histogram_name = "TrackPSpectrum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        hist_filler.book_histogram_fill(histogram_name,\
                                           calc_trkP,\
                                           selections = selections,\
                                           bins = p_bins,\
                                           xlabel ="Track P [GeV]",\
                                           ylabel = "Number of Tracks")

        ################################################################################
        histogram_name = "2dTrackPSpectrumVsEta"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        hist_filler.book_2dhistogram_fill(histogram_name,\
                                           calc_trkP,\
                                           calc_trkEta,\
                                           selections = selections,\
                                           bins_x = p_bins,\
                                           bins_y = list(np.linspace(-2.5, 2.5, 100)),\
                                           xlabel ="Track Pt [GeV]",\
                                           ylabel = "Track #eta")

        ################################################################################
        histogram_name = "TrackTruthPSpectrum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        hist_filler.book_histogram_fill(histogram_name,\
                                           calc_TruthMomentum,\
                                           selections = selections,\
                                           bins = p_bins,\
                                           xlabel ="Track Truth P [GeV]",\
                                           ylabel = "Number of Tracks")

        ################################################################################
        histogram_name = "2dTrackTruthPSpectrumVsEta"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        hist_filler.book_2dhistogram_fill(histogram_name,\
                                           calc_TruthMomentum,\
                                           calc_trkEta,\
                                           selections = selections,\
                                           bins_x = p_bins,\
                                           bins_y = list(np.linspace(-2.5, 2.5, 100)),\
                                           xlabel ="Track Truth P [GeV]",\
                                           ylabel = "Track #eta")
