
from HistogramFillingTools.HistogramFiller import HistogramFiller, get_p, get_bins, get_log_bins, create_selection_function
import ROOT

import os
from math import pi
import pickle
import numpy as np
from selections.selections import EtaBin, PBin, sel_SubleadingTrack
from variables.variables import calc_trkP, calc_trkPt, calc_TruthMomentum

def CreateTrackSpectrumPlots(hist_filler, base_selection, eta_ranges,  p_bins_for_eta_range, description):
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
        hist_filler.BookHistograms(histogram_name,\
                                           calc_trkPt,\
                                           list_selections = selections,\
                                           bins = ptbins,\
                                           xlabel ="Track P_{T} [GeV]",\
                                           ylabel = "Number of Tracks")

        ################################################################################
        histogram_name = "TrackPSpectrum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        hist_filler.BookHistograms(histogram_name,\
                                           calc_trkPt,\
                                           list_selections = selections,\
                                           bins = p_bins,\
                                           xlabel ="Track P [GeV]",\
                                           ylabel = "Number of Tracks")

        ################################################################################
        histogram_name = "TrackTruthPSpectrum"
        histogram_name = histogram_name + "_" + "_" + description + "_Eta_" + str(eta_count)
        hist_filler.BookHistograms(histogram_name,\
                                           calc_TruthMomentum,\
                                           list_selections = selections,\
                                           bins = p_bins,\
                                           xlabel ="Track Truth P [GeV]",\
                                           ylabel = "Number of Tracks")

