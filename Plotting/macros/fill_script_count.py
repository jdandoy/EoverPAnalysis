#!/usr/bin/env python
# coding: utf-8

from eop_plotting.histogram_filling import HistogramFiller, get_p, get_bins, get_log_bins, create_selection_function, create_inverse_selection_function
from eop_plotting.eop_histograms import put_binning_vectors_in_file, create_eop_histograms
from eop_plotting.track_spectrum_plots import create_spectrum_plots

import ROOT
import os
from math import pi
import pickle
import numpy as np
from selections import EtaBin, PBin, sel_SubleadingTrack, sel_Event, sel_hasHADExtrapolation, NTRTX

def write_histograms(histogram_dictionary, outFile):
    outFile.cd()
    for key in histogram_dictionary:
        if not outFile.cd(key):
            outFile.mkdir(key)
        outFile.cd(key)
        histogram_dictionary[key].Write()

#create the functions that select tracks in the different eta bins
eta_bin_edges = [0.0, 0.2, 0.7, 1.3, 1.8, 2.5]
eta_bin_tuples = [(eta_bin_edges[i], eta_bin_edges[i+1]) for i in range(0, len(eta_bin_edges)-1)]
eta_bin_descriptions = ["eta_extrapol00_02", "eta_extrapol02_07", "eta_extrapol07_13", "eta_extrapol13_18", "eta_extapol18_25"]
eta_bin_branches = ["trk_etaEMB2","trk_etaEME2","trk_phiEMB2", "trk_phiEME2"]
eta_bin_selections = [create_selection_function(EtaBin, eta_bin_branches, eta_bin_tuple[0], eta_bin_tuple[1]) for eta_bin_tuple in eta_bin_tuples]

from selections import sel_HardScatter, ParticlePDGID_ABS, ParticlePDGID, sel_TightIso
sel_Pion = create_selection_function(ParticlePDGID_ABS, ["trk_truthPdgId"], 211.0)
sel_PionPos = create_selection_function(ParticlePDGID, ["trk_truthPdgId"], 211.0)
sel_PionNeg = create_selection_function(ParticlePDGID, ["trk_truthPdgId"], -211.0)
sel_KaonPos = create_selection_function(ParticlePDGID, ["trk_truthPdgId"], 321.0)
sel_KaonNeg = create_selection_function(ParticlePDGID, ["trk_truthPdgId"], -321.0)
sel_ProtonPos = create_selection_function(ParticlePDGID, ["trk_truthPdgId"], 2212.0)
sel_ProtonNeg = create_selection_function(ParticlePDGID, ["trk_truthPdgId"], -2212.0)
sel_Other = create_inverse_selection_function([sel_PionPos, sel_PionNeg, sel_KaonPos, sel_KaonNeg, sel_ProtonPos, sel_ProtonNeg], name="sel_other")

pion_selections = [sel_Pion, sel_HardScatter]
pion_pos_selections = [sel_PionPos, sel_HardScatter]
pion_neg_selections = [sel_PionNeg, sel_HardScatter]
kaon_pos_selections = [sel_KaonPos, sel_HardScatter]
kaon_neg_selections = [sel_KaonNeg, sel_HardScatter]
proton_pos_selections = [sel_ProtonPos, sel_HardScatter]
proton_neg_selections = [sel_ProtonNeg, sel_HardScatter]
other_selections = [sel_Other, sel_HardScatter]
sel_Truth = [sel_HardScatter]
sel_TightIso = [sel_TightIso]

#This is a script that fills the histograms for
def fill_histograms(hist_filler, outputRootFileName):
    #import thje variables that we want to plot
    from variables import calc_trkNearestNeighbourEM2, calc_trkP, calc_EOP, calc_trkPt, calc_trkAverageMu, calc_trkEtaID, calc_trkEtaECAL, calc_trkNPV2, calc_trkCount, calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike, calc_TruthMomentum

    hist_filler.apply_selection_for_channel("LowMuDataTightIso", sel_TightIso) #Tighter isolation requirement
    hist_filler.apply_selection_for_channel("PythiaJetJetTightIso", sel_TightIso) #Tighter isolation requirement
    hist_filler.apply_selection_for_channel("PythiaJetJetHardScatter", sel_Truth) #Those tracks truth matched to pions
    hist_filler.apply_selection_for_channel("PythiaJetJetHardScatterTightIso", sel_Truth + sel_TightIso) #Tighter isolation requirement
    outFile = ROOT.TFile(outputRootFileName, "RECREATE")

    #count the number of tracks in each channel
    histogram_name = "trkCount"
    selections = []
    trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                         calc_trkCount,\
                                                         selections = selections,\
                                                         bins = 1,\
                                                         range_low = -0.5,\
                                                         range_high = +0.5,\
                                                         xlabel ='Always 0',\
                                                         ylabel = 'Number of Tracks')

    #count the number of events in each channel
    histogram_name = "EventCount"
    selections = [sel_Event]
    trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                         calc_trkCount,\
                                                         selections = selections,\
                                                         bins = 1,\
                                                         range_low = -0.5,\
                                                         range_high = +0.5,\
                                                         xlabel ='Always 0',\
                                                         ylabel = 'Number Events')

    histograms = hist_filler.DumpHistograms()
    for histogram_name in histograms:
        write_histograms(histograms[histogram_name], outFile)



    print("THEJOBFINISHED!")
