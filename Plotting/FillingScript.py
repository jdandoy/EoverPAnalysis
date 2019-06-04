#!/usr/bin/env python
# coding: utf-8

from HistogramFillingTools.HistogramFiller import HistogramFiller, get_p, get_bins, get_log_bins, create_selection_function
from HistogramFillingTools.EOPHistograms import PutBinningVectorsInFile, CreateEOPBinnedHistograms
from HistogramFillingTools.TrackSpectrumPlots import CreateTrackSpectrumPlots

import ROOT
import os
from math import pi
import pickle
import numpy as np
from selections.selections import EtaBin, PBin, sel_SubleadingTrack, sel_Event, sel_hasHADExtrapolation, NTRTX

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

from selections.selections import sel_HardScatter, ParticlePDGID_ABS
sel_Pion = create_selection_function(ParticlePDGID_ABS, ["trk_truthPdgId"], 211.0)
pion_selections = [sel_Pion, sel_HardScatter]

#This is a script that fills the histograms for
def fill_histograms(hist_filler, outputRootFileName):
    #import thje variables that we want to plot
    from variables.variables import calc_trkNearestNeighbourEM2, calc_trkP, calc_EOP, calc_trkPt, calc_trkAverageMu, calc_trkEtaID, calc_trkEtaECAL, calc_trkNPV2, calc_trkCount, calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike

    hist_filler.ApplySelectionsForChannel("PythiaJetJetPionsReweighted", pion_selections)

    #reweight the event count in MC to match the one from data
    event_count_reweight_file = ROOT.TFile("ReweightingHistograms/EventCountPythiaJetJetToData.root", "READ")
    hist = event_count_reweight_file.Get("EventCountPythiaJetJetToData")
    hist_filler.weightCalculator.addReweightHistogram("PythiaJetJet", calc_trkCount, hist, selection=[]) 

    trk_count_reweight_file = ROOT.TFile("ReweightingHistograms/PythiaJetJetPionsOnlyTrackCountReweightedToData.root", "READ")
    hist = trk_count_reweight_file.Get("PythiaJetJetPionsOnlyTrackCountReweightedToData")
    hist_filler.weightCalculator.addReweightHistogram("PythiaJetJetPionsReweighted", calc_trkCount, hist, selection=[])

    for i, eta_bin_selection in enumerate(eta_bin_selections):
        event_count_reweight_file = ROOT.TFile("ReweightingHistograms/PtSpectrumReweightLowMuDataOverPythiaJetJet_Eta"+str(i)+".root", "READ")
        hist = event_count_reweight_file.Get("PtSpectrumReweightLowMuDataOverPythiaJetJet_Eta"+str(i))
        hist_filler.weightCalculator.addReweightHistogram("PythiaJetJet", calc_trkPt, hist, selection=[eta_bin_selection]) 

        event_count_reweight_file = ROOT.TFile("ReweightingHistograms/PtSpectrumReweightLowMuDataOverSinglePion_Eta"+str(i)+".root", "READ")
        hist = event_count_reweight_file.Get("PtSpectrumReweightLowMuDataOverSinglePion_Eta"+str(i))
        hist_filler.weightCalculator.addReweightHistogram("SinglePion", calc_trkPt, hist, selection=[eta_bin_selection]) 

        event_count_reweight_file = ROOT.TFile("ReweightingHistograms/PtSpectrumReweightLowMuDataOverPythiaJetJetPionsReweighted_Eta"+str(i)+".root", "READ")
        hist = event_count_reweight_file.Get("PtSpectrumReweightLowMuDataOverPythiaJetJetPionsReweighted_Eta"+str(i))
        hist_filler.weightCalculator.addReweightHistogram("PythiaJetJetPionsReweighted", calc_trkPt, hist, selection=[eta_bin_selection]) 

    #import the selections that we want to plot
    outFile = ROOT.TFile(outputRootFileName, "RECREATE")

    ##just count the number of tracks in each histogram
    histogram_name = "trkCount"
    selections = []
    trkCountHist = hist_filler.BookHistograms(histogram_name,\
                                                         calc_trkCount,\
                                                         list_selections = selections,\
                                                         bins = 1,\
                                                         range_low = -0.5,\
                                                         range_high = +0.5,\
                                                         xlabel ='Always 0',\
                                                         ylabel = 'Number of Tracks')

    ##just count the number of events 
    histogram_name = "EventCount"
    selections = [sel_Event]
    trkCountHist = hist_filler.BookHistograms(histogram_name,\
                                                         calc_trkCount,\
                                                         list_selections = selections,\
                                                         bins = 1,\
                                                         range_low = -0.5,\
                                                         range_high = +0.5,\
                                                         xlabel ='Always 0',\
                                                         ylabel = 'Number Events')

    ################################################################################
    #plot the first set of variables that we're interested in
    #plot the trk avaerage mu histogram
    histogram_name = "trkAverageMu"
    selections = []
    trkAverageMuHist = hist_filler.BookHistograms(histogram_name,\
                                                        calc_trkAverageMu,\
                                                        list_selections = selections,\
                                                        bins = 10,\
                                                        range_low = 0.0,\
                                                        range_high = 10.0,\
                                                        xlabel ='Average #mu of Event',\
                                                        ylabel = 'Number of Tracks')

   ################################################################################
   #plot a histogram of the average event NPV
    histogram_name = "trkNPV2"
    trkNPV2Hist = hist_filler.BookHistograms(histogram_name,\
                                       calc_trkNPV2,\
                                       list_selections = [],\
                                       bins = 13,\
                                       range_low = -0.5,\
                                       range_high = 12.5,\
                                       xlabel ="NPV with 2 Tracks",\
                                       ylabel = "Number of Tracks")


#   ################################################################################
   #plot a histogram of the average event NPV
    histogram_name = "eventNPV2Hist"
    eventNPV2Hist = hist_filler.BookHistograms(histogram_name,\
                                       calc_trkNPV2,\
                                       list_selections = [sel_Event],\
                                       bins = 13,\
                                       range_low = -0.5,\
                                       range_high = 12.5,\
                                       xlabel ="NPV with 2 Tracks",\
                                       ylabel = "Number Events")

#   ################################################################################
    histogram_name = "eventAverageMu"
    selections = [sel_Event]
    eventAverageMuHist = hist_filler.BookHistograms(histogram_name,
                                          calc_trkAverageMu,\
                                          list_selections = selections,\
                                          bins = 10,\
                                          range_low = 0.0,\
                                          range_high = 10.0,\
                                          xlabel ='<#mu>',\
                                          ylabel = 'Number of Events')

#    ################################################################################
    #prepare the momentum bins
    binMax = 30.0
    binMin = 0.5
    nBins = 100
    p_bins = get_log_bins(binMin, binMax, nBins)
    p_bins_reference = p_bins
    histogram_name = "trkPtHist"
    hist_filler.BookHistograms(histogram_name,\
                                       calc_trkPt,\
                                       list_selections = [],\
                                       bins = p_bins,\
                                       xlabel ="Track P_{T} [GeV]",\
                                       ylabel = "Number of Tracks")

    hist_filler.BookHistograms(histogram_name + "HasExtrapolation",\
                                       calc_trkPt,\
                                       list_selections = [sel_hasHADExtrapolation],\
                                       bins = p_bins,\
                                       xlabel ="Track P_{T} [GeV]",\
                                       ylabel = "Number of Tracks")

#    ################################################################################
    #prepare the momentum bins
    binMax = 30.0
    binMin = 0.5
    nBins = 50
    p_bins = get_log_bins(binMin, binMax, nBins)
    p_bins_reference = p_bins
    histogram_name = "LeadingPtTrkHist"
    trkPtHistZoom = hist_filler.BookHistograms(histogram_name,\
                                       calc_trkPt,\
                                       list_selections = [sel_Event],\
                                       bins = p_bins,\
                                       xlabel ="Track P_{T} [GeV]",\
                                       ylabel = "Number of Tracks")

#    ################################################################################
    #prepare the momentum bins
    binMax = 30.0
    binMin = 0.5
    nBins = 50
    p_bins = get_log_bins(binMin, binMax, nBins)
    p_bins_reference = p_bins
    histogram_name = "SubleadingPtTrkHist"
    trkPtHistZoom = hist_filler.BookHistograms(histogram_name,\
                                       calc_trkPt,\
                                       list_selections = [sel_SubleadingTrack],\
                                       bins = p_bins,\
                                       xlabel ="Track P_{T} [GeV]",\
                                       ylabel = "Number of Tracks")

  
    from variables.variables import calc_trkEtaECAL, calc_trkPhiECAL
    from selections.selections import PBin, sel_NonZeroEnergy
    from calculation.calculation import calculation
    p_bin_selection = create_selection_function(PBin, ["trk_p"], 3., 4.)
    histogram_name = "TrkEtaPhiEMCal_MomentumBetween3And4GeV_Denomenator"
    hist_filler.Book2DHistograms(histogram_name,\
                            calc_trkEtaECAL,\
                            calc_trkPhiECAL,\
                            list_selections=[p_bin_selection],\
                            bins_x = 200,\
                            bins_y = 200,\
                            range_low_y = -3.14,\
                            range_high_y = +3.14,\
                            range_low_x = -2.5,\
                            range_high_x = +2.5,\
                            xlabel = "#eta_{EMCal}",\
                            ylabel = "#phi_{EMCal}",\
                            zlabel = "Number of Tracks")

    histogram_name = "TrkEtaPhiEMCal_MomentumBetween3And4GeV_Numerator"
    hist_filler.Book2DHistograms(histogram_name,\
                            calc_trkEtaECAL,\
                            calc_trkPhiECAL,\
                            list_selections=[p_bin_selection, sel_NonZeroEnergy],\
                            bins_x = 200,\
                            bins_y = 200,\
                            range_low_y = -3.14,\
                            range_high_y = +3.14,\
                            range_low_x = -2.5,\
                            range_high_x = +2.5,\
                            xlabel = "#eta_{EMCal}",\
                            ylabel = "#phi_{EMCal}",\
                            zlabel = "Number of Tracks")


    #try between 2 and 3 GeV

    p_bin_selection = create_selection_function(PBin, ["trk_p"], 2., 3.)
    histogram_name = "TrkEtaPhiEMCal_MomentumBetween2And3GeV_Denomenator"
    hist_filler.Book2DHistograms(histogram_name,\
                            calc_trkEtaECAL,\
                            calc_trkPhiECAL,\
                            list_selections=[p_bin_selection],\
                            bins_x = 200,\
                            bins_y = 200,\
                            range_low_y = -3.14,\
                            range_high_y = +3.14,\
                            range_low_x = -2.5,\
                            range_high_x = +2.5,\
                            xlabel = "#eta_{EMCal}",\
                            ylabel = "#phi_{EMCal}",\
                            zlabel = "Number of Tracks")

    histogram_name = "TrkEtaPhiEMCal_MomentumBetween2And3GeV_Numerator"
    hist_filler.Book2DHistograms(histogram_name,\
                            calc_trkEtaECAL,\
                            calc_trkPhiECAL,\
                            list_selections=[p_bin_selection, sel_NonZeroEnergy],\
                            bins_x = 200,\
                            bins_y = 200,\
                            range_low_y = -3.14,\
                            range_high_y = +3.14,\
                            range_low_x = -2.5,\
                            range_high_x = +2.5,\
                            xlabel = "#eta_{EMCal}",\
                            ylabel = "#phi_{EMCal}",\
                            zlabel = "Number of Tracks")


    ################################################################################
    histogramName = "TrackEtaID"
    hist_filler.BookHistograms(histogramName,\
                                       calc_trkEtaID,\
                                       list_selections = [],\
                                       bins = 100,\
                                       range_low = -5,\
                                       range_high = +5,\
                                       xlabel ="Track #eta ID",\
                                       ylabel = "Number of Tracks")

#   ################################################################################
    histogramName = "TwoDTrackPvsTrkEtaID"
    max_bin = 2.4
    min_bin = -2.4
    nBins = 48
    eta_bins = get_bins(min_bin, max_bin, nBins)
    hist_filler.Book2DHistograms(histogramName,\
                                             calc_trkEtaID,\
                                             calc_trkP,\
                                             list_selections=[],\
                                             bins_x=eta_bins,\
                                             xlabel="Track #eta ID",\
                                             bins_y=p_bins,\
                                             ylabel="Track P [GeV]",\
                                             zlabel="Number of Tracks",\
                                             )

#   ################################################################################
    histogramName = "TwoDTrackPtVsEtaHistogram"
    hist_filler.Book2DHistograms(histogramName,\
                                             calc_trkEtaID,\
                                             calc_trkPt,\
                                             list_selections=[],\
                                             bins_x=eta_bins,\
                                             xlabel="Track #eta ID",\
                                             bins_y=p_bins,\
                                             ylabel="Track P_{T} [GeV]",\
                                             zlabel="Number of Tracks",\
                                             )

#   ################################################################################
    histogramName = "TwoDTrackPtVsEtaHistogram_HasExtrapolation"
    hist_filler.Book2DHistograms(histogramName,\
                                             calc_trkEtaID,\
                                             calc_trkPt,\
                                             list_selections=[sel_hasHADExtrapolation],\
                                             bins_x=eta_bins,\
                                             xlabel="Track #eta ID",\
                                             bins_y=p_bins,\
                                             ylabel="Track P_{T} [GeV]",\
                                             zlabel="Number of Tracks",\
                                             )



#   ################################################################################
    histogramName = "trkEtaECALHist"
    hist_filler.BookHistograms(histogramName,\
                                          calc_trkEtaECAL,
                                          list_selections = [],
                                          bins = 100,
                                          range_low = -5,
                                          range_high = +5,
                                          xlabel ="Track #eta EM Layer 2",
                                          ylabel = "Number of Tracks",
                                       )

#   ################################################################################
    histogramName = "TwoDHistTrkPvsPhiInnerToExtrapolEM2"
    dPhi_bins = []
    min_bin = 0.0
    max_bin = pi
    NBins = 100
    dPhi_bins = get_bins(min_bin, max_bin, NBins)

    from variables.variables import calc_trkDPhi
    hist_filler.Book2DHistograms(histogramName,\
                                             calc_trkDPhi,\
                                             calc_trkPt,\
                                             list_selections=[],\
                                             bins_x=dPhi_bins,\
                                             xlabel="|#phi_{ID} - #phi_{EM2}|",\
                                             bins_y=p_bins,\
                                             ylabel="Track P_{T} [GeV]",\
                                             zlabel="Number of Tracks",\
                                             )

#   ################################################################################
    histogramName = "lowPTLess07_TwoDHistTrkEtavsDEtaInnerToExtrapolEM2"
    from variables.variables import calc_trkDEta
    from calculation.calculation import calculation
    def lowPT(trk):
        return trk["trk_pt"] < 0.7
    branches =["trk_pt"]
    sel_lowPT = calculation(lowPT, branches)
    hist_filler.Book2DHistograms(histogramName,\
                                             calc_trkEtaID,\
                                             calc_trkDEta,\
                                             list_selections=[sel_lowPT],\
                                             bins_x=50,\
                                             range_low_x=-2.5,\
                                             range_high_x=+2.5,\
                                             xlabel="Track #eta_{ID}",\
                                             bins_y=50,\
                                             range_low_y=0.0,\
                                             range_high_y=1.0,\
                                             ylabel="|#eta_{ID} - #eta_{EM2}|",\
                                             zlabel="Number of Tracks",\
                                             )

#   ################################################################################
    from calculation.calculation import calculation
    from selections.selections import EtaBin
    sel_Eta00_08 = create_selection_function(EtaBin, ["trk_etaID"], 0.0, 0.8)
    histogramName = "EtaLess08_TwoDHistTrkPvsPhiInnerToExtrapolEM2"
    hist_filler.Book2DHistograms(histogramName,\
                                             calc_trkDPhi,\
                                             calc_trkPt,\
                                             list_selections=[sel_Eta00_08],\
                                             bins_x=dPhi_bins,\
                                             xlabel="|#phi_{ID} - #phi_{EM2}|",\
                                             bins_y=p_bins,\
                                             ylabel="Track P_{T} [GeV]",\
                                             zlabel="Number of Tracks",\
                                             )

#    ################################################################################
    histogramName = "NearestDRHist"
    hist_filler.BookHistograms(histogramName,
                                    calc_trkNearestNeighbourEM2,
                                    list_selections = [],
                                    bins = 25,
                                     range_low = 0.0,
                                     range_high = 5,
                                     xlabel ="dR to Nearest Track",
                                    ylabel = "Number of Tracks",
                                    )

    from selections.selections import  sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum
    sel_NTRT20 = create_selection_function(NTRTX, ["trk_nTRT"], 20.0, 100000.0)
    MIP_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum]

    ################################################################################
    selections = []
    histogramName = "InclusiveEOP"
    hist_filler.BookHistograms(histogramName,\
                                     calc_EOP,
                                     list_selections = selections,
                                     bins = 50,
                                     range_low = -1,
                                     range_high = 5,
                                     xlabel ="E/p",
                                     ylabel = "Number of Tracks",
                                     )


    ################################################################################
    from selections.selections import sel_NonZeroEnergy
    selections = [sel_NonZeroEnergy]
    histogramName = "NonZeroEnergy_InclusiveEOP"
    hist_filler.BookHistograms(histogramName,\
                                     calc_EOP,
                                     list_selections = selections,
                                     bins = 50,
                                     range_low = -1,
                                     range_high = 5,
                                     xlabel ="E/p",
                                     ylabel = "Number of Tracks",
                                     )



    ################################################################################
    # This is figure 3a in the paper:

    selections = []
    binMax = 10.05
    binLow = 0.5
    nBins = 15
    bins = get_log_bins(binLow, binMax, nBins)

    histogramName = "InclusiveZeroFractionVsPDenomenator"
    hist_filler.BookHistograms(histogramName,\
                                          calc_trkP,\
                                          list_selections = selections,\
                                          bins = bins,\
                                          xlabel ="Track P [GeV]",\
                                          ylabel = "Number of Tracks",\
                                          )

    from selections.selections import sel_ELessEqual0
    histogramName = "InclusiveZeroFractionVsPNumerator"
    selections = [sel_ELessEqual0]
    hist_filler.BookHistograms(histogramName,\
                                                    calc_trkP,\
                                                    list_selections = selections,\
                                                    bins = bins,\
                                                    xlabel ="Track P [GeV]",\
                                                    ylabel = "N(E<=0)/N",\
                                                    )

    ################################################################################
    #This is figure 3b of the paper
    bins = [-2.3, -1.8, -1.5, -1.4, -1.1, -0.6, 0.0, 0.6, 1.1, 1.4, 1.5, 1.8, 2.3]
    selections = []
    histogramName = "InclusiveZeroFractionVsEtaDenomenator"
    hist_filler.BookHistograms(histogramName,\
                                              calc_trkEtaID,\
                                              list_selections = selections,\
                                              bins = bins,\
                                              xlabel ="Track |#eta|",\
                                              ylabel = "Number of Tracks",\
                                              )
    from selections.selections import sel_ELessEqual0
    histogramName = "InclusiveZeroFractionVsEtaNumerator"
    selections = [sel_ELessEqual0]
    hist_filler.BookHistograms(histogramName,\
                                                   calc_trkEtaID,\
                                                   list_selections = selections,\
                                                   bins = bins,\
                                                   xlabel ="Track |#eta|",\
                                                   ylabel = "N(E<=0)/N",\
                                                   )

    ################################################################################
    bins = [0.0, 0.6, 1.1, 1.4, 1.5, 1.8, 2.3]
    from variables.variables import calc_trkEta_ABS
    selections = []
    histogramName = "InclusiveZeroFractionVsAbsEtaDenomenator"
    hist_filler.BookHistograms(histogramName,\
                                              calc_trkEta_ABS,\
                                              list_selections = selections,\
                                              bins = bins,\
                                              xlabel ="Track |#eta|",\
                                              ylabel = "Number of Tracks",\
                                              )
    from selections.selections import sel_ELessEqual0
    histogramName = "InclusiveZeroFractionVsAbsEtaNumerator"
    selections = [sel_ELessEqual0]
    hist_filler.BookHistograms(histogramName,\
                                                   calc_trkEta_ABS,\
                                                   list_selections = selections,\
                                                   bins = bins,\
                                                   xlabel ="Track |#eta|",\
                                                   ylabel = "N(E<=0)/N",\
                                                   )

    ################################################################################
    from variables.variables import calc_trkEta_ABS
    from selections.selections import EtaBin

    canvases = []
    keep_histograms_alive = []
    for (eta_bin_selection, eta_bin_description, eta_bin_edge) in zip(eta_bin_selections, eta_bin_descriptions,eta_bin_tuples):
        center = (eta_bin_edge[0] + eta_bin_edge[1])/2.0
        binMax = 30.05
        binLow = get_p(0.5, center)
        nBins = 20
        bins = get_log_bins(binLow, binMax, nBins)

        #do the eta selection and count the inclusive number of tracks in the bin
        selections = [eta_bin_selection]
        histogramName = "ZeroFractionVsP" + eta_bin_description + "Denomenator"
        hist_filler.BookHistograms(histogramName,\
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = bins,\
                                                  xlabel ="Track P [GeV]",\
                                                  ylabel = "Number of tracks",\
                                                  )

        #do the eta selections and count the number of tracks with an energy deposity less than or equal to 0.0.
        from selections.selections import sel_ELessEqual0
        selections = [sel_ELessEqual0] + [eta_bin_selection]
        histogramName = "ZeroFractionVsP" + eta_bin_description + "Numerator"
        hist_filler.BookHistograms(histogramName,\
                                                       calc_trkP,\
                                                       list_selections = selections,\
                                                       bins = bins,\
                                                       xlabel ="Track P [GeV]",\
                                                       ylabel = "N(E<=0)/N",\
                                                       )

    ################################################################################
    from variables.variables import calc_trkEta_ABS
    base_description = ["N_{TRT hits} >= 20"]

    for (eta_bin_selection, eta_bin_description, eta_bin_edge) in zip(eta_bin_selections, eta_bin_descriptions, eta_bin_tuples):
        center = (eta_bin_edge[0] + eta_bin_edge[1])/2.0
        binMax = 30.05
        binLow = get_p(0.5, center)
        nBins = 20
        bins = get_log_bins(binLow, binMax, nBins)

        #do the eta selection and count the inclusive number of tracks in the bin
        selections = [eta_bin_selection] + [sel_NTRT20]
        histogramName = "NTRT20ZeroFractionVsP" + eta_bin_description + "Denomenator"
        hist_filler.BookHistograms(histogramName,\
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = bins,\
                                                  xlabel ="Track P [GeV]",\
                                                  ylabel = "Number of tracks",\
                                                  )

        #do the eta selections and count the number of tracks with an energy deposity less than or equal to 0.0.
        from selections.selections import sel_ELessEqual0
        selections = [sel_ELessEqual0] + [eta_bin_selection] + [sel_NTRT20]
        histogramName = "NTRT20ZeroFractionVsP" + eta_bin_description + "Numerator"
        hist_filler.BookHistograms(histogramName,\
                                                       calc_trkP,\
                                                       list_selections = selections,\
                                                       bins = bins,\
                                                       xlabel ="Track P [GeV]",\
                                                       ylabel = "N(E<=0)/N",\
                                                       )

    ################################################################################
    ##Create a set of p and eta bins for the measurement ##########################

    from calculation.calculation import calculation
    from selections.selections import EtaBin, PBin
    from variables.variables import calc_EOPBkg, calc_EnergyAnulus
    from selections.selections import  sel_NonZeroEnergy, sel_HardScatter

    ##select inclusive distributions
    ##Create a set of binned EOP response histograms 
    eta_ranges = eta_bin_tuples
    base_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum]
    p_bins_for_eta_range = []
    for eta_range in eta_bin_tuples:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
    description = "MIPSelectionBetween30and90OfMomentum"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    from selections.selections import sel_EHadFracAbove70, sel_Lar1_1GeV
    eta_ranges = eta_bin_tuples
    base_selection = [sel_EHadFracAbove70, sel_NTRT20, sel_Lar1_1GeV]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
    description = "MIPSelectionHadFracAbove70"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = eta_bin_tuples
    base_selection = [sel_NTRT20, sel_NonZeroEnergy]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
    description = "20TRTHitsNonZeroEnergy"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 300)
        p_bins_for_eta_range.append(p_bins)
    CreateTrackSpectrumPlots(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = eta_bin_tuples
    base_selection = [sel_NTRT20]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
    description = "20TRTHits"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 300)
        p_bins_for_eta_range.append(p_bins)
    CreateTrackSpectrumPlots(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = eta_bin_tuples
    base_selection = [sel_NonZeroEnergy]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
    description = "NonZeroEnergy"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 300)
        p_bins_for_eta_range.append(p_bins)
    CreateTrackSpectrumPlots(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = eta_bin_tuples
    base_selection = []
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
    description = "Inclusive"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 300)
        p_bins_for_eta_range.append(p_bins)
    CreateTrackSpectrumPlots(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    # Crete histograms for the hard scatter histograms
    eta_ranges = eta_bin_tuples
    base_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum, sel_HardScatter]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
    description = "MIPSelectionBetween30and90OfMomentumHardScatter"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = eta_bin_tuples
    from selections.selections import sel_EHadFracAbove70, sel_Lar1_1GeV
    base_selection = [sel_EHadFracAbove70, sel_NTRT20, sel_Lar1_1GeV, sel_HardScatter]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
   
    description = "MIPSelectionHadFracAbove70HardScatter"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    base_selection = [sel_NTRT20, sel_NonZeroEnergy, sel_HardScatter]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
    description = "20TRTHitsNonZeroEnergyHardScatter"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = eta_bin_tuples
    base_selection = [sel_NonZeroEnergy, sel_HardScatter]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
    description = "NonZeroEnergyHardScatter"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = eta_bin_tuples
    base_selection = [sel_HardScatter]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = get_p(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = get_log_bins(p_bins_min, 30.05, 15)
        p_bins_for_eta_range.append(p_bins)
    description = "InclusiveHardScatter"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    histograms = hist_filler.DumpHistograms()
    for histogram_name in histograms:
        write_histograms(histograms[histogram_name], outFile)

    print("THEJOBFINISHED!")
