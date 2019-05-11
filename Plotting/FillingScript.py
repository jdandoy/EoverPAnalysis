#!/usr/bin/env python
# coding: utf-8

from HistogramFillingTools.HistogramFiller import HistogramFiller, getP, getBins, getLogBins
from HistogramFillingTools.EOPHistograms import PutBinningVectorsInFile, CreateEOPBinnedHistograms

import ROOT
import os
from math import pi
import pickle
import numpy as np
from calculation.calculation import calculation
from selections.selections import EtaBin, PBin, sel_SubleadingTrack

def WriteToFile(histogram_dictionary, outFile):
    outFile.cd()
    for key in histogram_dictionary:
        if not outFile.cd(key):
            outFile.mkdir(key)
        outFile.cd(key)
        histogram_dictionary[key].Write()



#This is a script that fills the histograms for
def FillingScript(hist_filler, outputRootFileName):
    #import thje variables that we want to plot
    from variables.variables import calc_trkNearestNeighbourEM2, calc_trkP, calc_EOP, calc_trkPt, calc_trkAverageMu, calc_trkEtaID, calc_trkEtaECAL, calc_trkNPV2, calc_trkCount, calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike
    #import the selections that we want to plot
    from selections.selections import sel_NTRT20, sel_NTRT25, sel_NTRT30, sel_ECALEta0_6, sel_PGreater1 ,sel_PGreater1_5, sel_PGreater2, sel_PGreater2_5, sel_Z0SinThetaLess1_5, sel_d0Less1_5, sel_Event, sel_hasHADExtrapolation
    #impot the ID selections that we want to plot
    from selections.selections import sel_IDEta00_06,sel_IDEta06_11,sel_IDEta11_14,sel_IDEta14_15, sel_IDEta15_18, sel_IDEta18_23, EtaBin

    outFile = ROOT.TFile(outputRootFileName, "RECREATE")

    histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
    xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
    variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]
    selections = []

    for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
        hist_filler.BookHistograms(histogram_name,\
                               variable,\
                               list_selections = selections,\
                               bins = 10,\
                               range_low = -0.5,\
                               range_high = 9.5,\
                               xlabel=xlabel,\
                               ylabel="Number of Tracks")

    histogram_name = "trkTRTHits"
    from variables.variables import calc_nTRT
    hist_filler.BookHistograms(histogram_name,\
                           calc_nTRT,\
                           list_selections = selections,\
                           range_low = -0.5,\
                           range_high = 59.5,\
                           bins = 60,\
                           xlabel = "Number of TRT Hits",\
                           ylabel = "Number of Tracks")

    histogram_name = "trkEMDR100"
    from variables.variables import calc_EnergyEMDR100
    hist_filler.BookHistograms(histogram_name,\
                           calc_EnergyEMDR100,\
                           list_selections = selections,\
                           range_low = -2.0,\
                           range_high = + 10.0,\
                           bins = 48,\
                           xlabel = "E_{EM}^{#DeltaR<0.1}[GeV]",\
                           ylabel = "Number of Tracks")

    histogram_name = "MomentumHadFrac"
    from variables.variables import calc_MomentumHadFrac
    hist_filler.BookHistograms(histogram_name,\
                           calc_MomentumHadFrac,\
                           list_selections = selections,\
                           range_low = -1.0,\
                           range_high = + 5.0,\
                           bins = 48,\
                           xlabel = "E^{HAD}/P",\
                           ylabel = "Number of Tracks")

    histogram_name = "HadFrac"
    from variables.variables import calc_HadFrac
    hist_filler.BookHistograms(histogram_name,\
                           calc_HadFrac,\
                           list_selections = selections,\
                           range_low = -1.0,\
                           range_high = + 2.0,\
                           bins = 48,\
                           xlabel = "E^{HAD}/E^{Total}",\
                           ylabel = "Number of Tracks")


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
    p_bins = getLogBins(binMin, binMax, nBins)
    p_bins_reference = p_bins 
    histogram_name = "trkPtHist"
    trkPtHistZoom = hist_filler.BookHistograms(histogram_name,\
                                       calc_trkPt,\
                                       list_selections = [],\
                                       bins = p_bins,\
                                       xlabel ="Track P_{T} [GeV]",\
                                       ylabel = "Number of Tracks")

    trkPtHistZoom = hist_filler.BookHistograms(histogram_name + "HasExtrapolation",\
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
    p_bins = getLogBins(binMin, binMax, nBins)
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
    p_bins = getLogBins(binMin, binMax, nBins)
    p_bins_reference = p_bins
    histogram_name = "SubleadingPtTrkHist"
    trkPtHistZoom = hist_filler.BookHistograms(histogram_name,\
                                       calc_trkPt,\
                                       list_selections = [sel_SubleadingTrack],\
                                       bins = p_bins,\
                                       xlabel ="Track P_{T} [GeV]",\
                                       ylabel = "Number of Tracks")

#           ################################################################################
#           ## Look in different bins of pseudorapidity
    base_description = []
    etaSelections = [sel_IDEta00_06,\
                    sel_IDEta06_11,\
                    sel_IDEta11_14,\
                    sel_IDEta14_15,\
                    sel_IDEta15_18,\
                    sel_IDEta18_23]

    eta_selectionDescriptions = [\
                              "|#eta_{ID}|<0.6",\
                              "0.6<|#eta_{ID}|<1.1",\
                              "1.1<|#eta_{ID}|<1.4",\
                              "1.4<|#eta_{ID}|<1.5",\
                              "1.5<|#eta_{ID}|<1.8",\
                              "1.8<|#eta_{ID}|<2.3"\
                              ]

    file_descriptions = ["eta06", "eta06_11", "eta11_14", "eta14_15", "eta15_18", "eta18_23"]
  
    from variables.variables import calc_trkEtaECAL, calc_trkPhiECAL
    from selections.selections import PBin, sel_NonZeroEnergy
    from calculation.calculation import calculation
    PBinFunction = lambda x, y=3, z=4: PBin(x, y,z)
    PBinFunction.__name__ = "SelMomentumRange3_4"
    PBinSelection = calculation(PBinFunction, ["trk_p"])
    histogram_name = "TrkEtaPhiEMCal_MomentumBetween3And4GeV_Denomenator"
    hist_filler.Book2DHistograms(histogram_name,\
                            calc_trkEtaECAL,\
                            calc_trkPhiECAL,\
                            list_selections=[PBinSelection],\
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
                            list_selections=[PBinSelection, sel_NonZeroEnergy],\
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

    PBinFunction = lambda x, y=2, z=3: PBin(x, y,z)
    PBinFunction.__name__ = "SelMomentumRange2_3"
    PBinSelection = calculation(PBinFunction, ["trk_p"])
    histogram_name = "TrkEtaPhiEMCal_MomentumBetween2And3GeV_Denomenator"
    hist_filler.Book2DHistograms(histogram_name,\
                            calc_trkEtaECAL,\
                            calc_trkPhiECAL,\
                            list_selections=[PBinSelection],\
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
                            list_selections=[PBinSelection, sel_NonZeroEnergy],\
                            bins_x = 200,\
                            bins_y = 200,\
                            range_low_y = -3.14,\
                            range_high_y = +3.14,\
                            range_low_x = -2.5,\
                            range_high_x = +2.5,\
                            xlabel = "#eta_{EMCal}",\
                            ylabel = "#phi_{EMCal}",\
                            zlabel = "Number of Tracks")

    for (etaSelection, eta_selectionDescription, file_description) in zip(etaSelections, eta_selectionDescriptions, file_descriptions):
        #do the eta selection and count the inclusive number of tracks in the bin
        selections = [etaSelection]
        hist_filler.BookHistograms("TrkPtHist"+file_description,\
                                                  calc_trkPt,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="Track P_{T} [GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  )



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
    eta_bins = getBins(min_bin, max_bin, nBins)
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
    dPhi_bins = getBins(min_bin, max_bin, NBins)

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
    Eta00_08 = lambda x, y=0.0,z=0.8 : EtaBin(x, y,z)
    sel_Eta00_08 = calculation(Eta00_08, ["trk_etaID"])
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

    from selections.selections import sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum, sel_PGreater2, sel_PGreater2_5, sel_PGreater3
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
    selections = [sel_PGreater1, sel_ECALEta0_6]
    histogramName = "EtaID0_6_PGreater1_0_EOPHist"
    hist_filler.BookHistograms(histogramName,\
                                               calc_EOP,
                                               list_selections = selections,
                                               bins = 50,
                                               range_low = -1,
                                               range_high = 5,
                                               xlabel ="E/p",
                                               ylabel = "Number of Tracks")

    ################################################################################
    selections = [sel_PGreater1_5, sel_IDEta00_06]
    histogramName = "EtaID0_6_PGreater1_5_EOPHist"
    hist_filler.BookHistograms(histogramName,
                                                calc_EOP,
                                                list_selections = selections,
                                                 bins = 50,
                                                 range_low = -1,
                                                 range_high = 5,
                                                 xlabel ="E/p",
                                                 ylabel = "Number of Tracks")

    ################################################################################
    selections = [sel_PGreater2, sel_IDEta00_06]
    histogramName = "EtaID0_6_PGreater2_0_EOPHist"
    hist_filler.BookHistograms(histogramName,
                                              calc_EOP,
                                              list_selections = selections,
                                              bins = 50,
                                               range_low = -1,
                                              range_high = 5,
                                              xlabel ="E/p",
                                               ylabel = "Number of Tracks",
                                              )


    ################################################################################
    from selections.selections import sel_IDEta19_23, sel_IDEta00_06, sel_PBetween12_18, sel_PBetween22_28, sel_PBetween28_36
    # This is figure 2a and 2d in the paper:
    selections = [sel_PBetween12_18, sel_IDEta00_06]
    histogramName = "EtaID0_6_PBetween12_18_EOPHist"
    hist_filler.BookHistograms(histogramName,
                                              calc_EOP,
                                              list_selections = selections,
                                              bins = 50,
                                               range_low = -1,
                                              range_high = 5,
                                              xlabel ="E/p",
                                               ylabel = "Number of Tracks",
                                              )

    ################################################################################
    from selections.selections import sel_PBetween22_28
    histogramName = "EtaID0_6_PBetween22_28_EOPHist"
    selections = [sel_PBetween22_28, sel_IDEta00_06]
    hist_filler.BookHistograms(histogramName,\
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.75,\
                                          range_high = 4,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                         )

    ################################################################################
    # This is figure 2c in the paper:
    selections = [sel_PBetween28_36, sel_IDEta19_23]
    histogramName = "EtaIDBetween19_23_PBetween28_36_EOPHist"
    hist_filler.BookHistograms(histogramName,\
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.75,\
                                          range_high = 4,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                          )

    ################################################################################
    # This is figure 2c in the paper:
    selections = [sel_PBetween28_36, sel_IDEta00_06]
    histogramName = "EtaID0_6_PBetween28_36_EOPHist"
    hist_filler.BookHistograms(histogramName,
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.5,\
                                          range_high = 3,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                          )

    ################################################################################
    from selections.selections import sel_NonZeroEnergy
    # This is figure 2a and 2d in the paper:
    histogramName = "NonZero_EtaID0_6_PBetween12_18_EOPHist"
    selections = [sel_PBetween12_18, sel_IDEta00_06, sel_NonZeroEnergy]
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
    from selections.selections import sel_PBetween22_28
    selections = [sel_PBetween22_28, sel_IDEta00_06, sel_NonZeroEnergy]
    histogramName = "NonZero_EtaID0_6_PBetween22_28_EOPHist"
    hist_filler.BookHistograms(histogramName,\
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.75,\
                                          range_high = 4,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                         )

    ################################################################################
    # This is figure 2c in the paper:
    selections = [sel_PBetween28_36, sel_IDEta19_23, sel_NonZeroEnergy]
    histogramName = "NonZero_EtaIDBetween19_23_PBetween22_28_EOPHist"
    hist_filler.BookHistograms(histogramName,
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.75,\
                                          range_high = 4,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                          )

    ################################################################################
    # This is figure 2c in the paper:
    selections = [sel_PBetween28_36, sel_IDEta00_06, sel_NonZeroEnergy]
    histogramName = "NonZero_EtaID0_6_PBetween28_36_EOPHist"
    hist_filler.BookHistograms(histogramName,
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.5,\
                                          range_high = 3,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                          )


    ################################################################################
    # This is figure 3a in the paper:

    selections = []
    binMax = 10.05
    binLow = 0.5
    nBins = 15
    bins = getLogBins(binLow, binMax, nBins)

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
    from selections.selections import sel_IDEta00_02, sel_IDEta02_04, sel_IDEta04_06, sel_IDEta00_06

    etaSelections = [sel_IDEta00_02,\
                     sel_IDEta02_04,\
                     sel_IDEta04_06,\
                     sel_IDEta00_06,\
                     sel_IDEta06_11,\
                     sel_IDEta11_14,\
                     sel_IDEta14_15,\
                     sel_IDEta15_18,\
                     sel_IDEta18_23]


    canvases = []
    keep_histograms_alive = []

    file_descriptions = ["etaID00_02", "etaID02_04", "etaID04_06", "etaID00_06", "etaID06_11", "etaID11_14", "etaID14_15", "etaID15_18", "etaID18_23"]
    centers = [0.2, 0.4, 0.6, 0.6, 1.1, 1.4, 1.5, 1.8, 2.3]

    for (etaSelection, eta_selectionDescription, file_description, center) in zip(etaSelections, eta_selectionDescriptions, file_descriptions, centers):
        binMax = 15.05
        binLow = getP(0.5, center)
        nBins = 20
        bins = getLogBins(binLow, binMax, nBins)

        #do the eta selection and count the inclusive number of tracks in the bin
        selections = [etaSelection]
        histogramName = "ZeroFractionVsP" + file_description + "Denomenator"
        hist_filler.BookHistograms(histogramName,\
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = bins,\
                                                  xlabel ="Track P [GeV]",\
                                                  ylabel = "Number of tracks",\
                                                  )

        #do the eta selections and count the number of tracks with an energy deposity less than or equal to 0.0.
        from selections.selections import sel_ELessEqual0
        selections = [sel_ELessEqual0] + [etaSelection]
        histogramName = "ZeroFractionVsP" + file_description + "Numerator"
        hist_filler.BookHistograms(histogramName,\
                                                       calc_trkP,\
                                                       list_selections = selections,\
                                                       bins = bins,\
                                                       xlabel ="Track P [GeV]",\
                                                       ylabel = "N(E<=0)/N",\
                                                       )

    ################################################################################
    from variables.variables import calc_trkEta_ABS
    from selections.selections import sel_NTRT20
    base_description = ["N_{TRT hits} >= 20"]

    for (etaSelection, eta_selectionDescription, file_description, center) in zip(etaSelections, eta_selectionDescriptions, file_descriptions, centers):
        binMax = 15.05
        binLow = getP(0.5, center)
        nBins = 20
        bins = getLogBins(binLow, binMax, nBins)

        #do the eta selection and count the inclusive number of tracks in the bin
        selections = [etaSelection] + [sel_NTRT20]
        histogramName = "NTRT20ZeroFractionVsP" + file_description + "Denomenator"
        hist_filler.BookHistograms(histogramName,\
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = bins,\
                                                  xlabel ="Track P [GeV]",\
                                                  ylabel = "Number of tracks",\
                                                  )

        #do the eta selections and count the number of tracks with an energy deposity less than or equal to 0.0.
        from selections.selections import sel_ELessEqual0
        selections = [sel_ELessEqual0] + [etaSelection] + [sel_NTRT20]
        histogramName = "NTRT20ZeroFractionVsP" + file_description + "Numerator"
        hist_filler.BookHistograms(histogramName,\
                                                       calc_trkP,\
                                                       list_selections = selections,\
                                                       bins = bins,\
                                                       xlabel ="Track P [GeV]",\
                                                       ylabel = "N(E<=0)/N",\
                                                       )

    ################################################################################
    from selections.selections import sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum, sel_PGreater2, sel_PGreater2_5, sel_PGreater3
    MIP_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum]
    selections = [] + MIP_selection
    histogramName = "MIPSelection_HadBetween30And90OfMomentum_EOP"
    hist_filler.BookHistograms(histogramName,
                                     calc_EOP,
                                     list_selections = selections,
                                     bins = 50,
                                     range_low = -1,
                                     range_high = 5,
                                     xlabel ="E/p",
                                     ylabel = "Number of Tracks",
                                     )


    ################################################################################
    selections = [sel_ECALEta0_6] + MIP_selection
    histogramName = "MIPSelection_HadBetween30And90OfMomentum_ECALEta00_06_EOP"
    hist_filler.BookHistograms(histogramName,\
                                               calc_EOP,
                                               list_selections = selections,
                                               bins = 50,
                                               range_low = -1,
                                               range_high = 5,
                                               xlabel ="E/p",
                                               ylabel = "Number of Tracks")


    ################################################################################
    ##Create a set of p and eta bins for the measurement ##########################

    from calculation.calculation import calculation
    from selections.selections import EtaBin, PBin
    from variables.variables import calc_EOPBkg, calc_EnergyAnulus
    from selections.selections import sel_NTRT20, sel_NonZeroEnergy, sel_HardScatter, sel_AnyPion

    ##select inclusive distributions
    ##Create a set of binned EOP response histograms 
    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    base_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "MIPSelectionBetween30and90OfMomentum"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    from selections.selections import sel_EHadFracAbove70, sel_NTRT20, sel_Lar1_1GeV
    base_selection = [sel_EHadFracAbove70, sel_NTRT20, sel_Lar1_1GeV]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
   
    description = "MIPSelectionHadFracAbove70"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    base_selection = [sel_NTRT20, sel_NonZeroEnergy]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "20TRTHitsNonZeroEnergy"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    base_selection = [sel_NonZeroEnergy]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "NonZeroEnergy"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    base_selection = []
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "Inclusive"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    # Crete histograms for the hard scatter histograms
    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    base_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum, sel_HardScatter]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "MIPSelectionBetween30and90OfMomentumHardScatter"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    from selections.selections import sel_EHadFracAbove70, sel_NTRT20, sel_Lar1_1GeV
    base_selection = [sel_EHadFracAbove70, sel_NTRT20, sel_Lar1_1GeV, sel_HardScatter]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
   
    description = "MIPSelectionHadFracAbove70HardScatter"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    base_selection = [sel_NTRT20, sel_NonZeroEnergy, sel_HardScatter]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "20TRTHitsNonZeroEnergyHardScatter"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    base_selection = [sel_NonZeroEnergy, sel_HardScatter]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "NonZeroEnergyHardScatter"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    base_selection = [sel_HardScatter]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "InclusiveHardScatter"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    #create plots for only pions
    from selections.selections import sel_EHadFracAbove70, sel_NTRT20, sel_Lar1_1GeV
    base_selection = [sel_EHadFracAbove70, sel_NTRT20, sel_Lar1_1GeV, sel_HardScatter, sel_AnyPion]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
   
    description = "MIPSelectionHadFracAbove70HardScatterOnlyPion"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    base_selection = [sel_NTRT20, sel_NonZeroEnergy, sel_HardScatter, sel_AnyPion]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "20TRTHitsNonZeroEnergyHardScatterOnlyPion"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    base_selection = [sel_NonZeroEnergy, sel_HardScatter, sel_AnyPion]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "NonZeroEnergyHardScatterOnlyPion"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]
    base_selection = [sel_HardScatter, sel_AnyPion]
    p_bins_for_eta_range = []
    for eta_range in eta_ranges:
        p_bins_min = getP(0.5, (eta_range[0] + eta_range[1]) / 2.0)
        p_bins = getLogBins(p_bins_min, 15.05, 20)
        p_bins_for_eta_range.append(p_bins)
    description = "InclusiveHardScatterOnlyPion"
    PutBinningVectorsInFile(outFile, eta_ranges, p_bins_for_eta_range, description)
    #CreateEOPBinnedHistograms(hist_filler, base_selection, eta_ranges, p_bins_for_eta_range, description) 

    histograms = hist_filler.DumpHistograms()
    for histogram_name in histograms:
        WriteToFile(histograms[histogram_name], outFile)

    print("THEJOBFINISHED!")
