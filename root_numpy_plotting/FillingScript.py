#!/usr/bin/env python
# coding: utf-8
from PlottingTools.Plotter import Plotter, DrawDataVsMC, DivideHistograms,Draw2DHistogramOnCanvas, getP, getBins, getLogBins, GetBinsFromHistogram
import ROOT
#from variables.variables import calc_weight
#from inputs.samples import INPUT
import os
from math import pi
import pickle
import numpy as np

def WriteToFile(histogram_dictionary, outFile):
    outFile.cd()
    for key in histogram_dictionary:
        if not outFile.cd(key):
            outFile.mkdir(key)
        outFile.cd(key)
        print("Writing histogram " + histogram_dictionary[key].GetName())
        histogram_dictionary[key].Write()

#This is a script that fills the histograms for
def FillingScript(plotter, outputRootFileName):
    #import thje variables that we want to plot
    from variables.variables import calc_trkNearestNeighbourEM2, calc_trkP, calc_EOP, calc_trkPt, calc_trkAverageMu, calc_trkEtaID, calc_trkEtaECAL, calc_trkNPV2, calc_trkCount, calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike
    #import the selections that we want to plot
    from selections.selections import sel_NTRT20, sel_NTRT25, sel_NTRT30, sel_ECALEta0_6, sel_PGreater1 ,sel_PGreater1_5, sel_PGreater2, sel_PGreater2_5, sel_Z0SinThetaLess1_5, sel_d0Less1_5, sel_Event, sel_hasHADExtrapolation
    #impot the ID selections that we want to plot
    from selections.selections import sel_IDEta00_06,sel_IDEta06_11,sel_IDEta11_14,sel_IDEta14_15, sel_IDEta15_18, sel_IDEta18_23

    outFile = ROOT.TFile(outputRootFileName, "RECREATE")

    histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
    xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
    variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]
    selections = []

    for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
        plotter.BookHistograms(histogram_name,\
                               variable,\
                               list_selections = selections,\
                               bins = 10,\
                               range_low = -0.5,\
                               range_high = 9.5,\
                               xlabel=xlabel,\
                               ylabel="Number of Tracks")

    histogram_name = "trkTRTHits"
    from variables.variables import calc_nTRT
    plotter.BookHistograms(histogram_name,\
                           calc_nTRT,\
                           list_selections = selections,\
                           range_low = -0.5,\
                           range_high = 59.5,\
                           bins = 60,\
                           xlabel = "Number of TRT Hits",\
                           ylabel = "Number of Tracks")

    histogram_name = "trkEMDR100"
    from variables.variables import calc_EnergyEMDR100
    plotter.BookHistograms(histogram_name,\
                           calc_EnergyEMDR100,\
                           list_selections = selections,\
                           range_low = -2.0,\
                           range_high = + 10.0,\
                           bins = 48,\
                           xlabel = "E_{EM}^{#DeltaR<0.1}[GeV]",\
                           ylabel = "Number of Tracks")

    histogram_name = "MomentumHadFrac"
    from variables.variables import calc_MomentumHadFrac
    plotter.BookHistograms(histogram_name,\
                           calc_MomentumHadFrac,\
                           list_selections = selections,\
                           range_low = -1.0,\
                           range_high = + 5.0,\
                           bins = 48,\
                           xlabel = "E^{HAD}/P",\
                           ylabel = "Number of Tracks")

    histogram_name = "HadFrac"
    from variables.variables import calc_HadFrac
    plotter.BookHistograms(histogram_name,\
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
    trkCountHist = plotter.BookHistograms(histogram_name,\
                                                         calc_trkCount,\
                                                         list_selections = selections,\
                                                         bins = 1,\
                                                         range_low = -0.5,\
                                                         range_high = +0.5,\
                                                         xlabel ='Always 0',\
                                                         ylabel = 'Number of Tracks')

    ################################################################################
    #plot the first set of variables that we're interested in
    #plot the trk avaerage mu histogram
    histogram_name = "trkAverageMu"
    selections = []
    trkAverageMuHist = plotter.BookHistograms(histogram_name,\
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
    trkNPV2Hist = plotter.BookHistograms(histogram_name,\
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
    eventNPV2Hist = plotter.BookHistograms(histogram_name,\
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
    eventAverageMuHist = plotter.BookHistograms(histogram_name,
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
    p_bins_reference = p_bins # a set of bins to use later to measure the efficiency of certain selections
    histogram_name = "trkPtHist"
    trkPtHistZoom = plotter.BookHistograms(histogram_name,\
                                       calc_trkPt,\
                                       list_selections = [],\
                                       bins = p_bins,\
                                       xlabel ="Track P_{T} [GeV]",\
                                       ylabel = "Number of Tracks")

    trkPtHistZoom = plotter.BookHistograms(histogram_name + "HasExtrapolation",\
                                       calc_trkPt,\
                                       list_selections = [sel_hasHADExtrapolation],\
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
    plotter.Book2DHistograms(histogram_name,\
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
    plotter.Book2DHistograms(histogram_name,\
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
        trkMultiplicity_Eta = plotter.BookHistograms("TrkPtHist"+file_description,\
                                                  calc_trkPt,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="Track P_{T} [GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  )

        trkMultiplicity_Eta = plotter.BookHistograms("TrkPtHist"+file_description + "HasExtrapolation",\
                                                  calc_trkPt,\
                                                  list_selections = selections + [sel_hasHADExtrapolation],\
                                                  bins = p_bins,\
                                                  xlabel ="Track P_{T} [GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  )

    ################################################################################
    histogramName = "TrackEtaID"
    trkEtaIDHist = plotter.BookHistograms(histogramName,\
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
    TwoDtrkPvstrkEta = plotter.Book2DHistograms(histogramName,\
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
    TwoDtrkPtvstrkEta = plotter.Book2DHistograms(histogramName,\
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
    TwoDtrkPtvstrkEta = plotter.Book2DHistograms(histogramName,\
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
    trkEtaECALHist = plotter.BookHistograms(histogramName,\
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
    TwoDtrkPvstrkDPhi = plotter.Book2DHistograms(histogramName,\
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
    TwoDtrkEtavstrkDEta = plotter.Book2DHistograms(histogramName,\
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
    CentalTwoDtrkPvstrkDPhi = plotter.Book2DHistograms(histogramName,\
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
    trkNearestDRHist = plotter.BookHistograms(histogramName,
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
    trkEOPHist = plotter.BookHistograms(histogramName,\
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
    trkEOPHist = plotter.BookHistograms(histogramName,\
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
    trkEOPHistPGreater1 = plotter.BookHistograms(histogramName,\
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
    trkEOPHistPGreater1_5 = plotter.BookHistograms(histogramName,
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
    trkEOPHistPGreater2 = plotter.BookHistograms(histogramName,
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
    trkEOPHistFig2a = plotter.BookHistograms(histogramName,
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
    trkEOPHistFig2b = plotter.BookHistograms(histogramName,\
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
    trkEOPHistFig2c = plotter.BookHistograms(histogramName,\
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
    trkEOPHistFig2c = plotter.BookHistograms(histogramName,
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
    trkEOPHistFig2a = plotter.BookHistograms(histogramName,\
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
    trkEOPHistFig2b = plotter.BookHistograms(histogramName,\
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
    trkEOPHistFig2c = plotter.BookHistograms(histogramName,
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
    trkEOPHistFig2c = plotter.BookHistograms(histogramName,
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
    trkMultiplicity = plotter.BookHistograms(histogramName,\
                                          calc_trkP,\
                                          list_selections = selections,\
                                          bins = bins,\
                                          xlabel ="Track P [GeV]",\
                                          ylabel = "Number of Tracks",\
                                          )

    from selections.selections import sel_ELessEqual0
    histogramName = "InclusiveZeroFractionVsPNumerator"
    selections = [sel_ELessEqual0]
    trkMultiplicity_ELessZero = plotter.BookHistograms(histogramName,\
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
    trkMultiplicity_Eta = plotter.BookHistograms(histogramName,\
                                              calc_trkEtaID,\
                                              list_selections = selections,\
                                              bins = bins,\
                                              xlabel ="Track |#eta|",\
                                              ylabel = "Number of Tracks",\
                                              )
    from selections.selections import sel_ELessEqual0
    histogramName = "InclusiveZeroFractionVsEtaNumerator"
    selections = [sel_ELessEqual0]
    trkMultiplicity_Eta_Zero = plotter.BookHistograms(histogramName,\
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
    trkMultiplicity_AbsEta = plotter.BookHistograms(histogramName,\
                                              calc_trkEta_ABS,\
                                              list_selections = selections,\
                                              bins = bins,\
                                              xlabel ="Track |#eta|",\
                                              ylabel = "Number of Tracks",\
                                              )
    from selections.selections import sel_ELessEqual0
    histogramName = "InclusiveZeroFractionVsAbsEtaNumerator"
    selections = [sel_ELessEqual0]
    trkMultiplicity_AbsEta_Zero = plotter.BookHistograms(histogramName,\
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
        trkMultiplicity_Eta = plotter.BookHistograms(histogramName,\
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
        trkMultiplicity_Eta_Zero = plotter.BookHistograms(histogramName,\
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
        trkMultiplicity_Eta = plotter.BookHistograms(histogramName,\
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
        trkMultiplicity_Eta_Zero = plotter.BookHistograms(histogramName,\
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
    trkEOPHist = plotter.BookHistograms(histogramName,
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
    trkEOPHistPGreater1 = plotter.BookHistograms(histogramName,\
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

    #create a set of strings that could describe the eta or momentum selections
    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]

    #go and get the average E/P for MIP particles in each of the eta bins.
    for eta_range in eta_ranges:
        #get the function that selectts tracks in that bin
        EtaBinFunction = lambda x, y = eta_range[0], z=eta_range[1]: EtaBin(x, y,z)
        EtaBinFunction.__name__ = "EtaRangeSelection"+str(eta_range[0]) + "_" + str(eta_range[1])
        eta_binSelection = calculation(EtaBinFunction, ["trk_etaID"])

        #calculate the lowest momentum track that can end up in that bin
        #create nbins where there are at least 10,000 entries per bin
        eta = eta_range[1]
        p_bins_max = 15.05
        p_bins_min = getP(0.5, eta)
        nBins = 20
        p_bins = getLogBins(p_bins_min, p_bins_max, nBins)
        p_bins_fine = getLogBins(p_bins_min, p_bins_max, 10000)

        eop_bins_min = -1
        eop_bins_max = 5
        nBins = 120
        eop_bins = getBins(eop_bins_min, eop_bins_max, nBins)

        print("+" * 50)
        print("for Eta " + str(eta))
        print(eop_bins)
        print(p_bins)

        MIP_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum]
        selections = MIP_selection + [eta_binSelection]

        histogram_name = "trkPtHist"
        histogram_name = histogram_name + "_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkPtHistZoom = plotter.BookHistograms(histogram_name,\
                                           calc_trkPt,\
                                           list_selections = [eta_binSelection],\
                                           bins = p_bins_reference,\
                                           xlabel ="Track P_{T} [GeV]",\
                                           ylabel = "Number of Tracks")

   #    ################################################################################
        #prepare the momentum bins
        histogram_name = "trkPtHist"
        histogram_name = histogram_name + "_MIPSelection_HADBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkPtHistZoom = plotter.BookHistograms(histogram_name,\
                                           calc_trkPt,\
                                           list_selections = selections,\
                                           bins = p_bins_reference,\
                                           xlabel ="Track P_{T} [GeV]",\
                                           ylabel = "Number of Tracks")

        histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
        xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
        variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]

        for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
            histogram_name = histogram_name + "_MIPSelection_HADBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
            plotter.BookHistograms(histogram_name,\
                                   variable,\
                                   list_selections = selections,\
                                   bins = 10,\
                                   range_low = -0.5,\
                                   range_high = 9.5,\
                                   xlabel=xlabel,\
                                   ylabel="Number of Tracks")

        histogramName = "TrkMultiplicityVsP_Inclusive_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = [eta_binSelection],\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  )

        histogramName = "UnweightedTrkMultiplicityVsP_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  useWeights=False
                                                  )


        histogramName = "TrkMultiplicityVsP_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  )

        #get the binning histograms for the histograms to be used for fitting
        binningHistogram = plotter.binningHistograms[histogramName]
        FourThousandTracks_pbins = GetBinsFromHistogram(binningHistogram, 4000.0)[0]

        print FourThousandTracks_pbins

        histogramName = "EOPProfileVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.BookTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )

        histogramName = "EOPProfileVsMomentum_FourThousandTracks_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.BookTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = FourThousandTracks_pbins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )

        histogramName = "2DHist_EOPVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.Book2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p",\
                                                  )

        histogramName = "2DHist_EOPVsMomentum_FourThousandTracks_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.Book2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = FourThousandTracks_pbins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p",\
                                                  )

        histogramName = "EnergyAnulusProfileVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )

        histogramName = "EnergyAnulusProfileVsMomentum_FourThousandTracks_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = FourThousandTracks_pbins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )

        histogramName = "2DHist_EnergyAnulusVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )

        histogramName = "2DHist_EnergyAnulusVsMomentum_FourThousandTracks_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = FourThousandTracks_pbins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )

        histogramName = "EnergyBkgProfileVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )


        histogramName = "EnergyBkgProfileVsMomentum_FourThousandTracks_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = FourThousandTracks_pbins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogramName = "2DHist_EnergyBkgVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p BKG",\
                                                  )

        histogramName = "2DHist_EnergyBkgVsMomentum_FourThousandTracks_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = FourThousandTracks_pbins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p BKG",\
                                                  )

        #go and get the E/p distribution in each of the E/p bins
        p_ranges = [ (p_bins[i], p_bins[i+1])  for i in range(0, len(p_bins)-1) ]
        for p_range in p_ranges:
            print("The prange is " + str(p_range))
            PBinFunction = lambda x, y=p_range[0], z=p_range[1]: PBin(x, y,z)
            PBinFunction.__name__ = "SelMomentumRange"+str(p_range[0]) + "_" + str(p_range[1])
            sel_PBin = calculation(PBinFunction, ["trk_p"])
            selections = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum, eta_binSelection, sel_PBin]
            histogramName = "EOPDistribution_MIPSelection_HadBetween30And90OfMomentum_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOP,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )

            histogramName = "EOPBkgDistribution_MIPSelection_HadBetween30And90OfMomentum_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPBkgDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOPBkg,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p Bkg",\
                                                      )


        #go and get the E/p distribution in each of the E/p bins
        p_ranges = [ (FourThousandTracks_pbins[i], FourThousandTracks_pbins[i+1])  for i in range(0, len(FourThousandTracks_pbins)-1) ]
        for p_range in p_ranges:
            print("The prange is " + str(p_range))
            PBinFunction = lambda x,y=p_range[0],z=p_range[1]: PBin(x, y,z)
            PBinFunction.__name__ = "SelMomentumRange"+str(p_range[0]) + "_" + str(p_range[1])
            sel_PBin = calculation(PBinFunction, ["trk_p"])
            selections = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum, eta_binSelection, sel_PBin]
            histogramName = "EOPDistribution_FourThousandTracks_MIPSelection_HadBetween30And90OfMomentum_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOP,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )

            histogramName = "EOPBkgDistribution_FourThousandTracks_MIPSelection_HadBetween30And90OfMomentum_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPBkgDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOPBkg,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p Bkg",\
                                                      )

        from selections.selections import sel_EHadFracAbove70
        MIP_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadFracAbove70]
        selections = MIP_selection + [eta_binSelection]

    #    ################################################################################
        #prepare the momentum bins
        histogram_name = "trkPtHist"
        histogram_name = histogram_name + "_MIPSelection_HADFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkPtHistZoom = plotter.BookHistograms(histogram_name,\
                                           calc_trkPt,\
                                           list_selections = selections,\
                                           bins = p_bins_reference,\
                                           xlabel ="Track P_{T} [GeV]",\
                                           ylabel = "Number of Tracks")

        histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
        xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
        variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]

        for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
            histogram_name = histogram_name + "_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
            plotter.BookHistograms(histogram_name,\
                                   variable,\
                                   list_selections = selections,\
                                   bins = 10,\
                                   range_low = -0.5,\
                                   range_high = 9.5,\
                                   xlabel=xlabel,\
                                   ylabel="Number of Tracks")

        histogramName = "UnweightedTrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  useWeights=False
                                                  )

        histogramName = "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P [GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  )

        #get the binning histograms for the histograms to be used for fitting
        binningHistogram = plotter.binningHistograms[histogramName]
        FourThousandTracks_pbins = GetBinsFromHistogram(binningHistogram, 4000.0)[0]

        histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.BookTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )

        histogramName = "EOPProfileVsMomentum_FourThousandTracks_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.BookTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = FourThousandTracks_pbins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )

        histogramName = "2DHist_EOPVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.Book2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P [GeV]",\
                                                  ylabel = "E/p",\
                                                  )

        histogramName = "2DHist_EOPVsMomentum_FourThousandTracks_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.Book2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = FourThousandTracks_pbins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P [GeV]",\
                                                  ylabel = "E/p",\
                                                  )

        histogramName = "EnergyAnulusProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )

        histogramName = "EnergyAnulusProfileVsMomentum_FourThousandTracks_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = FourThousandTracks_pbins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )


        histogramName = "2DHist_EnergyAnulusVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )

        histogramName = "2DHist_EnergyAnulusVsMomentum_FourThousandTracks_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = FourThousandTracks_pbins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )

        histogramName = "EnergyBkgProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogramName = "EnergyBkgProfileVsMomentum_FourThousandTracks_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = FourThousandTracks_pbins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogramName = "2DHist_EnergyBkgVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p BKG",\
                                                  )

        histogramName = "2DHist_EnergyBkgVsMomentum_FourThousandTracks_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = FourThousandTracks_pbins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p BKG",\
                                                  )

        #go and get the E/p distribution in each of the E/p bins
        p_ranges = [ (p_bins[i], p_bins[i+1] ) for i in range(0, len(p_bins)-1) ]
        for p_range in p_ranges:
            PBinFunction = lambda x,y=p_range[0],z=p_range[1]: PBin(x,y,z)
            PBinFunction.__name__ = "SelMomentumRange"+str(p_range[0]) + "_" + str(p_range[1])
            sel_PBin = calculation(PBinFunction, ["trk_p"])
            selections = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadFracAbove70, eta_binSelection, sel_PBin]
            histogramName = "EOPDistribution_MIPSelection_HadFracAbove70_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOP,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )

            histogramName = "EOPBkgDistribution_MIPSelection_HadFracAbove70_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            AverageEOPBkg  =  plotter.BookHistograms(histogramName,
                                                      calc_EOPBkg,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p Bkg",\
                                                      )
            histogram_name = "trkTRTHits"
            histogram_name = histogram_name + "_MIPSelection_HadFracAbove70_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_nTRT
            plotter.BookHistograms(histogram_name,\
                                   calc_nTRT,\
                                   list_selections = selections,\
                                   range_low = -0.5,\
                                   range_high = 59.5,\
                                   bins = 60,\
                                   xlabel = "Number of TRT Hits",\
                                   ylabel = "Number of Tracks")

            histogram_name = "trkEMDR100"
            histogram_name = histogram_name + "_MIPSelection_HadFracAbove70_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_EnergyEMDR100
            plotter.BookHistograms(histogram_name,\
                                   calc_EnergyEMDR100,\
                                   list_selections = selections,\
                                   range_low = -2.0,\
                                   range_high = + 10.0,\
                                   bins = 48,\
                                   xlabel = "E_{EM}^{#DeltaR<0.1}[GeV]",\
                                   ylabel = "Number of Tracks")

            histogram_name = "MomentumHadFrac"
            histogram_name = histogram_name + "_MIPSelection_HadFracAbove70_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_MomentumHadFrac
            plotter.BookHistograms(histogram_name,\
                                   calc_MomentumHadFrac,\
                                   list_selections = selections,\
                                   range_low = -1.0,\
                                   range_high = + 5.0,\
                                   bins = 48,\
                                   xlabel = "E^{HAD}/P",\
                                   ylabel = "Number of Tracks")

            histogram_name = "HadFrac"
            histogram_name = histogram_name + "_MIPSelection_HadFracAbove70_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_HadFrac
            plotter.BookHistograms(histogram_name,\
                                   calc_HadFrac,\
                                   list_selections = selections,\
                                   range_low = -1.0,\
                                   range_high = + 2.0,\
                                   bins = 48,\
                                   xlabel = "E^{HAD}/E^{Total}",\
                                   ylabel = "Number of Tracks")

            histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
            xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
            variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]

            for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
                histogram_name = histogram_name + "_MIPSelection_HadFracAbove70_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
                plotter.BookHistograms(histogram_name,\
                                       variable,\
                                       list_selections = selections,\
                                       bins = 10,\
                                       range_low = -0.5,\
                                       range_high = 9.5,\
                                       xlabel=xlabel,\
                                       ylabel="Number of Tracks")

        #go and get the E/p distribution in each of the E/p bins
        p_ranges = [ (FourThousandTracks_pbins[i], FourThousandTracks_pbins[i+1] ) for i in range(0, len(FourThousandTracks_pbins)-1) ]
        for p_range in p_ranges:
            PBinFunction = lambda x,y=p_range[0],z=p_range[1]: PBin(x,y,z)
            PBinFunction.__name__ = "SelMomentumRange"+str(p_range[0]) + "_" + str(p_range[1])
            sel_PBin = calculation(PBinFunction, ["trk_p"])
            selections = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadFracAbove70, eta_binSelection, sel_PBin]
            histogramName = "EOPDistribution_FourThousandTracks_MIPSelection_HadFracAbove70_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOP,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )

            histogramName = "EOPBkgDistribution_FourThousandTracks_MIPSelection_HadFracAbove70_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            AverageEOPBkg  =  plotter.BookHistograms(histogramName,
                                                      calc_EOPBkg,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p Bkg",\
                                                      )


    ################################################################################
    ##Create a set of p and eta bins for the measurement ##########################
    from calculation.calculation import calculation
    from selections.selections import EtaBin, PBin
    from variables.variables import calc_EOPBkg, calc_EnergyAnulus

    #create a set of strings that could describe the eta or momentum selections
    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]

    #go and get the average E/P for MIP particles in each of the eta bins.
    for eta_range in eta_ranges:
        #get the function that selectts tracks in that bin
        EtaBinFunction = lambda x,y=eta_range[0],z=eta_range[1]: EtaBin(x,y,z)
        EtaBinFunction.__name__ = "EtaRangeSelection"+str(eta_range[0]) + "_" + str(eta_range[1])
        eta_binSelection = calculation(EtaBinFunction, ["trk_etaID"])

        #calculate the lowest momentum track that can end up in that bin
        #create nbins where there are at least 10,000 entries per bin
        eta = eta_range[1]
        p_bins_max = 15.05
        p_bins_min = getP(0.5, eta)
        nBins = 20
        p_bins = getLogBins(p_bins_min, p_bins_max, nBins)
        p_bins_fine = getLogBins(p_bins_min, p_bins_max, 10000)

        eop_bins_min = -1
        eop_bins_max = 5
        nBins = 120
        eop_bins = getBins(eop_bins_min, eop_bins_max, nBins)

        print("+" * 50)
        print("for Eta " + str(eta))
        print(eop_bins)
        print(p_bins)

        selections = [eta_binSelection]

        histogramName = "UnweightedTrkMultiplicityVsP_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  useWeights=False
                                                  )

        histogramName = "TrkMultiplicityVsP_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  )


        histogramName = "EOPProfileVsMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.BookTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )

        histogramName = "2DHist_EOPVsMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.Book2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p",\
                                                  )

        histogramName = "EnergyAnulusProfileVsMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )

        histogramName = "2DHist_EnergyAnulusVsMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )

        histogramName = "EnergyBkgProfileVsMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogramName = "2DHist_EnergyBkgVsMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p BKG",\
                                                  )

        #go and get the E/p distribution in each of the E/p bins
        p_ranges = [ (p_bins[i], p_bins[i+1])  for i in range(0, len(p_bins)-1) ]
        for p_range in p_ranges:
            print("The prange is " + str(p_range))
            PBinFunction = lambda x,y=p_range[0],z=p_range[1]: PBin(x, y,z)
            PBinFunction.__name__ = "SelMomentumRange"+str(p_range[0]) + "_" + str(p_range[1])
            sel_PBin = calculation(PBinFunction, ["trk_p"])
            selections = [eta_binSelection, sel_PBin]
            histogramName = "EOPDistribution_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOP,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )

            histogramName = "EOPBkgDistribution_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            AverageEOPBkg  =  plotter.BookHistograms(histogramName,
                                                      calc_EOPBkg,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p Bkg",\
                                                      )

    #Select only those tracks with a cluster
    ################################################################################
    ##Create a set of p and eta bins for the measurement ##########################
    from calculation.calculation import calculation
    from selections.selections import EtaBin, PBin, sel_NonZeroEnergy
    from variables.variables import calc_EOPBkg, calc_EnergyAnulus

    #create a set of strings that could describe the eta or momentum selections
    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]

    #go and get the average E/P for MIP particles in each of the eta bins.
    for eta_range in eta_ranges:
        #get the function that selectts tracks in that bin
        EtaBinFunction = lambda x,y=eta_range[0],z=eta_range[1]: EtaBin(x,y,z)
        EtaBinFunction.__name__ = "EtaRangeSelection"+str(eta_range[0]) + "_" + str(eta_range[1])
        eta_binSelection = calculation(EtaBinFunction, ["trk_etaID"])

        #calculate the lowest momentum track that can end up in that bin
        #create nbins where there are at least 10,000 entries per bin
        eta = eta_range[1]
        p_bins_max = 15.05
        p_bins_min = getP(0.5, eta)
        nBins = 20
        p_bins = getLogBins(p_bins_min, p_bins_max, nBins)
        p_bins_fine = getLogBins(p_bins_min, p_bins_max, 10000)

        eop_bins_min = -1
        eop_bins_max = 5
        nBins = 120
        eop_bins = getBins(eop_bins_min, eop_bins_max, nBins)

        print("+" * 50)
        print("for Eta " + str(eta))
        print(eop_bins)
        selections = [sel_NonZeroEnergy, eta_binSelection]

    #    ################################################################################
        #prepare the momentum bins
        histogram_name = "trkPtHist"
        histogram_name = histogram_name + "_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkPtHistZoom = plotter.BookHistograms(histogram_name,\
                                           calc_trkPt,\
                                           list_selections = selections,\
                                           bins = p_bins_reference,\
                                           xlabel ="Track P_{T} [GeV]",\
                                           ylabel = "Number of Tracks")

        histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
        xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
        variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]

        for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
            histogram_name = histogram_name + "_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
            plotter.BookHistograms(histogram_name,\
                                   variable,\
                                   list_selections = selections,\
                                   bins = 10,\
                                   range_low = -0.5,\
                                   range_high = 9.5,\
                                   xlabel=xlabel,\
                                   ylabel="Number of Tracks")

        histogram_name = "trkTRTHits"
        histogram_name = histogram_name + "_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        from variables.variables import calc_nTRT
        plotter.BookHistograms(histogram_name,\
                               calc_nTRT,\
                               list_selections = selections,\
                               range_low = -0.5,\
                               range_high = 59.5,\
                               bins = 60,\
                               xlabel = "Number of TRT Hits",\
                               ylabel = "Number of Tracks")

        histogram_name = "trkEMDR100"
        histogram_name = histogram_name + "_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        from variables.variables import calc_EnergyEMDR100
        plotter.BookHistograms(histogram_name,\
                               calc_EnergyEMDR100,\
                               list_selections = selections,\
                               range_low = -2.0,\
                               range_high = + 10.0,\
                               bins = 48,\
                               xlabel = "E_{EM}^{#DeltaR<0.1}[GeV]",\
                               ylabel = "Number of Tracks")

        histogram_name = "MomentumHadFrac"
        histogram_name = histogram_name + "_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        from variables.variables import calc_MomentumHadFrac
        plotter.BookHistograms(histogram_name,\
                               calc_MomentumHadFrac,\
                               list_selections = selections,\
                               range_low = -1.0,\
                               range_high = + 5.0,\
                               bins = 48,\
                               xlabel = "E^{HAD}/P",\
                               ylabel = "Number of Tracks")

        histogram_name = "HadFrac"
        histogram_name = histogram_name + "_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        from variables.variables import calc_HadFrac
        plotter.BookHistograms(histogram_name,\
                               calc_HadFrac,\
                               list_selections = selections,\
                               range_low = -1.0,\
                               range_high = + 2.0,\
                               bins = 48,\
                               xlabel = "E^{HAD}/E^{Total}",\
                               ylabel = "Number of Tracks")

        histogramName = "UnweightedTrkMultiplicityVsP_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  useWeights=False\
                                                  )

        histogramName = "TrkMultiplicityVsP_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  )

        #get the custom bins for the EOP Measurement
        hist = plotter.binningHistograms[histogramName]
        FourThousandTracks_pbins = GetBinsFromHistogram(hist, 4000.0)[0]

        histogramName = "EOPProfileVsMomentum_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.BookTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )

        histogramName = "EOPProfileVsMomentum_FourThousandTracks_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.BookTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = FourThousandTracks_pbins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )

        histogramName = "2DHist_EOPVsMomentum_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.Book2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p",\
                                                  )

        histogramName = "2DHist_FourThousandTracks_EOPVsMomentum_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.Book2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = FourThousandTracks_pbins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p",\
                                                  )

        histogramName = "EnergyAnulusProfileVsMomentum_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )

        histogramName = "EnergyAnulusProfileVsMomentum_FourThousandTracks_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = FourThousandTracks_pbins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )

        histogramName = "2DHist_EnergyAnulusVsMomentum_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )

        histogramName = "2DHist_EnergyAnulusVsMomentum_FourThousandTracks_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = FourThousandTracks_pbins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )

        histogramName = "EnergyBkgProfileVsMomentum_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogramName = "EnergyBkgProfileVsMomentum_FourThousandTracks_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = FourThousandTracks_pbins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )

        histogramName = "2DHist_EnergyBkgVsMomentum_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p BKG",\
                                                  )

        histogramName = "2DHist_EnergyBkgVsMomentum_FourThousandTracks_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = FourThousandTracks_pbins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p BKG",\
                                                  )

        #go and get the E/p distribution in each of the E/p bins
        p_ranges = [ (p_bins[i], p_bins[i+1])  for i in range(0, len(p_bins)-1) ]
        for p_range in p_ranges:
            print("The prange is " + str(p_range))
            PBinFunction = lambda x, y=p_range[0], z=p_range[1]: PBin(x, y, z)
            PBinFunction.__name__ = "SelMomentumRange"+str(p_range[0]) + "_" + str(p_range[1])
            sel_PBin = calculation(PBinFunction, ["trk_p"])
            selections = [sel_NonZeroEnergy,eta_binSelection, sel_PBin]
            histogramName = "EOPDistribution_NonZeroE_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOP,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )

            histogramName = "EOPBkgDistribution"
            histogramName = histogramName + "_NonZeroE_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPBkgDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOPBkg,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p Bkg",\
                                                      )
            histogram_name = "trkTRTHits"
            histogram_name = histogram_name + "_NonZeroE_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_nTRT
            plotter.BookHistograms(histogram_name,\
                                   calc_nTRT,\
                                   list_selections = selections,\
                                   range_low = -0.5,\
                                   range_high = 59.5,\
                                   bins = 60,\
                                   xlabel = "Number of TRT Hits",\
                                   ylabel = "Number of Tracks")

            histogram_name = "trkEMDR100"
            histogram_name = histogram_name + "_NonZeroE_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_EnergyEMDR100
            plotter.BookHistograms(histogram_name,\
                                   calc_EnergyEMDR100,\
                                   list_selections = selections,\
                                   range_low = -2.0,\
                                   range_high = + 10.0,\
                                   bins = 48,\
                                   xlabel = "E_{EM}^{#DeltaR<0.1}[GeV]",\
                                   ylabel = "Number of Tracks")

            histogram_name = "MomentumHadFrac"
            histogram_name = histogram_name + "_NonZeroE_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_MomentumHadFrac
            plotter.BookHistograms(histogram_name,\
                                   calc_MomentumHadFrac,\
                                   list_selections = selections,\
                                   range_low = -1.0,\
                                   range_high = + 5.0,\
                                   bins = 48,\
                                   xlabel = "E^{HAD}/P",\
                                   ylabel = "Number of Tracks")

            histogram_name = "HadFrac"
            histogram_name = histogram_name + "_NonZeroE_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_HadFrac
            plotter.BookHistograms(histogram_name,\
                                   calc_HadFrac,\
                                   list_selections = selections,\
                                   range_low = -1.0,\
                                   range_high = + 2.0,\
                                   bins = 48,\
                                   xlabel = "E^{HAD}/E^{Total}",\
                                   ylabel = "Number of Tracks")

            histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
            xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
            variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]

            for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
                histogram_name = histogram_name + "_NonZeroE_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
                plotter.BookHistograms(histogram_name,\
                                       variable,\
                                       list_selections = selections,\
                                       bins = 10,\
                                       range_low = -0.5,\
                                       range_high = 9.5,\
                                       xlabel=xlabel,\
                                       ylabel="Number of Tracks")


    from selections.selections import sel_NTRT20
    for eta_range in eta_ranges:
        #get the function that selectts tracks in that bin
        EtaBinFunction = lambda x,y=eta_range[0],z=eta_range[1]: EtaBin(x,y,z)
        EtaBinFunction.__name__ = "EtaRangeSelection"+str(eta_range[0]) + "_" + str(eta_range[1])
        eta_binSelection = calculation(EtaBinFunction, ["trk_etaID"])

        #calculate the lowest momentum track that can end up in that bin
        #create nbins where there are at least 10,000 entries per bin
        eta = eta_range[1]
        p_bins_max = 15.05
        p_bins_min = getP(0.5, eta)
        nBins = 20
        p_bins = getLogBins(p_bins_min, p_bins_max, nBins)
        p_bins_fine = getLogBins(p_bins_min, p_bins_max, 10000)

        eop_bins_min = -1
        eop_bins_max = 5
        nBins = 120
        eop_bins = getBins(eop_bins_min, eop_bins_max, nBins)

        print("+" * 50)
        print("for Eta " + str(eta))
        print(eop_bins)
        selections = [sel_NonZeroEnergy, eta_binSelection, sel_NTRT20]

    #    ################################################################################
        #prepare the momentum bins
        histogram_name = "trkPtHist"
        histogram_name = histogram_name + "_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkPtHistZoom = plotter.BookHistograms(histogram_name,\
                                           calc_trkPt,\
                                           list_selections = selections,\
                                           bins = p_bins_reference,\
                                           xlabel ="Track P_{T} [GeV]",\
                                           ylabel = "Number of Tracks")

        histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
        xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
        variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]

        for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
            histogram_name = histogram_name + "_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
            plotter.BookHistograms(histogram_name,\
                                   variable,\
                                   list_selections = selections,\
                                   bins = 10,\
                                   range_low = -0.5,\
                                   range_high = 9.5,\
                                   xlabel=xlabel,\
                                   ylabel="Number of Tracks")

        histogram_name = "trkTRTHits"
        histogram_name = histogram_name + "_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        from variables.variables import calc_nTRT
        plotter.BookHistograms(histogram_name,\
                               calc_nTRT,\
                               list_selections = selections,\
                               range_low = -0.5,\
                               range_high = 59.5,\
                               bins = 60,\
                               xlabel = "Number of TRT Hits",\
                               ylabel = "Number of Tracks")

        histogram_name = "trkEMDR100"
        histogram_name = histogram_name + "_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        from variables.variables import calc_EnergyEMDR100
        plotter.BookHistograms(histogram_name,\
                               calc_EnergyEMDR100,\
                               list_selections = selections,\
                               range_low = -2.0,\
                               range_high = + 10.0,\
                               bins = 48,\
                               xlabel = "E_{EM}^{#DeltaR<0.1}[GeV]",\
                               ylabel = "Number of Tracks")

        histogram_name = "MomentumHadFrac"
        histogram_name = histogram_name + "_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        from variables.variables import calc_MomentumHadFrac
        plotter.BookHistograms(histogram_name,\
                               calc_MomentumHadFrac,\
                               list_selections = selections,\
                               range_low = -1.0,\
                               range_high = + 5.0,\
                               bins = 48,\
                               xlabel = "E^{HAD}/P",\
                               ylabel = "Number of Tracks")

        histogram_name = "HadFrac"
        histogram_name = histogram_name + "_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        from variables.variables import calc_HadFrac
        plotter.BookHistograms(histogram_name,\
                               calc_HadFrac,\
                               list_selections = selections,\
                               range_low = -1.0,\
                               range_high = +2.0,\
                               bins = 48,\
                               xlabel = "E^{HAD}/E^{Total}",\
                               ylabel = "Number of Tracks")

        histogramName = "UnweightedTrkMultiplicityVsP_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  useWeights=False\
                                                  )

        histogramName = "TrkMultiplicityVsP_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        trkMultiplicity =  plotter.BookHistograms(histogramName,
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = p_bins_fine,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  )

        histogramName = "EOPProfileVsMomentum_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.BookTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )

        histogramName = "2DHist_EOPVsMomentum_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.Book2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p",\
                                                  )

        histogramName = "EnergyAnulusProfileVsMomentum_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )


        histogramName = "2DHist_EnergyAnulusVsMomentum_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )


        histogramName = "EnergyBkgProfileVsMomentum_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.BookTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )


        histogramName = "2DHist_EnergyBkgVsMomentum_NonZeroE_20TRT_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Book2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p BKG",\
                                                  )


        #go and get the E/p distribution in each of the E/p bins
        p_ranges = [ (p_bins[i], p_bins[i+1])  for i in range(0, len(p_bins)-1) ]
        for p_range in p_ranges:
            print("The prange is " + str(p_range))
            PBinFunction = lambda x, y=p_range[0], z=p_range[1]: PBin(x, y, z)
            PBinFunction.__name__ = "SelMomentumRange"+str(p_range[0]) + "_" + str(p_range[1])
            sel_PBin = calculation(PBinFunction, ["trk_p"])
            selections = [sel_NonZeroEnergy,eta_binSelection, sel_PBin, sel_NTRT20]
            histogramName = "EOPDistribution_NonZeroE_20TRT_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOP,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )

            histogramName = "EOPBkgDistribution"
            histogramName = histogramName + "_NonZeroE_20TRT_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            EOPBkgDist  =  plotter.BookHistograms(histogramName,
                                                      calc_EOPBkg,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p Bkg",\
                                                      )
            histogram_name = "trkTRTHits"
            histogram_name = histogram_name + "_NonZeroE_20TRT_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_nTRT
            plotter.BookHistograms(histogram_name,\
                                   calc_nTRT,\
                                   list_selections = selections,\
                                   range_low = -0.5,\
                                   range_high = 59.5,\
                                   bins = 60,\
                                   xlabel = "Number of TRT Hits",\
                                   ylabel = "Number of Tracks")

            histogram_name = "trkEMDR100"
            histogram_name = histogram_name + "_NonZeroE_20TRT_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_EnergyEMDR100
            plotter.BookHistograms(histogram_name,\
                                   calc_EnergyEMDR100,\
                                   list_selections = selections,\
                                   range_low = -2.0,\
                                   range_high = + 10.0,\
                                   bins = 48,\
                                   xlabel = "E_{EM}^{#DeltaR<0.1}[GeV]",\
                                   ylabel = "Number of Tracks")

            histogram_name = "MomentumHadFrac"
            histogram_name = histogram_name + "_NonZeroE_20TRT_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_MomentumHadFrac
            plotter.BookHistograms(histogram_name,\
                                   calc_MomentumHadFrac,\
                                   list_selections = selections,\
                                   range_low = -1.0,\
                                   range_high = + 5.0,\
                                   bins = 48,\
                                   xlabel = "E^{HAD}/P",\
                                   ylabel = "Number of Tracks")

            histogram_name = "HadFrac"
            histogram_name = histogram_name + "_NonZeroE_20TRT_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            from variables.variables import calc_HadFrac
            plotter.BookHistograms(histogram_name,\
                                   calc_HadFrac,\
                                   list_selections = selections,\
                                   range_low = -1.0,\
                                   range_high = + 2.0,\
                                   bins = 48,\
                                   xlabel = "E^{HAD}/E^{Total}",\
                                   ylabel = "Number of Tracks")

            histogram_names = ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]
            xlabels = ["Number of Clusters","Number of Clusters in EM Calorimeter","Number of Clusters in HAD Calorimeter","Number of Clusters with EM Prob > 0.5","Number of Clusters with EM Prob < 0.5"]
            variables = [calc_trkNClusters, calc_trkNClusters_EM, calc_trkNClusters_HAD, calc_trkNClusters_emlike, calc_trkNClusters_hadlike]

            for histogram_name, variable, xlabel in zip(histogram_names, variables, xlabels):
                histogram_name = histogram_name + "_NonZeroE_20TRT_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
                plotter.BookHistograms(histogram_name,\
                                       variable,\
                                       list_selections = selections,\
                                       bins = 10,\
                                       range_low = -0.5,\
                                       range_high = 9.5,\
                                       xlabel=xlabel,\
                                       ylabel="Number of Tracks")




#   #go and get the average eneryg in the EM anulus for MIP particles in each of the eta bins.
#   from variables.variables import calc_EnergyAnulus
#   for eta_range, eta_descriptor, eta_binSelection in zip(eta_ranges, eta_descriptors, eta_binSelections):
#       selections = MIP_selection + [eta_binSelection]
#       AverageAnulus =  plotter.BookTProfileHistograms(calc_trkP,\
#                                                  calc_EnergyAnulus,\
#                                                 list_selections = selections,\
#                                                 bins = p_bins,\
#                                                 xlabel ="P[GeV]",\
#                                                 ylabel = "<E_{EM Anulus}>[GeV]",\
#                                                 )
#       description = ["MIP Selection"] + [eta_descriptor] + ["P_{T} Reweighted"]
#       trkAnulus_canvas = DrawDataVsMC(     AverageAnulus,\
#                                         plotter.channelLabels,\
#                                         MCKey='PythiaJetJet',\
#                                         DataKey='LowMuData',\
#                                         extra_description = description,\
#                                         scale_factor = scale_factor,\
#                                         ratio_min = 0.6,\
#                                         ratio_max = 1.4,\
#                                         doLogx = True,\
#                                         doLogy = False\
#                                         )
#       trkAnulus_canvas.Print(plotter_directory+"/AverageEnergyInAnulusVsMomentumInEtaBin" + str(eta_range[0]) + "_" + str(eta_range[1]) + ".png")

#   #go and get the average eneryg in the EM anulus for MIP particles in each of the eta bins.
#   from variables.variables import calc_EOPBkg
#   for eta_range, eta_descriptor, eta_binSelection in zip(eta_ranges, eta_descriptors, eta_binSelections):
#       selections = MIP_selection + [eta_binSelection]
#       AverageAnulus =  plotter.BookTProfileHistograms(calc_trkP,\
#                                                  calc_EnergyAnulus,\
#                                                 list_selections = selections,\
#                                                 bins = p_bins,\
#                                                 xlabel ="P[GeV]",\
#                                                 ylabel = "<E/p>_{BKG}",\
#                                                 )
#       description = ["MIP Selection"] + [eta_descriptor] + ["P_{T} Reweighted"]
#       trkAnulus_canvas = DrawDataVsMC(     AverageAnulus,\
#                                         plotter.channelLabels,\
#                                         MCKey='PythiaJetJet',\
#                                         DataKey='LowMuData',\
#                                         extra_description = description,\
#                                         scale_factor = scale_factor,\
#                                         ratio_min = 0.6,\
#                                         ratio_max = 1.4,\
#                                         doLogx = True,\
#                                         doLogy = False\
#                                         )
#       trkAnulus_canvas.Print(plotter_directory+"/AverageEOPBkgVsMomentumInEtaBin" + str(eta_range[0]) + "_" + str(eta_range[1]) + ".png")
    histograms = plotter.DumpHistograms()
    for histogram_name in histograms:
        WriteToFile(histograms[histogram_name], outFile)
    print("THEJOBFINISHED!")
