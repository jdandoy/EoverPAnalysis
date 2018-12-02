#!/usr/bin/env python
# coding: utf-8
from PlottingTools.Plotter import Plotter, DrawDataVsMC, DivideHistograms,Draw2DHistogramOnCanvas
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
    from variables.variables import calc_trkNearestNeighbourEM2, calc_trkP, calc_EOP, calc_trkPt, calc_trkAverageMu, calc_trkEtaID, calc_trkEtaECAL, calc_trkNPV2, calc_trkCount
    #import the selections that we want to plot
    from selections.selections import sel_NTRT20, sel_NTRT25, sel_NTRT30, sel_ECALEta0_6, sel_PGreater1 ,sel_PGreater1_5, sel_PGreater2, sel_PGreater2_5, sel_Z0SinThetaLess1_5, sel_d0Less1_5, sel_Event
    #impot the ID selections that we want to plot
    from selections.selections import sel_IDEta00_06,sel_IDEta06_11,sel_IDEta11_14,sel_IDEta14_15, sel_IDEta15_18, sel_IDEta18_23

    outFile = ROOT.TFile(outputRootFileName, "RECREATE")

    ##just count the number of tracks in each histogram
    histogram_name = "trkCount"
    selections = []
    trkCountHist = plotter.GetHistograms(histogram_name,\
                                                         calc_trkCount,\
                                                         list_selections = selections,\
                                                         bins = 1,\
                                                         range_low = -0.5,\
                                                         range_high = +0.5,\
                                                         xlabel ='Always 0',\
                                                         ylabel = 'Number of Tracks')
    WriteToFile(trkCountHist, outFile)

    ################################################################################
    #plot the first set of variables that we're interested in
    #plot the trk avaerage mu histogram
    histogram_name = "trkAverageMu"
    selections = []
    trkAverageMuHist = plotter.GetHistograms(histogram_name,\
                                                        calc_trkAverageMu,\
                                                        list_selections = selections,\
                                                        bins = 10,\
                                                        range_low = 0.0,\
                                                        range_high = 10.0,\
                                                        xlabel ='Average #mu of Event',\
                                                        ylabel = 'Number of Tracks')
    WriteToFile(trkAverageMuHist, outFile)
   ################################################################################
   #plot a histogram of the average event NPV
    histogram_name = "trkNPV2"
    trkNPV2Hist = plotter.GetHistograms(histogram_name,\
                                       calc_trkNPV2,\
                                       list_selections = [],\
                                       bins = 13,\
                                       range_low = -0.5,\
                                       range_high = 12.5,\
                                       xlabel ="NPV with 2 Tracks",\
                                       ylabel = "Number of Tracks")
    WriteToFile(trkNPV2Hist, outFile)
#   ################################################################################yy
   #plot a histogram of the average event NPV
    histogram_name = "eventNPV2Hist"
    eventNPV2Hist = plotter.GetHistograms(histogram_name,\
                                       calc_trkNPV2,\
                                       list_selections = [sel_Event],\
                                       bins = 13,\
                                       range_low = -0.5,\
                                       range_high = 12.5,\
                                       xlabel ="NPV with 2 Tracks",\
                                       ylabel = "Number Events")
    WriteToFile(eventNPV2Hist, outFile)
#   ################################################################################yy
    histogram_name = "eventAverageMu"
    selections = [sel_Event]
    eventAverageMuHist = plotter.GetHistograms(histogram_name,
                                          calc_trkAverageMu,\
                                          list_selections = selections,\
                                          bins = 10,\
                                          range_low = 0.0,\
                                          range_high = 10.0,\
                                          xlabel ='<#mu>',\
                                          ylabel = 'Number of Events')
    WriteToFile(eventAverageMuHist, outFile)
#    ################################################################################yy
    #prepare the momentum bins
    binMax = 30.0
    binLow = 0.5
    nBins = 100
    base = (binMax/binLow) ** (1./float(nBins))
    p_bins = []
    min_p = []
    for i in range(0, nBins + 1):
        p_bins.append(0.5 * (base) ** i )
    histogram_name = "trkPtHist"
    trkPtHistZoom = plotter.GetHistograms(histogram_name,\
                                       calc_trkPt,\
                                       list_selections = [],\
                                       bins = p_bins,\
                                       xlabel ="Track P_{T} [GeV]",\
                                       ylabel = "Number of Tracks")

    WriteToFile(trkPtHistZoom, outFile)
#           ################################################################################yy
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

    for (etaSelection, eta_selectionDescription, file_description) in zip(etaSelections, eta_selectionDescriptions, file_descriptions):
        #do the eta selection and count the inclusive number of tracks in the bin
        selections = [etaSelection]
        trkMultiplicity_Eta = plotter.GetHistograms("TrkPtHist"+file_description,\
                                                  calc_trkPt,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="Track P_{T} [GeV]",\
                                                  ylabel = "Number of Tracks",\
                                                  )
        WriteToFile(trkMultiplicity_Eta, outFile)

    ################################################################################
    histogramName = "TrackEtaID"
    trkEtaIDHist = plotter.GetHistograms(histogramName,\
                                       calc_trkEtaID,\
                                       list_selections = [],\
                                       bins = 100,\
                                       range_low = -5,\
                                       range_high = +5,\
                                       xlabel ="Track #eta ID",\
                                       ylabel = "Number of Tracks")
    WriteToFile(trkEtaIDHist, outFile)

#   ################################################################################yy
    histogramName = "TwoDTrackPvsTrkEtaID"
    bin_size = 0.1
    max_bin = 2.4
    min_bin = -2.4
    eta_bins = []
    eta_bins.append(min_bin)
    while abs(eta_bins[-1] - max_bin) > 0.0001:
        eta_bins.append(eta_bins[-1] + bin_size)
    TwoDtrkPvstrkEta = plotter.Get2DHistograms(histogramName,\
                                             calc_trkEtaID,\
                                             calc_trkP,\
                                             list_selections=[],\
                                             bins_x=eta_bins,\
                                             xlabel="Track #eta ID",\
                                             bins_y=p_bins,\
                                             ylabel="Track P [GeV]",\
                                             zlabel="Number of Tracks",\
                                             )
    WriteToFile(TwoDtrkPvstrkEta, outFile)

#   ################################################################################yy
    histogramName = "TwoDTrackPtVsEtaHistogram"
    TwoDtrkPtvstrkEta = plotter.Get2DHistograms(histogramName,\
                                             calc_trkEtaID,\
                                             calc_trkPt,\
                                             list_selections=[],\
                                             bins_x=eta_bins,\
                                             xlabel="Track #eta ID",\
                                             bins_y=p_bins,\
                                             ylabel="Track P_{T} [GeV]",\
                                             zlabel="Number of Tracks",\
                                             )

    WriteToFile(TwoDtrkPtvstrkEta, outFile)

#   ################################################################################
    histogramName = "trkEtaECALHist"
    trkEtaECALHist = plotter.GetHistograms(histogramName,\
                                          calc_trkEtaECAL,
                                          list_selections = [],
                                          bins = 100,
                                          range_low = -5,
                                          range_high = +5,
                                          xlabel ="Track #eta EM Layer 2",
                                          ylabel = "Number of Tracks",
                                       )
    WriteToFile(trkEtaECALHist, outFile)

#   ################################################################################
    histogramName = "TwoDHistTrkPvsPhiInnerToExtrapolEM2"
    dPhi_bins = []
    min_bin = 0.0
    max_bin = pi
    NBins = 100.0
    bin_size = (max_bin-min_bin)/NBins
    dPhi_bins = []
    dPhi_bins.append(min_bin)
    while abs(dPhi_bins[-1] - max_bin) > 0.0001:
        dPhi_bins.append(dPhi_bins[-1] + bin_size)

    from variables.variables import calc_trkDPhi
    TwoDtrkPvstrkDPhi = plotter.Get2DHistograms(histogramName,\
                                             calc_trkDPhi,\
                                             calc_trkPt,\
                                             list_selections=[],\
                                             bins_x=dPhi_bins,\
                                             xlabel="|#phi_{ID} - #phi_{EM2}|",\
                                             bins_y=p_bins,\
                                             ylabel="Track P_{T} [GeV]",\
                                             zlabel="Number of Tracks",\
                                             )
    WriteToFile(TwoDtrkPvstrkDPhi, outFile)

#   ################################################################################
    histogramName = "lowPTLess07_TwoDHistTrkEtavsDEtaInnerToExtrapolEM2"
    from variables.variables import calc_trkDEta
    from calculation.calculation import calculation
    def lowPT(trk):
        return trk["trk_pt"] < 0.7
    branches =["trk_pt"]
    sel_lowPT = calculation(lowPT, branches)
    TwoDtrkEtavstrkDEta = plotter.Get2DHistograms(histogramName,\
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
    WriteToFile(TwoDtrkEtavstrkDEta, outFile)

#   ################################################################################
    from calculation.calculation import calculation
    from selections.selections import EtaBin
    Eta00_08 = lambda x : EtaBin(x, 0.0, 0.8)
    sel_Eta00_08 = calculation(Eta00_08, ["trk_etaID"])
    histogramName = "EtaLess08_TwoDHistTrkPvsPhiInnerToExtrapolEM2"
    CentalTwoDtrkPvstrkDPhi = plotter.Get2DHistograms(histogramName,\
                                             calc_trkDPhi,\
                                             calc_trkPt,\
                                             list_selections=[sel_Eta00_08],\
                                             bins_x=dPhi_bins,\
                                             xlabel="|#phi_{ID} - #phi_{EM2}|",\
                                             bins_y=p_bins,\
                                             ylabel="Track P_{T} [GeV]",\
                                             zlabel="Number of Tracks",\
                                             )
    WriteToFile(CentalTwoDtrkPvstrkDPhi, outFile)

#    ################################################################################
    histogramName = "NearestDRHist"
    trkNearestDRHist = plotter.GetHistograms(histogramName,
                                    calc_trkNearestNeighbourEM2,
                                    list_selections = [],
                                    bins = 25,
                                     range_low = 0.0,
                                     range_high = 5,
                                     xlabel ="dR to Nearest Track",
                                    ylabel = "Number of Tracks",
                                    )
    WriteToFile(trkNearestDRHist, outFile)

    from selections.selections import sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum, sel_PGreater2, sel_PGreater2_5, sel_PGreater3
    MIP_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum]

#  ################################################################################yy
#  #count the prevalence of certain hadrons in the measurement
#  from selections.selections import sel_Pion, sel_AntiPion, sel_Proton, sel_AntiProton, sel_Kaon, sel_AntiKaon, sel_Fake, sel_Pileup, sel_HardScatter, sel_Muon, sel_AntiMuon, sel_Electron, sel_AntiElectron

#  #The total number of tracks
#  selections = []
#  total_tracks = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_Fake]
#  total_fakes = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_Pileup]
#  total_pileup = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_HardScatter]
#  total_hardScatter = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  #What different species of hadrons are present??
#  selections = [sel_Pion]
#  total_Pion = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_AntiPion]
#  total_AntiPion = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_Proton]
#  total_Proton = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_AntiProton]
#  total_AntiProton = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_Kaon]
#  total_Kaon = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_AntiKaon]
#  total_AntiKaon = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_Electron]
#  total_Electron = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_AntiElectron]
#  total_AntiElectron = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_Muon]
#  total_Muon = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  selections = [sel_AntiMuon]
#  total_AntiMuon = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#  print(10 * "\n")
#  print(50 * "=")
#  print("The total number of tracks was " + str(total_tracks))
#  print("The total number of fakes was " + str(total_fakes))
#  print("The total number of pileup tracks was " + str(total_pileup))
#  print("The total number of truth-matched tracks was " + str(total_hardScatter))
#  print("The fraction of truth-matched tracks was " + str(total_hardScatter/total_tracks))
#  print("The fraction of pileup tracks was " + str(total_pileup/total_tracks))
#  print("The fraction of fake tracks was " + str(total_fakes/total_tracks))
#  print("\n")
#  print("The following fractions are w.r.t the total number of truth matched tracks")
#  print("The fraction of pi + tracks was " + str(total_Pion/total_hardScatter))
#  print("The fraction of pi - tracks was " + str(total_AntiPion/total_hardScatter))
#  print("The fraction of k + tracks was " + str(total_Kaon/total_hardScatter))
#  print("The fraction of k - tracks was " + str(total_AntiKaon/total_hardScatter))
#  print("The fraction of p + tracks was " + str(total_Proton/total_hardScatter))
#  print("The fraction of p - tracks was " + str(total_AntiProton/total_hardScatter))
#  print("The fraction of e - tracks was " + str(total_Electron/total_hardScatter))
#  print("The fraction of e + tracks was " + str(total_AntiElectron/total_hardScatter))
#  print("The fraction of mu - tracks was " + str(total_Muon/total_hardScatter))
#  print("The fraction of mu + tracks was " + str(total_AntiMuon/total_hardScatter))
#  print(5 * "\n")
#  print(50 * "=")

    ################################################################################yy
    selections = []
    histogramName = "InclusiveEOP"
    trkEOPHist = plotter.GetHistograms(histogramName,\
                                     calc_EOP,
                                     list_selections = selections,
                                     bins = 50,
                                     range_low = -1,
                                     range_high = 5,
                                     xlabel ="E/p",
                                     ylabel = "Number of Tracks",
                                     )
    WriteToFile(trkEOPHist, outFile)


    ################################################################################yy
    from selections.selections import sel_NonZeroEnergy
    selections = [sel_NonZeroEnergy]
    histogramName = "NonZeroEnergy_InclusiveEOP"
    trkEOPHist = plotter.GetHistograms(histogramName,\
                                     calc_EOP,
                                     list_selections = selections,
                                     bins = 50,
                                     range_low = -1,
                                     range_high = 5,
                                     xlabel ="E/p",
                                     ylabel = "Number of Tracks",
                                     )
    WriteToFile(trkEOPHist, outFile)

    ################################################################################yy
    selections = [sel_PGreater1, sel_ECALEta0_6]
    histogramName = "EtaID0_6_PGreater1_0_EOPHist"
    trkEOPHistPGreater1 = plotter.GetHistograms(histogramName,\
                                               calc_EOP,
                                               list_selections = selections,
                                               bins = 50,
                                               range_low = -1,
                                               range_high = 5,
                                               xlabel ="E/p",
                                               ylabel = "Number of Tracks")
    WriteToFile(trkEOPHistPGreater1, outFile)

    ################################################################################yy
    selections = [sel_PGreater1_5, sel_IDEta00_06]
    histogramName = "EtaID0_6_PGreater1_5_EOPHist"
    trkEOPHistPGreater1_5 = plotter.GetHistograms(histogramName,
                                                calc_EOP,
                                                list_selections = selections,
                                                 bins = 50,
                                                 range_low = -1,
                                                 range_high = 5,
                                                 xlabel ="E/p",
                                                 ylabel = "Number of Tracks")
    WriteToFile(trkEOPHistPGreater1_5, outFile)

    ################################################################################yy
    selections = [sel_PGreater2, sel_IDEta00_06]
    histogramName = "EtaID0_6_PGreater2_0_EOPHist"
    trkEOPHistPGreater2 = plotter.GetHistograms(histogramName,
                                              calc_EOP,
                                              list_selections = selections,
                                              bins = 50,
                                               range_low = -1,
                                              range_high = 5,
                                              xlabel ="E/p",
                                               ylabel = "Number of Tracks",
                                              )
    WriteToFile(trkEOPHistPGreater2, outFile)


    ################################################################################yy
    from selections.selections import sel_IDEta19_23, sel_IDEta00_06, sel_PBetween12_18, sel_PBetween22_28, sel_PBetween28_36
    # This is figure 2a and 2d in the paper:
    selections = [sel_PBetween12_18, sel_IDEta00_06]
    histogramName = "EtaID0_6_PBetween12_18_EOPHist"
    trkEOPHistFig2a = plotter.GetHistograms(histogramName,
                                              calc_EOP,
                                              list_selections = selections,
                                              bins = 50,
                                               range_low = -1,
                                              range_high = 5,
                                              xlabel ="E/p",
                                               ylabel = "Number of Tracks",
                                              )
    WriteToFile(trkEOPHistFig2a, outFile)

    ################################################################################yy
    from selections.selections import sel_PBetween22_28
    histogramName = "EtaID0_6_PBetween22_28_EOPHist"
    selections = [sel_PBetween22_28, sel_IDEta00_06]
    trkEOPHistFig2b = plotter.GetHistograms(histogramName,\
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.75,\
                                          range_high = 4,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                         )
    WriteToFile(trkEOPHistFig2b, outFile)

    ################################################################################yy
    # This is figure 2c in the paper:
    selections = [sel_PBetween28_36, sel_IDEta19_23]
    histogramName = "EtaIDBetween19_23_PBetween28_36_EOPHist"
    trkEOPHistFig2c = plotter.GetHistograms(histogramName,\
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.75,\
                                          range_high = 4,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                          )
    WriteToFile(trkEOPHistFig2c, outFile)

    ################################################################################yy
    # This is figure 2c in the paper:
    selections = [sel_PBetween28_36, sel_IDEta00_06]
    histogramName = "EtaID0_6_PBetween28_36_EOPHist"
    trkEOPHistFig2c = plotter.GetHistograms(histogramName,
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.5,\
                                          range_high = 3,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                          )
    WriteToFile(trkEOPHistFig2c, outFile)

    ################################################################################
    from selections.selections import sel_NonZeroEnergy
    # This is figure 2a and 2d in the paper:
    histogramName = "NonZero_EtaID0_6_PBetween12_18_EOPHist"
    selections = [sel_PBetween12_18, sel_IDEta00_06, sel_NonZeroEnergy]
    trkEOPHistFig2a = plotter.GetHistograms(histogramName,\
                                              calc_EOP,
                                              list_selections = selections,
                                              bins = 50,
                                               range_low = -1,
                                              range_high = 5,
                                              xlabel ="E/p",
                                               ylabel = "Number of Tracks",
                                              )
    WriteToFile(trkEOPHistFig2a, outFile)

    ################################################################################yy
    from selections.selections import sel_PBetween22_28
    selections = [sel_PBetween22_28, sel_IDEta00_06, sel_NonZeroEnergy]
    histogramName = "NonZero_EtaID0_6_PBetween22_28_EOPHist"
    trkEOPHistFig2b = plotter.GetHistograms(histogramName,\
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.75,\
                                          range_high = 4,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                         )
    WriteToFile(trkEOPHistFig2b, outFile)

    ################################################################################yy
    # This is figure 2c in the paper:
    selections = [sel_PBetween28_36, sel_IDEta19_23, sel_NonZeroEnergy]
    histogramName = "NonZero_EtaIDBetween19_23_PBetween22_28_EOPHist"
    trkEOPHistFig2c = plotter.GetHistograms(histogramName,
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.75,\
                                          range_high = 4,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                          )
    WriteToFile(trkEOPHistFig2c, outFile)

    ################################################################################yy
    # This is figure 2c in the paper:
    selections = [sel_PBetween28_36, sel_IDEta00_06, sel_NonZeroEnergy]
    histogramName = "NonZero_EtaID0_6_PBetween28_36_EOPHist"
    trkEOPHistFig2c = plotter.GetHistograms(histogramName,
                                          calc_EOP,\
                                          list_selections = selections,\
                                          bins = 50,\
                                          range_low = -0.5,\
                                          range_high = 3,\
                                          xlabel ="E/p",\
                                          ylabel = "Number of Tracks",\
                                          )

    WriteToFile(trkEOPHistFig2c, outFile)

    ################################################################################
    # This is figure 3a in the paper:

    selections = []
    binMax = 10.05
    binLow = 0.5
    nBins = 15
    base = (binMax/binLow) ** (1./float(nBins))
    bins = []
    min_p = []
    for i in range(0, nBins + 1):
        bins.append(binLow * (base) ** i )
    histogramName = "InclusiveZeroFractionVsPDenomenator"
    trkMultiplicity = plotter.GetHistograms(histogramName,\
                                          calc_trkP,\
                                          list_selections = selections,\
                                          bins = bins,\
                                          xlabel ="Track P [GeV]",\
                                          ylabel = "Number of Tracks",\
                                          )
    WriteToFile(trkMultiplicity, outFile)
    from selections.selections import sel_ELessEqual0
    histogramName = "InclusiveZeroFractionVsPNumerator"
    selections = [sel_ELessEqual0]
    trkMultiplicity_ELessZero = plotter.GetHistograms(histogramName,\
                                                    calc_trkP,\
                                                    list_selections = selections,\
                                                    bins = bins,\
                                                    xlabel ="Track P [GeV]",\
                                                    ylabel = "N(E<=0)/N",\
                                                    )
    WriteToFile(trkMultiplicity_ELessZero, outFile)
#  ratio_histogram = DivideHistograms(trkMultiplicity_ELessZero, trkMultiplicity)
#  description = ["Inclusive Selection",\
#                  "Track P_{T} Reweighted"]
#  scale_factor = 5.0
#  DataVsMCTrackLess0 = DrawDataVsMC(ratio_histogram,\
#                                     plotter.channelLabels,\
#                                     MCKey='PythiaJetJet',\
#                                     DataKey='LowMuData',\
#                                     extra_description = description,\
#                                     scale_factor = scale_factor,\
#                                     ratio_min = 0.8,\
#                                     ratio_max = 1.2,\
#                                     doLogx = True,\
#                                     doLogy = False,\
#                                     xTicksNumber = 510\
#                                     )
#  DataVsMCTrackLess0.Draw()
#  DataVsMCTrackLess0.Print(plotter_directory+"/EOPAcceptanceVsPInclusive_reweighted.png")

    ################################################################################yy
    #This is figure 3b of the paper
    bins = [-2.3, -1.8, -1.5, -1.4, -1.1, -0.6, 0.0, 0.6, 1.1, 1.4, 1.5, 1.8, 2.3]
    selections = []
    histogramName = "InclusiveZeroFractionVsEtaDenomenator"
    trkMultiplicity_Eta = plotter.GetHistograms(histogramName,\
                                              calc_trkEtaID,\
                                              list_selections = selections,\
                                              bins = bins,\
                                              xlabel ="Track |#eta|",\
                                              ylabel = "Number of Tracks",\
                                              )
    WriteToFile(trkMultiplicity_Eta, outFile)
    from selections.selections import sel_ELessEqual0
    histogramName = "InclusiveZeroFractionVsEtaNumerator"
    selections = [sel_ELessEqual0]
    trkMultiplicity_Eta_Zero = plotter.GetHistograms(histogramName,\
                                                   calc_trkEtaID,\
                                                   list_selections = selections,\
                                                   bins = bins,\
                                                   xlabel ="Track |#eta|",\
                                                   ylabel = "N(E<=0)/N",\
                                                   )
    WriteToFile(trkMultiplicity_Eta_Zero, outFile)
#  DataVsMCTrackLess0 = DrawDataVsMC(ratio_histogram,\
#                                     plotter.channelLabels,\
#                                     MCKey='PythiaJetJet',\
#                                     DataKey='LowMuData',\
#                                     extra_description = description,\
#                                     scale_factor = scale_factor,\
#                                     ratio_min = 0.9,\
#                                     ratio_max = 1.1,\
#                                     doLogx = False,\
#                                     doLogy = False\
#                                     )
#  DataVsMCTrackLess0.Draw()
#  DataVsMCTrackLess0.Print(plotter_directory+"/EOPAcceptanceVsEtaInclusive_reweighted.png")


    ################################################################################yy
    bins = [0.0, 0.6, 1.1, 1.4, 1.5, 1.8, 2.3]
    from variables.variables import calc_trkEta_ABS
    selections = []
    histogramName = "InclusiveZeroFractionVsAbsEtaDenomenator"
    trkMultiplicity_AbsEta = plotter.GetHistograms(histogramName,\
                                              calc_trkEta_ABS,\
                                              list_selections = selections,\
                                              bins = bins,\
                                              xlabel ="Track |#eta|",\
                                              ylabel = "Number of Tracks",\
                                              )
    WriteToFile(trkMultiplicity_AbsEta, outFile)
    from selections.selections import sel_ELessEqual0
    histogramName = "InclusiveZeroFractionVsAbsEtaNumerator"
    selections = [sel_ELessEqual0]
    trkMultiplicity_AbsEta_Zero = plotter.GetHistograms(histogramName,\
                                                   calc_trkEta_ABS,\
                                                   list_selections = selections,\
                                                   bins = bins,\
                                                   xlabel ="Track |#eta|",\
                                                   ylabel = "N(E<=0)/N",\
                                                   )
    WriteToFile(trkMultiplicity_AbsEta_Zero, outFile)

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

#    eta_selectionDescriptions = [\
#                               "|#eta_{ID}|<0.6",\
#                               "0.6<|#eta_{ID}|<1.1",\
#                               "1.1<|#eta_{ID}|<1.4",\
#                               "1.4<|#eta_{ID}|<1.5",\
#                               "1.5<|#eta_{ID}|<1.8",\
#                               "1.8<|#eta_{ID}|<2.3"\
#                               ]



    canvases = []
    keep_histograms_alive = []

    file_descriptions = ["etaID00_02", "etaID02_04", "etaID04_06", "etaID00_06", "etaID06_11", "etaID11_14", "etaID14_15", "etaID15_18", "etaID18_23"]
    centers = [0.2, 0.4, 0.6, 0.6, 1.1, 1.4, 1.5, 1.8, 2.3]

    for (etaSelection, eta_selectionDescription, file_description, center) in zip(etaSelections, eta_selectionDescriptions, file_descriptions, centers):
        binMax = 10.05
        binLow = 0.5 / np.cos(2 * np.arctan(np.exp(-1.0 * center))) ## get the lower bin right for each plot
        nBins = 15
        base = (binMax/binLow) ** (1./float(nBins))
        bins = []
        for i in range(0, nBins + 1):
           bins.append(binLow * (base) ** i )

        #do the eta selection and count the inclusive number of tracks in the bin
        selections = [etaSelection]
        histogramName = "ZeroFractionVsP" + file_description + "Denomenator"
        trkMultiplicity_Eta = plotter.GetHistograms(histogramName,\
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = bins,\
                                                  xlabel ="Track P [GeV]",\
                                                  ylabel = "Number of tracks",\
                                                  )
        WriteToFile(trkMultiplicity_Eta, outFile)

        #do the eta selections and count the number of tracks with an energy deposity less than or equal to 0.0.
        from selections.selections import sel_ELessEqual0
        selections = [sel_ELessEqual0] + [etaSelection]
        histogramName = "ZeroFractionVsP" + file_description + "Numerator"
        trkMultiplicity_Eta_Zero = plotter.GetHistograms(histogramName,\
                                                       calc_trkP,\
                                                       list_selections = selections,\
                                                       bins = bins,\
                                                       xlabel ="Track P [GeV]",\
                                                       ylabel = "N(E<=0)/N",\
                                                       )
        WriteToFile(trkMultiplicity_Eta_Zero, outFile)

    ################################################################################
    from variables.variables import calc_trkEta_ABS
    from selections.selections import sel_NTRT20
    base_description = ["N_{TRT hits} >= 20"]

    for (etaSelection, eta_selectionDescription, file_description, center) in zip(etaSelections, eta_selectionDescriptions, file_descriptions, centers):
        binMax = 10.05
        binLow = 0.5 / np.cos(2 * np.arctan(np.exp(-1.0 * center))) ## get the lower bin right for each plot
        nBins = 15
        base = (binMax/binLow) ** (1./float(nBins))
        bins = []
        for i in range(0, nBins + 1):
           bins.append(binLow * (base) ** i )

        #do the eta selection and count the inclusive number of tracks in the bin
        selections = [etaSelection] + [sel_NTRT20]
        histogramName = "NTRT20ZeroFractionVsP" + file_description + "Denomenator"
        trkMultiplicity_Eta = plotter.GetHistograms(histogramName,\
                                                  calc_trkP,\
                                                  list_selections = selections,\
                                                  bins = bins,\
                                                  xlabel ="Track P [GeV]",\
                                                  ylabel = "Number of tracks",\
                                                  )
        WriteToFile(trkMultiplicity_Eta, outFile)

        #do the eta selections and count the number of tracks with an energy deposity less than or equal to 0.0.
        from selections.selections import sel_ELessEqual0
        selections = [sel_ELessEqual0] + [etaSelection] + [sel_NTRT20]
        histogramName = "NTRT20ZeroFractionVsP" + file_description + "Numerator"
        trkMultiplicity_Eta_Zero = plotter.GetHistograms(histogramName,\
                                                       calc_trkP,\
                                                       list_selections = selections,\
                                                       bins = bins,\
                                                       xlabel ="Track P [GeV]",\
                                                       ylabel = "N(E<=0)/N",\
                                                       )
        WriteToFile(trkMultiplicity_Eta_Zero, outFile)

    ################################################################################
    from selections.selections import sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum, sel_PGreater2, sel_PGreater2_5, sel_PGreater3
    MIP_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum]
    selections = [] + MIP_selection
    histogramName = "MIPSelection_HadBetween30And90OfMomentum_EOP"
    trkEOPHist = plotter.GetHistograms(histogramName,
                                     calc_EOP,
                                     list_selections = selections,
                                     bins = 50,
                                     range_low = -1,
                                     range_high = 5,
                                     xlabel ="E/p",
                                     ylabel = "Number of Tracks",
                                     )


    ################################################################################yy
    selections = [sel_ECALEta0_6] + MIP_selection
    histogramName = "MIPSelection_HadBetween30And90OfMomentum_ECALEta00_06_EOP"
    trkEOPHistPGreater1 = plotter.GetHistograms(histogramName,\
                                               calc_EOP,
                                               list_selections = selections,
                                               bins = 50,
                                               range_low = -1,
                                               range_high = 5,
                                               xlabel ="E/p",
                                               ylabel = "Number of Tracks")

#  ################################################################################yy
#  selections = [sel_PGreater1_5, sel_IDEta00_06] + MIP_selection
#  trkEOPHistPGreater1_5 = plotter.GetHistograms(calc_EOP,
#                                              list_selections = selections,
#                                               bins = 50,
#                                               range_low = -1,
#                                               range_high = 5,
#                                               xlabel ="E/p",
#                                               ylabel = "Number of Tracks")
#  description = ["MIP Selection",\
#                 "P[GeV]>1.5",\
#                 "|#eta_{ID}|<0.6",\
#                 "P_{T} Reweighted"]
#  DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater1_5,
#                              plotter.channelLabels,
#                              MCKey='PythiaJetJet',
#                             DataKey='LowMuData',
#                             ratio_min = 0.5,\
#                             ratio_max = 1.5,\
#                             extra_description = description)
#  DataVsMCEOP.Draw()
#  DataVsMCEOP.Print(plotter_directory+"/EOPP1_5Eta0_6_reweighted_MIP.png")
#  CloseCanvas(DataVsMCEOP)

#  ################################################################################yy
#  from calculation.calculation import calculation
#  from selections.selections import EtaBin, PBin
#  PBinFunction = lambda x: PBin(x, 0.5, 1.0)
#  sel_Pbin = calculation(PBinFunction, ["trk_p"])
#  selections = [sel_Pbin, sel_IDEta00_06] + MIP_selection
#  trkEOP                 = plotter.GetHistograms(calc_EOP,
#                                              list_selections = selections,
#                                               bins = 50,
#                                               range_low = -1,
#                                               range_high = 3,
#                                               xlabel ="E/p",
#                                               ylabel = "Number of Tracks")
#  description = ["MIP Selection",\
#                 "0.5<p[GeV]<1.0",\
#                 "|#eta_{ID}|<0.6",\
#                 "P_{T} Reweighted"]
#  DataVsMCEOP = DrawDataVsMC(trkEOPHist,\
#                              plotter.channelLabels,
#                              MCKey='PythiaJetJet',
#                             DataKey='LowMuData',
#                             ratio_min = 0.6,\
#                             ratio_max = 1.4,\
#                             extra_description = description)
#  DataVsMCEOP.Draw()
#  DataVsMCEOP.Print(plotter_directory+"/EOPP05_08Eta06_reweighted_MIP.png")
#  CloseCanvas(DataVsMCEOP)



#  ################################################################################yy
#  selections = [sel_PGreater2, sel_IDEta00_06] + MIP_selection
#  trkEOPHistPGreater2 = plotter.GetHistograms(calc_EOP,
#                                            list_selections = selections,
#                                            bins = 50,
#                                             range_low = -1,
#                                            range_high = 5,
#                                            xlabel ="E/p",
#                                             ylabel = "Number of Tracks",
#                                            )
#  description = ["MIP Selection",\
#                 "P[GeV]>2",\
#                 "|#eta_{ID}|<0.6",\
#                 "P_{T} Reweighted"]
#  DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater2,
#                              plotter.channelLabels,
#                              MCKey='PythiaJetJet',
#                              DataKey='LowMuData',
#                             extra_description = description)
#  DataVsMCEOP.Draw()
#  DataVsMCEOP.Print(plotter_directory+"/EOPP2Eta0_6_reweighted_MIP.png")
#  CloseCanvas(DataVsMCEOP)


#  ################################################################################yy
#  selections = [sel_PGreater2_5, sel_IDEta00_06] + MIP_selection
#  trkEOPHistPGreater2 = plotter.GetHistograms(calc_EOP,
#                                            list_selections = selections,
#                                            bins = 50,
#                                             range_low = -1,
#                                            range_high = 5,
#                                            xlabel ="E/p",
#                                             ylabel = "Number of Tracks",
#                                            )
#  description = ["MIP Selection",\
#                 "P[GeV]>2.5",\
#                 "|#eta_{ID}|<0.6",\
#                 "P_{T} Reweighted"]
#  DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater2,
#                              plotter.channelLabels,
#                              MCKey='PythiaJetJet',
#                              DataKey='LowMuData',
#                             extra_description = description)
#  DataVsMCEOP.Draw()
#  DataVsMCEOP.Print(plotter_directory+"/EOPP2_5Eta0_6_reweighted_MIP.png")
#  CloseCanvas(DataVsMCEOP)

#  ################################################################################yy
#  selections = [sel_PGreater3, sel_IDEta00_06] + MIP_selection
#  trkEOPHistPGreater2 = plotter.GetHistograms(calc_EOP,
#                                            list_selections = selections,
#                                            bins = 50,
#                                             range_low = -1,
#                                            range_high = 5,
#                                            xlabel ="E/p",
#                                             ylabel = "Number of Tracks",
#                                            )
#  description = ["MIP Selection",\
#                 "P[GeV]>3",\
#                 "|#eta_{ID}|<0.6",\
#                 "P_{T} Reweighted"]
#  DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater2,
#                              plotter.channelLabels,
#                              MCKey='PythiaJetJet',
#                              DataKey='LowMuData',
#                             extra_description = description)
#  DataVsMCEOP.Draw()
#  DataVsMCEOP.Print(plotter_directory+"/EOPP3Eta0_6_reweighted_MIP.png")
#  CloseCanvas(DataVsMCEOP)


#  ################################################################################yy
#  from selections.selections import sel_IDEta19_23, sel_IDEta00_06, sel_PBetween12_18, sel_PBetween22_28, sel_PBetween28_36
#  selections = [sel_PBetween12_18, sel_IDEta00_06] + MIP_selection
#  trkEOPHistFig2a = plotter.GetHistograms(calc_EOP,
#                                            list_selections = selections,
#                                            bins = 50,
#                                             range_low = -1,
#                                            range_high = 5,
#                                            xlabel ="E/p",
#                                             ylabel = "Number of Tracks",
#                                            )
#  description = ["MIP Selection",\
#                 "1.2<P[GeV]<1.8",\
#                 "|#eta_{ID}|<0.6",\
#                 "P_{T} Reweighted"]
#  DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2a,\
#                             plotter.channelLabels,\
#                             MCKey='PythiaJetJet',\
#                             DataKey='LowMuData',\
#                             ratio_min = 0.5,\
#                             ratio_max = 1.5,\
#                             extra_description = description)
#  DataVsMCEOP.Draw()
#  DataVsMCEOP.Print(plotter_directory+"/EOPP12_18_Eta0_6_reweighted_MIP.png")

#  ################################################################################yy
#  from selections.selections import sel_PBetween22_28
#  selections = [sel_PBetween22_28, sel_IDEta00_06] + MIP_selection
#  trkEOPHistFig2b = plotter.GetHistograms(calc_EOP,\
#                                        list_selections = selections,\
#                                        bins = 50,\
#                                        range_low = -0.75,\
#                                        range_high = 4,\
#                                        xlabel ="E/p",\
#                                        ylabel = "Number of Tracks",\
#                                       )
#  description = ["MIP Selection",\
#                 "2.2<P[GeV]<2.8",\
#                 "|#eta_{ID}|<0.6",\
#                 "P_{T} Reweighted"]
#  DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2b,\
#                              plotter.channelLabels,\
#                              MCKey='PythiaJetJet',\
#                              DataKey='LowMuData',\
#                             ratio_min = 0.5,\
#                             ratio_max = 1.5,\
#                              extra_description = description)
#  DataVsMCEOP.Draw()
#  DataVsMCEOP.Print(plotter_directory+"/EOPP22_28_Eta0_6_reweighted_MIP.png")

#  ################################################################################yy
#  selections = [sel_PBetween28_36, sel_IDEta19_23] + MIP_selection
#  trkEOPHistFig2c = plotter.GetHistograms(calc_EOP,\
#                                        list_selections = selections,\
#                                        bins = 50,\
#                                        range_low = -0.75,\
#                                        range_high = 4,\
#                                        xlabel ="E/p",\
#                                        ylabel = "Number of Tracks",\
#                                        )
#  description = ["MIP Selection",\
#                 "2.8<P[GeV]<3.6",\
#                 "1.9<|#eta_{ID}|<2.3",\
#                 "P_{T} Reweighted"]
#  DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2c,\
#                             plotter.channelLabels,\
#                             MCKey='PythiaJetJet',\
#                             DataKey='LowMuData',\
#                             ratio_min = 0.5,\
#                             ratio_max = 1.5,\
#                             extra_description = description\
#                             )
#  DataVsMCEOP.Draw()
#  DataVsMCEOP.Print(plotter_directory+"/EOPP28_36_Eta19_23_reweighted_MIP.png")

#  ################################################################################yy
#  selections = [sel_PBetween28_36, sel_IDEta00_06] + MIP_selection
#  trkEOPHistFig2c = plotter.GetHistograms(calc_EOP,\
#                                        list_selections = selections,\
#                                        bins = 50,\
#                                        range_low = -0.5,\
#                                        range_high = 3,\
#                                        xlabel ="E/p",\
#                                        ylabel = "Number of Tracks",\
#                                        )
#  description = ["MIP Selection",\
#                 "2.8<P[GeV]<3.6",\
#                 "|#eta_{ID}|<0.6",\
#                 "P_{T} Reweighted"]
#  DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2c,\
#                             plotter.channelLabels,\
#                             MCKey='PythiaJetJet',\
#                             DataKey='LowMuData',\
#                             ratio_min = 0.5,\
#                             ratio_max = 1.5,\
#                             extra_description = description\
#                             )
#  DataVsMCEOP.Draw()
#  DataVsMCEOP.Print(plotter_directory+"/EOPP28_36_Eta0_6_reweighted_MIP.png")

#   ################################################################################yy
#   ## Look in different bins of pseudorapidity
#   from variables.variables import calc_trkEta_ABS
#   base_description = []
#   etaSelections = [sel_IDEta00_06,\
#                    sel_IDEta06_11,\
#                    sel_IDEta11_14,\
#                    sel_IDEta14_15,\
#                    sel_IDEta15_18,\
#                    sel_IDEta18_23]

#   eta_selectionDescriptions = [\
#                              "|#eta_{ID}|<0.6",\
#                              "0.6<|#eta_{ID}|<1.1",\
#                              "1.1<|#eta_{ID}|<1.4",\
#                              "1.4<|#eta_{ID}|<1.5",\
#                              "1.5<|#eta_{ID}|<1.8",\
#                              "1.8<|#eta_{ID}|<2.3"\
#                              ]

#   file_descriptions = ["eta06", "eta06_11", "eta11_14", "eta14_15", "eta15_18", "eta18_23"]

#           for (etaSelection, eta_selectionDescription, file_description) in zip(etaSelections, eta_selectionDescriptions, file_descriptions):
#               #do the eta selection and count the inclusive number of tracks in the bin
#               selections = [etaSelection] + MIP_selection
#               trkMultiplicity_Eta = plotter.GetHistograms(calc_trkPt,\
#                                                         list_selections = selections,\
#                                                         bins = 80,\
#                                                         range_low = 0.5,\
#                                                         range_high = 5,\
#                                                         xlabel ="Track P_{T} [GeV]",\
#                                                         ylabel = "Number of Tracks",\
#                                                         )

#               description = ["MIP Selection", eta_selectionDescription] + ["P_{T} Reweighted"]
#               trkMultiplicity_canvas = DrawDataVsMC(trkMultiplicity_Eta,\
#                                                 plotter.channelLabels,\
#                                                 MCKey='PythiaJetJet',\
#                                                 DataKey='LowMuData',\
#                                                 extra_description = description,\
#                                                 scale_factor = scale_factor,\
#                                                 ratio_min = 0.8,\
#                                                 ratio_max = 1.2,\
#                                                 doLogx = False,\
#                                                 doLogy = False\
#                                                 )
#               trkMultiplicity_canvas.Print(plotter_directory+"/TrkPtReweighted"+file_description+"MIP.png")

#           ################################################################################
#           from variables.variables import calc_trkEMFraction, calc_trkHADFraction
#           from variables.variables import calc_trkEta_ABS
#           from selections.selections import sel_NonZeroEnergy

#           base_description = []
#           etaSelections = [sel_IDEta00_06,\
#                            sel_IDEta06_11,\
#                            sel_IDEta11_14,\
#                            sel_IDEta14_15,\
#                            sel_IDEta15_18,\
#                            sel_IDEta18_23\
#                            ]

#           eta_selectionDescriptions = [\
#                                      "|#eta_{ID}|<0.6",\
#                                      "0.6<|#eta_{ID}|<1.1",\
#                                      "1.1<|#eta_{ID}|<1.4",\
#                                      "1.4<|#eta_{ID}|<1.5",\
#                                      "1.5<|#eta_{ID}|<1.8",\
#                                      "1.8<|#eta_{ID}|<2.3"\
#                                      ]

#           file_descriptions = ["eta06", "eta06_11", "eta11_14", "eta14_15", "eta15_18", "eta18_23"]

#           for (etaSelection, eta_selectionDescription, file_description) in zip(etaSelections, eta_selectionDescriptions, file_descriptions):
#               #do the eta selection and count the inclusive number of tracks in the bin
#               selections = [etaSelection] + [sel_NonZeroEnergy]
#               trkMultiplicity_Eta = plotter.GetHistograms(calc_trkEMFraction,\
#                                                         list_selections = selections,\
#                                                         bins = 22,\
#                                                         range_low = -1.0,\
#                                                         range_high = 1.2,\
#                                                         xlabel ="E_{EM}/E_{Total}",\
#                                                         ylabel = "Number of Tracks",\
#                                                         )

#               description = ["E_{Total} != 0", eta_selectionDescription] + ["P_{T} Reweighted"]
#               trkMultiplicity_canvas = DrawDataVsMC(trkMultiplicity_Eta,\
#                                                 plotter.channelLabels,\
#                                                 MCKey='PythiaJetJet',\
#                                                 DataKey='LowMuData',\
#                                                 extra_description = description,\
#                                                 scale_factor = scale_factor,\
#                                                 ratio_min = 0.6,\
#                                                 ratio_max = 1.4,\
#                                                 doLogx = False,\
#                                                 doLogy = False\
#                                                 )
#               trkMultiplicity_canvas.Print(plotter_directory+"/TrkEMFractionReweighted"+file_description+".png")

#           ################################################################################
#           from variables.variables import calc_trkEMFraction, calc_trkHADFraction
#           from variables.variables import calc_trkEta_ABS
#           from selections.selections import sel_NonZeroEnergy

#           base_description = []
#           etaSelections = [sel_IDEta00_06,\
#                            sel_IDEta06_11,\
#                            sel_IDEta11_14,\
#                            sel_IDEta14_15,\
#                            sel_IDEta15_18,\
#                            sel_IDEta18_23\
#                            ]

#           eta_selectionDescriptions = [\
#                                      "|#eta_{ID}|<0.6",\
#                                      "0.6<|#eta_{ID}|<1.1",\
#                                      "1.1<|#eta_{ID}|<1.4",\
#                                      "1.4<|#eta_{ID}|<1.5",\
#                                      "1.5<|#eta_{ID}|<1.8",\
#                                      "1.8<|#eta_{ID}|<2.3"\
#                                      ]

#           file_descriptions = ["eta06", "eta06_11", "eta11_14", "eta14_15", "eta15_18", "eta18_23"]

#           for (etaSelection, eta_selectionDescription, file_description) in zip(etaSelections, eta_selectionDescriptions, file_descriptions):
#               #do the eta selection and count the inclusive number of tracks in the bin
#               selections = [etaSelection] + [sel_NonZeroEnergy]
#               trkMultiplicity_Eta = plotter.GetHistograms(calc_trkHADFraction,\
#                                                         list_selections = selections,\
#                                                         bins = 22,\
#                                                         range_low = -1.0,\
#                                                         range_high = 1.2,\
#                                                         xlabel ="E_{HAD}/E_{Total}",\
#                                                         ylabel = "Number of Tracks",\
#                                                         )

#               description = ["E_{Total} != 0", eta_selectionDescription] + ["P_{T} Reweighted"]
#               trkMultiplicity_canvas = DrawDataVsMC(trkMultiplicity_Eta,\
#                                                 plotter.channelLabels,\
#                                                 MCKey='PythiaJetJet',\
#                                                 DataKey='LowMuData',\
#                                                 extra_description = description,\
#                                                 scale_factor = scale_factor,\
#                                                 ratio_min = 0.6,\
#                                                 ratio_max = 1.4,\
#                                                 doLogx = False,\
#                                                 doLogy = False\
#                                                 )
#               trkMultiplicity_canvas.Print(plotter_directory+"/TrkHADFractionReweighted"+file_description+".png")

    ################################################################################
    ##Create a set of p and eta bins for the measurement ##########################
    from calculation.calculation import calculation
    from selections.selections import EtaBin, PBin
    from variables.variables import calc_EOPBkg, calc_EnergyAnulus

    #create a set of strings that could describe the eta or momentum selections
    eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]

    #go and get the average E/P for MIP particles in each of the eta bins.
    for eta_range, eta_binSelection in zip(eta_ranges, eta_binSelections):
        #get the function that selectts tracks in that bin
        EtaBinFunction = lambda x: EtaBin(x, eta_range[0], eta_range[1])
        sel_EtaBin = calculation(EtaBinFunction, ["trk_etaID"])

        #calculate the lowest momentum track that can end up in that bin
        center = eta_range[1]
        binMax = 10.05
        binLow = 0.5 / np.cos(2 * np.arctan(np.exp(-1.0 * center))) ## get the lower bin right for each plot
        nBins = 15
        base = (binMax/binLow) ** (1./float(nBins))
        p_bins = []
        for i in range(0, nBins + 1):
           p_bins.append(binLow * (base) ** i )

        #get the bins of the EOP variable
        eop_bins_min = -1
        eop_bins_max = 5
        nbins = 120
        step = float(eop_bins_max-eop_bins_min)/float(nbins)
        eop_bins = []
        for i in range(0, nBins + 1):
            eop_bins.append(eop_bins_min + i * step)


        MIP_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum]
        selections = MIP_selection + [eta_binSelection]
        histogramName = "EOPProfileVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.GetTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )
        WriteToFile(AverageEOP, outFile)

        histogramName = "2DHist_EOPVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.Get2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p",\
                                                  )
        WriteToFile(AverageEOP, outFile)

        histogramName = "EnergyAnulusProfileVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.GetTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )
        WriteToFile(AverageAnulus, outFile)

        histogramName = "2DHist_EnergyAnulusVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Get2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )
        WriteToFile(AverageAnulus, outFile)

        histogramName = "EnergyBkgProfileVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.GetTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )
        WriteToFile(AverageAnulus, outFile)

        histogramName = "2DHist_EnergyBkgVsMomentum_MIPSelection_HadBetween30And90OfMomentum_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Get2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{BKG}/p",\
                                                  )
        WriteToFile(AverageAnulus, outFile)

        #go and get the E/p distribution in each of the E/p bins
        p_ranges = [(x,y) for x,y in zip(p_bins[0:-1], p_bins[1:])]
        for p_range in p_ranges:
            PBinFunction = lambda x: PBin(x, p_range[0], p_range[1])
            sel_PBin = calculation(PBinFunction, ["trk_p"])
            selections = MIP_selection + [eta_binSelection, sel_PBin]
            histogramName = "EOPDistribution_MIPSelection_HadBetween30And90OfMomentum_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            AverageEOP  =  plotter.GetHistograms(histogramName,
                                                      calc_EOP,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )
            WriteToFile(AverageEOP, outFile)

        from selections.selections import sel_EHadFracAbove70

        MIP_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadFracAbove70]
        selections = MIP_selection + [eta_binSelection]
        histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.GetTProfileHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>",\
                                                  )
        WriteToFile(AverageEOP, outFile)

        histogramName = "2DHist_EOPVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageEOP  =  plotter.Get2DHistograms(histogramName,
                                                  calc_trkP,\
                                                  calc_EOP,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E/p",\
                                                  )
        WriteToFile(AverageEOP, outFile)

        histogramName = "EnergyAnulusProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.GetTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E_{EM Anulus}>[GeV]",\
                                                  )
        WriteToFile(AverageAnulus, outFile)

        histogramName = "2DHist_EnergyAnulusVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Get2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EnergyAnulus,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{EM Anulus} [GeV]",\
                                                  )
        WriteToFile(AverageAnulus, outFile)

        histogramName = "EnergyBkgProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.GetTProfileHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins = p_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "<E/p>_{BKG}",\
                                                  )
        WriteToFile(AverageAnulus, outFile)

        histogramName = "2DHist_EnergyBkgVsMomentum_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
        AverageAnulus =  plotter.Get2DHistograms(histogramName,\
                                                  calc_trkP,\
                                                  calc_EOPBkg,\
                                                  list_selections = selections,\
                                                  bins_x = p_bins,\
                                                  bins_y = eop_bins,\
                                                  xlabel ="P[GeV]",\
                                                  ylabel = "E_{BKG}/p",\
                                                  )
        WriteToFile(AverageAnulus, outFile)

        #go and get the E/p distribution in each of the E/p bins
        p_ranges = [(x,y) for x,y in zip(p_bins[0:-2], p_bins[1:-1])]
        for p_range in p_ranges:
            PBinFunction = lambda x: PBin(x, p_range[0], p_range[1])
            sel_PBin = calculation(PBinFunction, ["trk_p"])
            selections = MIP_selection + [eta_binSelection, sel_PBin]
            histogramName = "EOPDistribution_MIPSelection_HadFracAbove70_InEtaBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
            AverageEOP  =  plotter.GetHistograms(histogramName,
                                                      calc_EOP,\
                                                      list_selections = selections,\
                                                      bins = eop_bins,\
                                                      xlabel ="E/p",\
                                                      )
            WriteToFile(AverageEOP, outFile)


#   #go and get the average eneryg in the EM anulus for MIP particles in each of the eta bins.
#   from variables.variables import calc_EnergyAnulus
#   for eta_range, eta_descriptor, eta_binSelection in zip(eta_ranges, eta_descriptors, eta_binSelections):
#       selections = MIP_selection + [eta_binSelection]
#       AverageAnulus =  plotter.GetTProfileHistograms(calc_trkP,\
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
#       AverageAnulus =  plotter.GetTProfileHistograms(calc_trkP,\
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
    print("THEJOBFINISHED!")
