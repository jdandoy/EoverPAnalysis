#!/usr/bin/env python
# coding: utf-8

from PlottingTools.Plotter import Plotter, DrawDataVsMC, DivideHistograms,Draw2DHistogramOnCanvas
import ROOT
#from variables.variables import calc_weight
#from inputs.samples import INPUT
import os
from math import pi
import pickle

def CloseCanvas(canv):
    canv.Close()
    ROOT.gSystem.ProcessEvents()
    del canv

def WriteToFile(histogram_dictionary, outFile):
    outFile.cd()
    for key in histogram_dictionary:
        if not outFile.cd(key):
            outFile.mkdir(key)
        outFile.cd(key)
        print("Writing histogram " + key)
        histogram_dictionary[key].Write()

#This is a script that fills the histograms for
def FillingScript(plotter, outputRootFileName):
    #import thje variables that we want to plot
    from variables.variables import calc_trkNearestNeighbourEM2, calc_trkP, calc_EOP, calc_trkPt, calc_trkAverageMu, calc_trkEtaID, calc_trkEtaECAL, calc_trkNPV2, calc_trkCount
    #import the selections that we want to plot
    from selections.selections import sel_NTRT20, sel_NTRT25, sel_NTRT30, sel_ECALEta0_6, sel_PGreater1 ,sel_PGreater1_5, sel_PGreater2, sel_PGreater2_5, sel_Z0SinThetaLess1_5, sel_d0Less1_5, sel_Event
    #impot the ID selections that we want to plot
    from selections.selections import sel_IDEta0_6,sel_IDEta06_11,sel_IDEta11_14,sel_IDEta14_15, sel_IDEta15_18, sel_IDEta18_23

    outFile = ROOT.TFile(outputRootFileName, "RECREATE")

    ##just count the number of tracks in each histogram
    histogram_name = "trkCount"
    selections = []
    trkCountHist = plotter.GetHistograms(histogram_name,\
                                                         calc_trkCount,\
                                                         list_selections = selections,\
                                                         bins = 1,\
                                                         range_low = -0.5,\
                                                         range_high = +1.5,\
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
##    description = "Inclusive Selection"
##    DataVsMC = DrawDataVsMC(trkAverageMuHist,\
##                            plotter.channelLabels,\
##                            MCKey='PythiaJetJet',\
##                            DataKey='LowMuData',\
##                            extra_description = description)
##    DataVsMC.Draw()
##    DataVsMC.Print(plotter_directory+"/TrkAverageMu_unweighted.png")
##    CloseCanvas(DataVsMC)


   ################################################################################
   #plot a histogram of the average event NPV
    histogram_name = "trkNPV2"
    eventNPV2Hist = plotter.GetHistograms(histogram_name,\
                                       calc_trkNPV2,\
                                       list_selections = [],\
                                       bins = 13,\
                                       range_low = -0.5,\
                                       range_high = 12.5,\
                                       xlabel ="NPV with 2 Tracks",\
                                       ylabel = "Number of Events")
    WriteToFile(eventNPV2Hist, outFile)
#    description = "Inclusive Selection"
#    eventNPV2HistCanvas   = DrawDataVsMC(eventNPV2Hist,\
#                                       plotter.channelLabels,\
#                                       MCKey='PythiaJetJet',\
#                                       DataKey='LowMuData',\
#                                       extra_description = description)
#    eventNPV2HistCanvas.Draw()
#    eventNPV2HistCanvas.Print(plotter_directory+"/EventNPV2Tracks_beforeNPVReweight.png")
#    CloseCanvas(DataVsMC)


#    plotter.UseVariableAndHistogramToNormalize(calc_trkNPV2,
#                                             eventNPV2Hist,
#                                             "PythiaJetJet",
#                                             "LowMuData",
#                                            )

#    plotter_directory = base_plotter_directory + "/NPVReweighted"
#    if not os.path.exists(plotter_directory):
#        os.mkdir(plotter_directory)

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
##   description = "Inclusive Selection"
##   eventNPV2HistCanvas = DrawDataVsMC(eventNPV2Hist,\
##                                      plotter.channelLabels,\
##                                      MCKey='PythiaJetJet',\
##                                      DataKey='LowMuData',\
##                                      extra_description = description)
##   eventNPV2HistCanvas.Draw()
##   eventNPV2HistCanvas.Print(plotter_directory+"/EventNPV2Tracks_afterNPVReweight.png")


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
#    description = "Inclusive Selection"
#    DataVsMC = DrawDataVsMC(trkAverageMuHist,\
#                            plotter.channelLabels,\
#                            MCKey='PythiaJetJet',\
#                            DataKey='LowMuData',\
#                            extra_description = description)
#    DataVsMC.Draw()
#    DataVsMC.Print(plotter_directory+"/TrkAverageMu_NPVReweighted.png")
#    CloseCanvas(DataVsMC)


#    ################################################################################yy
#    #plot a histogram of the average trk NPV
#    trkNPV2Hist = plotter.GetHistograms(calc_trkNPV2,\
#                                        list_selections = [],\
#                                        bins = 13,\
#                                        range_low = -0.5,\
#                                       range_high = 12.5,\
#                                       xlabel ="NPV with 2 Tracks",\
#                                       ylabel = "Number of tracks")
   
#   description = "Inclusive Selection"
#   trkNPV2HistCanvas = DrawDataVsMC(trkNPV2Hist,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',\
#                                      ratio_min = 0.8,\
#                                      ratio_max = 1.2,\
#                                      extra_description = description)
#   trkNPV2HistCanvas.Draw()
#   trkNPV2HistCanvas.Print(plotter_directory+"/TrkNPV2Tracks_afterNPVReweight.png")


    #prepare the momentum bins
    binMax = 30.0
    binLow = 0.5
    nBins = 100
    base = (binMax/binLow) ** (1./float(nBins))
    p_bins = []
    min_p = []
    for i in range(0, nBins + 1):
        p_bins.append(0.5 * (base) ** i )

    ################################################################################yy
    #plot a histogram of the average track pt
    histogram_name = "trkPtHist"
    trkPtHistZoom = plotter.GetHistograms(histogram_name,\
                                        calc_trkPt,\
                                        list_selections = [],\
                                        bins = p_bins,\
                                        xlabel ="Track P_{T} [GeV]",\
                                        ylabel = "Number of Tracks")

    WriteToFile(trkPtHistZoom, outFile)
#   description = "Inclusive Selection"
#   trkPtHistZoomCanvas = DrawDataVsMC(trkPtHistZoom,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',\
#                                      ratio_min = 0.8,\
#                                      ratio_max = 1.2,\
#                                      doLogy = True,
#                                      doLogx = True,
#                                      extra_description = description)
#   trkPtHistZoomCanvas.Draw()
#   trkPtHistZoomCanvas.Print(plotter_directory+"/TrkPtZoomHistogram_unweighted.png")


#           ################################################################################yy
#           ## Look in different bins of pseudorapidity
    base_description = []
    etaSelections = [sel_IDEta0_6,\
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


#               description = [eta_selectionDescription]
#               trkMultiplicity_canvas = DrawDataVsMC(trkMultiplicity_Eta,\
#                                                 plotter.channelLabels,\
#                                                 MCKey='PythiaJetJet',\
#                                                 DataKey='LowMuData',\
#                                                 extra_description = description,\
#                                                 scale_factor = scale_factor,\
#                                                 ratio_min = 0.6,\
#                                                 ratio_max = 1.4,\
#                                                 doLogx = True,\
#                                                 doLogy = True\
#                                                 )
#               trkMultiplicity_canvas.Print(plotter_directory+"/TrkPtUnweighted"+file_description+".png")

#   ################################################################################
#   trkEtaIDHist = plotter.GetHistograms(calc_trkEtaID,\
#                                      list_selections = [],\
#                                      bins = 100,\
#                                      range_low = -5,\
#                                      range_high = +5,\
#                                      xlabel ="Track #eta ID",\
#                                      ylabel = "Number of Tracks")
#   description = "Inclusive Selection"
#   DataVsMCHistCanvas = DrawDataVsMC(trkEtaIDHist,\
#                                     plotter.channelLabels,\
#                                     MCKey='PythiaJetJet',\
#                                     DataKey='LowMuData',\
#                                     ratio_min = 0.9,\
#                                     ratio_max = 1.1,\
#                                     extra_description = description)
#   DataVsMCHistCanvas.Draw()
#   DataVsMCHistCanvas.Print(plotter_directory+"/TrkEtaInInnerDetector_unweighted.png")
#   CloseCanvas(DataVsMCHistCanvas)


#   ################################################################################yy
#   bin_size = 0.1
#   max_bin = 2.4
#   min_bin = -2.4
#   eta_bins = []
#   eta_bins.append(min_bin)
#   while abs(eta_bins[-1] - max_bin) > 0.0001:
#       eta_bins.append(eta_bins[-1] + bin_size)
#   TwoDtrkPvstrkEta = plotter.Get2DHistograms(calc_trkEtaID,\
#                                            calc_trkP,\
#                                            list_selections=[],\
#                                            bins_x=eta_bins,\
#                                            xlabel="Track #eta ID",\
#                                            bins_y=p_bins,\
#                                            ylabel="Track P [GeV]",\
#                                            zlabel="Number of Tracks",\
#                                            normalize=False)
#   can = Draw2DHistogramOnCanvas(TwoDtrkPvstrkEta["PythiaJetJet"], doLogx = False, doLogy = True)
#   can.Draw()
#   can.Print(plotter_directory+"/PythiaJetJetTwoDHistogramEtaTrackP.png")

#   can = Draw2DHistogramOnCanvas(TwoDtrkPvstrkEta["LowMuData"], doLogx = False, doLogy = True)
#   can.Draw()
#   can.Print(plotter_directory+"/LowMuDataTwoDHistogramEtaTrackP.png")

#   ################################################################################yy

#   TwoDtrkPvstrkEta = plotter.Get2DHistograms(calc_trkEtaID,\
#                                            calc_trkPt,\
#                                            list_selections=[],\
#                                            bins_x=eta_bins,\
#                                            xlabel="Track #eta ID",\
#                                            bins_y=p_bins,\
#                                            ylabel="Track P_{T} [GeV]",\
#                                            zlabel="Number of Tracks",\
#                                            normalize=False)

#   can = Draw2DHistogramOnCanvas(TwoDtrkPvstrkEta["PythiaJetJet"], doLogx = False, doLogy = True)
#   can.Draw()
#   can.Print(plotter_directory+"/PythiaJetJetTwoDHistogramEtaTrackPt.png")

#   can = Draw2DHistogramOnCanvas(TwoDtrkPvstrkEta["LowMuData"], doLogx = False, doLogy = True)
#   can.Draw()
#   can.Print(plotter_directory+"/LowMuDataTwoDHistogramEtaTrackPt.png")


#   ################################################################################
#   trkEtaECALHist = plotter.GetHistograms(calc_trkEtaECAL,
#                                         list_selections = [],
#                                         bins = 100,
#                                         range_low = -5,
#                                         range_high = +5,
#                                         xlabel ="Track #eta EM Layer 2",
#                                         ylabel = "Number of Tracks",
#                                      )
#   description = "Inclusive Selection"
#   DataVsMCHistECALCanvas = DrawDataVsMC(trkEtaECALHist,
#                                          plotter.channelLabels,
#                                          MCKey='PythiaJetJet',
#                                         DataKey='LowMuData',
#                                         ratio_min = 0.5,
#                                         ratio_max = 1.5,
#                                         extra_description = description)
#   DataVsMCHistECALCanvas.Draw()
#   DataVsMCHistECALCanvas.Print(plotter_directory+"/TrackEtaInTheEndcaCalorimeter_unweighted.png")
#   CloseCanvas(DataVsMCHistECALCanvas)


#   ################################################################################
#   dPhi_bins = []
#   min_bin = 0.0
#   max_bin = pi
#   NBins = 100.0
#   bin_size = (max_bin-min_bin)/NBins
#   dPhi_bins = []
#   dPhi_bins.append(min_bin)
#   while abs(dPhi_bins[-1] - max_bin) > 0.0001:
#       dPhi_bins.append(dPhi_bins[-1] + bin_size)

#   from variables.variables import calc_trkDPhi
#   TwoDtrkPvstrkDPhi = plotter.Get2DHistograms(calc_trkDPhi,\
#                                            calc_trkPt,\
#                                            list_selections=[],\
#                                            bins_x=dPhi_bins,\
#                                            xlabel="|#phi_{ID} - #phi_{EM2}|",\
#                                            bins_y=p_bins,\
#                                            ylabel="Track P_{T} [GeV]",\
#                                            zlabel="Number of Tracks",\
#                                            normalize=False)

#   can = Draw2DHistogramOnCanvas(TwoDtrkPvstrkDPhi["PythiaJetJet"], doLogx = False, doLogy = True)
#   can.Draw()
#   can.Print(plotter_directory+"/PythiaJetJetTwoDHistogramDPhiTrackPt.png")

#   can = Draw2DHistogramOnCanvas(TwoDtrkPvstrkDPhi["LowMuData"], doLogx = False, doLogy = True)
#   can.Draw()
#   can.Print(plotter_directory+"/LowMuDataTwoDHistogramDPhiTrackPt.png")

#   ################################################################################

#   from calculation.calculation import calculation
#   from selections.selections import EtaBin
#   Eta00_08 = lambda x : EtaBin(x, 0.0, 0.8)
#   sel_Eta00_08 = calculation(Eta00_08, ["trk_etaID"])

#   from variables.variables import calc_trkDPhi
#   TwoDtrkPvstrkDPhi = plotter.Get2DHistograms(calc_trkDPhi,\
#                                            calc_trkPt,\
#                                            list_selections=[sel_Eta00_08],\
#                                            bins_x=dPhi_bins,\
#                                            xlabel="|#phi_{ID} - #phi_{EM2}|",\
#                                            bins_y=p_bins,\
#                                            ylabel="Track P_{T} [GeV]",\
#                                            zlabel="Number of Tracks",\
#                                            normalize=False)

#   can = Draw2DHistogramOnCanvas(TwoDtrkPvstrkDPhi["PythiaJetJet"], doLogx = False, doLogy = True)
#   can.Draw()
#   can.Print(plotter_directory+"/Central00To08PythiaJetJetTwoDHistogramDPhiTrackPt.png")

#   can = Draw2DHistogramOnCanvas(TwoDtrkPvstrkDPhi["LowMuData"], doLogx = False, doLogy = True)
#   can.Draw()
#   can.Print(plotter_directory+"/Cental00To08LowMuDataTwoDHistogramDPhiTrackPt.png")

#   ################################################################################
#   plotter.UseVariableAndHistogramToNormalize(calc_trkPt,
#                                            trkPtHist,
#                                             "PythiaJetJet",
#                                             "LowMuData",
#                                           )

#   plotter_directory = base_plotter_directory + "/NPVandPTReweighted"
#   if not os.path.exists(plotter_directory):
#       os.mkdir(plotter_directory)

#   trkNearestDRHist = plotter.GetHistograms(calc_trkNearestNeighbourEM2,
#                                   list_selections = [],
#                                   bins = 25,
#                                    range_low = 0.0,
#                                    range_high = 5,
#                                    xlabel ="dR to Nearest Track",
#                                   ylabel = "Number of Tracks",
#                                   )

#   description = ["Inclusive Selection" ,"P_{T} Reweighted"]
#   DataVsMC = DrawDataVsMC(trkNearestDRHist,
#                           plotter.channelLabels,
#                           MCKey='PythiaJetJet',
#                           DataKey='LowMuData',
#                           ratio_min = 0.2,
#                           ratio_max = 1.8,
#                           extra_description = description)
#   DataVsMC.Draw()
#   DataVsMC.Print(plotter_directory+"/dRToNearestNeighbour_reweighted.png")
#   CloseCanvas(DataVsMC)

#   from selections.selections import sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum, sel_PGreater2, sel_PGreater2_5, sel_PGreater3
#   MIP_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum]

#   if True:
#           ################################################################################yy
#           #count the prevalence of certain hadrons in the measurement
#           from selections.selections import sel_Pion, sel_AntiPion, sel_Proton, sel_AntiProton, sel_Kaon, sel_AntiKaon, sel_Fake, sel_Pileup, sel_HardScatter, sel_Muon, sel_AntiMuon, sel_Electron, sel_AntiElectron

#           #The total number of tracks
#           selections = []
#           total_tracks = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_Fake]
#           total_fakes = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_Pileup]
#           total_pileup = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_HardScatter]
#           total_hardScatter = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           #What different species of hadrons are present??
#           selections = [sel_Pion]
#           total_Pion = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_AntiPion]
#           total_AntiPion = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_Proton]
#           total_Proton = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_AntiProton]
#           total_AntiProton = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_Kaon]
#           total_Kaon = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_AntiKaon]
#           total_AntiKaon = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_Electron]
#           total_Electron = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_AntiElectron]
#           total_AntiElectron = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_Muon]
#           total_Muon = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           selections = [sel_AntiMuon]
#           total_AntiMuon = plotter.GetNumberOfTracks("PythiaJetJet", list_selections = selections)

#           print(10 * "\n")
#           print(50 * "=")
#           print("The total number of tracks was " + str(total_tracks))
#           print("The total number of fakes was " + str(total_fakes))
#           print("The total number of pileup tracks was " + str(total_pileup))
#           print("The total number of truth-matched tracks was " + str(total_hardScatter))
#           print("The fraction of truth-matched tracks was " + str(total_hardScatter/total_tracks))
#           print("The fraction of pileup tracks was " + str(total_pileup/total_tracks))
#           print("The fraction of fake tracks was " + str(total_fakes/total_tracks))
#           print("\n")
#           print("The following fractions are w.r.t the total number of truth matched tracks")
#           print("The fraction of pi + tracks was " + str(total_Pion/total_hardScatter))
#           print("The fraction of pi - tracks was " + str(total_AntiPion/total_hardScatter))
#           print("The fraction of k + tracks was " + str(total_Kaon/total_hardScatter))
#           print("The fraction of k - tracks was " + str(total_AntiKaon/total_hardScatter))
#           print("The fraction of p + tracks was " + str(total_Proton/total_hardScatter))
#           print("The fraction of p - tracks was " + str(total_AntiProton/total_hardScatter))
#           print("The fraction of e - tracks was " + str(total_Electron/total_hardScatter))
#           print("The fraction of e + tracks was " + str(total_AntiElectron/total_hardScatter))
#           print("The fraction of mu - tracks was " + str(total_Muon/total_hardScatter))
#           print("The fraction of mu + tracks was " + str(total_AntiMuon/total_hardScatter))
#           print(5 * "\n")
#           print(50 * "=")

#           ################################################################################yy
#           trkPtReweightedHist = plotter.GetHistograms(calc_trkPt,
#                                           list_selections = [],
#                                           bins = p_bins,
#                                            xlabel ="Track P_{T} [GeV]",
#                                           ylabel = "Number of Tracks",
#                                           )
#           description = ["Inclusive Selection" , "P_{T} Reweighted"]
#           DataVsMC = DrawDataVsMC(trkPtReweightedHist,
#                                    plotter.channelLabels,
#                                   MCKey='PythiaJetJet',
#                                    DataKey='LowMuData',
#                                    doLogx=True,
#                                    doLogy=True,
#                                   extra_description = description)
#           DataVsMC.Draw()
#           DataVsMC.Print(plotter_directory+"/trkMultiplicityVsPt_reweighted.png")
#           CloseCanvas(DataVsMC)

#           ################################################################################yy
#           ## Look in different bins of pseudorapidity
#           from variables.variables import calc_trkEta_ABS
#           base_description = []
#           etaSelections = [sel_IDEta0_6,\
#                            sel_IDEta06_11,\
#                            sel_IDEta11_14,\
#                            sel_IDEta14_15,\
#                            sel_IDEta15_18,\
#                            sel_IDEta18_23]

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
#               selections = [etaSelection]
#               trkMultiplicity_Eta = plotter.GetHistograms(calc_trkPt,\
#                                                         list_selections = selections,\
#                                                         bins = p_bins,\
#                                                         xlabel ="Track P_{T} [GeV]",\
#                                                         ylabel = "Number of Tracks",\
#                                                         normalize = False,\
#                                                         )

#               description = [eta_selectionDescription] + ["P_{T} Reweighted"]
#               trkMultiplicity_canvas = DrawDataVsMC(trkMultiplicity_Eta,\
#                                                 plotter.channelLabels,\
#                                                 MCKey='PythiaJetJet',\
#                                                 DataKey='LowMuData',\
#                                                 extra_description = description,\
#                                                 scale_factor = scale_factor,\
#                                                 ratio_min = 0.8,\
#                                                 ratio_max = 1.2,\
#                                                 doLogx = True,\
#                                                 doLogy = True\
#                                                 )
#               trkMultiplicity_canvas.Print(plotter_directory+"/TrkPtReweighted"+file_description+".png")

#           ################################################################################yy
#           selections = []
#           trkEOPHist = plotter.GetHistograms(calc_EOP,
#                                            list_selections = selections,
#                                            bins = 50,
#                                            range_low = -1,
#                                            range_high = 5,
#                                            xlabel ="E/p",
#                                            ylabel = "Number of Tracks",
#                                            )
#           description  = ["Inclusive Selection",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHist,
#                                       plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                       DataKey='LowMuData',
#                                       ratio_min = 0.6,
#                                      ratio_max = 1.4,
#                                      extra_description = description,
#                                      )
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/InclusiveEOPDisitrbution.png")
#           CloseCanvas(DataVsMCEOP)


#           ################################################################################yy
#           from selections.selections import sel_NonZeroEnergy
#           selections = [sel_NonZeroEnergy]
#           trkEOPHist = plotter.GetHistograms(calc_EOP,
#                                            list_selections = selections,
#                                            bins = 50,
#                                            range_low = -1,
#                                            range_high = 5,
#                                            xlabel ="E/p",
#                                            ylabel = "Normalized Distribution",
#                                            normalize = True,
#                                            )
#           description  = ["E_{Total} != 0",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHist,
#                                       plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                       DataKey='LowMuData',
#                                       ratio_min = 0.5,
#                                       ratio_max = 1.5,
#                                       extra_description = description,
#                                       )
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/InclusiveEOPDisitrbution_NonZeroEnergy.png")
#           CloseCanvas(DataVsMCEOP)

#           ################################################################################yy
#           selections = [sel_PGreater1, sel_ECALEta0_6]
#           trkEOPHistPGreater1 = plotter.GetHistograms(calc_EOP,
#                                                      list_selections = selections,
#                                                      bins = 50,
#                                                      range_low = -1,
#                                                      range_high = 5,
#                                                      xlabel ="E/p",
#                                                      ylabel = "Number of Tracks")
#           description = ["P[GeV]>1",                "|#eta_{EM2}|<0.6",                "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater1,
#                                      plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                      DataKey='LowMuData',
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/Eta06EOPDistribution.png")
#           CloseCanvas(DataVsMCEOP)

#           ################################################################################yy
#           selections = [sel_PGreater1_5, sel_IDEta0_6]
#           trkEOPHistPGreater1_5 = plotter.GetHistograms(calc_EOP,
#                                                       list_selections = selections,
#                                                        bins = 50,
#                                                        range_low = -1,
#                                                        range_high = 5,
#                                                        xlabel ="E/p",
#                                                        ylabel = "Number of Tracks")
#           description = ["P[GeV]>1.5",               "|#eta_{ID}|<0.6",               "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater1_5,
#                                       plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                      DataKey='LowMuData',
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP1_5Eta0_6_reweighted.png")
#           CloseCanvas(DataVsMCEOP)

#           ################################################################################yy
#           selections = [sel_PGreater2, sel_IDEta0_6]
#           trkEOPHistPGreater2 = plotter.GetHistograms(calc_EOP,
#                                                     list_selections = selections,
#                                                     bins = 50,
#                                                      range_low = -1,
#                                                     range_high = 5,
#                                                     xlabel ="E/p",
#                                                      ylabel = "Number of Tracks",
#                                                     )
#           description = ["P[GeV]>2",                "|#eta_{ID}|<0.6",               "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater2,
#                                       plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                       DataKey='LowMuData',
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP2Eta0_6_reweighted.png")
#           CloseCanvas(DataVsMCEOP)


#           ################################################################################yy
#           from selections.selections import sel_IDEta19_23, sel_IDEta0_6, sel_PBetween12_18, sel_PBetween22_28, sel_PBetween28_36
#           # This is figure 2a and 2d in the paper:
#           selections = [sel_PBetween12_18, sel_IDEta0_6]
#           trkEOPHistFig2a = plotter.GetHistograms(calc_EOP,
#                                                     list_selections = selections,
#                                                     bins = 50,
#                                                      range_low = -1,
#                                                     range_high = 5,
#                                                     xlabel ="E/p",
#                                                      ylabel = "Number of Tracks",
#                                                     )
#           description = ["1.2<P[GeV]<1.8",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2a,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP12_18_Eta0_6_reweighted.png")

#           ################################################################################yy
#           from selections.selections import sel_PBetween22_28
#           selections = [sel_PBetween22_28, sel_IDEta0_6]
#           trkEOPHistFig2b = plotter.GetHistograms(calc_EOP,\
#                                                 list_selections = selections,\
#                                                 bins = 50,\
#                                                 range_low = -0.75,\
#                                                 range_high = 4,\
#                                                 xlabel ="E/p",\
#                                                 ylabel = "Number of Tracks",\
#                                                )
#           description = ["2.2<P[GeV]<2.8",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2b,\
#                                       plotter.channelLabels,\
#                                       MCKey='PythiaJetJet',\
#                                       DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                       extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP22_28_Eta0_6_reweighted.png")

#           ################################################################################yy
#           # This is figure 2c in the paper:
#           selections = [sel_PBetween28_36, sel_IDEta19_23]
#           trkEOPHistFig2c = plotter.GetHistograms(calc_EOP,\
#                                                 list_selections = selections,\
#                                                 bins = 50,\
#                                                 range_low = -0.75,\
#                                                 range_high = 4,\
#                                                 xlabel ="E/p",\
#                                                 ylabel = "Number of Tracks",\
#                                                 )
#           description = ["2.8<P[GeV]<3.6",\
#                          "1.9<|#eta_{ID}|<2.3",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2c,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description\
#                                      )
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP28_36_Eta19_23_reweighted.png")

#           ################################################################################yy
#           # This is figure 2c in the paper:
#           selections = [sel_PBetween28_36, sel_IDEta0_6]
#           trkEOPHistFig2c = plotter.GetHistograms(calc_EOP,\
#                                                 list_selections = selections,\
#                                                 bins = 50,\
#                                                 range_low = -0.5,\
#                                                 range_high = 3,\
#                                                 xlabel ="E/p",\
#                                                 ylabel = "Number of Tracks",\
#                                                 )
#           description = ["2.8<P[GeV]<3.6",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2c,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description\
#                                      )
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP28_36_Eta0_6_reweighted.png")

#           ################################################################################
#           from selections.selections import sel_NonZeroEnergy
#           # This is figure 2a and 2d in the paper:
#           selections = [sel_PBetween12_18, sel_IDEta0_6, sel_NonZeroEnergy]
#           trkEOPHistFig2a = plotter.GetHistograms(calc_EOP,
#                                                     list_selections = selections,
#                                                     bins = 50,
#                                                      range_low = -1,
#                                                     range_high = 5,
#                                                     xlabel ="E/p",
#                                                     normalize = True,
#                                                      ylabel = "Normalized Distribution",
#                                                     )
#           description = ["E_{Total} != 0",\
#                          "1.2<P[GeV]<1.8",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2a,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP12_18_Eta0_6_reweighted_nonzeroE.png")

#           ################################################################################yy
#           from selections.selections import sel_PBetween22_28
#           selections = [sel_PBetween22_28, sel_IDEta0_6, sel_NonZeroEnergy]
#           trkEOPHistFig2b = plotter.GetHistograms(calc_EOP,\
#                                                 list_selections = selections,\
#                                                 bins = 50,\
#                                                 range_low = -0.75,\
#                                                 range_high = 4,\
#                                                 xlabel ="E/p",\
#                                                 normalize = True,\
#                                                 ylabel = "Normalized Distribution",\
#                                                )
#           description = ["E_{Total} != 0",\
#                          "2.2<P[GeV]<2.8",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2b,\
#                                       plotter.channelLabels,\
#                                       MCKey='PythiaJetJet',\
#                                       DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                       extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP22_28_Eta0_6_reweighted_nonzeroE.png")

#           ################################################################################yy
#           # This is figure 2c in the paper:
#           selections = [sel_PBetween28_36, sel_IDEta19_23, sel_NonZeroEnergy]
#           trkEOPHistFig2c = plotter.GetHistograms(calc_EOP,\
#                                                 list_selections = selections,\
#                                                 bins = 50,\
#                                                 range_low = -0.75,\
#                                                 range_high = 4,\
#                                                 xlabel ="E/p",\
#                                                 normalize = True,\
#                                                 ylabel = "Normalized Distribution",\
#                                                 )
#           description = ["E_{Total} != 0",\
#                          "2.8<P[GeV]<3.6",\
#                          "1.9<|#eta_{ID}|<2.3",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2c,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description\
#                                      )
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP28_36_Eta19_23_reweighted_nonzeroE.png")

#           ################################################################################yy
#           # This is figure 2c in the paper:
#           selections = [sel_PBetween28_36, sel_IDEta0_6, sel_NonZeroEnergy]
#           trkEOPHistFig2c = plotter.GetHistograms(calc_EOP,\
#                                                 list_selections = selections,\
#                                                 bins = 50,\
#                                                 range_low = -0.5,\
#                                                 range_high = 3,\
#                                                 xlabel ="E/p",\
#                                                 normalize = True,\
#                                                 ylabel = "Normalized Distribution",\
#                                                 )
#           description = ["E_{Total} != 0",\
#                          "2.8<P[GeV]<3.6",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2c,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description\
#                                      )
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP28_36_Eta0_6_reweighted_nonzeroE.png")

#           ################################################################################
#           # This is figure 3a in the paper:

#           selections = []
#           binMax = 10.05
#           binLow = 0.5
#           nBins = 15
#           base = (binMax/binLow) ** (1./float(nBins))
#           bins = []
#           min_p = []
#           for i in range(0, nBins + 1):
#               bins.append(0.5 * (base) ** i )
#           trkMultiplicity = plotter.GetHistograms(calc_trkP,\
#                                                 list_selections = selections,\
#                                                 bins = bins,\
#                                                 xlabel ="Track P [GeV]",\
#                                                 ylabel = "Number of Tracks",\
#                                                 )
#           from selections.selections import sel_ELessEqual0
#           selections = [sel_ELessEqual0]
#           trkMultiplicity_ELessZero = plotter.GetHistograms(calc_trkP,\
#                                                           list_selections = selections,\
#                                                           bins = bins,\
#                                                           xlabel ="Track P [GeV]",\
#                                                           ylabel = "N(E<=0)/N",\
#                                                           )
#           ratio_histogram = DivideHistograms(trkMultiplicity_ELessZero, trkMultiplicity)
#           description = ["Inclusive Selection",\
#                          "Track P_{T} Reweighted"]
#           scale_factor = 5.0
#           DataVsMCTrackLess0 = DrawDataVsMC(ratio_histogram,\
#                                             plotter.channelLabels,\
#                                             MCKey='PythiaJetJet',\
#                                             DataKey='LowMuData',\
#                                             extra_description = description,\
#                                             scale_factor = scale_factor,\
#                                             ratio_min = 0.8,\
#                                             ratio_max = 1.2,\
#                                             doLogx = True,\
#                                             doLogy = False,\
#                                             xTicksNumber = 510\
#                                             )
#           DataVsMCTrackLess0.Draw()
#           DataVsMCTrackLess0.Print(plotter_directory+"/EOPAcceptanceVsPInclusive_reweighted.png")

#           ################################################################################yy
#           #This is figure 3b of the paper
#           bins = [-2.3, -1.8, -1.5, -1.4, -1.1, -0.6, 0.0, 0.6, 1.1, 1.4, 1.5, 1.8, 2.3]
#           selections = []
#           trkMultiplicity_Eta = plotter.GetHistograms(calc_trkEtaID,\
#                                                     list_selections = selections,\
#                                                     bins = bins,\
#                                                     xlabel ="Track |#eta|",\
#                                                     ylabel = "Number of Tracks",\
#                                                     )
#           from selections.selections import sel_ELessEqual0
#           selections = [sel_ELessEqual0]
#           trkMultiplicity_Eta_Zero = plotter.GetHistograms(calc_trkEtaID,\
#                                                          list_selections = selections,\
#                                                          bins = bins,\
#                                                          xlabel ="Track |#eta|",\
#                                                          ylabel = "N(E<=0)/N",\
#                                                          )
#           ratio_histogram = DivideHistograms(trkMultiplicity_Eta_Zero, trkMultiplicity_Eta)
#           description = ["Inclusive Selection",\
#                          "Track P_{T} Reweighted"]
#           scale_factor = 5.0
#           DataVsMCTrackLess0 = DrawDataVsMC(ratio_histogram,\
#                                             plotter.channelLabels,\
#                                             MCKey='PythiaJetJet',\
#                                             DataKey='LowMuData',\
#                                             extra_description = description,\
#                                             scale_factor = scale_factor,\
#                                             ratio_min = 0.9,\
#                                             ratio_max = 1.1,\
#                                             doLogx = False,\
#                                             doLogy = False\
#                                             )
#           DataVsMCTrackLess0.Draw()
#           DataVsMCTrackLess0.Print(plotter_directory+"/EOPAcceptanceVsEtaInclusive_reweighted.png")


#           ################################################################################yy
#           bins = [0.0, 0.6, 1.1, 1.4, 1.5, 1.8, 2.3]
#           selections = []
#           from variables.variables import calc_trkEta_ABS
#           trkMultiplicity_Eta = plotter.GetHistograms(calc_trkEta_ABS,\
#                                                     list_selections = selections,\
#                                                     bins = bins,\
#                                                     xlabel ="Track |#eta|",\
#                                                     ylabel = "Number of Tracks",\
#                                                     )
#           from selections.selections import sel_ELessEqual0
#           selections = [sel_ELessEqual0]
#           trkMultiplicity_Eta_Zero = plotter.GetHistograms(calc_trkEta_ABS,\
#                                                          list_selections = selections,\
#                                                          bins = bins,\
#                                                          xlabel ="Track |#eta|",\
#                                                          ylabel = "N(E<=0)/N",\
#                                                          )
#           ratio_histogram = DivideHistograms(trkMultiplicity_Eta_Zero, trkMultiplicity_Eta)
#           description = ["Inclusive Selection",\
#                          "Track P_{T} Reweighted"]
#           scale_factor = 5.0
#           DataVsMCTrackLess0 = DrawDataVsMC(ratio_histogram,\
#                                             plotter.channelLabels,\
#                                             MCKey='PythiaJetJet',\
#                                             DataKey='LowMuData',\
#                                             extra_description = description,\
#                                             scale_factor = scale_factor,\
#                                             ratio_min = 0.9,\
#                                             ratio_max = 1.1,\
#                                             doLogx = False,\
#                                             doLogy = False\
#                                             )
#           DataVsMCTrackLess0.Draw()
#           DataVsMCTrackLess0.Print(plotter_directory+"/EOPAcceptanceVsEtaAbs_reweighted.png")


#           ################################################################################
#           from variables.variables import calc_trkEta_ABS
#           base_description = ["Track P_{T} Reweighted"]

#           etaSelections = [sel_IDEta0_6,\
#                            sel_IDEta06_11,\
#                            sel_IDEta11_14,\
#                            sel_IDEta14_15,\
#                            sel_IDEta15_18,\
#                            sel_IDEta18_23]

#           eta_selectionDescriptions = [\
#                                      "|#eta_{ID}|<0.6",\
#                                      "0.6<|#eta_{ID}|<1.1",\
#                                      "1.1<|#eta_{ID}|<1.4",\
#                                      "1.4<|#eta_{ID}|<1.5",\
#                                      "1.5<|#eta_{ID}|<1.8",\
#                                      "1.8<|#eta_{ID}|<2.3"\
#                                      ]

#           binMax = 10.05
#           binLow = 0.5
#           nBins = 15
#           base = (binMax/binLow) ** (1./float(nBins))
#           bins = []
#           for i in range(0, nBins + 1):
#               bins.append(0.5 * (base) ** i )

#           canvases = []
#           keep_histograms_alive = []

#           file_descriptions = ["eta06", "eta06_11", "eta11_14", "eta14_15", "eta15_18", "eta18_23"]
#           file

#           for (etaSelection, eta_selectionDescription, file_description) in zip(etaSelections, eta_selectionDescriptions, file_descriptions):
#               #do the eta selection and count the inclusive number of tracks in the bin
#               selections = [etaSelection]
#               trkMultiplicity_Eta = plotter.GetHistograms(calc_trkP,\
#                                                         list_selections = selections,\
#                                                         bins = bins,\
#                                                         xlabel ="Track P [GeV]",\
#                                                         ylabel = "Number of tracks",\
#                                                         normalize = False,\
#                                                         )

#               #do the eta selections and count the number of tracks with an energy deposity less than or equal to 0.0.
#               from selections.selections import sel_ELessEqual0
#               selections = [sel_ELessEqual0] + [etaSelection]
#               trkMultiplicity_Eta_Zero = plotter.GetHistograms(calc_trkP,\
#                                                              list_selections = selections,\
#                                                              bins = bins,\
#                                                              xlabel ="Track P [GeV]",\
#                                                              ylabel = "N(E<=0)/N",\
#                                                              normalize = False,\
#                                                              )
#               ratio_histogram = DivideHistograms(trkMultiplicity_Eta_Zero, trkMultiplicity_Eta)
#               keep_histograms_alive.append(ratio_histogram)
#               description = [eta_selectionDescription] + base_description
#               scale_factor = 5.0
#               DataVsMCTrackLess0 = DrawDataVsMC(ratio_histogram,\
#                                                 plotter.channelLabels,\
#                                                 MCKey='PythiaJetJet',\
#                                                 DataKey='LowMuData',\
#                                                 extra_description = description,\
#                                                 scale_factor = scale_factor,\
#                                                 ratio_min = 0.6,\
#                                                 ratio_max = 1.4,\
#                                                 doLogx = True,\
#                                                 doLogy = False\
#                                                 )
#               canvases.append(DataVsMCTrackLess0)
#               DataVsMCTrackLess0.Print(plotter_directory+"/EOPAcceptanceVsPtEtaBin"+file_description+".png")

#           ################################################################################

#           ################################################################################
#           from variables.variables import calc_trkEta_ABS
#           from selections.selections import sel_NTRT20
#           base_description = ["N_{TRT hits} >= 20"]

#           etaSelections = [sel_IDEta0_6,\
#                            sel_IDEta06_11,\
#                            sel_IDEta11_14,\
#                            sel_IDEta14_15,\
#                            sel_IDEta15_18,\
#                            sel_IDEta18_23]

#           eta_selectionDescriptions = [\
#                                      "|#eta_{ID}|<0.6",\
#                                      "0.6<|#eta_{ID}|<1.1",\
#                                      "1.1<|#eta_{ID}|<1.4",\
#                                      "1.4<|#eta_{ID}|<1.5",\
#                                      "1.5<|#eta_{ID}|<1.8",\
#                                      "1.8<|#eta_{ID}|<2.3"\
#                                      ]

#           binMax = 10.05
#           binLow = 0.5
#           nBins = 15
#           base = (binMax/binLow) ** (1./float(nBins))
#           bins = []
#           for i in range(0, nBins + 1):
#               bins.append(0.5 * (base) ** i )

#           canvases = []
#           keep_histograms_alive = []

#           file_descriptions = ["eta06", "eta06_11", "eta11_14", "eta14_15", "eta15_18", "eta18_23"]

#           for (etaSelection, eta_selectionDescription, file_description) in zip(etaSelections, eta_selectionDescriptions, file_descriptions):
#               #do the eta selection and count the inclusive number of tracks in the bin
#               selections = [etaSelection] + [sel_NTRT20]
#               trkMultiplicity_Eta = plotter.GetHistograms(calc_trkP,\
#                                                         list_selections = selections,\
#                                                         bins = bins,\
#                                                         xlabel ="Track P [GeV]",\
#                                                         ylabel = "Number of tracks",\
#                                                         normalize = False,\
#                                                         )

#               #do the eta selections and count the number of tracks with an energy deposity less than or equal to 0.0.
#               from selections.selections import sel_ELessEqual0
#               selections = [sel_ELessEqual0] + [etaSelection] + [sel_NTRT20]
#               trkMultiplicity_Eta_Zero = plotter.GetHistograms(calc_trkP,\
#                                                              list_selections = selections,\
#                                                              bins = bins,\
#                                                              xlabel ="Track P [GeV]",\
#                                                              ylabel = "N(E<=0)/N",\
#                                                              normalize = False,\
#                                                              )
#               ratio_histogram = DivideHistograms(trkMultiplicity_Eta_Zero, trkMultiplicity_Eta)
#               keep_histograms_alive.append(ratio_histogram)
#               description = [eta_selectionDescription] + base_description + ["Track P_{T} Reweighted"]
#               scale_factor = 5.0
#               DataVsMCTrackLess0 = DrawDataVsMC(ratio_histogram,\
#                                                 plotter.channelLabels,\
#                                                 MCKey='PythiaJetJet',\
#                                                 DataKey='LowMuData',\
#                                                 extra_description = description,\
#                                                 scale_factor = scale_factor,\
#                                                 ratio_min = 0.6,\
#                                                 ratio_max = 1.4,\
#                                                 doLogx = True,\
#                                                 doLogy = False\
#                                                 )
#               canvases.append(DataVsMCTrackLess0)
#               DataVsMCTrackLess0.Print(plotter_directory+"/EOPAcceptanceVsPtEtaBinTRTHits20"+file_description+".png")

#           ################################################################################
#           from selections.selections import sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum, sel_PGreater2, sel_PGreater2_5, sel_PGreater3
#           MIP_selection = [sel_NTRT20, sel_Lar1_1GeV, sel_EHadBetween30And90OfMomentum]
#           selections = [] + MIP_selection
#           trkEOPHist = plotter.GetHistograms(calc_EOP,
#                                            list_selections = selections,
#                                            bins = 50,
#                                            range_low = -1,
#                                            range_high = 5,
#                                            xlabel ="E/p",
#                                            ylabel = "Number of Tracks",
#                                            )
#           description  = ["MIP Selection",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHist,
#                                       plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                       DataKey='LowMuData',
#                                       ratio_min = 0.6,
#                                      ratio_max = 1.4,
#                                      extra_description = description,
#                                      )
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/InclusiveEOPDisitrbution_MIP.png")
#           CloseCanvas(DataVsMCEOP)


#           ################################################################################yy
#           selections = [sel_PGreater1, sel_ECALEta0_6] + MIP_selection
#           trkEOPHistPGreater1 = plotter.GetHistograms(calc_EOP,
#                                                      list_selections = selections,
#                                                      bins = 50,
#                                                      range_low = -1,
#                                                      range_high = 5,
#                                                      xlabel ="E/p",
#                                                      ylabel = "Number of Tracks")
#           description = ["MIP Selection",\
#                          "P[GeV]>1",\
#                          "|#eta_{EM2}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater1,
#                                      plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                      DataKey='LowMuData',
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/Eta06EOPDistribution_MIP.png")
#           CloseCanvas(DataVsMCEOP)

#           ################################################################################yy
#           selections = [sel_PGreater1_5, sel_IDEta0_6] + MIP_selection
#           trkEOPHistPGreater1_5 = plotter.GetHistograms(calc_EOP,
#                                                       list_selections = selections,
#                                                        bins = 50,
#                                                        range_low = -1,
#                                                        range_high = 5,
#                                                        xlabel ="E/p",
#                                                        ylabel = "Number of Tracks")
#           description = ["MIP Selection",\
#                          "P[GeV]>1.5",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater1_5,
#                                       plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                      DataKey='LowMuData',
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP1_5Eta0_6_reweighted_MIP.png")
#           CloseCanvas(DataVsMCEOP)

#           ################################################################################yy
#           from calculation.calculation import calculation
#           from selections.selections import EtaBin, PBin
#           PBinFunction = lambda x: PBin(x, 0.5, 1.0)
#           sel_Pbin = calculation(PBinFunction, ["trk_p"])
#           selections = [sel_Pbin, sel_IDEta0_6] + MIP_selection
#           trkEOP                 = plotter.GetHistograms(calc_EOP,
#                                                       list_selections = selections,
#                                                        bins = 50,
#                                                        range_low = -1,
#                                                        range_high = 3,
#                                                        xlabel ="E/p",
#                                                        ylabel = "Number of Tracks")
#           description = ["MIP Selection",\
#                          "0.5<p[GeV]<1.0",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHist,\
#                                       plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                      DataKey='LowMuData',
#                                      ratio_min = 0.6,\
#                                      ratio_max = 1.4,\
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP05_08Eta06_reweighted_MIP.png")
#           CloseCanvas(DataVsMCEOP)



#           ################################################################################yy
#           selections = [sel_PGreater2, sel_IDEta0_6] + MIP_selection
#           trkEOPHistPGreater2 = plotter.GetHistograms(calc_EOP,
#                                                     list_selections = selections,
#                                                     bins = 50,
#                                                      range_low = -1,
#                                                     range_high = 5,
#                                                     xlabel ="E/p",
#                                                      ylabel = "Number of Tracks",
#                                                     )
#           description = ["MIP Selection",\
#                          "P[GeV]>2",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater2,
#                                       plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                       DataKey='LowMuData',
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP2Eta0_6_reweighted_MIP.png")
#           CloseCanvas(DataVsMCEOP)


#           ################################################################################yy
#           selections = [sel_PGreater2_5, sel_IDEta0_6] + MIP_selection
#           trkEOPHistPGreater2 = plotter.GetHistograms(calc_EOP,
#                                                     list_selections = selections,
#                                                     bins = 50,
#                                                      range_low = -1,
#                                                     range_high = 5,
#                                                     xlabel ="E/p",
#                                                      ylabel = "Number of Tracks",
#                                                     )
#           description = ["MIP Selection",\
#                          "P[GeV]>2.5",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater2,
#                                       plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                       DataKey='LowMuData',
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP2_5Eta0_6_reweighted_MIP.png")
#           CloseCanvas(DataVsMCEOP)

#           ################################################################################yy
#           selections = [sel_PGreater3, sel_IDEta0_6] + MIP_selection
#           trkEOPHistPGreater2 = plotter.GetHistograms(calc_EOP,
#                                                     list_selections = selections,
#                                                     bins = 50,
#                                                      range_low = -1,
#                                                     range_high = 5,
#                                                     xlabel ="E/p",
#                                                      ylabel = "Number of Tracks",
#                                                     )
#           description = ["MIP Selection",\
#                          "P[GeV]>3",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistPGreater2,
#                                       plotter.channelLabels,
#                                       MCKey='PythiaJetJet',
#                                       DataKey='LowMuData',
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP3Eta0_6_reweighted_MIP.png")
#           CloseCanvas(DataVsMCEOP)


#           ################################################################################yy
#           from selections.selections import sel_IDEta19_23, sel_IDEta0_6, sel_PBetween12_18, sel_PBetween22_28, sel_PBetween28_36
#           selections = [sel_PBetween12_18, sel_IDEta0_6] + MIP_selection
#           trkEOPHistFig2a = plotter.GetHistograms(calc_EOP,
#                                                     list_selections = selections,
#                                                     bins = 50,
#                                                      range_low = -1,
#                                                     range_high = 5,
#                                                     xlabel ="E/p",
#                                                      ylabel = "Number of Tracks",
#                                                     )
#           description = ["MIP Selection",\
#                          "1.2<P[GeV]<1.8",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2a,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP12_18_Eta0_6_reweighted_MIP.png")

#           ################################################################################yy
#           from selections.selections import sel_PBetween22_28
#           selections = [sel_PBetween22_28, sel_IDEta0_6] + MIP_selection
#           trkEOPHistFig2b = plotter.GetHistograms(calc_EOP,\
#                                                 list_selections = selections,\
#                                                 bins = 50,\
#                                                 range_low = -0.75,\
#                                                 range_high = 4,\
#                                                 xlabel ="E/p",\
#                                                 ylabel = "Number of Tracks",\
#                                                )
#           description = ["MIP Selection",\
#                          "2.2<P[GeV]<2.8",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2b,\
#                                       plotter.channelLabels,\
#                                       MCKey='PythiaJetJet',\
#                                       DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                       extra_description = description)
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP22_28_Eta0_6_reweighted_MIP.png")

#           ################################################################################yy
#           selections = [sel_PBetween28_36, sel_IDEta19_23] + MIP_selection
#           trkEOPHistFig2c = plotter.GetHistograms(calc_EOP,\
#                                                 list_selections = selections,\
#                                                 bins = 50,\
#                                                 range_low = -0.75,\
#                                                 range_high = 4,\
#                                                 xlabel ="E/p",\
#                                                 ylabel = "Number of Tracks",\
#                                                 )
#           description = ["MIP Selection",\
#                          "2.8<P[GeV]<3.6",\
#                          "1.9<|#eta_{ID}|<2.3",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2c,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description\
#                                      )
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP28_36_Eta19_23_reweighted_MIP.png")

#           ################################################################################yy
#           selections = [sel_PBetween28_36, sel_IDEta0_6] + MIP_selection
#           trkEOPHistFig2c = plotter.GetHistograms(calc_EOP,\
#                                                 list_selections = selections,\
#                                                 bins = 50,\
#                                                 range_low = -0.5,\
#                                                 range_high = 3,\
#                                                 xlabel ="E/p",\
#                                                 ylabel = "Number of Tracks",\
#                                                 )
#           description = ["MIP Selection",\
#                          "2.8<P[GeV]<3.6",\
#                          "|#eta_{ID}|<0.6",\
#                          "P_{T} Reweighted"]
#           DataVsMCEOP = DrawDataVsMC(trkEOPHistFig2c,\
#                                      plotter.channelLabels,\
#                                      MCKey='PythiaJetJet',\
#                                      DataKey='LowMuData',\
#                                      ratio_min = 0.5,\
#                                      ratio_max = 1.5,\
#                                      extra_description = description\
#                                      )
#           DataVsMCEOP.Draw()
#           DataVsMCEOP.Print(plotter_directory+"/EOPP28_36_Eta0_6_reweighted_MIP.png")

#   if True:

#           ################################################################################yy
#           ## Look in different bins of pseudorapidity
#           from variables.variables import calc_trkEta_ABS
#           base_description = []
#           etaSelections = [sel_IDEta0_6,\
#                            sel_IDEta06_11,\
#                            sel_IDEta11_14,\
#                            sel_IDEta14_15,\
#                            sel_IDEta15_18,\
#                            sel_IDEta18_23]

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
#               selections = [etaSelection] + MIP_selection
#               trkMultiplicity_Eta = plotter.GetHistograms(calc_trkPt,\
#                                                         list_selections = selections,\
#                                                         bins = 80,\
#                                                         range_low = 0.5,\
#                                                         range_high = 5,\
#                                                         xlabel ="Track P_{T} [GeV]",\
#                                                         ylabel = "Number of Tracks",\
#                                                         normalize = False,\
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
#           etaSelections = [sel_IDEta0_6,\
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
#                                                         normalize = False,\
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
#           etaSelections = [sel_IDEta0_6,\
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
#                                                         normalize = False,\
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

#   ################################################################################
#   ##Create a set of p and eta bins for the measurement ##########################
#   from calculation.calculation import calculation
#   from selections.selections import EtaBin, PBin

#   #prepare the momentum bins
#   binMax = 10.05
#   binLow = 0.5
#   nBins = 10
#   base = (binMax/binLow) ** (1./float(nBins))
#   p_bins = []
#   min_p = []
#   for i in range(0, nBins + 1):
#       p_bins.append(0.5 * (base) ** i )

#   #These are the pbins that we will use for the measurement
#   p_ranges = [(x, y) for x, y in zip(p_bins[0:-2], p_bins[1:-1])]
#   p_descriptors = []
#   p_binSelection = []
#   for p_range in p_ranges:
#       PBinFunction = lambda x: PBin(x, p_range[0], p_range[1])
#       sel_PBin = calculation(PBinFunction, ["trk_p"])
#       p_binSelection.append(sel_PBin)
#       p_descriptors.append('{0:.2f}'.format(p_range[0]) + '<P[GeV]<' + '{0:.2f}'.format(p_range[1]))

#   #create a set of strings that could describe the eta or momentum selections
#   eta_bins = [0.0, 0.6, 1.1, 1.7, 2.3]
#   eta_ranges = [(0.0, 0.6),(0.6,1.1),(1.1, 1.7),(1.7,2.3)]
#   eta_descriptors = []
#   eta_binSelections = []
#   for eta_range in eta_ranges:
#       EtaBinFunction = lambda x: EtaBin(x, eta_range[0], eta_range[1])
#       sel_EtaBin = calculation(EtaBinFunction, ["trk_etaID"])
#       eta_binSelections.append(sel_EtaBin)
#       eta_descriptors.append('{0:.1f}'.format(eta_range[0]) + "<|#eta_{ID}|<" + '{0:.1f}'.format(eta_range[1]))

#   #go and get the average E/P for MIP particles in each of the eta bins.
#   for eta_range, eta_descriptor, eta_binSelection in zip(eta_ranges, eta_descriptors, eta_binSelections):
#       selections = MIP_selection + [eta_binSelection]
#       AverageEOP  =  plotter.GetTProfileHistograms(calc_trkP,\
#                                                  calc_EOP,\
#                                                 list_selections = selections,\
#                                                 bins = p_bins,\
#                                                 xlabel ="P[GeV]",\
#                                                 ylabel = "<E/p>",\
#                                                 normalize = False,\
#                                                 )
#       description = ["MIP Selection"] + [eta_descriptor] + ["P_{T} Reweighted"]
#       trkEOP_canvas = DrawDataVsMC(     AverageEOP,\
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
#       trkEOP_canvas.Print(plotter_directory+"/EOPProfileVsMomentumInEtaBin" + str(eta_range[0]) + "_" + str(eta_range[1]) + ".png")

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
#                                                 normalize = False,\
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
#                                                 normalize = False,\
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
