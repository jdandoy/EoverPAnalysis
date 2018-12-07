from PlottingTools.HistogramManager import HistogramManager
import ROOT
from PlottingTools.Plotter import *
import os
import time

def CloseCanvas(canv):
    canv.Close()
    ROOT.gSystem.ProcessEvents()
    del canv

filename = "SunDec3.root"

HM = HistogramManager(filename)
HM.listHistograms()

base_description = ["P_{T} Reweighted"]
channelLabels = {"PythiaJetJet" : "Pythia8 MinBias and Dijet", "LowMuData": "2017 Data, Low-<#mu> Run 341294"}
plotter_directory = (filename.split("/")[-1]).replace(".root","") + "plots"

if not os.path.exists("Plots"):
    os.makedirs("Plots")

if not os.path.exists("Plots/" + plotter_directory):
    os.makedirs("Plots/" + plotter_directory)

plotter_directory = "Plots/" + plotter_directory

histogramName = "trkAverageMu"
hist = HM.getHistograms(histogramName)
description = base_description + ["Inclusive Selection"]
DataVsMC1 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        extra_description = description)
DataVsMC1.Draw()
DataVsMC1.Print(plotter_directory + "/" + histogramName + ".png")

histogramName = "trkCount"
hist = HM.getHistograms(histogramName)
description = base_description + ["Inclusive Selection"]
DataVsMC2 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        extra_description = description)
DataVsMC2.Draw()
DataVsMC2.Print(plotter_directory + "/" + histogramName + ".png")

histogramName = "trkNPV2"
hist = HM.getHistograms(histogramName)
description = base_description + ["Inclusive Selection"]
DataVsMC3 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        extra_description = description)
DataVsMC3.Draw()
DataVsMC3.Print(plotter_directory + "/" + histogramName + ".png")

histogramName = "trkPtHist"
hist = HM.getHistograms(histogramName)
description = base_description + ["Inclusive Selection"]
DataVsMC4 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx = True,\
                        doLogy = True,\
                        ratio_min = 0.8,\
                        ratio_max = 1.2,\
                        extra_description = description)
DataVsMC4.Draw()
DataVsMC4.Print(plotter_directory + "/" + histogramName + ".png")

histogramName =  "eventNPV2Hist"
hist = HM.getHistograms(histogramName)
description = base_description + ["Inclusive Selection"]
DataVsMC5 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        extra_description = description)
DataVsMC5.Draw()
DataVsMC5.Print(plotter_directory + "/" + histogramName + ".png")

histogramName =  "eventAverageMu"
hist = HM.getHistograms(histogramName)
description = base_description + ["Inclusive Selection"]
DataVsMC6 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        extra_description = description)
DataVsMC6.Draw()
DataVsMC6.Print(plotter_directory + "/" + histogramName + ".png")

histogramName =  "InclusiveEOP"
hist = HM.getHistograms(histogramName)
description = base_description + ["Inclusive Selection"]
DataVsMC7 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        extra_description = description)
DataVsMC7.Draw()
DataVsMC7.Print(plotter_directory + "/" + histogramName + ".png")


histogramName = "TrkPtHisteta06"
hist = HM.getHistograms(histogramName)
description = base_description + ["|#eta_{ID}|<0.6"]
DataVsMC8 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC8.Draw()
DataVsMC8.Print(plotter_directory + "/" + histogramName + ".png")

histogramName = "TrkPtHisteta06_11"
hist = HM.getHistograms(histogramName)
description = base_description + ["0.6<|#eta_{ID}|<1.1"]
DataVsMC8 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC8.Draw()
DataVsMC8.Print(plotter_directory + "/" + histogramName + ".png")

histogramName = "TrkPtHisteta11_14"
hist = HM.getHistograms(histogramName)
description = base_description + ["1.1<|#eta_{ID}|<1.4"]
DataVsMC8 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC8.Draw()
DataVsMC8.Print(plotter_directory + "/" + histogramName + ".png")

histogramName = "TrkPtHisteta14_15"
hist = HM.getHistograms(histogramName)
description = base_description + ["1.4<|#eta_{ID}|<1.5"]
DataVsMC8 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC8.Draw()
DataVsMC8.Print(plotter_directory + "/" + histogramName + ".png")

histogramName =  "TrkPtHisteta15_18"
hist = HM.getHistograms(histogramName)
description = base_description + ["1.5<|#eta_{ID}|<1.8"]
DataVsMC9 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC9.Draw()
DataVsMC9.Print(plotter_directory + "/" + histogramName + ".png")

histogramName =  "TrkPtHisteta18_23"
hist = HM.getHistograms(histogramName)
description = base_description + ["1.8<|#eta_{ID}|<2.3"]
DataVsMC9 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC9.Draw()
DataVsMC9.Print(plotter_directory + "/" + histogramName + ".png")

histogramName_num =  "InclusiveZeroFractionVsPNumerator"
hist_num = HM.getHistograms(histogramName_num)

histogramName_den = "InclusiveZeroFractionVsPDenomenator"
hist_den = HM.getHistograms(histogramName_den)
ratio_hist = DivideHistograms(hist_num, hist_den)

description = base_description + ["Inclusive Selection"]
DataVsMC9 = DrawDataVsMC(ratio_hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        doLogy=False,
                        ratio_min = 0.8,\
                        ratio_max = 1.2,\
                        extra_description = description)
DataVsMC9.Draw()
DataVsMC9.Print(plotter_directory + "/" + histogramName_num.replace("Numerator", "") + ".png")
description = base_description + ["Inclusive Selection", "0.2<|#eta_{ID}|<0.4"]
histogramName_num =  "ZeroFractionVsPetaID02_04Numerator"
hist_num = HM.getHistograms(histogramName_num)

histogramName_den = "ZeroFractionVsPetaID02_04Denomenator"
hist_den = HM.getHistograms(histogramName_den)
ratio_hist = DivideHistograms(hist_num, hist_den)

DataVsMC9 = DrawDataVsMC(ratio_hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        doLogy=False,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC9.Draw()
DataVsMC9.Print(plotter_directory + "/" + histogramName_num.replace("Numerator", "") + ".png")

histogramName_num =  "ZeroFractionVsPetaID00_02Numerator"
hist_num = HM.getHistograms(histogramName_num)

histogramName_den = "ZeroFractionVsPetaID00_02Denomenator"
hist_den = HM.getHistograms(histogramName_den)
ratio_hist = DivideHistograms(hist_num, hist_den)

DataVsMC9 = DrawDataVsMC(ratio_hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        doLogy=False,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC9.Draw()
DataVsMC9.Print(plotter_directory + "/" + histogramName_num.replace("Numerator", "") + ".png")

description = base_description + ["0.4<|#eta_{ID}|<0.6"]
histogramName_num =  "ZeroFractionVsPetaID04_06Numerator"
hist_num = HM.getHistograms(histogramName_num)

histogramName_den = "ZeroFractionVsPetaID04_06Denomenator"
hist_den = HM.getHistograms(histogramName_den)
ratio_hist = DivideHistograms(hist_num, hist_den)

DataVsMC9 = DrawDataVsMC(ratio_hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        doLogy=False,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC9.Draw()
DataVsMC9.Print(plotter_directory + "/" + histogramName_num.replace("Numerator", "") + ".png")

histogramName = "EtaLess08_TwoDHistTrkPvsPhiInnerToExtrapolEM2"
description = base_description + ["|#eta_{ID}|<0.8"]
hist = HM.getHistograms(histogramName)
DataVsMC10 = Draw2DHistogramOnCanvas(hist["PythiaJetJet"], doLogx = False, doLogy = True)
DataVsMC10.Draw()
DataVsMC10.Print(plotter_directory + "/" + histogramName.replace("Numerator", "") + "PythiaJetJet" + ".png")

histogramName = "EtaLess08_TwoDHistTrkPvsPhiInnerToExtrapolEM2"
description = base_description + ["|#eta_{ID}|<0.8"]
hist = HM.getHistograms(histogramName)
DataVsMC10 = Draw2DHistogramOnCanvas(hist["LowMuData"], doLogx = False, doLogy = True)
DataVsMC10.Draw()
DataVsMC10.Print(plotter_directory + "/" + histogramName.replace("Numerator", "") + "LowMuData" + ".png")

histogramName = "EtaID0_6_PGreater2_0_EOPHist"
hist = HM.getHistograms(histogramName)
DataVsMC9 = DrawDataVsMC(hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        doLogy=False,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC9.Draw()
DataVsMC9.Print(plotter_directory + "/" + histogramName_num.replace("Numerator", "") + ".png")


description = base_description + ["0.0<|#eta_{ID}|<0.6"]
histogramName_num =  "ZeroFractionVsPetaID00_06Numerator"
hist_num = HM.getHistograms(histogramName_num)

histogramName_den = "ZeroFractionVsPetaID00_06Denomenator"
hist_den = HM.getHistograms(histogramName_den)
ratio_hist = DivideHistograms(hist_num, hist_den)

DataVsMC9 = DrawDataVsMC(ratio_hist,\
                        channelLabels,\
                        MCKey='PythiaJetJet',\
                        DataKey='LowMuData',\
                        doLogx=True,
                        doLogy=False,
                        ratio_min = 0.6,\
                        ratio_max = 1.4,\
                        extra_description = description)
DataVsMC9.Draw()
DataVsMC9.Print(plotter_directory + "/" + histogramName_num.replace("Numerator", "") + ".png")

description = base_description + ["0.0<|#eta_{ID}|<0.4", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_0_4"
histograms = HM.getHistograms(histogramName)
for key in histograms:
    histograms[key] = histograms[key].ProjectionX(histogramName + "_px", "E")

DataVSMC10 = DrawDataVsMC(histograms,\
                          channelLabels,\
                          MCKey = "PythiaJetJet",\
                          DataKey = "LowMuData",\
                          doLogx = True,\
                          doLogy = False,
                          ratio_min = 0.9,\
                          ratio_max = 1.1,\
                          extra_description = description)
DataVSMC10.Draw()
DataVSMC10.Print(plotter_directory + "/" + histogramName + ".png")

description = base_description + ["0.4<|#eta_{ID}|<0.8", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_4_8"
histograms = HM.getHistograms(histogramName)
for key in histograms:
    histograms[key] = histograms[key].ProjectionX(histogramName + "_px", "E")

DataVSMC10 = DrawDataVsMC(histograms,\
                          channelLabels,\
                          MCKey = "PythiaJetJet",\
                          DataKey = "LowMuData",\
                          doLogx = True,\
                          doLogy = False,
                          ratio_min = 0.9,\
                          ratio_max = 1.1,\
                          extra_description = description)
DataVSMC10.Draw()
DataVSMC10.Print(plotter_directory + "/" + histogramName + ".png")


description = base_description + ["0.8<|#eta_{ID}|<1.2", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_8_12"
histograms = HM.getHistograms(histogramName)
for key in histograms:
    histograms[key] = histograms[key].ProjectionX(histogramName + "_px", "E")

DataVSMC10 = DrawDataVsMC(histograms,\
                          channelLabels,\
                          MCKey = "PythiaJetJet",\
                          DataKey = "LowMuData",\
                          doLogx = True,\
                          doLogy = False,
                          ratio_min = 0.9,\
                          ratio_max = 1.1,\
                          extra_description = description)
DataVSMC10.Draw()
DataVSMC10.Print(plotter_directory + "/" + histogramName + ".png")


description = base_description + ["1.2<|#eta_{ID}|<1.6", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_12_16"
histograms = HM.getHistograms(histogramName)
for key in histograms:
    histograms[key] = histograms[key].ProjectionX(histogramName + "_px", "E")

DataVSMC10 = DrawDataVsMC(histograms,\
                          channelLabels,\
                          MCKey = "PythiaJetJet",\
                          DataKey = "LowMuData",\
                          doLogx = True,\
                          doLogy = False,
                          ratio_min = 0.9,\
                          ratio_max = 1.1,\
                          extra_description = description)
DataVSMC10.Draw()
DataVSMC10.Print(plotter_directory + "/" + histogramName + ".png")

description = base_description + ["1.6<|#eta_{ID}|<2.0", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_16_20"
histograms = HM.getHistograms(histogramName)
for key in histograms:
    histograms[key] = histograms[key].ProjectionX(histogramName + "_px", "E")

DataVSMC10 = DrawDataVsMC(histograms,\
                          channelLabels,\
                          MCKey = "PythiaJetJet",\
                          DataKey = "LowMuData",\
                          doLogx = True,\
                          doLogy = False,
                          ratio_min = 0.9,\
                          ratio_max = 1.1,\
                          extra_description = description)
DataVSMC10.Draw()
DataVSMC10.Print(plotter_directory + "/" + histogramName + ".png")

description = base_description + ["2.0<|#eta_{ID}|<2.4", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_20_24"
histograms = HM.getHistograms(histogramName)
for key in histograms:
    histograms[key] = histograms[key].ProjectionX(histogramName + "_px", "E")

DataVSMC10 = DrawDataVsMC(histograms,\
                          channelLabels,\
                          MCKey = "PythiaJetJet",\
                          DataKey = "LowMuData",\
                          doLogx = True,\
                          doLogy = False,
                          ratio_min = 0.9,\
                          ratio_max = 1.1,\
                          extra_description = description)
DataVSMC10.Draw()
DataVSMC10.Print(plotter_directory + "/" + histogramName + ".png")

description = base_description + ["0.0<|#eta_{ID}|<0.4", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "2DHist_EOPVsMomentum_MIPSelection_HadFracAbove70_InBin_0_4"
histograms = HM.getHistograms(histogramName)

MCCanvas = Draw2DHistogramOnCanvas(histograms["PythiaJetJet"], doLogx = True, doLogy = False)
MCCanvas.Draw()
MCCanvas.Print(plotter_directory + "/" + histogramName + "PythiaJetJet.png")

DataCanvas = Draw2DHistogramOnCanvas(histograms["LowMuData"], doLogx = True, doLogy = False)
DataCanvas.Draw()
DataCanvas.Print(plotter_directory + "/" + histogramName + "LowMuData.png")

description = base_description + ["0.4<|#eta_{ID}|<0.8", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "2DHist_EOPVsMomentum_MIPSelection_HadFracAbove70_InBin_4_8"
histograms = HM.getHistograms(histogramName)

MCCanvas = Draw2DHistogramOnCanvas(histograms["PythiaJetJet"], doLogx = True, doLogy = False)
MCCanvas.Draw()
MCCanvas.Print(plotter_directory + "/" + histogramName + "PythiaJetJet.png")

DataCanvas = Draw2DHistogramOnCanvas(histograms["LowMuData"], doLogx = True, doLogy = False)
DataCanvas.Draw()
DataCanvas.Print(plotter_directory + "/" + histogramName + "LowMuData.png")

description = base_description + ["0.8<|#eta_{ID}|<1.2", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "2DHist_EOPVsMomentum_MIPSelection_HadFracAbove70_InBin_8_12"
histograms = HM.getHistograms(histogramName)

MCCanvas = Draw2DHistogramOnCanvas(histograms["PythiaJetJet"], doLogx = True, doLogy = False)
MCCanvas.Draw()
MCCanvas.Print(plotter_directory + "/" + histogramName + "PythiaJetJet.png")

DataCanvas = Draw2DHistogramOnCanvas(histograms["LowMuData"], doLogx = True, doLogy = False)
DataCanvas.Draw()
DataCanvas.Print(plotter_directory + "/" + histogramName + "LowMuData.png")

description = base_description + ["1.2<|#eta_{ID}|<1.6", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "2DHist_EOPVsMomentum_MIPSelection_HadFracAbove70_InBin_12_16"
histograms = HM.getHistograms(histogramName)

MCCanvas = Draw2DHistogramOnCanvas(histograms["PythiaJetJet"], doLogx = True, doLogy = False)
MCCanvas.Draw()
MCCanvas.Print(plotter_directory + "/" + histogramName + "PythiaJetJet.png")

DataCanvas = Draw2DHistogramOnCanvas(histograms["LowMuData"], doLogx = True, doLogy = False)
DataCanvas.Draw()
DataCanvas.Print(plotter_directory + "/" + histogramName + "LowMuData.png")

description = base_description + ["1.6<|#eta_{ID}|<2.0", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "2DHist_EOPVsMomentum_MIPSelection_HadFracAbove70_InBin_16_20"
histograms = HM.getHistograms(histogramName)

MCCanvas = Draw2DHistogramOnCanvas(histograms["PythiaJetJet"], doLogx = True, doLogy = False)
MCCanvas.Draw()
MCCanvas.Print(plotter_directory + "/" + histogramName + "PythiaJetJet.png")

DataCanvas = Draw2DHistogramOnCanvas(histograms["LowMuData"], doLogx = True, doLogy = False)
DataCanvas.Draw()
DataCanvas.Print(plotter_directory + "/" + histogramName + "LowMuData.png")

description = base_description + ["2.0<|#eta_{ID}|<2.4", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
histogramName = "2DHist_EOPVsMomentum_MIPSelection_HadFracAbove70_InBin_20_24"
histograms = HM.getHistograms(histogramName)

MCCanvas = Draw2DHistogramOnCanvas(histograms["PythiaJetJet"], doLogx = True, doLogy = False)
MCCanvas.Draw()
MCCanvas.Print(plotter_directory + "/" + histogramName + "PythiaJetJet.png")

DataCanvas = Draw2DHistogramOnCanvas(histograms["LowMuData"], doLogx = True, doLogy = False)
DataCanvas.Draw()
DataCanvas.Print(plotter_directory + "/" + histogramName + "LowMuData.png")

histogramNames = ["TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_0_4" , "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_12_16", "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_16_20", "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_20_24", "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_4_8", "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_8_12"]

descriptors = ["0.0<|#eta_{ID}|<0.4", "1.2<|#eta_{ID}|<1.6", "1.6<|#eta_{ID}|<2.0", "2.0<|#eta_{ID}|<2.4", "0.4<|#eta_{ID}|<0.8", "0.8<|#eta_{ID}|<1.2"]

for histogramName, descriptor in zip(histogramNames, descriptors):
    description = base_description + [descriptor, "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
    histograms = HM.getHistograms(histogramName)

    for key in histograms:
        histograms[key].Rebin(20)

    DataVSMC10 = DrawDataVsMC(histograms,\
                          channelLabels,\
                          MCKey = "PythiaJetJet",\
                          DataKey = "LowMuData",\
                          doLogx = True,\
                          doLogy = False,
                          ratio_min = 0.6,\
                          ratio_max = 1.4,\
                          extra_description = description)
    DataVSMC10.Draw()
    DataVSMC10.Print(plotter_directory + "/" + histogramName + ".png")

raw_input()
#['NTRT20ZeroFractionVsPetaID00_06Denomenator', 'TrkPtHisteta15_18', 'trkNPV2', 'eventAverageMu', 'trkCount', 'trkAverageMu', 'NTRT20ZeroFractionVsPetaID06_11Denomenator', 'NTRT20ZeroFractionVsPetaID02_04Numerator', 'EtaID0_6_PGreater2_0_EOPHist', 'NonZero_EtaID0_6_PBetween12_18_EOPHist', 'NTRT20ZeroFractionVsPetaID00_02Denomenator', 'NonZero_EtaIDBetween19_23_PBetween22_28_EOPHist', 'EtaID0_6_PGreater1_5_EOPHist', 'EtaLess08_TwoDHistTrkPvsPhiInnerToExtrapolEM2', 'InclusiveZeroFractionVsPNumerator', 'TrkPtHisteta06', 'NTRT20ZeroFractionVsPetaID11_14Numerator', 'NTRT20ZeroFractionVsPetaID00_06Numerator', 'NonZeroEnergy_InclusiveEOP', 'TrkPtHisteta18_23', 'InclusiveEOP', 'ZeroFractionVsPetaID06_11Denomenator', 'ZeroFractionVsPetaID02_04Denomenator', 'NTRT20ZeroFractionVsPetaID04_06Numerator', 'trkPtHist', 'NTRT20ZeroFractionVsPetaID00_02Numerator', 'TrkPtHisteta06_11', 'InclusiveZeroFractionVsAbsEtaNumerator', 'ZeroFractionVsPetaID04_06Denomenator', 'NearestDRHist', 'TrackEtaID', 'NonZero_EtaID0_6_PBetween22_28_EOPHist', 'NTRT20ZeroFractionVsPetaID06_11Numerator', 'InclusiveZeroFractionVsEtaNumerator', 'InclusiveZeroFractionVsPDenomenator', 'TwoDTrackPtVsEtaHistogram', 'ZeroFractionVsPetaID11_14Denomenator', 'ZeroFractionVsPetaID00_02Denomenator', 'eventNPV2Hist', 'trkEtaECALHist', 'TrkPtHisteta14_15', 'InclusiveZeroFractionVsEtaDenomenator', 'EtaID0_6_PBetween22_28_EOPHist', 'EtaID0_6_PGreater1_0_EOPHist', 'EtaIDBetween19_23_PBetween28_36_EOPHist', 'NonZero_EtaID0_6_PBetween28_36_EOPHist', 'TrkPtHisteta11_14', 'lowPTLess07_TwoDHistTrkEtavsDEtaInnerToExtrapolEM2', 'NTRT20ZeroFractionVsPetaID04_06Denomenator', 'ZeroFractionVsPetaID00_06Denomenator', 'TwoDHistTrkPvsPhiInnerToExtrapolEM2', 'EtaID0_6_PBetween12_18_EOPHist', 'NTRT20ZeroFractionVsPetaID11_14Denomenator', 'InclusiveZeroFractionVsAbsEtaDenomenator', 'TwoDTrackPvsTrkEtaID', 'NTRT20ZeroFractionVsPetaID02_04Denomenator', 'EtaID0_6_PBetween28_36_EOPHist']



