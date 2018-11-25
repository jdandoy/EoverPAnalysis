from PlottingTools.HistogramManager import HistogramManager
import ROOT
from PlottingTools.Plotter import *
import os

def CloseCanvas(canv):
    canv.Close()
    ROOT.gSystem.ProcessEvents()
    del canv

def CreatePlots(filename):
    HM = HistogramManager(filename)
    print HM.channels
    print HM.histograms

    base_description = []
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
    DataVsMC = DrawDataVsMC(hist,\
                            channelLabels,\
                            MCKey='PythiaJetJet',\
                            DataKey='LowMuData',\
                            extra_description = description)
    DataVsMC.Draw()
    DataVsMC.Print(plotter_directory + "/" + histogramName + ".png")
    CloseCanvas(DataVsMC)

    histogramName = "trkCount"
    hist = HM.getHistograms(histogramName)
    description = base_description + ["Inclusive Selection"]
    DataVsMC = DrawDataVsMC(hist,\
                            channelLabels,\
                            MCKey='PythiaJetJet',\
                            DataKey='LowMuData',\
                            extra_description = description)
    DataVsMC.Draw()
    DataVsMC.Print(plotter_directory + "/" + histogramName + ".png")
    CloseCanvas(DataVsMC)

    histogramName = "trkNPV2"
    hist = HM.getHistograms(histogramName)
    description = base_description + ["Inclusive Selection"]
    DataVsMC = DrawDataVsMC(hist,\
                            channelLabels,\
                            MCKey='PythiaJetJet',\
                            DataKey='LowMuData',\
                            extra_description = description)
    DataVsMC.Draw()
    DataVsMC.Print(plotter_directory + "/" + histogramName + ".png")
    CloseCanvas(DataVsMC)

    histogramName = "trkPtHist"
    hist = HM.getHistograms(histogramName)
    description = base_description + ["Inclusive Selection"]
    DataVsMC = DrawDataVsMC(hist,\
                            channelLabels,\
                            MCKey='PythiaJetJet',\
                            DataKey='LowMuData',\
                            doLogx = True,\
                            doLogy = True,\
                            ratio_min = 0.8,\
                            ratio_max = 1.2,\
                            extra_description = description)
    DataVsMC.Draw()
    DataVsMC.Print(plotter_directory + "/" + histogramName + ".png")
    CloseCanvas(DataVsMC)

    histogramName =  "eventNPV2Hist"
    hist = HM.getHistograms(histogramName)
    description = base_description + ["Inclusive Selection"]
    DataVsMC = DrawDataVsMC(hist,\
                            channelLabels,\
                            MCKey='PythiaJetJet',\
                            DataKey='LowMuData',\
                            extra_description = description)
    DataVsMC.Draw()
    DataVsMC.Print(plotter_directory + "/" + histogramName + ".png")
    CloseCanvas(DataVsMC)

    histogramName =  "eventAverageMu"
    hist = HM.getHistograms(histogramName)
    description = base_description + ["Inclusive Selection"]
    DataVsMC = DrawDataVsMC(hist,\
                            channelLabels,\
                            MCKey='PythiaJetJet',\
                            DataKey='LowMuData',\
                            extra_description = description)
    DataVsMC.Draw()
    DataVsMC.Print(plotter_directory + "/" + histogramName + ".png")
    CloseCanvas(DataVsMC)
