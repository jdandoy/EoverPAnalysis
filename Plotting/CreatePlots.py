from PlottingTools.HistogramManager import HistogramManager
import ROOT
from PlottingTools.Plotter import *
from array import array
import os
import time

ROOT.gROOT.SetBatch(ROOT.kTRUE)

def CloseCanvas(canv):
    canv.Close()
    ROOT.gSystem.ProcessEvents()
    del canv

filename = "PTSpectrumReweighted.root"
#filename = "CountReweightedPlots.root"

HM = HistogramManager(filename)
HM.listHistograms()

if "Count" in filename:
    base_description = []
if "Spectrum" in filename:
    base_description = ["P_{T} Reweighted"]
channelLabels = {"SinglePion": "Single Pion", "PythiaJetJet" : "Pythia8 MinBias and Dijet", "LowMuData": "2017 Low-<#mu> Data", "PythiaJetJetPionsReweighted":"Pythia8 MB+DJ Pions Only"}
plotter_directory = (filename.split("/")[-1]).replace(".root","") + "plots"

if not os.path.exists("Plots"):
    os.makedirs("Plots")

if not os.path.exists("Plots/" + plotter_directory):
    os.makedirs("Plots/" + plotter_directory)

plotter_directory = "Plots/" + plotter_directory

def fitHistograms(histograms, fit_function, histogramName, channels=[], eta_low=-1,eta_high=-1,p_low=-1,p_high=-1, refit=False):
    '''
    Fit function fit_function to the histograms
    fit_function can be gaus, landau or convolution
    '''

    low_value = 0.0
    max_x = {}
    max_x_err = {}
    chisq = {}

    low_fit = 1000.0
    high_fit = -1000.0
    low_rms = 0.0
    high_rms = 1.2
    means = {}
    sigmas = {}
    if not refit:
       for channel in channels:
           histogram = histograms[channel]
           low_bin = histogram.FindBin(low_rms)
           high_bin = histogram.FindBin(high_rms)
           histogram.GetXaxis().SetRange(low_bin, high_bin)
           mean = histogram.GetMean()
           rms = histogram.GetRMS()
           sigmas[channel]=rms
           #gaus fit:
           #new_low_fit = mean - rms * 1.0
           #new_high_fit = mean + rms * 1.0
           #convolution fit:
           new_low_fit = mean - rms * 1.4
           new_high_fit = mean + rms * 1.2
           means[channel]=mean
           if new_low_fit < low_fit:
               low_fit = new_low_fit
           if new_high_fit > high_fit:
               high_fit = new_high_fit
    if refit:
       histogram = histograms["LowMuData"]
       histogram_name = histogram.GetName()
       fit_function_string = fit_function + histogram_name
       fit = histogram.GetFunction(fit_function_string)
       low_fit = fit.GetParameter(0) - fit.GetParameter(1)
       high_fit = fit.GetParameter(0) + fit.GetParameter(1)
       for channel in channels:
           means[channel] = fit.GetParameter(0)
           sigmas[channel] = fit.GetParameter(1)

    for channel in channels:
       histogram = histograms[channel]

       histogram_name = histogram.GetName()
       #this is the landau distribution that will be fit to the histogramsograms
       landau = ROOT.TF1("landau_" + channel +histogram_name, "[2]*TMath::Landau(x, [0], [1])", -1.0, 5.0)
       landau.SetParName(0, "mpv")
       landau.SetParameter(0, means[channel])
       landau.SetParLimits(0, 0.3, 1.1)
       landau.SetParName(1, "sigma")
       landau.SetName("landau" +  histogram_name)
       landau.SetParameter(1, sigmas[channel]/4.0)
       landau.SetParLimits(1, sigmas[channel]/100.0, sigmas[channel]*2.0)
       landau.SetParName(2, "Norm")
       landau.SetParameter(2, histogram.Integral())

       #this is the gaus distribution that will be fit to the histogramsograms
       gaus = ROOT.TF1("gaus_" + channel +histogram_name, "[2]*TMath::Gaus(x, [0], [1])", -1.0, 5.0)
       gaus.SetParName(0, "mu")
       gaus.SetParameter(0, means[channel])
       gaus.SetParLimits(0, 0.45, 0.95)
       gaus.SetParName(1, "sigma")
       gaus.SetName("gaus" + histogram_name)
       gaus.SetParameter(1, sigmas[channel])
       gaus.SetParLimits(1, sigmas[channel]/3.0, 1.2*sigmas[channel])
       gaus.SetParName(2, "Norm")
       gaus.SetParameter(2, histogram.Integral())

       #Create a gaussian convoluted with a landau histogramsogram
       gaus_forconvolution = ROOT.TF1("gaus_forconvolution_" + channel +histogram_name, "TMath::Gaus(x, 0.0, [0])", -10.0, +10.0)

       landau_forconvolution = ROOT.TF1("landau_forconvolution_" + channel +histogram_name, "[2]*TMath::Landau(x, [0], [1])", -1.0, 5.0)

       convolution = ROOT.TF1Convolution(gaus_forconvolution, landau_forconvolution,-1,6,True)
       convolution.SetRange(-1.,5.)
       convolution.SetNofPointsFFT(10000)
       convolution_tofit = ROOT.TF1("f",convolution, -1.0, 5., convolution.GetNpar())
       convolution_tofit.SetName("convolution" +histogram_name)

       convolution_tofit.SetParName(0, "SigmaSmear")
       convolution_tofit.SetParLimits(0, 0.0, 100.0)
       convolution_tofit.SetParameter(0, 0.5)

       convolution_tofit.SetParName(1, "mpv")
       convolution_tofit.SetParameter(1, means[channel])
       convolution_tofit.SetParLimits(1, 0.3, 1.1)
       convolution_tofit.SetParName(2, "sigma")
       convolution_tofit.SetParameter(2, sigmas[channel]/4.0)
       convolution_tofit.SetParLimits(2, sigmas[channel]/100.0, sigmas[channel]*2.0)
       convolution_tofit.SetParName(3, "Norm")
       convolution_tofit.SetParameter(3, histogram.Integral())
       print("Created convolution function")
       convolution_tofit.Print()

       #Good settings for a landau x gaus
       fit_function_string = fit_function + histogram_name
       print("convolution" +histogram_name)
       print(fit_function_string)
       histogram.GetXaxis().SetRange(histogram.FindBin(0.15), histogram.FindBin(1.3))
       histogram.GetXaxis().SetRange()

       #need to figure out a good range to do the fit in. For now use this, but need to experiment with it
       #low=0.1
       #high=0.8
       #if p_high > 2:
       #    extra = ((28.0 - (p_high-2.0))/28.0) * 0.15
       #    high += extra
       #else:
       #    extra = ((2.0 - p_high)/2.0) * 0.8
       #    high += extra


       histogram.Fit(fit_function_string, "", "", low_fit, high_fit)

       fit = histogram.GetFunction(fit_function_string)
       if  fit:
           fit.SetLineColor(histogram.GetLineColor())
           if fit_function == "convolution":
 #              max_x_tmp, max_val_tmp = FindMostProbableValue(fit, low_fit, high_fit)
               max_x_tmp=0.0
               max_val_tmp=0.0
               max_x_tmp_err = 0.0
               chisq[channel]=fit.GetChisquare()
           elif fit_function == "gaus":
               max_x_tmp = fit.GetParameter(0)
               max_val_tmp = -1.0
               max_x_tmp_err = fit.GetParError(0)
               chisq[channel]=fit.GetChisquare()
       else:
           max_x_tmp = -1.0
           max_val_tmp = -1.0
           max_x_tmp_err = 0.0
           chisq[channel] = -1.0

       max_x[channel]=max_x_tmp
       max_x_err[channel]=max_x_tmp_err

       print "The maximum value in data was at eop " + str(max_x_tmp)

    return max_x, max_x_err, chisq

def FindMostProbableValue(fit_function, low, high, ndivisions=10000):
    x = np.linspace(low, high, ndivisions)
    vals = np.zeros(len(x))
    max_val = 0
    max_x = -999999.0
    for i in range(0, len(x)):
        vals[i] = fit_function.Eval(x[i])
        if vals[i] > max_val:
            max_val = vals[i]
            max_x = x[i]

    return max_x, max_val

def CreateZeroFractionPlotsFromSelection(HM, numerator_selection_name, denomenator_selection_name, filename, base_description=[]):
    #get the binning vectors
    f = ROOT.TFile(filename, "READ")
    tree = f.Get(numerator_selection_name + "BinningTree")
    for bins in tree:
        break
    eta_bins_low = getattr(bins, numerator_selection_name+"EtaBinsLow")
    eta_bins_high = getattr(bins, numerator_selection_name+"EtaBinsHigh")
    p_bins_low_for_eta_bin = []
    p_bins_high_for_eta_bin = []

    #get all of the binning information that we need
    for i in range(0, eta_bins_low.size()):
        p_bins_low_for_eta_bin.append(getattr(bins, numerator_selection_name+"PBinsLow_Eta"+str(i)))
        p_bins_high_for_eta_bin.append(getattr(bins, numerator_selection_name+"PBinsHigh_Eta"+str(i)))

    for i, eta_low, eta_high in zip(list(range(0, eta_bins_low.size())), eta_bins_low, eta_bins_high):
        for histogram in ["TrackPtSpectrum", "TrackPSpectrum"]:
            selection_name = numerator_selection_name
            histogram_name = histogram + "__" + selection_name + "_Eta_" + str(i)
            hist_numerator = HM.getHistograms(histogram_name, rebin = 10)

            selection_name = denomenator_selection_name
            histogram_name = histogram + "__" + selection_name + "_Eta_" + str(i)
            hist_denomenator = HM.getHistograms(histogram_name, rebin = 10)

            #divide the numerator by the denomenator
            ZeroFraction = DivideHistograms(hist_numerator, hist_denomenator, efficiency_error=True)
            description = base_description + [str(round(eta_low, 2)) + " < |#eta| < " + str(round(eta_high, 2))]
            histogram_name = "NonZeroFraction" +numerator_selection_name + denomenator_selection_name + histogram + "_" + str(i)

            DataVsMC1 = DrawDataVsMC(ZeroFraction,\
                                    channelLabels,\
                                    MCKeys = ['PythiaJetJet'],\
                                    #MCKeys = ['PythiaJetJet'],\
                                    DataKey='LowMuData',\
                                    doLogx=True,\
                                    doLogy=False,\
                                    ylabel="N(E!=0)/N",\
                                    ratio_min=0.6,\
                                    ratio_max=1.4,\
                                    extra_description = description)

            DataVsMC1[0].Draw()
            DataVsMC1[0].Print(plotter_directory + "/" + histogram_name + ".png")
            DataVsMC1[0].Close()

def CreatePlotsFromSelection(HM, selection_name, filename, base_description = [], doFit = False, fitfunction="gaus", refit=False):
    #get the binning vectors
    f = ROOT.TFile(filename, "READ")

    tree = f.Get(selection_name + "BinningTree")
    for bins in tree:
        break

    eta_bins_low = getattr(bins, selection_name+"EtaBinsLow")
    eta_bins_high = getattr(bins, selection_name+"EtaBinsHigh")
    p_bins_low_for_eta_bin = []
    p_bins_high_for_eta_bin = []

    #get all of the binning information that we need
    for i in range(0, eta_bins_low.size()):
        p_bins_low_for_eta_bin.append(getattr(bins, selection_name+"PBinsLow_Eta"+str(i)))
        p_bins_high_for_eta_bin.append(getattr(bins, selection_name+"PBinsHigh_Eta"+str(i)))

    histograms_in_eta_bins = ["TrackPtSpectrum",\
                                   "TrackPSpectrum",\
                                   "EOPProfileVsMomentum",\
                                   "EnergyAnulusProfileVsMomentum",\
                                   "EnergyBkgProfileVsMomentum"]

    histograms_in_momentum_bins = ["EOPDistribution",\
                                   #"EOPBkgDistribution",\
                                   #"trkTRTHits",\
                                   #"trkEMDR100",\
                                   #"MomentumHadFrac",\
                                   #"HadFrac",\
                                   #"NClusters",\
                                   #"NClusters_EM",\
                                   #"NClusters_HAD",\
                                   #"NClusters_emlike",\
                                   #"NClusters_hadlike",\
                                   ]

    #OK now lets make all of the plots in all of the bins!

    for i, eta_low, eta_high in zip(list(range(0, eta_bins_low.size())), eta_bins_low, eta_bins_high):
        #the plots binned in eta
        for histogram in histograms_in_eta_bins:
            if "MIP" in selection_name and "Spectrum" in histogram:
                continue
            if "pectrum" in histogram:
                rebin = 10
            else:
                rebin = 1
            xAxis_range = None
            print(eta_low, eta_high)

            histogram_name = histogram + "__" + selection_name + "_Eta_" + str(i)
            hist = HM.getHistograms(histogram_name, rebin = rebin)
            #shouldILogy=True
            to_plot = ['PythiaJetJet']
            if "Profile" in histogram_name and not "EnergyBkg" in histogram_name:
                shouldILogy = False
                ratio_min = 0.9
                ratio_max = 1.1

            else:
                shouldILogy = False
                ratio_min = 0.8
                ratio_max = 1.2

            if "Profile" in histogram_name:
                hist = ProjectProfiles(hist)

            if "EnergyBkgProfileVsMomentum" == hist:
                ratio_min =0.5
                ratio_min =1.5

            description = base_description + [str(round(eta_low, 2)) + " < |#eta| < " + str(round(eta_high, 2))]

            DataVsMC1 = DrawDataVsMC(hist,\
                                    channelLabels,\
                                    MCKeys = to_plot,\
                                    DataKey='LowMuData',\
                                    doLogy=shouldILogy,\
                                    doLogx=True,\
                                    ratio_min=ratio_min,\
                                    ratio_max=ratio_max,\
                                    extra_description = description)

            DataVsMC1[0].Draw()
            DataVsMC1[0].Print(plotter_directory + "/" + histogram_name + ".png")
            DataVsMC1[0].Close()
        p_bins_low = p_bins_low_for_eta_bin[i]
        p_bins_high = p_bins_high_for_eta_bin[i]

        #the plots binned in eta and momentum
        mpv_eop = []
        chisq=[]
        p_centers = []

        chisq_hists = {}
        eop_hists = {}
        p_bins = list(p_bins_low) + [p_bins_high[-1]]
        for channel in HM.channels:
            chisq_hists[channel] =  ROOT.TH1D(channel + "chisq_hist_Eta_" +  str(i),"hist", len(p_bins) -1, array('d', p_bins))
            eop_hists[channel] = ROOT.TH1D(channel + "eop_hist_Eta_" +  str(i),"hist" , len(p_bins)-1, array('d', p_bins))

        for histogram in histograms_in_momentum_bins:

            for j, p_low, p_high in zip(list(range(0, p_bins_low.size())), p_bins_low, p_bins_high):
                histattributes = []
                histogram_name = histogram + "_" + selection_name + "_Eta_" + str(i) + "_Momentum_" + str(j)
                hist = HM.getHistograms(histogram_name)
                description = base_description + [str(round(eta_low, 2)) + " < |#eta| < " + str(round(eta_high, 2))]
                description += [str(round(p_low, 3)) + " < P/GeV < " + str(round(p_high, 3))]

                if doFit and histogram == "EOPDistribution":
                    print(fitHistograms(hist, fitfunction, histogram_name, HM.channels, eta_low=eta_low, eta_high=eta_high, p_low=p_low, p_high=p_high))
                    max_eop, max_eop_err, chisq_fit = fitHistograms(hist, fitfunction, histogram_name, HM.channels, eta_low=eta_low, eta_high=eta_high, p_low=p_low, p_high=p_high)
                    if refit:
                        max_eop, max_eop_err, chisq_fit = fitHistograms(hist, fitfunction, histogram_name, HM.channels, eta_low=eta_low, eta_high=eta_high, p_low=p_low, p_high=p_high, refit=refit)
                    for channel in HM.channels:
                        chisq_hists[channel].SetBinContent(j+1, chisq_fit[channel])
                        chisq_hists[channel].SetBinError(j+1, 0.0)
                        eop_hists[channel].SetBinContent(j+1, max_eop[channel])
                        eop_hists[channel].SetBinError(j+1, max_eop_err[channel])
                    p_centers.append((p_low + p_high)/2.0)

                xAxis_range = None
                if "istribution" in histogram:
                    #rebin = 2
                    xAxis_range = (-0.5, 2.0)

                DataVsMC1 = DrawDataVsMC(hist,\
                                        channelLabels,\
                                        MCKeys = ['PythiaJetJet'],\
                                        #MCKeys = ['PythiaJetJet'],\
                                        xAxis_range=xAxis_range,\
                                        rebin=rebin,\
                                        doLogy=False,\
                                        DataKey='LowMuData',\
                                        extra_description = description)

                DataVsMC1[0].Draw()
                DataVsMC1[0].Print(plotter_directory + "/" + histogram_name + ".png")
                DataVsMC1[0].Close()

        if doFit:
           for hists, descr in zip([chisq_hists, eop_hists], ["chisq", "EOP"]):
              eta_low_str = str(eta_low)
              eta_high_str = str(eta_high)
              hists['LowMuData'].SetMarkerStyle(0)
              description = base_description + [eta_low_str + " < |#eta| < " + eta_high_str]
              DataVsMC = DrawDataVsMC(hists,\
                                    channelLabels,\
                                    MCKeys = ['PythiaJetJet'],\
                                    DataKey = "LowMuData",\
                                    doLogx = True,\
                                    doLogy = False,
                                    ratio_min = 0.8,\
                                    ratio_max = 1.2,\
                                    extra_description = description)

              DataVsMC[0].Draw()
              DataVsMC[0].Print(descr + "".join(base_description) + "DataFitResults" + eta_low_str + "_" + eta_high_str + ".png")
              raw_input("How do the fits look?")


#CreateZeroFractionPlotsFromSelection(HM, "NonZeroEnergy", "Inclusive", filename, base_description= base_description + [])
#CreateZeroFractionPlotsFromSelection(HM, "20TRTHitsNonZeroEnergy", "20TRTHits", filename, base_description= base_description + ["N_{TRT} >= 20"])

#test the plot creation
CreatePlotsFromSelection(HM,"20TRTHitsNonZeroEnergy", filename, base_description = base_description + ["N_{TRT} >= 20", "E_{TOTAL} != 0.0"], doFit = True, fitfunction="convolution", refit=False)
#CreatePlotsFromSelection(HM,"MIPSelectionHadFracAbove70", filename, base_description = base_description + ["MIP Selection"],doFit = True ,fitfunction="gaus", refit=True)
#CreatePlotsFromSelection(HM,"NonZeroEnergy", filename, base_description = base_description + ["E_{TOTAL} != 0.0"],doFit = True , fitfunction="convolution")
#CreatePlotsFromSelection(HM,"Inclusive", filename, base_description = base_description + [],doFit = True,fitfunction="convolution")

#CreatePlotsFromSelection(HM,"20TRTHitsNonZeroEnergyHardScatter", filename, base_description = ["N_{TRT} >= 20", "E_{TOTAL} != 0.0"], doFit = False)
#CreatePlotsFromSelection(HM,"MIPSelectionHadFracAbove70HardScatter", filename, base_description = ["MIP Selection"],doFit=False)
#CreatePlotsFromSelection(HM,"NonZeroEnergyHardScatter", filename, base_description = ["E_{TOTAL} != 0.0"],doFit=False)
#CreatePlotsFromSelection(HM,"InclusiveHardScatter", filename, base_description = [],doFit=False)

#CreatePlotsFromSelection(HM,"20TRTHitsNonZeroEnergyHardScatterOnlyPion", filename, base_description = ["N_{TRT} >= 20", "E_{TOTAL} != 0.0"], doFit = False)
#CreatePlotsFromSelection(HM,"MIPSelectionHadFracAbove70HardScatterOnlyPion", filename, base_description = ["MIP Selection"],doFit=False)
#CreatePlotsFromSelection(HM,"NonZeroEnergyHardScatterOnlyPion", filename, base_description = ["E_{TOTAL} != 0.0"],doFit=False)
#CreatePlotsFromSelection(HM,"InclusiveHardScatterOnlyPion", filename, base_description = [],doFit=False)


if True:
        histogramName = "TwoDTrackPtVsEtaHistogram_HasExtrapolation"
        hist = HM.getHistograms(histogramName)
        description = base_description + ["Inclusive Selection"]
        DataVsMC = Draw2DHistogramOnCanvas(hist["PythiaJetJet"], doLogx = False, doLogy = True)
        DataVsMC.Print(plotter_directory + "/PythiaJetJet" + histogramName + ".png")
        DataVsMC.Close()

        histogramName = "TrkEtaPhiEMCal_MomentumBetween3And4GeV_Denomenator"
        hist = HM.getHistograms(histogramName)
        DataVsMC = Draw2DHistogramOnCanvas(hist["PythiaJetJet"], doLogx = False, doLogy = False, x_range = (-1.0, +1.0), zlabel = "N(E!=0)")
        DataVsMC.Print(plotter_directory + "/PythiaJetJet" + histogramName + ".png")
        DataVsMC.Draw()
        DataVsMC.Close()
        DataVsMC = Draw2DHistogramOnCanvas(hist["LowMuData"], doLogx = False, doLogy = False, x_range = (-1.0, +1.0), zlabel = "N(E!=0)")
        DataVsMC.Print(plotter_directory + "/LowMuData" + histogramName + ".png")
        DataVsMC.Draw()
        DataVsMC.Close()
        hist_den = hist

        histogramName = "TrkEtaPhiEMCal_MomentumBetween3And4GeV_Numerator"
        hist = HM.getHistograms(histogramName)
        DataVsMC = Draw2DHistogramOnCanvas(hist["PythiaJetJet"], doLogx = False, doLogy = False, x_range = (-1.0, +1.0), zlabel = "N(E!=0)")
        DataVsMC.Print(plotter_directory + "/PythiaJetJet" + histogramName + ".png")
        DataVsMC.Draw()
        DataVsMC.Close()
        DataVsMC = Draw2DHistogramOnCanvas(hist["LowMuData"], doLogx = False, doLogy = False, x_range = (-1.0, +1.0), zlabel = "N(E!=0)")
        DataVsMC.Print(plotter_directory + "/LowMuData" + histogramName + ".png")
        DataVsMC.Draw()
        DataVsMC.Close()
        hist_num = hist

        ratio = DivideHistograms(hist_num, hist_den)
        DataVsMC = Draw2DHistogramOnCanvas(ratio["PythiaJetJet"], doLogx = False, doLogy = False, x_range = (-1.0, +1.0), zlabel = "N(E!=0)/N(Inclusive)")
        DataVsMC.Draw()
        DataVsMC.Print(plotter_directory + "/PythiaJetJet" + "TrkEtaPhiEMCal_MomentumBetween3And4GeV_ZeroFraction" + ".png")
        DataVsMC.Close()

        DataVsMC = Draw2DHistogramOnCanvas(ratio["LowMuData"], doLogx = False, doLogy = False, x_range = (-1.0, +1.0), zlabel = "N(E!=0)/N(Inclusive)")
        DataVsMC.Draw()
        DataVsMC.Print(plotter_directory + "/LowMuData" + "TrkEtaPhiEMCal_MomentumBetween3And4GeV_ZeroFraction" + ".png")
        DataVsMC.Close()

        histogramName = "trkAverageMu"
        hist = HM.getHistograms(histogramName)
        description = base_description + ["Inclusive Selection"]
        DataVsMC1 = DrawDataVsMC(hist,\
                                channelLabels,\
                                MCKeys = ['PythiaJetJet'],\
                                DataKey='LowMuData',\
                                ratio_min=0.6,\
                                ratio_max=1.4,\
                                extra_description = description)
        DataVsMC1[0].Draw()
        DataVsMC1[0].Print(plotter_directory + "/" + histogramName + ".png")
        DataVsMC1[0].Close()

        histogramName = "LeadingPtTrkHist"
        hist = HM.getHistograms(histogramName)
        description = base_description + ["Inclusive Selection"]
        DataVsMC1 = DrawDataVsMC(hist,\
                                channelLabels,\
                                MCKeys = ['PythiaJetJet'],\
                                DataKey='LowMuData',\
                                ratio_min=0.6,\
                                ratio_max=1.4,\
                                doLogx=True,\
                                xlabel="Leading Track P_{T} [GeV]",\
                                ylabel="Number of Events",\
                                extra_description = description)
        DataVsMC1[0].Draw()
        DataVsMC1[0].Print(plotter_directory + "/" + histogramName + ".png")
        DataVsMC1[0].Close()

        histogramName = "SubleadingPtTrkHist"
        hist = HM.getHistograms(histogramName)
        description = base_description + ["Inclusive Selection"]
        DataVsMC1 = DrawDataVsMC(hist,\
                                channelLabels,\
                                MCKeys = ['PythiaJetJet'],\
                                DataKey='LowMuData',\
                                ratio_min=0.2,\
                                ratio_max=1.8,\
                                doLogx=True,\
                                xlabel="Subleading Track P_{T} [GeV]",\
                                ylabel="Number of Events",\
                                extra_description = description)
        DataVsMC1[0].Draw()
        DataVsMC1[0].Print(plotter_directory + "/" + histogramName + ".png")
        DataVsMC1[0].Close()


        histogramName = "trkNPV2"
        hist = HM.getHistograms(histogramName)
        description = base_description + ["Inclusive Selection"]
        DataVsMC3 = DrawDataVsMC(hist,\
                                channelLabels,\
                                MCKeys = ['PythiaJetJet'],\
                                DataKey='LowMuData',\
                                extra_description = description)
        DataVsMC3[0].Draw()
        DataVsMC3[0].Print(plotter_directory + "/" + histogramName + ".png")
        DataVsMC3[0].Close()

        histogramName =  "eventNPV2Hist"
        hist = HM.getHistograms(histogramName)
        description = base_description + ["Inclusive Selection"]
        DataVsMC5 = DrawDataVsMC(hist,\
                                channelLabels,\
                                MCKeys = ['PythiaJetJet'],\
                                DataKey='LowMuData',\
                                ratio_min=0.2,\
                                ratio_max=1.8,\
                                extra_description = description)
        DataVsMC5[0].Draw()
        DataVsMC5[0].Print(plotter_directory + "/" + histogramName + ".png")

        histogramName =  "eventAverageMu"
        hist = HM.getHistograms(histogramName)
        description = base_description + ["Inclusive Selection"]
        DataVsMC6 = DrawDataVsMC(hist,\
                                channelLabels,\
                                MCKeys = ['PythiaJetJet'],\
                                DataKey='LowMuData',\
                                ratio_min=0.2,\
                                ratio_max=1.8,\
                                extra_description = description)
        DataVsMC6[0].Draw()
        DataVsMC6[0].Print(plotter_directory + "/" + histogramName + ".png")

        histogramName =  "InclusiveEOP"
        hist = HM.getHistograms(histogramName)
        description = base_description + ["Inclusive Selection"]
        DataVsMC7 = DrawDataVsMC(hist,\
                                channelLabels,\
                                MCKeys = ['PythiaJetJet'],\
                                DataKey='LowMuData',\
                                extra_description = description)
        DataVsMC7[0].Draw()
        DataVsMC7[0].Print(plotter_directory + "/" + histogramName + ".png")


        for extraString in ["", "HasExtrapolation"]:
            histogramName = "trkPtHist" + extraString
            hist = HM.getHistograms(histogramName)
            description = base_description + ["Inclusive Selection"]
            DataVsMC4 = DrawDataVsMC(hist,\
                                    channelLabels,\
                                    MCKeys = ['PythiaJetJet'],\
                                    DataKey='LowMuData',\
                                    doLogx = True,\
                                    doLogy = True,\
                                    ratio_min = 0.8,\
                                    ratio_max = 1.2,\
                                    extra_description = description)
            DataVsMC4[0].Draw()
            DataVsMC4[0].Print(plotter_directory + "/" + histogramName + ".png")



        histogramName_num =  "InclusiveZeroFractionVsPNumerator"
        hist_num = HM.getHistograms(histogramName_num)

        histogramName_den = "InclusiveZeroFractionVsPDenomenator"
        hist_den = HM.getHistograms(histogramName_den)
        ratio_hist = DivideHistograms(hist_num, hist_den)

        description = base_description + ["Inclusive Selection"]
        DataVsMC9 = DrawDataVsMC(ratio_hist,\
                                channelLabels,\
                                MCKeys = ['PythiaJetJet'],\
                                DataKey='LowMuData',\
                                doLogx=True,
                                doLogy=False,
                                ratio_min = 0.8,\
                                ratio_max = 1.2,\
                                extra_description = description)
        DataVsMC9[0].Draw()
        DataVsMC9[0].Print(plotter_directory + "/" + histogramName_num.replace("Numerator", "") + ".png")

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

        DataVSMC10 = DrawDataVsMC(histograms,\
                                  channelLabels,\
                                  MCKeys = ['PythiaJetJet'],\
                                  DataKey = "LowMuData",\
                                  doLogx = True,\
                                  doLogy = False,
                                  ratio_min = 0.9,\
                                  ratio_max = 1.1,\
                                  extra_description = description)
        DataVSMC10[0].Draw()
        DataVSMC10[0].Print(plotter_directory + "/" + histogramName + ".png")


        description = base_description + ["0.8<|#eta_{ID}|<1.2", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
        histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_8_12"
        histograms = HM.getHistograms(histogramName)
        for key in histograms:
            histograms[key] = histograms[key].ProjectionX(histogramName + "_px", "E")

        DataVSMC10 = DrawDataVsMC(histograms,\
                                  channelLabels,\
                                  MCKeys = ['PythiaJetJet'],\
                                  DataKey = "LowMuData",\
                                  doLogx = True,\
                                  doLogy = False,
                                  ratio_min = 0.9,\
                                  ratio_max = 1.1,\
                                  extra_description = description)
        DataVSMC10[0].Draw()
        DataVSMC10[0].Print(plotter_directory + "/" + histogramName + ".png")


        description = base_description + ["1.2<|#eta_{ID}|<1.6", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
        histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_12_16"
        histograms = HM.getHistograms(histogramName)
        for key in histograms:
            histograms[key] = histograms[key].ProjectionX(histogramName + "_px", "E")

        DataVSMC10 = DrawDataVsMC(histograms,\
                                  channelLabels,\
                                  MCKeys = ['PythiaJetJet'],\
                                  DataKey = "LowMuData",\
                                  doLogx = True,\
                                  doLogy = False,
                                  ratio_min = 0.9,\
                                  ratio_max = 1.1,\
                                  extra_description = description)
        DataVSMC10[0].Draw()
        DataVSMC10[0].Print(plotter_directory + "/" + histogramName + ".png")

        description = base_description + ["1.6<|#eta_{ID}|<2.0", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
        histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_16_20"
        histograms = HM.getHistograms(histogramName)
        for key in histograms:
            histograms[key] = histograms[key].ProjectionX(histogramName + "_px", "E")

        DataVSMC10 = DrawDataVsMC(histograms,\
                                  channelLabels,\
                                  MCKeys = ['PythiaJetJet'],\
                                  DataKey = "LowMuData",\
                                  doLogx = True,\
                                  doLogy = False,
                                  ratio_min = 0.9,\
                                  ratio_max = 1.1,\
                                  extra_description = description)
        DataVSMC10[0].Draw()
        DataVSMC10[0].Print(plotter_directory + "/" + histogramName + ".png")

        description = base_description + ["2.0<|#eta_{ID}|<2.4", "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "MIP Selection"]
        histogramName = "EOPProfileVsMomentum_MIPSelection_HadFracAbove70_InBin_20_24"
        histograms = HM.getHistograms(histogramName)
        for key in histograms:
            histograms[key] = histograms[key].ProjectionX(histogramName + "_px", "E")

        DataVSMC10 = DrawDataVsMC(histograms,\
                                  channelLabels,\
                                  MCKeys = ['PythiaJetJet'],\
                                  DataKey = "LowMuData",\
                                  doLogx = True,\
                                  doLogy = False,
                                  ratio_min = 0.9,\
                                  ratio_max = 1.1,\
                                  extra_description = description)
        DataVSMC10[0].Draw()
        DataVSMC10[0].Print(plotter_directory + "/" + histogramName + ".png")

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

        #histogramNames = ["TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_0_4" , "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_12_16", "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_16_20", "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_20_24", "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_4_8", "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_8_12"]
        #
        eta_descriptors = ["0.0<|#eta_{ID}|<0.4", "0.4<|#eta_{ID}|<0.8", "0.8<|#eta_{ID}|<1.2", "1.2<|#eta_{ID}|<1.6", "1.6<|#eta_{ID}|<2.0", "2.0<|#eta_{ID}|<2.4"]
        #
        #for histogramName, descriptor in zip(histogramNames, eta_descriptors):
        #    description = base_description + [descriptor, "#frac{E_{HAD}}{E_{TOTAL}} > 0.7", "E^{dR<0.1}_{EM} < 1.1 GeV"]
        #    histograms = HM.getHistograms(histogramName)
        #
        #    for key in histograms:
        #        histograms[key].Rebin(20)
        #
        #    DataVSMC10 = DrawDataVsMC(histograms,\
        #                          channelLabels,\
        #                          MCKeys = ['PythiaJetJet'],\
        #                          DataKey = "LowMuData",\
        #                          doLogx = True,\
        #                          doLogy = False,
        #                          ratio_min = 0.6,\
        #                          ratio_max = 1.4,\
        #                          extra_description = description)
        #    DataVSMC10.Draw()
        #    DataVSMC10.Print(plotter_directory + "/" + histogramName + ".png")
        eta_ranges = [(0.0, 0.4), (0.4, 0.8), (0.8, 1.2), (1.2, 1.6), (1.6, 2.0), (2.0, 2.4)]
        profileNames = ["EOPProfileVsMomentum", "EOPProfileVsMomentum_MIPSelection_HadBetween30And90OfMomentum", "EOPProfileVsMomentum_MIPSelection_HadFracAbove70", "EOPProfileVsMomentum_NonZeroE"]
        bkgProfileNames = ["EnergyBkgProfileVsMomentum", "EnergyBkgProfileVsMomentum_MIPSelection_HadBetween30And90OfMomentum", "EnergyBkgProfileVsMomentum_MIPSelection_HadFracAbove70", "EnergyBkgProfileVsMomentum_NonZeroE"]
        TwoDHistNames = ["2DHist_EOPVsMomentum", "2DHist_EOPVsMomentum_MIPSelection_HadBetween30And90OfMomentum", "2DHist_EOPVsMomentum_MIPSelection_HadFracAbove70", "2DHist_EOPVsMomentum_NonZeroE"]
        #bkgTwoDHistNames = ["2DHist_EOPBkgVsMomentum", "2DHist_EOPBkgVsMomentum_MIPSelection_HadBetween30And90OfMomentum", "2DHist_EOPBkgVsMomentum_MIPSelection_HadFracAbove70", "2DHist_EOPBkgVsMomentum_NonZeroE"]
        plotDescriptors = [ [], ["0.3 P < E_{HAD} < 0.9 P", "E^{dR<0.1}_{EM} < 1.1 GeV", "N_{TRT} >= 20"], ["E_{HAD}/E_{TOTAL} > 0.7", "E^{dR<0.1}_{EM} < 1.1 GeV", "N_{TRT} >= 20"], ["E_{TOTAL} != 0.0"]]

        for eta_range, eta_descriptor in zip(eta_ranges, eta_descriptors):
            print eta_range
            print eta_descriptor
            histogramName = "TrkMultiplicityVsP_NonZeroE" + "_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
            hist_NonZero = HM.getHistograms(histogramName, rebin = 100)
            histogramName = "TrkMultiplicityVsP" + "_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
            hist_Inclusive = HM.getHistograms(histogramName, rebin = 100)
            histogramName = "TrkMultiplicityVsP_MIPSelection_HadFracAbove70" + "_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) 
            hist_HadFracAbove70 = HM.getHistograms(histogramName, rebin = 100)

            frac_MIP_of_NonZero = DivideHistograms(hist_HadFracAbove70, hist_NonZero)
            frac_MIP_of_Inclusive = DivideHistograms(hist_HadFracAbove70, hist_Inclusive)

            DataVSMC10 = DrawDataVsMC(frac_MIP_of_NonZero,\
                                  channelLabels,\
                                  MCKeys = ['PythiaJetJet'],\
                                  DataKey = "LowMuData",\
                                  doLogx = True,\
                                  doLogy = False,
                                  ratio_min = 0.2,\
                                  ratio_max = 1.8,\
                                  ylabel="N(MIP)/N(E!=0)",\
                                  extra_description = base_description  + [eta_descriptor])
            DataVSMC10[0].Draw()
            DataVSMC10[0].Print(plotter_directory + "/" + "FracNonZero_MIP_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + ".png")
            DataVSMC10[0].Close()

            DataVSMC10 = DrawDataVsMC(frac_MIP_of_Inclusive,\
                                  channelLabels,\
                                  MCKeys = ['PythiaJetJet'],\
                                  DataKey = "LowMuData",\
                                  doLogx = True,\
                                  doLogy = False,
                                  ratio_min = 0.6,\
                                  ratio_max = 1.4,\
                                  ylabel="N(MIP)/N(Inclusive)",\
                                  extra_description = base_description +  [eta_descriptor])
            DataVSMC10[0].Draw()
            DataVSMC10[0].Print(plotter_directory + "/" + "FracInclusive_MIP_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + ".png")
            DataVSMC10[0].Close()




            for bkgProfileName, profileName, TwoDHistName, plotDescriptor in zip(bkgProfileNames, profileNames, TwoDHistNames, plotDescriptors):
                description = ["P_{T} Reweighted", eta_descriptor] + plotDescriptor

                histogramName = profileName + "_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
                histograms = HM.getHistograms(histogramName)
                histograms = ProjectProfiles(histograms)
                histogram_num = histograms
                DataVSMC10 = DrawDataVsMC(histograms,\
                                      channelLabels,\
                                      MCKeys = ['PythiaJetJet'],\
                                      DataKey = "LowMuData",\
                                      doLogx = True,\
                                      doLogy = False,
                                      ratio_min = 0.5,\
                                      ratio_max = 1.5,\
                                      extra_description = description)
                DataVSMC10[0].Draw()
                DataVSMC10[0].Print(plotter_directory + "/" + histogramName + ".png")
                DataVSMC10[0].Close()

                histogramName = bkgProfileName + "_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
                histograms = HM.getHistograms(histogramName)
                histograms = ProjectProfiles(histograms)
                histogram_den = histograms
                DataVSMC10 = DrawDataVsMC(histograms,\
                                      channelLabels,\
                                      MCKeys = ['PythiaJetJet'],\
                                      DataKey = "LowMuData",\
                                      doLogx = True,\
                                      doLogy = False,
                                      ratio_min = 0.5,\
                                      ratio_max = 1.5,\
                                      extra_description = description)
                DataVSMC10[0].Draw()
                DataVSMC10[0].Print(plotter_directory + "/" + histogramName + ".png")
                DataVSMC10[0].Close()

                EOPCorrHistograms = SubtractHistograms(histogram_num, histogram_den)
                histogramName = "EOPCorrHistogram_" + bkgProfileName + "_"  + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
                DataVSMC10 = DrawDataVsMC(EOPCorrHistograms,\
                                      channelLabels,\
                                      MCKeys = ['PythiaJetJet'],\
                                      DataKey = "LowMuData",\
                                      doLogx = True,\
                                      doLogy = False,\
                                      ylabel = "<E/p>_{Corr}",\
                                      ratio_min = 0.95,\
                                      ratio_max = 1.05,\
                                      extra_description = description)
                DataVSMC10[0].Draw()
                DataVSMC10[0].Print(plotter_directory + "/" + histogramName + ".png")
                DataVSMC10[0].Close()

                histogramName = TwoDHistName + "_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
                histograms = HM.getHistograms(histogramName)

                DataCanvas = Draw2DHistogramOnCanvas(histograms["LowMuData"], doLogx = True, doLogy = False)
                DataCanvas.Print(plotter_directory + "/" + histogramName + "LowMuData.png")

                MCCanvas = Draw2DHistogramOnCanvas(histograms["PythiaJetJet"], doLogx = True, doLogy = False)
                MCCanvas.Print(plotter_directory + "/" + histogramName + "PythiaJetJet.png")

            histogram_name =  "TrkMultiplicityVsP_NonZeroE_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
            histograms = HM.getHistograms(histogram_name)
            DataVsMC = DrawDataVsMC(histograms,\
                                    channelLabels,\
                                    MCKeys = ['PythiaJetJet'],\
                                    DataKey="LowMuData",\
                                    ratio_min=0.6,\
                                    ratio_max=1.4,\
                                    rebin=100,\
                                    extra_description = ["P_{T} Reweighted", "E_{Total} != 0.0", eta_descriptor])
            DataVsMC[0].Draw()
            DataVsMC[0].Print(plotter_directory + "/" + histogram_name + ".png")
            DataVsMC[0].Close()
        #plot the average energy in the anulus
        profileNames = ["EnergyAnulusProfileVsMomentum_MIPSelection_HadBetween30And90OfMomentum", "EnergyAnulusProfileVsMomentum_MIPSelection_HadFracAbove70"]
        TwoDHistNames = ["2DHist_EnergyAnulusVsMomentum_MIPSelection_HadBetween30And90OfMomentum", "2DHist_EnergyAnulusVsMomentum_MIPSelection_HadFracAbove70"]
        profileDescriptors = [ ["0.3 P < E_{HAD} < 0.9 P", "E^{dR<0.1}_{EM} < 1.1 GeV", "N_{TRT} >= 20"], ["E_{HAD}/E_{TOTAL} > 0.7", "E^{dR<0.1}_{EM} < 1.1 GeV", "N_{TRT} >= 20"], ["E_{TOTAL} != 0.0"] ]

        for eta_range, eta_descriptor in zip(eta_ranges, eta_descriptors):
            eta = eta_range[1]
            p_bins_max = 15.05
            p_bins_min = getP(0.5, eta)
            nBins = 20
            p_bins = getLogBins(p_bins_min, p_bins_max, nBins)
            for profileName, TwoDHistName, profileDescriptor in zip(profileNames, TwoDHistNames,  profileDescriptors):
                description = ["P_{T} Reweighted", eta_descriptor] + profileDescriptor

                histogramName = profileName + "_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
                histograms = HM.getHistograms(histogramName)
                histograms = ProjectProfiles(histograms)
                DataVSMC10 = DrawDataVsMC(histograms,\
                                      channelLabels,\
                                      MCKeys = ['PythiaJetJet'],\
                                      DataKey = "LowMuData",\
                                      doLogx = True,\
                                      doLogy = False,
                                      ratio_min = 0.4,\
                                      ratio_max = 1.6,\
                                      extra_description = description)
                DataVSMC10[0].Draw()
                DataVSMC10[0].Print(plotter_directory + "/" + histogramName + ".png")
                DataVSMC10[0].Close()

                histogramName = TwoDHistName + "_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
                histograms = HM.getHistograms(histogramName)
                DataCanvas = Draw2DHistogramOnCanvas(histograms["LowMuData"], doLogx = True, doLogy = False)
                DataCanvas.Print(plotter_directory + "/" + histogramName + "LowMuData.png")

                MCCanvas = Draw2DHistogramOnCanvas(histograms["PythiaJetJet"], doLogx = True, doLogy = False)
                MCCanvas.Print(plotter_directory + "/" + histogramName + "PythiaJetJet.png")

                p_ranges = [ (p_bins[i], p_bins[i+1])  for i in range(0, len(p_bins)-1) ]
                for p_range in p_ranges:
                    histogramNames = ["trkEMDR100", "MomentumHadFrac", "HadFrac", "trkTRTHits"]
                    p_high_str = "{:.2f}".format(p_range[1])
                    p_low_str= "{:.2f}".format(p_range[0])
                    for selection_type in ["_NonZeroE_20TRT_InEtaBin_", "_NonZeroE_InEtaBin_", "_MIPSelection_HadFracAbove70_InEtaBin_"]:
                        if "_NonZeroE_20TRT_InEtaBin_" == selection_type:
                            extra_stuff = ["E_{TOTAL} != 0", "N_{TRT Hits} >= 20"]
                        elif "_MIPSelection_HadFracAbove70_InEtaBin_" == selection_type:
                            extra_stuff = ["MIP Seleciton"]
                        else:
                            extra_stuff = ["E_{TOTAL} !=0"]
                        for histogramName in histogramNames:
                            if selection_type == "_MIPSelection_HadFracAbove70_":
                                break
                            histogramName = histogramName + selection_type + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" + str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
                            histograms = HM.getHistograms(histogramName)
                            DataVSMC = DrawDataVsMC(histograms,\
                                           channelLabels,\
                                           MCKeys = ['PythiaJetJet'],\
                                           DataKey = "LowMuData",\
                                           doLogx = False,\
                                           doLogy = False,
                                           ratio_min = 0.4,\
                                           ratio_max = 1.6,\
                                           extra_description = ["P_{T} Reweighted", eta_descriptor, p_low_str + " < |P/GeV| < " + p_high_str] + extra_stuff)
                            DataVSMC[0].Draw()
                            DataVSMC[0].Print(plotter_directory + "/" + histogramName + ".png")
                            DataVSMC[0].Close()

                        for histogramName in ["NClusters","NClusters_EM","NClusters_HAD","NClusters_emlike","NClusters_hadlike"]:
                            histogram_name = histogramName + selection_type + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1])) + "_InPBin_" +  str(int(100*p_range[0])) + "_" + str(int(100*p_range[1]))
                            hist = HM.getHistograms(histogram_name)
                            DataVsMC1 = DrawDataVsMC(hist,\
                                                    channelLabels,\
                                                    MCKeys = ['PythiaJetJet'],\
                                                    DataKey='LowMuData',\
                                                    extra_description =  ["P_{T} Reweighted", eta_descriptor, p_low_str + " < |P/GeV| < " + p_high_str] + extra_stuff)
                            DataVsMC1[0].Draw()
                            DataVsMC1[0].Print(plotter_directory + "/" + histogram_name + ".png")
                            DataVsMC1[0].Close()
                        ROOT.gSystem.ProcessEvents()

raw_input()
#['NTRT20ZeroFractionVsPetaID00_06Denomenator', 'TrkPtHisteta15_18', 'trkNPV2', 'eventAverageMu', 'trkCount', 'trkAverageMu', 'NTRT20ZeroFractionVsPetaID06_11Denomenator', 'NTRT20ZeroFractionVsPetaID02_04Numerator', 'EtaID0_6_PGreater2_0_EOPHist', 'NonZero_EtaID0_6_PBetween12_18_EOPHist', 'NTRT20ZeroFractionVsPetaID00_02Denomenator', 'NonZero_EtaIDBetween19_23_PBetween22_28_EOPHist', 'EtaID0_6_PGreater1_5_EOPHist', 'EtaLess08_TwoDHistTrkPvsPhiInnerToExtrapolEM2', 'InclusiveZeroFractionVsPNumerator', 'TrkPtHisteta06', 'NTRT20ZeroFractionVsPetaID11_14Numerator', 'NTRT20ZeroFractionVsPetaID00_06Numerator', 'NonZeroEnergy_InclusiveEOP', 'TrkPtHisteta18_23', 'InclusiveEOP', 'ZeroFractionVsPetaID06_11Denomenator', 'ZeroFractionVsPetaID02_04Denomenator', 'NTRT20ZeroFractionVsPetaID04_06Numerator', 'trkPtHist', 'NTRT20ZeroFractionVsPetaID00_02Numerator', 'TrkPtHisteta06_11', 'InclusiveZeroFractionVsAbsEtaNumerator', 'ZeroFractionVsPetaID04_06Denomenator', 'NearestDRHist', 'TrackEtaID', 'NonZero_EtaID0_6_PBetween22_28_EOPHist', 'NTRT20ZeroFractionVsPetaID06_11Numerator', 'InclusiveZeroFractionVsEtaNumerator', 'InclusiveZeroFractionVsPDenomenator', 'TwoDTrackPtVsEtaHistogram', 'ZeroFractionVsPetaID11_14Denomenator', 'ZeroFractionVsPetaID00_02Denomenator', 'eventNPV2Hist', 'trkEtaECALHist', 'TrkPtHisteta14_15', 'InclusiveZeroFractionVsEtaDenomenator', 'EtaID0_6_PBetween22_28_EOPHist', 'EtaID0_6_PGreater1_0_EOPHist', 'EtaIDBetween19_23_PBetween28_36_EOPHist', 'NonZero_EtaID0_6_PBetween28_36_EOPHist', 'TrkPtHisteta11_14', 'lowPTLess07_TwoDHistTrkEtavsDEtaInnerToExtrapolEM2', 'NTRT20ZeroFractionVsPetaID04_06Denomenator', 'ZeroFractionVsPetaID00_06Denomenator', 'TwoDHistTrkPvsPhiInnerToExtrapolEM2', 'EtaID0_6_PBetween12_18_EOPHist', 'NTRT20ZeroFractionVsPetaID11_14Denomenator', 'InclusiveZeroFractionVsAbsEtaDenomenator', 'TwoDTrackPvsTrkEtaID', 'NTRT20ZeroFractionVsPetaID02_04Denomenator', 'EtaID0_6_PBetween28_36_EOPHist']

#Loop through the data and MC histograms and perform fits with landau distributions




#create a set of strings that could describe the eta or momentum selections
eta_ranges = [(0.0, 0.4),(0.4,0.8),(0.8,1.2),(1.2,1.6),(1.6,2.0),(2.0,2.4)]

#go and get the average E/P for MIP particles in each of the eta bins.
for eta_range in eta_ranges:

    eta_high_str = "{:.1f}".format(eta_range[1])
    eta_low_str= "{:.1f}".format(eta_range[0])

    #create nbins where there are at least 10,000 entries per bin
    #histogramName = "TrkMultiplicityVsP_MIPSelection_HadFracAbove70_InBin_" + str(int(10*eta_range[0])) + "_" + str(int(10*eta_range[1]))
    #binningHistogram = HM.getHistograms(histogramName)["LowMuData"]
    #FourThousandTracks_pbins = GetBinsFromHistogram(binningHistogram, 4000.0)[0]
    p_bins_max = 15.05
    p_bins_min = getP(0.5, eta_range[1])
    nBins = 20
    p_bins = getLogBins(p_bins_min, p_bins_max, nBins)
    if eta_range[0] < 0.2:
        p_bins = p_bins[1:]

    #loop through the different eta and p bins and perform the fits:
    #p_ranges = [ (FourThousandTracks_pbins[i], FourThousandTracks_pbins[i+1])  for i in range(0, len(FourThousandTracks_pbins)-1) ]

    #create a histogram with each of the bins set to those from the four thousand track histograms

    bin_edges = p_bins
    bin_mpv_data = []
    bin_mpv_error_data = []
    bin_mpv_MC = []
    bin_mpv_error_MC = []

    count = 0
    p_ranges = [ (p_bins[i], p_bins[i+1])  for i in range(0, len(p_bins)-1) ]

    if False:

        description = ["P_{T} Reweighted", "MIP Selection", eta_low_str + " < |#eta| < " + eta_high_str , p_low_str + " < P/GeV < " + p_high_str]

        DataVsMC = DrawDataVsMC(histograms,\
                              channelLabels,\
                              MCKeys = ['PythiaJetJet'],\
                              DataKey = "LowMuData",\
                              doLogx = False,\
                              doLogy = False,
                              ratio_min = 0.4,\
                              ratio_max = 1.6,\
                              xAxis_range = (0.0, 2.0),\
                              extra_description = description)
        DataVsMC[0].Draw()
        DataVsMC[0].cd()
        raw_input()
        top_pad = DataVsMC[1]
        top_pad.cd()
        fit_function_MC.Draw("Same")

        if fit_function_MC.GetNDF() <= 1:
            continue

        #draw a little text thing describing the fit result
        chisq_MC_str = "{:.3f}".format(fit_function_MC.GetChisquare()/fit_function_MC.GetNDF())
        prob_MC_str = "{:.3f}".format(fit_function_MC.GetProb())
        mpv_MC_str = "{:.3f}".format(fit_function_MC.GetParameter(1))
        mpvErr_MC_str = "{:.3f}".format(fit_function_MC.GetParError(1))
        sigma_MC_str = "{:.3f}".format(fit_function_MC.GetParameter(1))
        sigmaErr_MC_str = "{:.3f}".format(fit_function_MC.GetParError(1))

        chisq_data_str = "{:.3f}".format(fit_function_data.GetChisquare()/fit_function_data.GetNDF())
        prob_data_str = "{:.3f}".format(fit_function_data.GetProb())
        mpv_data_str = "{:.3f}".format(fit_function_data.GetParameter(1))
        mpvErr_data_str = "{:.3f}".format(fit_function_data.GetParError(1))
        sigma_data_str = "{:.3f}".format(fit_function_data.GetParameter(1))
        sigmaErr_data_str = "{:.3f}".format(fit_function_data.GetParError(1))

        chisq_MC = fit_function_MC.GetChisquare()/fit_function_MC.GetNDF()
        prob_MC = fit_function_MC.GetProb()

        if fit_function == "convolution":
            mpv_MC = fit_function_MC.GetParameter(1)
            mpvErr_MC = fit_function_MC.GetParError(1)

        else:
            mpv_MC = fit_function_MC.GetParameter(0)
            mpvErr_MC = fit_function_MC.GetParError(0)

        sigma_MC = fit_function_MC.GetParameter(1)
        sigmaErr_MC = fit_function_MC.GetParError(1)

        chisq_data = fit_function_data.GetChisquare()/fit_function_data.GetNDF()
        prob_data = fit_function_data.GetProb()
        if fit_function == "convolution":
            mpv_data = fit_function_data.GetParameter(1)
            mpvErr_data = fit_function_data.GetParError(1)
        else:
            mpv_data = fit_function_data.GetParameter(0)
            mpvErr_data = fit_function_data.GetParError(0)

        sigma_data = fit_function_data.GetParameter(1)
        sigmaErr_data = fit_function_data.GetParError(1)

        bin_mpv_data.append(mpv_data)
        bin_mpv_error_data.append(mpvErr_data)
        bin_mpv_MC.append(mpv_MC)
        bin_mpv_error_MC.append(mpvErr_MC)

        if fit_function != "convolution":
            dataResult = "Data: #mu=" + mpv_data_str + "#pm" + mpvErr_data_str + " #sigma=" + sigma_data_str + "#pm" + sigmaErr_data_str + " #Chi^2/NDOF=" + chisq_data_str + " Prob=" + prob_data_str
            MCResult = "MC: #mu=" + mpv_MC_str + "#pm" + mpvErr_MC_str + " #sigma=" + sigma_MC_str + "#pm" + sigmaErr_MC_str + " #Chi^2/NDOF=" + chisq_MC_str + " Prob=" + prob_MC_str
        else:
            dataResult = "Data: #mu=" + mpv_data_str + "#pm" + mpvErr_data_str + " #Chi^2/NDOF=" + chisq_data_str + " Prob=" + prob_data_str
            MCResult = "MC: #mu=" + mpv_MC_str + "#pm" + mpvErr_MC_str + " #Chi^2/NDOF=" + chisq_MC_str + " Prob=" + prob_MC_str

        DrawText(0.5, 0.55, dataResult, size=0.035)
        DrawText(0.5, 0.5, MCResult, size=0.035)

        print "Creating the histogram for"
        print plotter_directory + "/" + histogramName + ".png"

        DataVsMC[0].Modified()
        DataVsMC[0].Update()
        DataVsMC[0].Print(plotter_directory + "/" + histogramName + ".png")
        raw_input()
    #create an MC and data histogram with each bin set to the result from the fit

    data_hist = ROOT.TH1D("DataFitResults" + str(eta_range[0]) + "_" + str(eta_range[1]),"DataFitResults" + str(eta_range[0]) + "_" + str(eta_range[1]), len(p_bins[1:])-1, array('d', p_bins[1:]))
    MC_hist = ROOT.TH1D("MCFitResults" + str(eta_range[0]) + "_" + str(eta_range[1]),"DataFitResults" + str(eta_range[0]) + "_" + str(eta_range[1]), len(p_bins[1:])-1, array('d', p_bins[1:]))

    for i in range(1, data_hist.GetNbinsX() + 1):
        data_hist.SetBinContent(i, bin_mpv_data[i-1])
        data_hist.SetBinError(i, bin_mpv_error_MC[i-1])
        MC_hist.SetBinContent(i, bin_mpv_MC[i-1])
        MC_hist.SetBinError(i, bin_mpv_error_MC[i-1])

    histograms = {"PythiaJetJet":MC_hist, "LowMuData":data_hist}

    description = ["P_{T} Reweighted", "MIP Selection", eta_low_str + " < |#eta| < " + eta_high_str]
    DataVsMC = DrawDataVsMC(histograms,\
                          channelLabels,\
                          MCKeys = ['PythiaJetJet'],\
                          DataKey = "LowMuData",\
                          doLogx = True,\
                          doLogy = False,
                          ratio_min = 0.9,\
                          ratio_max = 1.1,\
                          extra_description = description)

    DataVsMC[0].Draw()
    DataVsMC[0].Print("DataFitResults" + str(eta_range[0]) + "_" + str(eta_range[1]) + ".png")
