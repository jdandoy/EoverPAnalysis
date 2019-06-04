import ROOT
from PlottingTools.Plotter import ProjectProfiles, DrawDataVsMC
from array import array
from FittingTools.FittingTools import fitHistograms


def CreateZeroFractionPlotsFromSelection(HM, numerator_selection_name, denomenator_selection_name, filename, base_description=[], channelLabels={}, plotting_directory=""):
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
            DataVsMC1[0].Print(plotting_directory + "/" + histogram_name + ".png")
            DataVsMC1[0].Close()

def CreatePlotsFromSelection(HM, selection_name, filename, base_description = [], doFit = False, fitfunction="gaus", refit=False,channelLabels={}, plotting_directory=""):
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
            DataVsMC1[0].Print(plotting_directory + "/" + histogram_name + ".png")
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
                DataVsMC1[0].Print(plotting_directory + "/" + histogram_name + ".png")
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
