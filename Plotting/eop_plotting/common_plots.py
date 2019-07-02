import ROOT
from plotting_tools import ProjectProfiles, DrawDataVsMC
from array import array
from fitting_tools import fitHistograms


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
        for histogram in ["TrackPtSpectrum", "TrackPSpectrum", "TrackTruthPSpectrum"]:
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
                                    MCKeys = MCKeys,\
                                    #MCKeys = MCKeys,\
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

def CreatePlotsFromSelection(HM, selection_name, filename, base_description = [], channelLabels={}, MCKeys=[], plotting_directory=""):
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
                                   "TrackTruthPSpectrum",\
                                   "EOPProfileVsMomentum",\
                                   "EnergyAnulusProfileVsMomentum",\
                                   "EnergyBkgProfileVsMomentum"]

    histograms_in_momentum_bins = ["EOPDistribution",\
                                    "EOPBkgDistribution",\
                                    "trkTRTHits",\
                                    "trkEMDR100",\
                                    "MomentumHadFrac",\
                                    "HadFrac",\
#                                   "NClusters",\
#                                   "NClusters_EM",\
#                                   "NClusters_HAD",\
#                                   "NClusters_emlike",\
#                                   "NClusters_hadlike",\
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
            to_plot = MCKeys
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
            DataVsMC1[0].Print(plotting_directory + "/" + histogram_name + "_"  + selection_name + ".png")
            DataVsMC1[0].Close()

        p_bins_low = p_bins_low_for_eta_bin[i]
        p_bins_high = p_bins_high_for_eta_bin[i]

        #the plots binned in eta and momentum
        for histogram in histograms_in_momentum_bins:
            print("Fetching histograms {}".format(histogram))
            if histogram == "EOPDistribution":
                rebin = 20
            if histogram == "HadFrac":
                rebin = 2
            else:
                rebin = None
            to_plot = []
            for j, p_low, p_high in zip(list(range(0, p_bins_low.size())), p_bins_low, p_bins_high):
                print("Got histogram in bin {} with momentum between {} and {}".format(j, p_low, p_high))
                histattributes = []
                histogram_name = histogram + "_" + selection_name + "_Eta_" + str(i) + "_Momentum_" + str(j)
                hist = HM.getHistograms(histogram_name)
                description = base_description + [str(round(eta_low, 2)) + " < |#eta| < " + str(round(eta_high, 2))]
                description += [str(round(p_low, 3)) + " < P/GeV < " + str(round(p_high, 3))]

                DataVsMC1 = DrawDataVsMC(hist,\
                                        channelLabels,\
                                        MCKeys = MCKeys,\
                                        xAxis_range=xAxis_range,\
                                        rebin=rebin,\
                                        doLogy=False,\
                                        DataKey='LowMuData',\
                                        extra_description = description)

                #DataVsMC1[0].Draw()
                #DataVsMC1[0].Print(plotting_directory + "/" + histogram_name + ".png")
                to_plot.append(DataVsMC1[0].Clone())
                #DataVsMC1[0].Close()

            #create the canvas of all of the plots
            total_canvas = ROOT.TCanvas(histogram, histogram, 3000, 3000)
            total_canvas.Draw()
            total_canvas.Divide(5, 3)
            for j, canvas in enumerate(to_plot):
                total_canvas.cd(j+1)
                canvas.DrawClonePad()
            total_canvas.Draw()
            canvas_name = plotting_directory + "/" + histogram + "_"  + selection_name + "_Eta_" + str(i) + ".png"
            print("Printing on png file {}".format(canvas_name))
            total_canvas.Print(canvas_name)
            total_canvas.Close()
