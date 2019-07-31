import ROOT
from plotting_tools import ProjectProfiles, DrawDataVsMC, DivideHistograms, cleanUpHistograms
from array import array
from fitting_tools import fitHistograms
from histogram_manager import HistogramManager


# a funciton to show the composition of tracks from the hard scatter as a function of track pT
def CreateCompositionPlot(HM, folder):
    ROOT.gROOT.SetBatch(ROOT.kFALSE)
    particles = ["Pion", "Kaon", "Proton"]
    charges = ["Pos", "Neg"]
    channels = ["PythiaJetJetHardScatter"]
    histogram_name_base = "trkPtHist"

    labels = []
    channel_names = []
    for p in particles:
        label_str = ""
        if p == "Kaon":
            label_str += "#\Kappa"
        if p == "Pion":
            label_str += "#\pi"
        if p == "Proton":
            label_str += "#\mathrm{p}"
        for c in charges:
            if p != "Proton":
                if c == "Pos":
                    label_str = label_str.rstrip("#") +  "^{+}"
                if c == "Neg":
                    label_str = label_str.rstrip("#") +  "^{-}"
            else:
                if c == "Neg":
                    label_str = "#\\bar{" + label_str.strip("#") + "}"
            for ch in channels:
                channel_names.append("{}{}{}".format(ch,p,c))
                labels.append(label_str)

    for ch in channels:
        channel_names.append("{}{}".format(ch,"Other"))
        labels.append("Other")

    CB_color_cycle = [r'#377eb8', r'#ff7f00',r'#dede00', r'#4daf4a',\
                  r'#f781bf', r'#a65628', r'#e41a1c', r'#984ea3',\
                  ]
    root_colors = [ROOT.TColor.GetColor(c) for c in CB_color_cycle]


    histograms = HM.getHistograms("trkPtHist")

    histogram_to_normalize_to = histograms["PythiaJetJetHardScatter"]
    stack = ROOT.THStack("Stack", "Stack")

    l=ROOT.TLegend()
    l.SetBorderSize(0)
    l.SetNColumns(2)

    c1 = ROOT.TCanvas()
    for channel, color, label in zip(channel_names, root_colors, labels):
        hist = histograms[channel]
        hist.Divide(histogram_to_normalize_to)
        hist.SetLineColor(color)
        hist.SetFillColor(color)
        hist.SetMarkerColor(color)
        l.AddEntry(hist,label)
        stack.Add(hist)
    c1.Draw()
    stack.Draw("Hist ][")
    stack.GetXaxis().SetTitle("Track P_{T} [GeV]")
    stack.GetYaxis().SetTitle("Fractional Composition")
    l.Draw("SAME")
    c1.SetLogx()
    c1.Update()
    c1.Modified()
    c1.Print(folder + "/FractionalComposition.png")
    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    count = -1
    for channel, color, label in zip(channel_names, root_colors, labels):
        count += 1
        hist = histograms[channel]
        hist = cleanUpHistograms(hist)
        hist.SetLineColor(color)
        hist.SetFillColor(color)
        hist.SetMarkerColor(color)
        hist.SetMinimum(0.0)
        if count == 0:
            hist.Draw("][ Hist")
            stack.GetXaxis().SetTitle("Track P_{T} [GeV]")
            stack.GetYaxis().SetTitle("Fractional Composition")
        else:
            hist.Draw("][ HIST SAME")
    c1.Draw()
    l.Draw("SAME")
    c1.SetLogx()
    c1.Update()
    c1.Modified()
    c1.Print(folder + "/FractionalComposition_nostack.png")
    ROOT.gROOT.SetBatch(ROOT.kTRUE)



def CreateZeroFractionPlotsFromSelection(HM, numerator_selection_name, denomenator_selection_name, filename, base_description=[], channelLabels={}, plotting_directory="",MCKeys=[]):
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

    for i, eta_low, eta_high in zip(list(range(0, eta_bins_low.size())), eta_bins_low, eta_bins_high):
        #the plots binned in eta
        for histogram in histograms_in_eta_bins:
            ylabel = None
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
            if "EOPProfileVs" in histogram_name and not "EnergyBkg" in histogram_name:
                shouldILogy = False
                ratio_min = 0.9
                ratio_max = 1.1
                ylabel = "<E/P>_{Raw}"

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
            ratio_min = 0.6
            ratio_max = 1.4
            if "EOPDistribution" in histogram:
                rebin = 20
                ratio_min = 0.6
                ratio_max = 1.4
            elif histogram == "HadFrac":
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
                                        ratio_min = ratio_min,\
                                        ratio_max = ratio_max,\
                                        doLogy=False,\
                                        marker_size = 1.0,\
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
