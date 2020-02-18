import histogram_manager
import ROOT
import argparse

def get_binning(tf, selection):
    tree = tf.Get(selection + "BinningTree")
    for bins in tree:
        break

    eta_bins_low = getattr(bins, selection+"EtaBinsLow")
    eta_bins_high = getattr(bins, selection+"EtaBinsHigh")

    p_bins_low_for_eta_bin = []
    p_bins_high_for_eta_bin = []

    for i in range(0, eta_bins_low.size()):
        p_bins_low_for_eta_bin.append(getattr(bins, selection+"PBinsLow_Eta"+str(i)))
        p_bins_high_for_eta_bin.append(getattr(bins, selection+"PBinsHigh_Eta"+str(i)))
    return {"eta":(eta_bins_low, eta_bins_high), "momentum":(p_bins_low_for_eta_bin, p_bins_high_for_eta_bin)}

def get_histograms(tf, binnings, hist_manager, base_name = "EOPDistribution", selection = ""):
    histograms = {}
    eta_bins = binnings["eta"]
    p_bins = binnings["momentum"]
    for eta_count, (eta_low, eta_high, p_bins_low, p_bins_high) in enumerate(zip(eta_bins[0], eta_bins[1], p_bins[0], p_bins[1])):
        for p_count, (p_low, p_high) in enumerate(zip(p_bins_low, p_bins_high)):
            histogram_name = base_name + "_" + selection + "_Eta_" + str(eta_count) + "_Momentum_" + str(p_count)
            print("retrieving {}".format(histogram_name))
            histograms[histogram_name] ={}
            histograms[histogram_name]["histograms"] = hist_manager.getHistograms(histogram_name)
            histograms[histogram_name]["eta_bins"] = (eta_low, eta_high)
            histograms[histogram_name]["p_bins"] = (p_bins_low, p_bins_high)
            histograms[histogram_name]["p_count"] = p_count
            histograms[histogram_name]["eta_count"] = eta_count
    return histograms

parser = argparse.ArgumentParser()
parser.add_argument("--root_file", dest="root_file", required=False, default="inclusive.root")
parser.add_argument("--selection", dest="selection", required=False, default="NonZeroEnergy")
args = parser.parse_args()
root_file = args.root_file
selection = args.selection

nsigma = 1.1
fitter = ROOT.JES_BalanceFitter(nsigma)
hm = histogram_manager.HistogramManager(args.root_file)
tf = ROOT.TFile(args.root_file, "READ")
binnings = get_binning(tf, selection)
histograms = get_histograms(tf, binnings, hm, "EOPDistribution", selection = args.selection)

for h in histograms:
    canvas = ROOT.TCanvas(h, h)
    needed_for_fit = histograms[h]
    these_histograms = needed_for_fit["histograms"]
    eta_bins = needed_for_fit["eta_bins"]
    p_bins = needed_for_fit["p_bins"]
    p_count = needed_for_fit["p_count"]
    eta_count = needed_for_fit["eta_count"]
    for process in these_histograms:
        to_fit = these_histograms[process]
        fitter.FitAndDraw(to_fit,0.0);
        canvas.Print("{}.eps".format(to_fit.GetName()))
