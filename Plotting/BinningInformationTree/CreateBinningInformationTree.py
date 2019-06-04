import ROOT
from os import path as os_path
from sys import path as sys_path
from array import array

Curr_DIR = os_path.expandvars('$EOPPlottingDir')
sys_path.insert(1, Curr_DIR)

from PlottingTools.HistogramManager import HistogramManager
import argparse
parser = argparse.ArgumentParser(description='Create a root file that has a histogram required to reweight events')
parser.add_argument('--file', '-f', dest="file", type=str, required=True, help='Create a histogram required for reweighting')
parser.add_argument('--selectionName', '-selname', dest="selname", type=str, required=True, help='The name of the selection')
parser.add_argument('--histogramName', '-hm', dest="hist_name", type=str, required=True, help='The name of the histogram to use for reweighting')
parser.add_argument('--outputFile', '-on', dest="output_name", type=str, required=True, help='the name of the output file. the tree inside will have the same name')

args = parser.parse_args()

filename = args.file
selname = args.selname
output_name = args.output_name
hist_name = args.hist_name

hm = HistogramManager(filename)

f=ROOT.TFile(filename, "READ")
tree = f.Get(selname + "BinningTree")
for bins in tree:
    break

eta_bins_low = getattr(bins, selname+"EtaBinsLow")
eta_bins_high = getattr(bins, selname+"EtaBinsHigh")
p_bins_low_for_eta_bin = []
p_bins_high_for_eta_bin = []

#get all of the binning information that we need
for i in range(0, eta_bins_low.size()):
    p_bins_low_for_eta_bin.append(getattr(bins, selname+"PBinsLow_Eta"+str(i)))
    p_bins_high_for_eta_bin.append(getattr(bins, selname+"PBinsHigh_Eta"+str(i)))


#create a new output file
if not ".root" in output_name:
    output_name += ".root"

out_file = ROOT.TFile("BinningInformationTree/" + output_name, "RECREATE")
for j, eta_bin in enumerate(eta_bins_low):
    for m, p_bin in enumerate(p_bins_low_for_eta_bin[j]):
        name = selname+"_Eta_"+str(j) + "_PBin_" + str(m)
        histogram_name = hist_name + "_" + selname + "_Eta_" + str(i) + "_Momentum_" + str(j)
        histograms = hm.getHistograms(histogram_name)
        for channel in histograms:
             out_tree = ROOT.TTree(channel + "_" + name, channel + "_" + name)
             hist = histograms[channel]
             n = hist.Integral()
             sigma = hist.GetRMS()

             n = array( 'f', [n] )
             sigma = array( 'f', [sigma] )

             out_tree.Branch("NTracks", n, "NTracks/F")
             out_tree.Branch("Sigma", sigma, "Sigma/F")
             out_tree.Fill()
             out_file.cd()
             out_tree.Write()
out_file.Close()

