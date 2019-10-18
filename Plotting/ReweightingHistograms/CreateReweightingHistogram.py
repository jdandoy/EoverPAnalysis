import ROOT
from os import path as os_path
from sys import path as sys_path

Curr_DIR = os_path.expandvars('$EOPPlottingDir')
sys_path.insert(1, Curr_DIR)

from histogram_manager import HistogramManager
import argparse
parser = argparse.ArgumentParser(description='Create a root file that has a histogram required to reweight events')
parser.add_argument('--file', '-f', dest="file", type=str, required=True, help='Create a histogram required for reweighting')
parser.add_argument('--numerator', '-num', dest="num", type=str, required=True, help='The channel to reweight to')
parser.add_argument('--denomenator', '-den', dest="den", type=str, required=True, help='The channel to be reweighted')
parser.add_argument('--rebin_x', '-rbx', dest="rbx", type=int, required=False, default=0, help='The rebinning to be applied to the x-axis')
parser.add_argument('--rebin_y', '-rby', dest="rby", type=int, required=False, default=0, help='The rebinning to be applied to the y-axis')
parser.add_argument('--histogramName', '-hm', dest="hist_name", type=str, required=True, help='The name of the histogram to use for reweighting')
parser.add_argument('--outputName', '-on', dest="output_name", type=str, required=True, help='the name of the root file and histogram.')

args = parser.parse_args()

filename = args.file
num_channel = args.num
den_channel = args.den
output_name = args.output_name
hist_name = args.hist_name

hm = HistogramManager(filename)
hm.listHistograms()

histograms = hm.getHistograms(hist_name)

if args.rbx != 0:
    for hist in histograms:
        histograms[hist].RebinX(args.rbx)
if args.rby != 0:
    for hist in histograms:
        histograms[hist].RebinY(args.rby)

#numerator:
ratio_hist = histograms[num_channel].Clone(output_name)
ratio_hist.Divide(histograms[den_channel])

#create a new output file
if not ".root" in output_name:
    output_name += ".root"

out_file = ROOT.TFile("ReweightingHistograms/" + output_name, "RECREATE")
ratio_hist.Write()
out_file.Close()

