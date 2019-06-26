from argparse import ArgumentParser
import os
parser = ArgumentParser(description='Tell me where you want them plots!')
parser.add_argument("--loc", type=str, dest="loc", required=True)
args = parser.parse_args()

plot_dir = os.path.join(os.getenv("EOPPlottingDir"), "Plots", "PTSpectrumReweightedplots")
plot_list = ["EtaLess08_TwoDHistTrkPvsPhiInnerToExtrapolEM2LowMuData.png"]

for plot in plot_list:
    plot = os.path.join(plot_dir, plot)
    command = "cp {} {}".format(plot, args.loc)
    print("calling {}".format(command))
    os.system(command)

