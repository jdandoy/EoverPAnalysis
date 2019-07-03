from argparse import ArgumentParser
import os
parser = ArgumentParser(description='Tell me where you want them plots!')
parser.add_argument("--loc", type=str, dest="loc", required=True)
args = parser.parse_args()

plot_dir = os.path.join(os.getenv("EOPPlottingDir"), "Plots", "PTSpectrumReweightedplots")
plot_list = ["EtaLess08_TwoDHistTrkPvsPhiInnerToExtrapolEM2LowMuData.png"]

for i in range(0, 5):
    plot_list += ["EOPDistribution_Inclusive_Eta_{}.png".format(i)]
    plot_list += ["EnergyBkgProfileVsMomentum__MIPSelectionHadFracAbove70_Eta_{}_MIPSelectionHadFracAbove70.png".format(i)]
    plot_list += ["EnergyProfileVsMomentum__MIPSelectionHadFracAbove70_Eta_{}_MIPSelectionHadFracAbove70.png".format(i)]
    plot_list += ["NonZeroFractionNonZeroEnergyInclusiveTrackPSpectrum_{}.png".format(i)]
    plot_list += ["NonZeroFraction20TRTHitsNonZeroEnergy20TRTHitsTrackPSpectrum_{}.png".format(i)]

for plot in plot_list:
    plot = os.path.join(plot_dir, plot)
    command = "cp {} {}".format(plot, args.loc)
    print("calling {}".format(command))
    os.system(command)

