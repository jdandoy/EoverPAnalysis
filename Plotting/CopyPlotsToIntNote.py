from argparse import ArgumentParser
import os
parser = ArgumentParser(description='Tell me where you want them plots!')
parser.add_argument("--loc", type=str, dest="loc", required=True)
parser.add_argument("--flavour", type=str, dest="flavour", required=True)
args = parser.parse_args()


if args.flavour == "reweighted":
    plot_dir = os.path.join(os.getenv("EOPPlottingDir"), "Plots", "PTSpectrumReweightedplots")
    plot_list = ["EtaLess08_TwoDHistTrkPvsPhiInnerToExtrapolEM2LowMuData.png"]

    for i in range(0, 5):
        plot_list += ["EOPDistribution_Inclusive_Eta_{}.png".format(i)]
        plot_list += ["EnergyBkgProfileVsMomentum__MIPSelectionHadFracAbove70_Eta_{}_MIPSelectionHadFracAbove70.png".format(i)]
        plot_list += ["EOPProfileVsMomentum__MIPSelectionHadFracAbove70_Eta_{}_MIPSelectionHadFracAbove70.png".format(i)]
        plot_list += ["NonZeroFractionNonZeroEnergyInclusiveTrackPSpectrum_{}.png".format(i)]
        plot_list += ["NonZeroFraction20TRTHitsNonZeroEnergy20TRTHitsTrackPSpectrum_{}.png".format(i)]

    for plot in plot_list:
        plot = os.path.join(plot_dir, plot)
        command = "cp {} {}".format(plot, args.loc)
        print("calling {}".format(command))
        os.system(command)

elif args.flavour == "event_count":
    plot_dir = os.path.join(os.getenv("EOPPlottingDir"), "Plots", "count_reweight_plotsplots")
    plot_list = ["trkPtHist.png", "LeadingPtTrkHist.png", "eventNPV2Hist.png"]
    base = "FractionalComposition_{}.png"
    for i in range(0, 5):
        plot_list.append(base.format(i))

    base = "FractionalComposition_nostack_{}.png"
    for i in range(0, 5):
        plot_list.append(base.format(i))

    base = "TrackPtSpectrum__Inclusive_Eta_{}_Inclusive.png"
    for i in range(0, 5):
        plot_list.append(base.format(i))

    for plot in plot_list:
        plot = os.path.join(plot_dir, plot)
        command = "cp {} {}".format(plot, args.loc)
        print("calling {}".format(command))
        os.system(command)

elif args.flavour == "npv_reweight":
    plot_dir = os.path.join(os.getenv("EOPPlottingDir"), "Plots", "event_npv2_reweightedplots")
    plot_list = ["trkPtHist.png", "LeadingPtTrkHist.png", "eventNPV2Hist.png"]
    base = "TrackPtSpectrum__Inclusive_Eta_{}_Inclusive.png"
    for i in range(0, 5):
        plot_list.append(base.format(i))
    for plot in plot_list:
        plot = os.path.join(plot_dir, plot)
        command = "cp {} {}".format(plot, os.path.join(args.loc,"{}_{}.png".format(plot.split("/")[-1].rstrip(".png"),"npv_reweight")))
        print("calling {}".format(command))
        os.system(command)

elif args.flavour == "pt_reweight":
    plot_dir = os.path.join(os.getenv("EOPPlottingDir"), "Plots", "pt_reweightedplots")
    plot_list = ["trkPtHist.png", "LeadingPtTrkHist.png", "eventNPV2Hist.png"]
    base = "TrackPtSpectrum__Inclusive_Eta_{}_Inclusive.png"
    for i in range(0, 5):
        plot_list.append(base.format(i))
    base = "TrackPtSpectrum__Inclusive_Eta_{}_Inclusive.png"
    for i in range(0, 5):
        plot_list.append(base.format(i))
    for plot in plot_list:
        plot = os.path.join(plot_dir, plot)
        command = "cp {} {}".format(plot, os.path.join(args.loc,"{}_{}.png".format(plot.split("/")[-1].rstrip(".png"),"pt_reweight")))
        print("calling {}".format(command))
        os.system(command)
