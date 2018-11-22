#First get the number of entries in each of the trees in inputs
from inputs.samples import INPUT
from ClusterPlots import FillingScript
from PlottingTools.Plotter import Plotter
from variables.variables import calc_weight
import ROOT

EntriesPerFile = {}
treeName = "EoverP_ClusterEnergyInDetTrackParticlesLooseIsolatedVertexAssociated_tree"

for key in INPUT:
    filenames = INPUT[key][1]
    for filename in filenames:
        file = ROOT.TFile(filename, "READ")
        tree = file.Get(treeName)
        EntriesPerFile[filename] = tree.GetEntries()

NPartitions = 3
print EntriesPerFile

def GeneratePartitions(Entries, NPartitions):
    return_list = []
    if Entries <= NPartitions * 10:
        raise ValueError("The number of events per partition is so small. What is the point?")
    step = int(Entries/NPartitions) - 10

    print(step)
    count = 1
    for i in range(0, NPartitions):
        return_list.append( ( (i * step), (i+1) * step) )
        count += 1
    return_list.append((return_list[-1][1],Entries))

    if return_list[-1][0] > return_list[-1][1]:
        raise ValueError("The order of the last partition doesn't makse sense")

    print(len(return_list))

    if len(return_list) != NPartitions + 1:
        raise ValueError("The number of partitions for this file was not NPartions = " + str(NPartitions))

    return return_list

def GenerateListOfPartitions(EntriesPerFile, NPartitions):
    return_list = []
    for i in range(0, NPartitions + 1):
        return_list.append({})

    for key in EntriesPerFile:
        partitions_for_file = GeneratePartitions(EntriesPerFile[key], NPartitions)
        for partition, dictionary in zip(partitions_for_file, return_list):
            dictionary[key] = partition

    return return_list

partitions = GenerateListOfPartitions(EntriesPerFile, NPartitions)

#Create a plotter for each partition:
for i in range(0, NPartitions):
    partition = partitions[i]
    plots = Plotter(INPUT, treeName, calc_weight, base_selections = "", partition_dictionary = partition)
    outfile = "testing" + str(i) + ".root"
    FillingScript(plots, outfile)

