#First get the number of entries in each of the trees in inputs
from sys import path as sys_path
from os import path as os_path
Curr_DIR = os_path.expandvars('$PWD')
sys_path.insert(1, Curr_DIR)
import os
from inputs.samples import INPUT
from ClusterPlots import FillingScript
from PlottingTools.Plotter import Plotter
from variables.variables import calc_weight
import ROOT
import pickle

import argparse
parser = argparse.ArgumentParser(description='Submit plotting batch jobs for the EoverPAnalysis plotting')
parser.add_argument('--treeName', '-tn', dest="treeName", type=str, required=True, help='the name of the tree to read from')
parser.add_argument('--NPartitions', '-np', dest="NPartitions", type=int, default='""', help='the number of plotting jobs to submit')
parser.add_argument('--jobName', '-jobName', dest="jobName", type=str, default='""', help='the name of the job to be submitted')

args = parser.parse_args()

#Create a pickle file and list for each submission
treeName = args.treeName
jobName = args.jobName
NPartitions = args.NPartitions

submission_list = []
EntriesPerFile = {}

for key in INPUT:
    filenames = INPUT[key][1]
    for filename in filenames:
        tfile = ROOT.TFile(filename, "READ")
        tree = tfile.Get(treeName)
        EntriesPerFile[filename] = tree.GetEntries()

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

#make the output directory
if not os.path.exists("Outputs"):
    os.makedirs("Outputs")

if not os.path.exists("Outputs/"+jobName):
    os.makedirs("Outputs/"+jobName)

outputDir = "Outputs/"+jobName

submission_pickle_file = outputDir+ "/" + jobName + ".pickle"
partitions = GenerateListOfPartitions(EntriesPerFile, NPartitions)
leading_script = file("submit_" + jobName + ".sh", "w")

cwd = os.getcwd()
#Create a plotter for each partition, and also a submission script:
for i in range(0, len(partitions)):
    partition = partitions[i]
    plots = Plotter(INPUT, treeName, calc_weight, base_selections = "", partition_dictionary = partition)
    submission_list.append(plots)
    file = open(outputDir + "/submission" + str(i) + ".sh", "w")
    file.write("cd " + cwd + "\n")
    file.write("source batchPlittingSubmission/env.sh\n")
    file.write("python batchPlottingSubmission/submit.py --num " + str(i) + " --picklefile " + submission_pickle_file + " --jobname " + jobName + "\n")
    #file.write("cp *.root " + cwd + "\n")
    leading_script.write("bsub -q 8nm -J sub" + str(i) +  " < " + outputDir + "/submission" + str(i) + ".sh" + "\n")

#create a pickle file for each submission
pickle.dump( submission_list, open(submission_pickle_file, "wb" ) )


