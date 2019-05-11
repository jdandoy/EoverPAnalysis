#First get the number of entries in each of the trees in inputs
from sys import path as sys_path
from os import path as os_path
Curr_DIR = os_path.expandvars('$EOPPlottingDir')
sys_path.insert(1, Curr_DIR)
import os
from inputs.samples import INPUT
from HistogramFillingTools.HistogramFiller import HistogramFiller
from variables.variables import calc_weight
import ROOT
import pickle
import glob

import argparse
parser = argparse.ArgumentParser(description='Submit plotting batch jobs for the EoverPAnalysis plotting')
parser.add_argument('--treeName', '-tn', dest="treeName", type=str, required=True, help='the name of the tree to read from')
parser.add_argument('--NPartitions', '-np', dest="NPartitions", type=int, default='0', help='the number of plotting jobs to submit')
parser.add_argument('--jobName', '-jobName', dest="jobName", type=str, default='""', help='the name of the job to be submitted')
parser.add_argument('--jobFlavour', '-jobFlavour', dest="jobFlavour", type=str, default='tomorrow', help='What condor queue should the jobs run on?')

args = parser.parse_args()

#Create a pickle file and list for each submission
treeName = args.treeName
job_name = args.jobName
NPartitions = args.NPartitions
flavour = args.jobFlavour

submission_list = []
EntriesPerFile = {}
tree_dict = {}

for key in INPUT:
    filenames = INPUT[key][1]
    for filename in filenames:
        print "Adding files in " + filename
        isFile = os.path.isfile(filename)
        tree = ROOT.TChain(treeName)
        if not isFile:
            if filename[-1] != "/":
                search_for = filename + "/*.root"
            else:
                search_for = filename + "*.root*"
            files =  glob.glob(search_for)
            files.sort()
            for f in files:
                if "eos" in f:
                     tfile = ROOT.TFile.Open("root://eosatlas/" + f, "READ")
                else:
                     tfile = ROOT.TFile.Open(f, "READ")
                print "Adding "+f
                tree.Add(f)
                tfile.Close()
        else:
            tree.Add("root://eosatlas/" + filename)

        EntriesPerFile[filename] = tree.GetEntries()
        f = tree.GetCurrentFile()
        tree.SetDirectory(0)
        f.Close()

        tree_dict[filename] = tree

def GeneratePartitions(Entries, NPartitions):
    return_list = []
    if Entries <= NPartitions * 10:
        raise ValueError("The number of events per partition is so small. What is the point?")
    step = int(Entries/NPartitions) - 10

    #print(step)
    count = 1
    for i in range(0, NPartitions):
        return_list.append( ( (i * step), (i+1) * step) )
        count += 1
    return_list.append((return_list[-1][1],Entries))

    if return_list[-1][0] > return_list[-1][1]:
        raise ValueError("The order of the last partition doesn't makse sense")

    #print(len(return_list))

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
leading_script = file("condor_" + job_name + ".sub", "w")
cwd = os.getcwd()

leading_script.write("Universe = vanilla\n")
leading_script.write("Executable = condorSubmission/plot.sh\n")

#create the output directories for then job
if not os.path.exists(job_name):
    os.makedirs(job_name)
if not os.path.exists(job_name+"/Output"):
    os.makedirs(job_name+"/Output")
if not os.path.exists(job_name+"/Error"):
    os.makedirs(job_name+"/Error")
if not os.path.exists(job_name+"/Log"):
    os.makedirs(job_name+"/Log")
if not os.path.exists(job_name+"/Submission"):
    os.makedirs(job_name+"/Submission")
submission_pickle_file = job_name + "/Submission/" + job_name + ".pickle"

leading_script.write("Error = " +job_name + "/Error/job.$(Process)\n")
leading_script.write("Output = " +job_name + "/Output/job.$(Process)\n")
leading_script.write("Log = "+job_name+"/Log/job.$(Process)\n")
leading_script.write("+ProjectName='atlas-eopplotting'\n")
leading_script.write('+JobFlavour = ' + flavour + '\n')
leading_script.write("should_transfer_files = YES\n")
leading_script.write("when_to_transfer_output = ON_Exit\n")
leading_script.write("transfer_output         = True\n")
leading_script.write("transfer_input_files    = CondorPythonLocal, variables, HistorgamFillingTools, selections, condorSubmission/submit.py ,calculation, FillingScript.py, " + submission_pickle_file + "\n")
leading_script.write("transfer_output_files   = " + job_name + "_$(Process).root\n")
leading_script.write("\n")

for i in range(0, len(partitions)):
    partition = partitions[i]
    hist_filler = HistogramFiller(INPUT, treeName, calc_weight, base_selections = "", partition_dictionary = partition)
    hist_filler.tree_dict = tree_dict
    submission_list.append(hist_filler)
    leading_script.write("Arguments = $(Process) "  +  submission_pickle_file.split("/")[-1] + " " + job_name + "\n")
    leading_script.write("Queue 1\n")
    leading_script.write("\n")

#create a pickle file for each submission
pickle.dump( submission_list, open(submission_pickle_file, "wb" ) )
