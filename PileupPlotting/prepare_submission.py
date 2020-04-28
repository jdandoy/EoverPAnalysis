#First get the number of entries in each of the trees in inputs
import SampleSorter

#import utils
import os
#import ROOT
#import pickle

import argparse
parser = argparse.ArgumentParser(description='Submit plotting batch jobs for the EoverPAnalysis plotting')
parser.add_argument('--tree_name', type=str, default='LA_EoverP_InDetTrackParticlesSortedLooseIsolatedVertexAssociated_tree', help='The name of the input TTree.')
parser.add_argument('--njobs', type=int, default=100, help='The number of plotting jobs to submit.')
parser.add_argument('--job_name', type=str, default='testdir_1', help='The name of the job to be submitted.')
parser.add_argument('--queue_flavour', '-queue_flavour', dest="queue_flavour", type=str, default='longlunch', help='The condor queue to use.')
args = parser.parse_args()

top_dir = os.getcwd()+'/condor/'+args.job_name
dir_list = ['error', 'log', 'output','scripts', 'results']
for this_dir in dir_list:
  if not os.path.exists(top_dir+'/'+this_dir):
    os.makedirs(top_dir+'/'+this_dir)

#Make file lists for all the jobs and place in a pickle file in the scripts directory
file_lists = SampleSorter.make_sample_lists(args.njobs, args.tree_name)
SampleSorter.make_file_pickles(file_name='{0}/scripts/job_file_list.pickle'.format(top_dir), files_per_job=file_lists)
njobs = len(file_lists) #In case a different number of jobs was needed for the file list
#njobs = args.njobs

#Copy needed scripts to the new condor scripts directory
os.system('cp EoP_simple_plotter.py {0}/scripts/'.format(top_dir))
os.system('cp PileupHistogrammer.py {0}/scripts/'.format(top_dir))
os.system('cp execute_job.sh {0}/scripts/'.format(top_dir))

condor_sub = \
"""
Universe = vanilla
Executable = {0}/scripts/execute_job.sh
Error =  {0}/error/job.$(Process)
Output = {0}/output/job.$(Process)
Log =    {0}/log/job.$(Process)
+JobFlavour    = "{1}"
Request_memory = 10000
should_transfer_files   = YES
when_to_transfer_output = ON_Exit
transfer_output         = True
initialdir              = {0}/results/

transfer_input_files = {0}/scripts/EoP_simple_plotter.py,{0}/scripts/job_file_list.pickle,{0}/scripts/PileupHistogrammer.py

Arguments = $(Process)
Queue {2}
""".format(top_dir, args.queue_flavour, njobs)

print(condor_sub)
with open( '{0}/scripts/{1}.sub'.format(top_dir,args.job_name), "w") as submission_script:
  submission_script.write(condor_sub)

##Create a pickle file and list for each submission
#tree_name = args.tree_name
#job_name = args.job_name
#n_jobs = args.n_jobs
#flavour = args.queue_flavour
#file_flavour = args.file_flavour
#filling_script = args.filling_script
#condor_directories = ["condor", args.job_name]
#
#
#files = utils.get_files(file_flavour)
#trees = utils.tchain_files_together(args.tree_name, files)
#
#condor_directory = project_dir
#for path in condor_directories:
#    if not os.path.exists(os.path.join(condor_directory, path)):
#        os.makedirs(os.path.join(condor_directory, path))
#    condor_directory = os.path.join(condor_directory, path)
#
#submission_script_dir = os.path.join(project_dir, condor_directory, "{}_{}".format(job_name, "scripts"))
#if not os.path.exists(submission_script_dir):
#    os.makedirs(submission_script_dir)
#
#
##create the executables for the condor jobs
#executable = os.path.join(submission_script_dir, "plot.sh")
#executable_local = os.path.join(submission_script_dir, "plot_local.sh")
#python_executable = os.path.join(submission_script_dir, "plot.py")
#
##create the python script that is needed for plotting
#plotting_instructions_python = []
#plotting_instructions_python.append("from histogram_filling import HistogramFiller")
#plotting_instructions_python.append("from {} import fill_histograms".format(filling_script.split("/")[-1].split(".")[0]))
#plotting_instructions_python.append("import pickle")
#plotting_instructions_python.append("import argparse")
#plotting_instructions_python.append("parser = argparse.ArgumentParser(description=\'Submit plotting batch jobs for the EoverPAnalysis plotting\')")
#plotting_instructions_python.append("parser.add_argument(\'--num\', '-n', dest=\"num\", type=int, required=True, help=\'Which submission number was this?\')")
#plotting_instructions_python.append("parser.add_argument('--picklefile' '-p', dest='picklefile', type=str, default=\"\", help='Where to get the plotter')")
#plotting_instructions_python.append("parser.add_argument('--jobName', '-j', dest=\"jobname\", type=str, default='\"\"', help='the names of the batch jobs')")
#plotting_instructions_python.append("args = parser.parse_args()")
#plotting_instructions_python.append("i = args.num")
#plotting_instructions_python.append("file = args.picklefile")
#plotting_instructions_python.append("name = args.jobname")
#plotting_instructions_python.append("plots = pickle.load(open(file, \"rb\"))[i]")
#plotting_instructions_python.append("fill_histograms(plots, name + \"_\" + str(i) + \".root\")")
#with open(python_executable, 'w') as f:
#    for line in plotting_instructions_python:
#        f.write(line+"\n")
#
##create the shell script to be executed
#plotting_instruction_script = []
#plotting_instruction_script.append("#!/bin/bash")
#activate_location = os.path.join(os.getenv("EOPPlottingDir"),"venv_EOPPlotting/bin/activate")
#plotting_instruction_script.append("source {}".format(activate_location))
#plotting_instruction_script.append("source ./setup_condor.sh")
#plotting_instruction_script.append("printf \"Start time: \"; /bin/date")
#plotting_instruction_script.append("printf \"Job is running on node: \"; /bin/hostname")
#plotting_instruction_script.append("printf \"Job running as user: \"; /usr/bin/id")
#plotting_instruction_script.append("printf \"Job is running in directory: \"; /bin/pwd")
#plotting_instruction_script.append("ls -al")
#plotting_instruction_script.append("python {} ".format("plot.py") + " --num ${1} --picklefile ${2} --jobName ${3}")
#with open(executable, 'w') as f:
#    for line in plotting_instruction_script:
#        f.write(line+"\n")
#
##create the shell script that runs jobs locally
#plotting_instruction_local_script = ["cp {} {}".format(filling_script, filling_script.split("/")[-1])]
#plotting_instruction_local_script += ["cp {} {}".format(python_executable, python_executable.split("/")[-1])]
#plotting_instruction_local_script +=  plotting_instruction_script 
#plotting_instruction_local_script += ["rm {}".format(filling_script.split("/")[-1])]
#plotting_instruction_local_script += ["rm {}".format(python_executable.split("/")[-1])]
#with open(executable_local, 'w') as f:
#    for line in plotting_instruction_local_script:
#        f.write(line+"\n")
#
#with open(os.path.join(submission_script_dir,"condor_{}.sub".format(job_name)),"w") as leading_script:
#    leading_script.write("Universe = vanilla\n")
#    leading_script.write("Executable = {}\n".format(executable))
#    #create the output directories for then job
#    if not os.path.exists(condor_directory+"/Output"):
#        os.makedirs(condor_directory+"/Output")
#    if not os.path.exists(condor_directory+"/Error"):
#        os.makedirs(condor_directory+"/Error")
#    if not os.path.exists(condor_directory+"/Log"):
#        os.makedirs(condor_directory+"/Log")
#    if not os.path.exists(condor_directory+"/Submission"):
#        os.makedirs(condor_directory+"/Submission")
#    submission_pickle_file = os.path.join(condor_directory, "Submission", "{}.pickle".format(job_name))
#
#    leading_script.write("Error = " +condor_directory + "/Error/job.$(Process)\n")
#    leading_script.write("Output = " +condor_directory + "/Output/job.$(Process)\n")
#    leading_script.write("Log = "+condor_directory+"/Log/job.$(Process)\n")
#    leading_script.write('+JobFlavour = "' + flavour + '"\n')
#    leading_script.write('Request_memory = 10000\n')
#    leading_script.write("should_transfer_files = YES\n")
#    leading_script.write("when_to_transfer_output = ON_Exit\n")
#    leading_script.write("transfer_output         = True\n")
#    leading_script.write("transfer_input_files    = {rw},{eop},{u},{p},{py},{setup},{fs},{bin},{cert}\n"\
#            .format(\
#            rw=os.path.join(project_dir, "ReweightingHistograms"),\
#            eop=os.path.join(project_dir,"eop_plotting"),\
#            u=os.path.join(project_dir,"utils"),\
#            p=os.path.join(submission_pickle_file),\
#            py=python_executable,\
#            setup=os.path.join(project_dir, "setup_condor.sh"),\
#            fs=os.path.abspath(args.filling_script),\
#            bin=os.path.join(project_dir,"bin"),\
#            cert = os.path.join(os.getenv("EOPPlottingDir"), "grid_proxy")))
#    leading_script.write("transfer_output_files   = " + job_name + "_$(Process).root\n")
#    leading_script.write('transfer_output_remaps = "{} = {}"\n'.format(job_name + "_$(Process).root" , os.path.join(condor_directory, job_name + "_$(Process).root") ) )
#    leading_script.write("\n")
#
#    #a list of plotters for each of the jobs
#    submission_list = []
#    for i in range(0, n_jobs):
#        partition = {}
#        for channel in partitions:
#            partition[channel] = {}
#            for f in partitions[channel]:
#                assert  len(partitions[channel][f]) == n_jobs
#                partition[channel][f] =  partitions[channel][f][i]
#        hist_filler = HistogramFiller(trees, tree_name, calc_weight, selection_string = "", partitions = partition)
#        submission_list.append(hist_filler)
#        leading_script.write("Arguments = $(Process) "  + submission_pickle_file.split("/")[-1] + " " + job_name + "\n")
#        leading_script.write("Queue 1\n")
#        leading_script.write("\n")
#
##create a pickle file for each submission
#pickle.dump( submission_list, open(submission_pickle_file, "wb" ) )
#print("Created the submission file. Ready to go!")
#os._exit(0)
