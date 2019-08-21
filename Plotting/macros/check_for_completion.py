##This is a script that checks for job completion
import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description='check that all of the plotting jobs have finished running')
parser.add_argument('--jobDir', '-jd', dest="jobDir", type=str, required=True, help='where to look for the jobs that ran')
parser.add_argument('--resub', '-r', dest="resub", required=False, action='store_true', help="resubmit the jobs")

args = parser.parse_args()
jobDir = args.jobDir

outputFiles = []
for root, dirs, files in os.walk(jobDir+"/Output/", topdown=False):
   for name in files:
      outputFiles.append(os.path.join(root, name))

#First check for the existence of LSF folders in the directory
not_finished_count = 0
#Root File Exists Count
root_file_count = 0
not_finished = []
allJobsDone = True
for outputFile in outputFiles:
    #get the STDOUT file

    found_finished = False
    with open(outputFile) as f:
        lines = f.readlines()
        for line in lines:
            if "FINISHED" in line:
                found_finished = True
            if "killed" in line:
                print(outputFile + " killed prematurely")

    if not found_finished:
        allJobsDone = False
        not_finished_count += 1
        filename = outputFile.split("/")[-1]
        job_id = int(filename.split('.')[-1])
        not_finished.append(job_id)
        print("the job with ID " + outputFile.split(".")[-1] + " didn't finish")

    job_number = outputFile.split(".")[-1]
    has_file = os.path.isfile(jobDir + "_" + job_number + ".root")
    if not has_file:
        root_file_count += 1

if allJobsDone:
    print("all plotting jobs have finished")
else:
    print("This many jobs didn't finish " + str(not_finished_count))
    print("This many root files don't exist " + str(root_file_count))

if args.resub:
    jobName = jobDir

    leading_script = file("condor_" + jobName + ".resub", "w")
    leading_script.write("Universe = vanilla\n")
    leading_script.write("Executable = condorSubmission/plot.sh\n")

    #create the output directories for then job
    if not os.path.exists(jobName):
        os.makedirs(jobName)
    if not os.path.exists(jobName+"/Output"):
        os.makedirs(jobName+"/Output")
    if not os.path.exists(jobName+"/Error"):
        os.makedirs(jobName+"/Error")
    if not os.path.exists(jobName+"/Log"):
        os.makedirs(jobName+"/Log")
    if not os.path.exists(jobName+"/Submission"):
        os.makedirs(jobName+"/Submission")
    submission_pickle_file = jobName + "/Submission/" + jobName + ".pickle"

    leading_script.write("+ProjectName='atlas-eopplotting'\n")
    leading_script.write('+JobFlavour = "testmatch"\n')
    leading_script.write("should_transfer_files = YES\n")
    leading_script.write("when_to_transfer_output = ON_Exit\n")
    leading_script.write("transfer_output         = True\n")
    leading_script.write("transfer_input_files    = CondorPythonLocal, PlottingTools,ReweightingHistograms, variables, HistogramFillingTools, selections, condorSubmission/submit.py ,calculation, FillingScript.py, " + submission_pickle_file + "\n")
    leading_script.write("\n")

    for jobID in not_finished:
        jobIDStr = str(jobID)
        leading_script.write("transfer_output_files   = " + jobName + "_"+jobIDStr+".root\n")
        leading_script.write("Error = " +jobName + "/Error/job."+jobIDStr+"\n")
        leading_script.write("Output = " +jobName + "/Output/job."+jobIDStr+"\n")
        leading_script.write("Log = "+jobName+"/Log/job."+jobIDStr+"\n")
        leading_script.write("Arguments = "+ jobIDStr +" "  +  submission_pickle_file.split("/")[-1] + " " + jobName + "\n")
        leading_script.write("Queue 1\n")
        leading_script.write("\n")


