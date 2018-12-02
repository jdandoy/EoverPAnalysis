##This is a script that checks for job completion
import os
import argparse

parser = argparse.ArgumentParser(description='check that all of the plotting jobs have finished running')
parser.add_argument('--logFile', '-tn', dest="logFile", type=str, required=True, help='the input log file from the submission')

args = parser.parse_args()
logFile = args.logFile

with open(logFile) as f:
    content = f.readlines()

content = [x.strip() for x in content] #get the lines
content = [x for x in content if x[0:3] == "Job"] #only get the lines that describe a job ID
split_content = [x.split(" ") for x in content] #split each line by spaces
jobIDs_with_carrots = [x[1] for x in split_content] #this is the job ID in the form <NUMBER>
jobIDs = [x.lstrip("<").rstrip(">") for x in jobIDs_with_carrots] #remove the carrots around the ID

print "checking for the completion of jobs " + str(jobIDs)

#First check for the existence of LSF folders in the directory
job_count = -1
allJobsDone = True
for jobID in jobIDs:
    exists = os.path.exists("LSFJOB_"+jobID)
    job_count += 1
    if not exists:
        print "Couldn't find a LSF submission folder for " + jobID + " corresponding to " + str(job_count)
        allJobsDone = False
        continue

    #get the STDOUT file
    found_finished = False
    with open("LSFJOB_" + jobID + "/STDOUT") as f:
        content = f.readlines()
        for line in content:
            if "Killed" in line:
                print("The job was killed prematurely")
                allJobsDone = False
            if "THEJOBFINISHED!" in line:
                found_finished = True

    if not found_finished:
        allJobsDone = False
        print("the job with ID " + jobID + " didn't finish! This corresponds to submission " + str(job_count))

if allJobsDone:
    print "all plotting jobs have finished"
