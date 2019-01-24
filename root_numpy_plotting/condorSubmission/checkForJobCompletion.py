##This is a script that checks for job completion
import os
import argparse

parser = argparse.ArgumentParser(description='check that all of the plotting jobs have finished running')
parser.add_argument('--jobDir', '-jd', dest="jobDir", type=str, required=True, help='where to look for the jobs that ran')

args = parser.parse_args()
jobDir = args.jobDir

outputFiles = []
for root, dirs, files in os.walk(jobDir+"/Output/", topdown=False):
   for name in files:
      outputFiles.append(os.path.join(root, name))

#First check for the existence of LSF folders in the directory
not_finished_count = 0
allJobsDone = True
for outputFile in outputFiles:
    #get the STDOUT file
    found_finished = False
    with open(outputFile) as f:
        content = f.readlines()
        for line in content[::-1]:
            if "THEJOBFINISHED!" in line:
                found_finished = True

    if not found_finished:
        allJobsDone = False
        not_finished_count += 1
        print("the job with ID " + outputFile.split(".")[-1] + " didn't finish")

if allJobsDone:
    print "all plotting jobs have finished"
else:
    print "This many jobs didn't finish " + str(not_finished_count)
