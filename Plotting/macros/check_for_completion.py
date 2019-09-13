##This is a script that checks for job completion
import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description='check that all of the plotting jobs have finished running')
parser.add_argument('--job_dir', '-jd', dest="job_dir", type=str, required=True, help='where to look for the jobs that ran')
parser.add_argument('--resub', '-r', dest="resub", required=False, action='store_true', help="resubmit the jobs")

args = parser.parse_args()
job_dir = args.job_dir

outputFiles = []
for root, dirs, files in os.walk(job_dir+"/Output/", topdown=False):
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
    has_file = os.path.isfile(job_dir + "_" + job_number + ".root")
    if not has_file:
        root_file_count += 1

if allJobsDone:
    print("all plotting jobs have finished")
else:
    print("This many jobs didn't finish " + str(not_finished_count))
    print("This many root files don't exist " + str(root_file_count))

if args.resub:
    jobname = job_dir.rstrip("/").split("/")[-1]
    #create the output directories for then job
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)
    if not os.path.exists(job_dir+"/Output"):
        os.makedirs(job_dir+"/Output")
    if not os.path.exists(job_dir+"/Error"):
        os.makedirs(job_dir+"/Error")
    if not os.path.exists(job_dir+"/Log"):
        os.makedirs(job_dir+"/Log")
    if not os.path.exists(job_dir+"/Submission"):
        os.makedirs(job_dir+"/Submission")
    original_condor_file = os.path.join(job_dir, "{}_scripts".format(jobname), "condor_{}.sub".format(jobname))
    new_condor_file = os.path.join(job_dir, "{}_scripts".format(jobname), "condor_{}.resub".format(jobname))

    lines_to_keep = []
    lines_to_rewrite = []
    first_sub = False
    with open(original_condor_file, "r") as f:
        print("opened: {}".format(original_condor_file))
        lines = f.readlines()
        for i,l in enumerate(lines):
            if "Arguments = " in l:
                break
            if "$(Process)" in l:
                lines_to_rewrite.append(l)
                continue
            lines_to_keep.append(l)

        one_argument_line = None
        for i,l in enumerate(lines):
            if "Argument" in l:
                one_argument_line = l
                break
        for i in not_finished:
            lines_to_keep.append("\n")
            lines_to_keep.append(one_argument_line.replace("$(Process)", str(i)))
            for extra_l in lines_to_rewrite:
                lines_to_keep.append(extra_l.replace("$(Process)", str(i)))
            lines_to_keep.append("Queue 1")
            lines_to_keep.append("\n")

    with open(new_condor_file, 'w') as f:
        for l in lines_to_keep:
            print("Writing: {}".format(l))
            f.write(l)

    print("Finished")
    os.system("condor_submit {}".format(new_condor_file))


