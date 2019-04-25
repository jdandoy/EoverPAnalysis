#!/bin/bash

#NEEDED
export HOME=$(pwd)
export PROOFANADIR=$(pwd)
#ROOT STUFF
export EOPPlottingDir=$(pwd)
export PYTHONPATH=$PWD/CondorPythonLocal/lib64/python2.7/site-packages:$PYTHONPATH

printf "Start time: "; /bin/date
printf "Job is running on node: "; /bin/hostname
printf "Job running as user: "; /usr/bin/id
printf "Job is running in directory: "; /bin/pwd

python submit.py --num ${1} --picklefile ${2} --jobName ${3}
