#!/bin/bash
source /afs/cern.ch/user/j/jdandoy/work/Ep/TreeMaker/source/EoverPAnalysis/Plotting/venv_EOPPlotting/bin/activate
#!!source ./setup_condor.sh
printf "Start time: "; /bin/date
printf "Job is running on node: "; /bin/hostname
printf "Job running as user: "; /usr/bin/id
printf "Job is running in directory: "; /bin/pwd
ls -al
python EoP_simple_plotter.py --job_index ${1}
printf "Finished!"
