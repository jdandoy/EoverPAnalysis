# EoverPPlotting

These plotting macros use root_numpy, cython and a module called atlas-plots. They ship jobs to condor for all of your plotting needs. The histograms are filled in the Fillingscript.py file.

## Setup
This creates local installations of root_numpy, cython and psutils. These packages are needed for these plotting macros.
```
source setup.sh
```

## When logging back in
```
source login.sh
```

## Prepare batch plotting jobs for submission
```
python condorSubmission/PrepareSubmission.py --treeName LA_EoverP_InDetTrackParticlesSortedLooseIsolatedVertexAssociated_tree --jobName Plots --NPartitions 200
```

## Test one of the jobs locally
```
python condorSubmission/submit.py --num 0 --picklefile Plots/Submission/Plots.pickle --jobName Plots
```

## Submit all of the batch jobs
```
condor_submit condor_testPlots.sub
```

## Upon job completion
```
hadd Plots_hadded.root Plots*.root
```

## Draw histograms
```
git clone git@github.com:joeycarter/atlas-plots.git
cd atlas-plots
pip install -e . --user
cd ..
python
from CreatePlots import CreatePlots
CreatePlots("Outputs/Plots/Plots_hadded.root")
```
