# EoverPPlotting

These plotting macros use root_numpy, cython and a module called atlas-plot. To keep the same root version and environment consistent, always run these macros after setting up AnalysisBase,21.2.23

## Setup
This creates local installations of root_numpy, cython and psutils. These packages are needed for these plotting macros.
```
source setup.sh
```

## When logging back in
```
asetup AnalysisBase,21.2.23
export PYTHONPATH=$HOME/CondorPythonLocal/lib/python2.7/site-packages:$PYTHONPATH
export EOPPlottingDir=$(pwd)
```

## Submit batch plotting macros
```
python condorSubmission/PrepareSubmission.py --treeName NAME --NPartitions 100 --jobName testPlotting
condor_submit condor_testPlotting.sub
```

## Upon job completion
```
hadd testPlotting_hadded.root testPlotting*.root
```

## Draw histograms
```
git clone git@github.com:joeycarter/atlas-plots.git
cd atlas-plots
pip install -e . --user
cd ..
python
from CreatePlots import CreatePlots
CreatePlots("Outputs/testPlotting/testPlotting_hadded.root")
```
