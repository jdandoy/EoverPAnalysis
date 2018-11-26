# EoverPPlotting

These plotting macros use root_numpy, cython and a module called atlas-plot. To keep the same root version and environment consistent, always run these macros after setting up AnalysisBase,21.2.23

## Setup
```
asetup AnalysisBase,21.2.23
pip --install cython --user
pip --install root_numpy --user
git clone git@github.com:joeycarter/atlas-plots.git
cd atlas-plots
pip install -e . --user
cd ..
print_env batchPlottingSubmission/env.sh
sed '/PWD/d' batchPlottingSubmission/env.sh
```

## When logging back in
```
asetup AnalysisBase,21.2.23
```

## Submit batch plotting macros
```
python batchPlottingSubmission/PrepareSubmission.py --treeName NAME --NPartitions 100 --jobName testPlotting
source submit_testPlotting.sh > testPlotting.log
```

## Check for job completion
```
python batchPlottingSubmission/checkForJobCompletion.py --logFile testPlotting.log
```

## Upon job completion
```
cd Output/testPlotting
hadd testPlotting_hadded.root *.root
```

#draw histograms
```
python
from CreatePlots import CreatePlots
CreatePlots("Outputs/testPlotting/testPlotting_hadded.root")
```
