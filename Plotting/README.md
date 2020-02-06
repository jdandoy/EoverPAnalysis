# EoverPPlotting

These macros ship jobs to condor for all of your plotting needs. The histograms are filled in the Fillingscript.py file.

## Setup
Create the environment needed for plotting. This package uses a python virtual environment (https://docs.python.org/3/tutorial/venv.html), in conjuction with an lcg virtual environment (https://ep-dep-sft.web.cern.ch/document/lcg-releases). You must be on a fresh login of lxplus (centos7) for everything to work. The setup will not work if an analysisbase release has be set up, for example. Running this package on any computing site with condor and cvmfs should (in principle) work. 

You must set the CMT config to match your operating system. "CMTCONFIG=x86_64-centos7-gcc8-opt" for centos7, and "CMTCONFIG=x86_64-slc6-gcc8-opt" for slc6
```
cd EoverPAnalysis/Plotting
git clone ssh://git@gitlab.cern.ch:7999/cburr/lcg_virtualenv.git ./lcg_virtualenv
CMTCONFIG=x86_64-centos7-gcc8-opt ./lcg_virtualenv/create_lcg_virtualenv venv_EOPPlotting LCG_94python3
source venv_EOPPlotting/bin/activate
```

Install all of the dependencies. This package uses uproot (https://uproot.readthedocs.io/en/latest/), along with python bindings for XRootD (https://github.com/xrootd/xrootd).
```
pip install --upgrade pip
git clone https://github.com/joeycarter/atlas-plots.git
cd atlas-plots
pip install -e .
cd ..
git clone https://github.com/xrootd/xrootd.git
cd xrootd/bindings/python/
python setup_standalone.py install
cd ../../../
pip install uproot
```

Finish the setup by running the setup script
```
source ./setup.sh
```

Every few months, you will need to generate a grid proxy to access files on eos. You can do this with the following lines:
```
voms-proxy-init --voms atlas --hours 10000
```

## When logging back in
```
cd $TestArea/EoverPAnalysis/Plotting
source ./setup.sh
```

## Package Philosophy
This script prepares a condor job. The job is defined in the macros/fill_script.py, which must have a function called fill_histograms. fill_histograms takes a HistogramFiller and books many histograms for plotting. Take a look inside of the script fill_scipt.py to get an idea of how this works.

All selections and variables to be plotted are defined as instances of the class "Caclulation". To define a new selection for the number of tracks, you could write the following function:
```
from calculation import Calculation
def TightIso(trk):
    return trk["trk_nearest_dR_EM"] > 0.55
sel_TightIso = Calculation(TightIso, ["trk_nearest_dR_EM"])
```
The initialization of Calculation takes a function that calculates the selection (trk_nearest_dR_EM > 0.55), and a list of branches needed to perform the calculation (trk_nearest_dR_EM).

Similarly, one can create a function that calculates E/P
```
from calculation import Calculation
def EOP(trk):
    return (trk["trk_ClusterEnergy_EM_200"] + trk["trk_ClusterEnergy_HAD_200"])/trk["trk_p"]
branches = ["trk_ClusterEnergy_EM_200", "trk_ClusterEnergy_HAD_200", "trk_p"]
calc_EOP = Calculation(EOP, branches)
```

To fill an EOP histogram for all tracks passing the TightIso selection, define a histogram filling script like macros/fill_script_test.py .
```
def fill_histograms(hist_filler, outputRootFileName):
    outFile = ROOT.TFile(outputRootFileName, "RECREATE")

    #count the number of tracks in each channel
    histogram_name = "EOP"
    selections = [sel_TightIso]
    trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                         calc_EOP,\
                                                         selections = selections,\
                                                         bins = 100,\
                                                         range_low = -0.5,\
                                                         range_high = +3.5,\
                                                         xlabel ='E/P',\
                                                         ylabel = 'Number of Tracks')
    histograms = hist_filler.DumpHistograms()
    for histogram_name in histograms:
        write_histograms(histograms[histogram_name], outFile)

    outFile.Close()
```

The following lines would prepare a batch job for submission. This job would use 100 condor jobs, and run over a test set of files. The condor jobqueue corresponds to the --queue_flavour flag and can be set to espresso (20mins), longlunch (2hr), workday (8hr), etc.
```
python macros/prepare_submission.py --tree_name LA_EoverP_InDetTrackParticlesSortedLooseIsolatedVertexAssociated_tree --n_jobs 100 --queue_flavour longlunch --file_flavour test --filling_script macros/fill_script.py --job_name test
```
You can change the filling script to any filling script that you have written.

These lines will prepare a batch job for plotting from the identified hadron trees. 
```
python macros/prepare_submission.py --tree_name Tree_Ks --n_jobs 100 --queue_flavour longlunch --file_flavour identified --filling_script macros/fill_script_identified.py --job_name identified_ks
```


## Test one of the jobs locally
```
source condor/test/test_scripts/plot_local.sh 1 condor/test/Submission/test.pickle test
```

## Submit all of the batch jobs
```
condor_submit condor/test/test_scripts/condor_test.sub
```

## Upon job completion
```
hadd Plots_hadded.root Plots*.root
```

## Create TTrees to use for binning. These contain information about the total # of tracks and std deviation.
```
python BinningInformationTree/CreateBinningInformationTree.py -selname MIPSelectionHadFracAbove70 --file PTSpectrumReweighted.root --histogramName EOPDistribution --outputFile MIPSelectonHadFracAbove70Binning
python BinningInformationTree/CreateBinningInformationTree.py -selname 20TRTHitsNonZeroEnergy --file PTSpectrumReweighted.root --histogramName EOPDistribution --outputFile 20TRTHitsNonZeroEnergyBinning
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
