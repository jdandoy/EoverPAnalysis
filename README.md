# EoverPAnalysis
[The complete documentation of this package is hosted on ReadTheDocs.org](http://eoverp.readthedocs.io/en/latest/).

This package is based on the [xAODAnaHelpers (xAH)](https://github.com/UCATLAS/xAODAnaHelpers) RootCore package, thus I strongly recommend you check out the [xAH documentation first](https://xaodanahelpers.readthedocs.io/en/latest/).

For questions please contact: joakim.olsson[at]cern.ch

Ported to release 21 by: Lukas Adamek (lukas.adamek[at]cern.ch)

## Setup in Release 21
```
mkdir myAnalysis; cd myAnalysis
mkdir source && mkdir run && mkdir build && mkdir run/results
cd source
git clone http://github.com/UCATLAS/xAODAnaHelpers xAODAnaHelpers
git clone http://github.com/luadamek/EoverPAnalysis
asetup AnalysisBase,21.2.22,here
cd ../build
cmake ../source && make
```

## Running locally in Release 21
```
cd ../build
source x86_64-slc6-gcc62-opt/setup.sh
cd ../run
mkdir results
rucio download user.luadamek.14704913.EXT1._000001.pool.root
xAH_run.py --files $TestArea/../run/user.luadamek/user.luadamek.14704913.EXT1._000001.pool.root --config $TestArea/EoverPAnalysis/scripts/config_eop_mc_lowmu_runII_general.py --submitDir $TestArea/../run/results/eop_mc_test_0 --force direct
cp results/eop_mc_test_0/hist-user.luadamek.root results/
python $TestArea/EoverPAnalysis/scripts/plotting/make_plots_1d_hists.py --selections EoverP_LoosePrimaryTrks_ClusterEnergy_Run1paper
```

## Submitting jobs to CONDOR in Release 21
First you need to make sure that you have a list of files that you want to run over in a txt file in the filelists folder. To generate the .txt file, use GenerateListOfFiles.py script.
python GenerateListOfFiles --grid_site GRID_SITE --dataset_name DATASET_NAME --output_file OUTPUTFILENAME

As an example, we can call this command with the following dataset:
```
cd /source/EoverPAnalysis/filelists/
python GenerateListOfFiles.py --grid_site CA-VICTORIA-WESTGRID-T2_LOCALGROUPDISK --dataset_name user.luadamek.mc16_13TeV.361022.jetjet.DAOD_EOP.e3668_s3170_r10572.20180724_alpha1_EXT1 --output_filename 361022_filelist.txt
```

## Setup in Release 20.7

```
mkdir myAnalysis; cd myAnalysis
git clone -b RootCore http://github.com/UCATLAS/xAODAnaHelpers xAODAnaHelpers # checkout R20.7 branch
git clone http://github.com/UCATLAS/EoverPAnalysis EoverPAnalysis
lsetup 'rcsetup Base,2.4.37' # or later version of (Ath)AnalysisBase
rc clean && rc find_packages && rc compile && rc make_par
```

## Running in Release 20.7

### Grid proxy

If your datasets are located on the grid (the ones in the file lists that comes with this package are), you need to have a valid grid proxy in order to access them.

```
voms-proxy-init -voms atlas
``` 

If you haven't done so already, you might want to add the following lines to your ~/.bash_profile:

```
alias grid="voms-proxy-init -voms atlas -out $HOME/.globus/gridproxy.cert -valid 1000:00"
export X509_USER_PROXY=$HOME/.globus/gridproxy.cert
```

The datasets that come with the default package are located on MWT2_UC_LOCALGROUPDISK, and are accessed via FAX. To set up fax, do:

```
lsetup fax; fax-get-best-redirector
```

NOTE: If you are submitting jobs to condor from lxplus, you'll need to put your gridproxy.cert in a location accessible by the CERN condor nodes:

```
cp /afs/cern.ch/user/j/jolsson/.globus/gridproxy.cert /eos/user/j/jolsson/
export X509_USER_PROXY=/eos/user/j/jolsson/gridproxy.cert
```


### First condor test run

```
source $ROOTCOREBIN/../EoverPAnalysis/scripts/run_condor_test_eop_lowmu.sh 0 # where '0' is a tag for the run
```

The output will then be located in 'results', e.g. $ROOTCOREBIN/../results/condor_test_eop_lowmu_{mc,data}_YYYYMMDD_0/

The condor output histograms and cutflows can easily be merged, jusrun the script below after your condor jobs have finished

```
source $ROOTCOREBIN/../EoverPAnalysis/scripts/merge_condor_eop.py $ROOTCOREBIN/../results/run_condor_eop_lowmu_latest.log
```

## Configuration

### config_* scripts

In 'scripts' you'll find files with names like 'config_*' (ex. 'config_data.py'). These files set the run options, i.e. what event and track selection to apply, what histograms to make, etc. Create your own as needed! 

### run_condor_* scripts

In 'scripts' you'll also find files with names such as 'run_condor_*' (ex. 'run_condor_test_eop_lowmu.sh'). These let you automate submission to condor (see 'First condor test run' above for instructions).
