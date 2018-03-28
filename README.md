# EoverPAnalysis
[The complete documentation of this package is hosted on ReadTheDocs.org](http://eoverp.readthedocs.io/en/latest/).

This package is based on the [xAODAnaHelpers (xAH)](https://github.com/UCATLAS/xAODAnaHelpers) RootCore package, thus I strongly recommend you check out the [xAH documentation first](https://xaodanahelpers.readthedocs.io/en/latest/).

For questions please contact: joakim.olsson[at]cern.ch

## Setup

```
mkdir myAnalysis; cd myAnalysis
git clone -b RootCore http://github.com/UCATLAS/xAODAnaHelpers xAODAnaHelpers # checkout R20.7 branch
git clone http://github.com/UCATLAS/EoverPAnalysis EoverPAnalysis
lsetup 'rcsetup Base,2.4.37' # or later version of (Ath)AnalysisBase
rc clean && rc find_packages && rc compile && rc make_par
```

## Running

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

### Local test run

```
mkdir results
xAH_run.py --files $ROOTCOREBIN/../EoverPAnalysis/filelists/data15_13TeV_lowmu_test1.txt --inputList --config $ROOTCOREBIN/../EoverPAnalysis/scripts/config_eop_data_lowmu.py --submitDir $ROOTCOREBIN/../results/eop_data_test_0 --verbose --force direct
```

### First condor test run

```
source $ROOTCOREBIN/../EoverPAnalysis/scripts/run_condor_test_eop_lowmu.sh 0 # where '0' is a tag for the run
```

The output will then be located in 'results', e.g. $ROOTCOREBIN/../results/condor_test_eop_lowmu_{mc,data}_YYYYMMDD_0/

The condor output histograms and cutflows can easily be merged, just run the script below after your condor jobs have finished

```
source $ROOTCOREBIN/../EoverPAnalysis/scripts/merge_condor_eop.py $ROOTCOREBIN/../results/run_condor_eop_lowmu_latest.log
```

## Configuration

### config_* scripts

In 'scripts' you'll find files with names like 'config_*' (ex. 'config_data.py'). These files set the run options, i.e. what event and track selection to apply, what histograms to make, etc. Create your own as needed! 

### run_condor_* scripts

In 'scripts' you'll also find files with names such as 'run_condor_*' (ex. 'run_condor_test_eop_lowmu.sh'). These let you automate submission to condor (see 'First condor test run' above for instructions).
