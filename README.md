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
git clone https://github.com/mattleblanc/IDTrackSel.git
asetup AnalysisBase,21.2.44,here
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
xAH_run.py --files user.luadamek/user.luadamek.14704913.EXT1._000001.pool.root --config ../source/EoverPAnalysis/scripts/config_eop_tree_dump.py --submitDir testing --force direct
```

## Submitting Grid Jobs in Release 21
Grid jobs are handled by a submission script located in $TestArea/EoverPAnalysis/scripts/. The grid job script takes four arguments as input: the submission directory, a txt file with all samples listed, a descriptor to label the output, and a configuration file. As an example, the following command will submit grid jobs to run over the 361022 jet jet MC sample.

```
cd ../run
python $TestArea/EoverPAnalysis/scripts/submit_grid.py --user luadamek --tag 21.2.44 --submitDir results --FileList $TestArea/EoverPAnalysis/filelists/test_list.txt --config $TestArea/EoverPAnalysis/scripts/config_eop_tree_dump.py --descriptor test
```


## Plotting
Plotting macros can be found the root_numpy_plotting folder. See the ipynb files for examples of how to create plots.
