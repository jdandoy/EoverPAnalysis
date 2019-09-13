
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_94python3 x86_64-centos7-gcc62-opt
source bin/activate
export EOPPlottingDir=$PWD
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_94python3 x86_64-centos7-gcc62-opt
export PYTHONPATH="${EOPPlottingDir}/eop_plotting:$PYTHONPATH"
export PYTHONPATH="${EOPPlottingDir}/utils:$PYTHONPATH"
export PYTHONPATH="${EOPPlottingDir}/macros:$PYTHONPATH"
export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2016/bin/x86_64-linux:$PATH

echo "READY TO GO!"
