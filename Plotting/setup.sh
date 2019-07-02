
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_94python3 x86_64-centos7-gcc62-opt
source bin/activate
export EOPPlottingDir=$PWD
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_94python3 x86_64-centos7-gcc62-opt
export PYTHONPATH="${EOPPlottingDir}/eop_plotting:$PYTHONPATH"
export PYTHONPATH="${EOPPlottingDir}/utils:$PYTHONPATH"

echo "READY TO GO!"
