
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
export EOPPlottingDir=$(pwd)
asetup AnalysisBase,21.2.23

pip uninstall cython
pip uninstall root_numpy
pip unistall psutil

pip install --install-option="--prefix=~/CondorPythonLocal" --ignore-installed cython
pip install --install-option="--prefix=~/CondorPythonLocal" --ignore-installed root_numpy
pip install --install-option="--prefix=~/CondorPythonLocal" --ignore-installed psutil

echo "READY TO GO!"
