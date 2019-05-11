
export EOPPlottingDir=$PWD
mkdir -p $PWD/CondorPythonLocal
export PYTHONPATH=$PWD/CondorPythonLocal/lib64/python2.7/site-packages:$PYTHONPATH
export BLAS=None LAPACK=None ATLAS=None

pip install --install-option="--prefix=$PWD/CondorPythonLocal" --ignore-installed cython
pip install --install-option="--prefix=$PWD/CondorPythonLocal" --ignore-installed root_numpy
pip install --install-option="--prefix=$PWD/CondorPythonLocal" --ignore-installed psutil 

echo "READY TO GO!"
