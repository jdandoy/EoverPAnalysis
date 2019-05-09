
export EOPPlottingDir=$PWD
mkdir -p $PWD/CondorPythonLocal
export PYTHONPATH=$PWD/CondorPythonLocal/lib64/python2.7/site-packages:$PYTHONPATH
export BLAS=None LAPACK=None ATLAS=None

echo "READY TO GO!"
