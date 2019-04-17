#!/bin/bash

# Test on a Pythia JZ0 file
xAH_run.py \
    --config src/EoverPAnalysis/scripts/xAH_EoverP.py \
    --files src/EoverPAnalysis/data/input_test.txt \
    --inputList \
    --submitDir localtest/ \
    -f \
    --nevents 1500 \
    direct
    
# Test on a data file
xAH_run.py \
    --config src/EoverPAnalysis/scripts/xAH_EoverP.py \
    --files src/EoverPAnalysis/data/input_test_data.txt \
    --inputList \
    --submitDir localtest_data/ \
    -f \
    --nevents 1500 \
    --extraOptions="--isData" \
    direct
