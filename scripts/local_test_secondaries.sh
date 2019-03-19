#!/bin/bash

xAH_run.py \
    --config src/EoverPAnalysis/scripts/xAH_EoverP_Secondaries.py \
    --files src/EoverPAnalysis/data/input_test.txt \
    --inputList \
    --submitDir localtest/ \
    -f \
    --nevents 50000 \
    direct
