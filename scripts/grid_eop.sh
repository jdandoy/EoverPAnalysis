#!/bin/bash

# Pythia JZ0 - JZ2
xAH_run.py \
    --config src/EoverPAnalysis/scripts/xAH_EoverP.py \
    --files src/EoverPAnalysis/filelists/mc16/mc16_JZ012_calibhit.txt \
    --inputList \
    --inputRucio \
    -f \
    prun \
    --optGridOutputSampleName="user.mleblanc.%in:name[2]%.%in:name[3]%.EOPv01-mc"

# 2017 low-mu data
xAH_run.py \
    --config src/EoverPAnalysis/scripts/xAH_EoverP.py \
    --files src/EoverPAnalysis/filelists/data17/data_lowmu.txt \
    --inputList \
    --inputRucio \
    --extraOptions="--isData" \
    -f \
    prun \
    --optGridOutputSampleName="user.mleblanc.%in:name[2]%.%in:name[3]%.EOPv01-data"
