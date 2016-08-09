#!/bin/bash
if [ $# -eq 0 ]

  then
    echo "Usage: source run_condor_all_eop_pileup.tag"

  else

    cd $ROOTCOREBIN/..

    tag=$1
    today=$(date +"%Y%m%d")
    files_data=EoverPAnalysis/filelists/data15_13TeV_pileup.txt
    files_JZ0W=EoverPAnalysis/filelists/mc15_13TeV_pileup_JZ0W_test200k.txt
    files_JZ1W=EoverPAnalysis/filelists/mc15_13TeV_pileup_JZ1W_all.txt

    mkdir -p results

    echo "---> Running JZxW pileup MC samples:"
    echo xAH_run.py --files ${files_data} --inputList --config EoverPAnalysis/scripts/config_eop_data.py --submitDir results/condor_all_eop_pileup_data_${today}_${tag} --verbose --force condor --optFilesPerWorker 10
    xAH_run.py --files ${files_data} --inputList --config EoverPAnalysis/scripts/config_eop_data.py --submitDir results/condor_all_eop_pileup_data_${today}_${tag} --verbose --force condor --optFilesPerWorker 10
    echo xAH_run.py --files ${files_JZ0W} --inputList --config EoverPAnalysis/scripts/config_eop_mc.py --submitDir results/condor_all_eop_pileup_JZ0W_${today}_${tag} --verbose --force condor --optFilesPerWorker 10
    xAH_run.py --files ${files_JZ0W} --inputList --config EoverPAnalysis/scripts/config_eop_mc.py --submitDir results/condor_all_eop_pileup_JZ0W_${today}_${tag} --verbose --force condor --optFilesPerWorker 10
    echo xAH_run.py --files ${files_JZ1W} --inputList --config EoverPAnalysis/scripts/config_eop_mc.py --submitDir results/condor_all_eop_pileup_JZ1W_${today}_${tag} --verbose --force condor --optFilesPerWorker 10
    xAH_run.py --files ${files_JZ1W} --inputList --config EoverPAnalysis/scripts/config_eop_mc.py --submitDir results/condor_all_eop_pileup_JZ1W_${today}_${tag} --verbose --force condor --optFilesPerWorker 10

    echo "---> Write to logfile:"
    echo "# ---> "$(date +"%Y-%m-%d:%H:%M:%S") >> results/run_condor_eop_pileup.log
    echo ${files_data} >> results/run_condor_eop_pileup.log
    echo results/condor_all_eop_pileup_data_${today}_${tag} >> results/run_condor_eop_pileup.log
    echo ${files_JZ0W} >> results/run_condor_eop_pileup.log
    echo results/condor_all_eop_pileup_JZ0W_${today}_${tag} >> results/run_condor_eop_pileup.log
    echo ${files_JZ1W} >> results/run_condor_eop_pileup.log
    echo results/condor_all_eop_pileup_JZ1W_${today}_${tag} >> results/run_condor_eop_pileup.log
    echo ${files_data} > results/run_condor_eop_pileup_latest.log
    echo results/condor_all_eop_pileup_data_${today}_${tag} >> results/run_condor_eop_pileup_latest.log
    echo ${files_JZ0W} >> results/run_condor_eop_pileup_latest.log
    echo results/condor_all_eop_pileup_JZ0W_${today}_${tag} >> results/run_condor_eop_pileup_latest.log
    echo ${files_JZ1W} >> results/run_condor_eop_pileup_latest.log
    echo results/condor_all_eop_pileup_JZ1W_${today}_${tag} >> results/run_condor_eop_pileup_latest.log

    echo "--> Jobs submitted!"
    echo "source $ROOTCOREBIN/../EoverPAnalysis/scripts/merge_condor_eop.py $ROOTCOREBIN/../results/run_condor_eop_pileup_latest.log # when condor jobs are finished to merge output files"

fi
