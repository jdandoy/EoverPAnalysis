#!/usr/bin/env python 
import os

import subprocess
try:
  __version__ = subprocess.check_output(['git','log', '--pretty=format:%h', '-n 1'], cwd=os.path.dirname(os.path.realpath(__file__))).strip()
except:
  print('git not available to extract current tag')
  __version__ = 'private'

import argparse
parser = argparse.ArgumentParser(description='Submit gridjobs for EoverPAnalysis')

parser.add_argument('--user', '-u', type=str, required=True, help='Your (CERN id) grid username')
parser.add_argument('--tag', dest='analysisTag', type=str, default=__version__, help='Tag for the analysis')
parser.add_argument('--submitDir', type=str, default='submitDir', help='dir to store the output')
parser.add_argument('--overwrite', '-w', action='store_true', default=True, help='overwrite previous submitDir')

args = parser.parse_args()

dataFileList='$ROOTCOREBIN/../EoverPAnalysis/filelists/data15_13TeV_lowmu_all_rucio.txt'
mcFileList='$ROOTCOREBIN/../EoverPAnalysis/filelists/mc15_13TeV_lowmu_all_rucio.txt'

samples = []
try:
  samples.extend([sample.rstrip('\n') for sample in open(os.path.expandvars(dataFileList))])
  samples.extend([sample.rstrip('\n') for sample in open(os.path.expandvars(mcFileList))])
except IOError, e:
  print e
  exit()

for sample in samples:

  sampleTag = '.'.join(sample.split('.')[2:4:])
  optGridOutputSampleName = 'user.{0:s}.{1:s}.{2:s}'.format(args.user, sampleTag, args.analysisTag)

  if 'data15_13TeV' in sample:
    config = os.path.expandvars('$ROOTCOREBIN/../EoverPAnalysis/scripts/config_eop_data_lowmu.py')
  elif 'mc15_13TeV' in sample:
    config = os.path.expandvars('$ROOTCOREBIN/../EoverPAnalysis/scripts/config_eop_mc_lowmu.py')
  else:
    'Invalid sample: {0:s}'.format(samples)
    continue

  command = 'xAH_run.py --files {0:s} --inputRucio --config {1:s} --submitDir {2:s}'.format(sample, config, args.submitDir)
  if args.overwrite:
    command += ' --force'
  command += ' prun --optGridOutputSampleName {}'.format(optGridOutputSampleName)

  print('submitting {0:s}'.format(sample))
  print('running: {0:s}'.format(command))
  subprocess.call(command, shell=True)
