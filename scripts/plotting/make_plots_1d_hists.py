#! /usr/bin/env python

# @file:    make_eop_plots.py
# @purpose: make plots for ATLAS E/p analysis
# @author:  Joakim Olsson <joakim.olsson@cern.ch>
# @date:    June 2015

import os, sys
import re

import argparse

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='E/p plotting')

  parser.add_argument('--data', required=False, type=str, dest='data', metavar='<datafile.root>', default='$TestArea/../run/results/hist-user.luadamek.root', help='Input data file')
  parser.add_argument('--mc', required=False, type=str, dest='mc', metavar='<mcfile.root>', default='$TestArea/../run/results/hist-user.luadamek.root', help='Input MC file')
  parser.add_argument('--output', '-o', required=False, type=str, dest='output', default='./plots', help='Directory for output files')
  parser.add_argument('--selections', '-f', type=str, nargs='+', dest='selections', default=['EoverP_LoosePrimaryTrks_ClusterEnergy_Tile_defaultCuts'], help='Folders where to look for plots in the input files')

  args = parser.parse_args()

  data = args.data
  mc = args.mc

  output_dir = args.output
  ensure_dir(output_dir)

  selections = args.selections

  # Default setup (low-mu samples) TODO: Make this easier to customize
  if re.search("(?<=_)\d{1,8}_.+(?=\.)", data):
    tag = re.search("(?<=_)\d{1,8}_.+(?=\.)", data).group()
  else:
    tag = "test"

  testArea = os.environ['TestArea']
  plot_input_list = testArea +"/EoverPAnalysis/scripts/plotting/plotlist_eop_1d_hists.txt"
  file_tag = "lowmu_mc"
  # sample_tag = "13 TeV 2015 data from period B1 low-#mu run"
  sample_tag = ""
  leg_labels = ["Data 2017, Low-#LT#font[50]{#mu}#GT", "Pythia MinBias"]
  lumi_label = "#sqrt{s} = 13 TeV, ? b^{-1}"
  options = ""
  #options = "public"


  from plot_maker import plot_1d_hists
  for selection in selections:

    selection_tag = "tilecuts"

    print "----> plot config: "
    print "selection:", selection
    print "data:", data
    print "mc:", mc
    print "tag:", tag
    print "plot_input_list:",plot_input_list
    print "output_dir:",output_dir
    print "file_tag:",file_tag
    print "sample_tag:",sample_tag

    print "\n----> running script to save plots: "
    print "selection: "+selection
    file_tag_new = file_tag+"_"+selection_tag

    plot_1d_hists(data, mc, selection, plot_input_list,
                  output_dir, file_tag_new, sample_tag, options, leg_labels, lumi_label)
