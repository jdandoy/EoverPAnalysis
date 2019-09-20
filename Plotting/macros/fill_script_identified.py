#!/usr/bin/env python
# coding: utf-8

from eop_plotting.histogram_filling import HistogramFiller, get_p, get_bins, get_log_bins, create_selection_function, create_inverse_selection_function
from eop_plotting.eop_histograms import put_binning_vectors_in_file, create_eop_histograms
from eop_plotting.track_spectrum_plots import create_spectrum_plots

import ROOT
import os
from math import pi
import pickle
import numpy as np

def write_histograms(histogram_dictionary, outFile):
    outFile.cd()
    for key in histogram_dictionary:
        if not outFile.cd(key):
            outFile.mkdir(key)
        outFile.cd(key)
        histogram_dictionary[key].Write()

#This is a script that fills the histograms for
def fill_histograms(hist_filler, outputRootFileName):

    from variables_identified import calc_vertex_mass, calc_vertex_Rxy, calc_weight, calc_vertex_count
    hist_filler.weight_calculator = calc_weight

    histogram_name = "VertexCount"
    selections = []
    trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                         calc_vertex_count,\
                                                         selections = selections,\
                                                         bins = 1,\
                                                         range_low = -0.5,\
                                                         range_high = +0.5,\
                                                         xlabel ='Always 0',\
                                                         ylabel = 'Number of vertices')

    histogram_name = "KsVertexMass"
    selections = []
    trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                         calc_vertex_mass,\
                                                         selections = selections,\
                                                         bins = 100,\
                                                         range_low = 480.0,\
                                                         range_high = 520.0,\
                                                         xlabel ='M [GeV]',\
                                                         ylabel = 'Number of vertices')

    histogram_name = "PhiVertexMass"
    selections = []
    trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                         calc_vertex_mass,\
                                                         selections = selections,\
                                                         bins = 100,\
                                                         range_low = 987.354,\
                                                         range_high = 1200.0,\
                                                         xlabel ='M [GeV]',\
                                                         ylabel = 'Number of vertices')

    histogram_name = "LambdaVertexMass"
    selections = []
    trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                         calc_vertex_mass,\
                                                         selections = selections,\
                                                         bins = 100,\
                                                         range_low = 1105.0,\
                                                         range_high = 1125.0,\
                                                         xlabel ='M [GeV]',\
                                                         ylabel = 'Number of vertices')

    histogram_name = "vertex_Rxy"
    selections = []
    trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                         calc_vertex_Rxy,\
                                                         selections = selections,\
                                                         bins = 110,\
                                                         range_low = 0.0,\
                                                         range_high = 40.0,\
                                                         xlabel ='rxy [mm]',\
                                                         ylabel = 'Number of vertices')

    histograms = hist_filler.DumpHistograms()
    outFile = ROOT.TFile(outputRootFileName, "RECREATE")
    for histogram_name in histograms:
        write_histograms(histograms[histogram_name], outFile)


    print("THEJOBFINISHED!")
