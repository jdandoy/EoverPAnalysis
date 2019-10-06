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

    from variables_identified import calc_vertex_mass, calc_vertex_Rxy, calc_weight, calc_vertex_count, calc_cos_theta, calc_vertex_pt
    hist_filler.weight_calculator = calc_weight

    from selections_identified import sel_tight_cos_theta, sel_chi_square_fifteen

    for selections in [[], [sel_tight_cos_theta, sel_chi_square_fifteen]]:
       descriptor = "_".join([sel.name for sel in selections])
       histogram_name = "VertexCount" + descriptor
       trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                            calc_vertex_count,\
                                                            selections = selections,\
                                                            bins = 1,\
                                                            range_low = -0.5,\
                                                            range_high = +0.5,\
                                                            xlabel ='Always 0',\
                                                            ylabel = 'Number of Vertices')

       histogram_name = "KsVertexMass" + descriptor
       trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                            calc_vertex_mass,\
                                                            selections = selections,\
                                                            bins = 100,\
                                                            range_low = 480.0,\
                                                            range_high = 520.0,\
                                                            xlabel ='M [MeV]',\
                                                            ylabel = 'Number of Vertices')

       histogram_name = "PhiVertexMass" + descriptor
       trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                            calc_vertex_mass,\
                                                            selections = selections,\
                                                            bins = 100,\
                                                            range_low = 987.354,\
                                                            range_high = 1200.0,\
                                                            xlabel ='M [MeV]',\
                                                            ylabel = 'Number of Vertices')

       histogram_name = "LambdaVertexMass" + descriptor
       trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                            calc_vertex_mass,\
                                                            selections = selections,\
                                                            bins = 100,\
                                                            range_low = 1105.0,\
                                                            range_high = 1125.0,\
                                                            xlabel ='M [MeV]',\
                                                            ylabel = 'Number of Vertices')

       histogram_name = "vertex_Rxy" + descriptor
       trkCountHist = hist_filler.book_histogram_fill(histogram_name,\
                                                            calc_vertex_Rxy,\
                                                            selections = selections,\
                                                            bins = 110,\
                                                            range_low = 0.0,\
                                                            range_high = 40.0,\
                                                            xlabel ='rxy [mm]',\
                                                            ylabel = 'Number of Vertices')

       histogram_name = "vetex_pt" + descriptor
       hist_filler.book_histogram_fill(histogram_name,\
                                      calc_vertex_pt,\
                                      selections = selections,\
                                      bins = 110,\
                                      range_low = 0.0,\
                                      range_high = 50.0,\
                                      xlabel = "Vertex P_{T} [GeV]",\
                                      ylabel = "Number of Vertices")

       histogram_name = "cos_theta" + descriptor
       hist_filler.book_histogram_fill(histogram_name,\
                                       calc_cos_theta,\
                                       selections = selections,\
                                       bins = 100,\
                                       range_low = -1.0,\
                                       range_high = 1.0,\
                                       xlabel = 'cos(#theta)',\
                                       ylabel = 'Number of Vertices')

       histograms = hist_filler.DumpHistograms()
       outFile = ROOT.TFile(outputRootFileName, "RECREATE")
       for histogram_name in histograms:
           write_histograms(histograms[histogram_name], outFile)

    print("THEJOBFINISHED!")
