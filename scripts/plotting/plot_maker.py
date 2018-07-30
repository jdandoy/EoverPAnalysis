#! /usr/bin/env python
import os
import sys
import re
import subprocess
import numpy as np
from ROOT import *
from eophelper.plotting import *
from eophelper.util import *

# Joakim Olsson (joakim.olsson@cern.ch)
# Created 2016-02-22

def plot_cutflows(datafile, mcfile, plot_input_list,
                  output_dir, sample_tag):

    data = TFile.Open(datafile, "READ")
    mc = TFile.Open(mcfile, "READ")

    # create output dir(s) if necessary
    ensure_dir(output_dir)

    # list of plots to make
    plot_names = plot_input_list

    # set basic plotting options
    options = "" # pass string options to PlotOptions
    do_auto_scale = True # scale y-axis appropriately
    do_data = False # plot data as solid black markers
    do_data_mc = True # plot data_mc_ratio
    invert_data_mc_ratio = True # invert data mc ratio order
    do_ratio = False # add a ratio plot under main plot
    do_norm = False # normalize all plots to unity
    do_norm_mc_to_data = False #
    do_auto_scale = False # scale y-axis appropriately
    draw_style = "hist][" # pass args to h.Draw("args"), NB: not applicable for do_data
    marker_type = "data_mc" #sets the markers and colors
    leg_labels = [""]
    atlas_label = "Internal"
    lumi_label = "#sqrt{s} = 13 TeV, 1.6 nb^{-1}"

    for plot_name in plot_names:
        hdata = data.Get(plot_name)
        hmc = mc.Get(plot_name)

        if sample_tag != "":
            plot_tags = [sample_tag]
        else:
            plot_tags = []


        if "TH1" in str(hdata):
            plot_opts = PlotOptions(plot_name,
                                    options,
                                    do_data,
                                    do_data_mc,
                                    invert_data_mc_ratio,
                                    do_ratio,
                                    do_norm,
                                    do_norm_mc_to_data,
                                    do_auto_scale,
                                    draw_style,
                                    marker_type,
                                    leg_labels,
                                    plot_tags,
                                    atlas_label,
                                    lumi_label)

            plot_output_path = output_dir+"/"+plot_name
            save_hist(hdata, plot_opts, plot_output_path)

def plot_1d_hists(datafile, mcfile, plot_input_dir, plot_input_list,
                  output_dir, file_tag, sample_tag, options="", leg_labels = "", lumi_label = ""):

    data = TFile.Open(datafile, "READ")
    mc = TFile.Open(mcfile, "READ")

    # create output dir(s) if necessary
    ensure_dir(output_dir+"/"+plot_input_dir+"/")

    # list of plots to make
    plot_names = [line.rstrip('\n') for line in open(plot_input_list) if "#" not in line]

    # set basic plotting options
    options = options # pass string options to PlotOptions
    do_data = True # plot data as solid black markers
    do_data_mc = True # plot data_mc_ratio
    invert_data_mc_ratio = True # invert data mc ratio order
    do_ratio = True # add a ratio plot under main plot
    do_norm = True # normalize all plots to unity
    do_norm_mc_to_data = False #
    do_auto_scale = True # scale y-axis appropriately
    draw_style = "hist][" # pass args to h.Draw("args"), NB: not applicable for do_data
    # marker_type = "data_mc" #sets the markers and colors
    marker_type = "data_mc" #sets the markers and colors
    if leg_labels == "":
        leg_labels = ["MC 2015", "MC15c"] #, "MC (Pythia JZ1W)"]
    #atlas_label = "Internal"
    atlas_label = "Internal"
    if lumi_label == "":
        lumi_label = "#sqrt{s} = 13 TeV, ? nb^{-1}"

    for plot_name in plot_names:
        hdata = data.Get(plot_input_dir+"/"+plot_name)
        hmc = mc.Get(plot_input_dir+"/"+plot_name)
        if plot_name == "mu_avg_0":
            hdata = data.Get(plot_input_dir+"/mu_avg")
            hmc = mc.Get(plot_input_dir+"/"+plot_name)


        if sample_tag != "":
            plot_tags = [sample_tag]
        else:
            plot_tags = []

        draw_avg = False
        if re.search("^eop", plot_name) or re.search("^p.*eop", plot_name):
            draw_avg = True

        if "TH1" in str(hdata):
            plot_opts = PlotOptions(plot_input_dir+"_"+plot_name,
                                    options,
                                    do_data,
                                    do_data_mc,
                                    invert_data_mc_ratio,
                                    do_ratio,
                                    do_norm,
                                    do_norm_mc_to_data,
                                    do_auto_scale,
                                    draw_style,
                                    marker_type,
                                    leg_labels,
                                    plot_tags,
                                    atlas_label,
                                    lumi_label,
                                    draw_avg)

            plot_output_path = output_dir+"/"+plot_input_dir+"/"+file_tag+"_"+plot_name

            save_hists_overlay([hdata, hmc], plot_opts, plot_output_path)

def plot_1d_hists_cf_mc(mcfiles, plot_input_dir, plot_input_list,
                  output_dir, file_tag, sample_tag, options = "", leg_labels = "", avg_ratio_labels = ""):


    mcfiles_open = []
    for mcfile in mcfiles:
        mcfiles_open.append(TFile.Open(mcfile, "READ"))

    # create output dir(s) if necessary
    ensure_dir(output_dir+"/"+plot_input_dir+"/")

    # list of plots to make
    plot_names = [line.rstrip('\n') for line in open(plot_input_list) if "#" not in line]

    # set basic plotting options
    options = options # pass string options to PlotOptions
    do_data = False # plot data as solid black markers
    do_data_mc = True # plot data_mc_ratio
    invert_data_mc_ratio = True # invert data mc ratio order
    do_ratio = True # add a ratio plot under main plot
    do_norm = True # normalize all plots to unity
    do_norm_mc_to_data = False #
    do_auto_scale = True # scale y-axis appropriately
    draw_style = "hist][" # pass args to h.Draw("args"), NB: not applicable for do_data
    # marker_type = "data_mc" #sets the markers and colors
    marker_type = "mc" #sets the markers and colors
    atlas_label = "Internal"
    lumi_label = "" #"#sqrt{s} = 13 TeV, 1.6 nb^{-1}"

    for plot_name in plot_names:
        hists = []
        for mc in mcfiles_open:
            hists.append(mc.Get(plot_input_dir+"/"+plot_name))

        if sample_tag != "":
            plot_tags = [sample_tag]
        else:
            plot_tags = []

        draw_avg = False
        # if (not "noavgline" in plot_name):
        #     if re.search("^eop", plot_name) or re.search("^p.*eop", plot_name):
        #         draw_avg = True
        draw_avg_ratio = False
        if re.search("^eop", plot_name) or re.search("^p.*eop", plot_name):
            draw_avg_ratio = True

        if "TH1" in str(hists[0]):
            plot_opts = PlotOptions(plot_input_dir+"_"+plot_name,
                                    options,
                                    do_data,
                                    do_data_mc,
                                    invert_data_mc_ratio,
                                    do_ratio,
                                    do_norm,
                                    do_norm_mc_to_data,
                                    do_auto_scale,
                                    draw_style,
                                    marker_type,
                                    leg_labels,
                                    plot_tags,
                                    atlas_label,
                                    lumi_label,
                                    draw_avg,
                                    draw_avg_ratio,
                                    avg_ratio_labels)

            plot_output_path = output_dir+"/"+plot_input_dir+"/"+file_tag+"_"+plot_name

            save_hists_overlay(hists, plot_opts, plot_output_path)

def plot_1d_hists_cf_data(datafiles, plot_input_dir, plot_input_list,
                          output_dir, file_tag, sample_tag, options = "", leg_labels = "", avg_ratio_labels = ""):


    datafiles_open = []
    for datafile in datafiles:
        datafiles_open.append(TFile.Open(datafile, "READ"))

    # create output dir(s) if necessary
    ensure_dir(output_dir+"/"+plot_input_dir+"/")

    # list of plots to make
    plot_names = [line.rstrip('\n') for line in open(plot_input_list) if "#" not in line]

    # set basic plotting options
    options = options # pass string options to PlotOptions
    do_data = True # plot data as solid black markers
    do_data_mc = True # plot data_mc_ratio
    invert_data_mc_ratio = True # invert data mc ratio order
    do_data_data = True # plot data vs. data comparisons
    do_ratio = True # add a ratio plot under main plot
    do_norm = True # normalize all plots to unity
    do_norm_mc_to_data = False #
    do_auto_scale = True # scale y-axis appropriately
    draw_style = "hist][" # pass args to h.Draw("args"), NB: not applicable for do_data
    marker_type = "data_data" #sets the markers and colors
    atlas_label = "Internal"
    lumi_label = "" #"#sqrt{s} = 13 TeV, 1.6 nb^{-1}"

    for plot_name in plot_names:
        hists = []
        for h in datafiles_open:
            hists.append(h.Get(plot_input_dir+"/"+plot_name))

        if sample_tag != "":
            plot_tags = [sample_tag]
        else:
            plot_tags = []

        draw_avg = False
        # if (not "noavgline" in plot_name):
        #     if re.search("^eop", plot_name) or re.search("^p.*eop", plot_name):
        #         draw_avg = True
        draw_avg_ratio = False
        if re.search("^eop", plot_name) or re.search("^p.*eop", plot_name):
            draw_avg_ratio = True

        if "TH1" in str(hists[0]):
            plot_opts = PlotOptions(plot_input_dir+"_"+plot_name,
                                    options,
                                    do_data,
                                    do_data_mc,
                                    invert_data_mc_ratio,
                                    do_ratio,
                                    do_norm,
                                    do_norm_mc_to_data,
                                    do_auto_scale,
                                    draw_style,
                                    marker_type,
                                    leg_labels,
                                    plot_tags,
                                    atlas_label,
                                    lumi_label,
                                    draw_avg,
                                    draw_avg_ratio,
                                    avg_ratio_labels,
                                    do_data_data)

            plot_output_path = output_dir+"/"+plot_input_dir+"/"+file_tag+"_"+plot_name

            save_hists_overlay(hists, plot_opts, plot_output_path)

def plot_2d_hists(inputfile, plot_input_dir, plot_input_list,
                  output_dir, file_tag, sample_tag):

    tfile = TFile.Open(inputfile, "READ")

    # create output dir(s) if necessary
    ensure_dir(output_dir+"/"+plot_input_dir+"/")

    # list of plots to make
    plot_names = [line.rstrip('\n') for line in open(plot_input_list) if "#" not in line]

    # set basic plotting options
    options = "2D" # pass string options to PlotOptions
    do_data = False
    do_data_mc = True # plot data_mc_ratio
    invert_data_mc_ratio = True # invert data mc ratio order
    do_ratio = False
    do_norm = False
    do_norm_mc_to_data = False
    do_auto_scale = False
    draw_style = "colz"
    marker_type = "" #sets the markers and colors
    leg_labels = []
    leg_labels = ["MC 2015", "MC15c"] #, "MC (Pythia JZ1W)"]
    atlas_label = "Internal"
    lumi_label = "#sqrt{s} = 13 TeV, 1.6 nb^{-1}"

    do_profile_x = False

    for plot_name in plot_names:
        hist = tfile.Get(plot_input_dir+"/"+plot_name)

        if sample_tag != "":
            plot_tags = [sample_tag]
        else:
            plot_tags = []

        if "TH2" in str(hist):
            plot_opts = PlotOptions(plot_input_dir+"_"+plot_name,
                                    options,
                                    do_data,
                                    do_data_mc,
                                    invert_data_mc_ratio,
                                    do_ratio,
                                    do_norm,
                                    do_norm_mc_to_data,
                                    do_auto_scale,
                                    draw_style,
                                    marker_type,
                                    leg_labels,
                                    plot_tags,
                                    atlas_label,
                                    lumi_label)

            plot_output_path = output_dir+"/"+plot_input_dir+"/"+file_tag+"_"+plot_name

            save_hist(hist, plot_opts, plot_output_path)

def plot_profilex(datafile, mcfile, plot_input_dir, plot_input_list,
                  output_dir, file_tag, sample_tag, options, leg_labels, lumi_label):

    data = TFile.Open(datafile, "READ")
    mc = TFile.Open(mcfile, "READ")

    # create output dir(s) if necessary
    ensure_dir(output_dir+"/"+plot_input_dir+"/")

    # list of plots to make
    plot_names = [line.rstrip('\n') for line in open(plot_input_list) if "#" not in line]

    # set basic plotting options
    options = options # pass string options to PlotOptions
    do_data = True # plot data as solid black markers
    do_data_mc = True # plot data_mc_ratio
    invert_data_mc_ratio = True # invert data mc ratio order
    do_ratio = True  # add a ratio plot under main plot
    do_norm = False # normalize all plots to unity
    do_norm_mc_to_data = False #
    do_auto_scale = False # scale y-axis appropriately
    draw_style = "hist" # pass args to h.Draw("args"), NB: not applicable for do_data
    marker_type = "data_data" #sets the markers and colors
    if leg_labels == "":
        leg_labels = ["MC 2015", "MC15c"] #, "MC (Pythia JZ1W)"]
    atlas_label = "Internal"
    if lumi_label == "":
        lumi_label = "#sqrt{s} = 13 TeV, 1.6 nb^{-1}"

    for plot_name in plot_names:
        # print "!!!!!!!"
        # print plot_input_dir
        # print plot_name
        hdata = data.Get(plot_input_dir+"/"+plot_name)
        hmc = mc.Get(plot_input_dir+"/"+plot_name)

        if sample_tag != "":
            plot_tags = [sample_tag]
        else:
            plot_tags = []

        draw_avg_ratio = False
        # if re.search("^eop.*_vs_", plot_name) or re.search("^p.*eop.*_vs_", plot_name):
        #     draw_avg_ratio = True

        plot_opts = PlotOptions(plot_input_dir+"_"+plot_name,
                                options,
                                do_data,
                                do_data_mc,
                                invert_data_mc_ratio,
                                do_ratio,
                                do_norm,
                                do_norm_mc_to_data,
                                do_auto_scale,
                                draw_style,
                                marker_type,
                                leg_labels,
                                plot_tags,
                                atlas_label,
                                lumi_label,
                                False,
                                draw_avg_ratio)

        if "TH2" in str(hdata):

            hdata.Rebin2D(plot_opts.rebinx, plot_opts.rebiny)
            hmc.Rebin2D(plot_opts.rebinx, plot_opts.rebiny)
            hdata_p = tprofile2hist(hdata.ProfileX(hdata.GetName()+"_profile_data"))
            hmc_p = tprofile2hist(hmc.ProfileX(hmc.GetName()+"_profile_mc"))
            # print plot_opts.rebinx
            # print plot_opts.rebiny
            # hdata_p = hdata.ProfileX("profile_data")
            # hmc_p = hmc.ProfileX("profile_mc")

            plot_output_path = output_dir+"/"+plot_input_dir+"/"+file_tag+"_"+plot_name
            save_hists_overlay([hdata_p, hmc_p], plot_opts, plot_output_path)

        elif "TProfile" in str(hdata):
            hdata_p = tprofile2hist(hdata)
            hmc_p = tprofile2hist(hmc)
            plot_output_path = output_dir+"/"+plot_input_dir+"/"+file_tag+"_"+plot_name
            save_hists_overlay([hdata_p, hmc_p], plot_opts, plot_output_path)

def plot_profilex_cf_data(datafiles, plot_input_dir, plot_input_list,
                  output_dir, file_tag, sample_tag, options, leg_labels, lumi_label = ""):

    datafiles_open = []
    for datafile in datafiles:
        datafiles_open.append(TFile.Open(datafile, "READ"))

    # create output dir(s) if necessary
    ensure_dir(output_dir+"/"+plot_input_dir+"/")

    # list of plots to make
    plot_names = [line.rstrip('\n') for line in open(plot_input_list) if "#" not in line]

    # set basic plotting options
    options = options # pass string options to PlotOptions
    do_data = True # plot data as solid black markers
    do_data_mc = True # plot data_mc_ratio
    invert_data_mc_ratio = True # invert data mc ratio order
    do_data_data = True # plot data vs. data comparisons
    do_ratio = True  # add a ratio plot under main plot
    do_norm = False # normalize all plots to unity
    do_norm_mc_to_data = False #
    do_auto_scale = False # scale y-axis appropriately
    draw_style = "hist" # pass args to h.Draw("args"), NB: not applicable for do_data
    marker_type = "data_mc" #sets the markers and colors
    atlas_label = "Internal"

    for plot_name in plot_names:

        hists = []
        for i,f in enumerate(datafiles_open):
            h = f.Get(plot_input_dir+"/"+plot_name)
            hists.append(h)

        if sample_tag != "":
            plot_tags = [sample_tag]
        else:
            plot_tags = []

        draw_avg_ratio = False
        if re.search("^eop.*_vs_", plot_name) or re.search("^p.*eop.*_vs_", plot_name):
            draw_avg_ratio = True

        if "TH2" in str(hists[0]):
            plot_opts = PlotOptions(plot_input_dir+"_"+plot_name,
                                    options,
                                    do_data,
                                    do_data_mc,
                                    invert_data_mc_ratio,
                                    do_ratio,
                                    do_norm,
                                    do_norm_mc_to_data,
                                    do_auto_scale,
                                    draw_style,
                                    marker_type,
                                    leg_labels,
                                    plot_tags,
                                    atlas_label,
                                    lumi_label,
                                    False,
                                    draw_avg_ratio,
                                    [],
                                    do_data_data)

            hists_p = []
            for i,h in enumerate(hists):
                h.Rebin2D(plot_opts.rebinx, plot_opts.rebiny)
                h_p = tprofile2hist(h.ProfileX(h.GetName()+"_profile_{}".format(i)))
                hists_p.append(h_p)

            plot_output_path = output_dir+"/"+plot_input_dir+"/"+file_tag+"_"+plot_name

            save_hists_overlay(hists_p, plot_opts, plot_output_path)

def plot_bg_subtr(datafile, mcfile, plot_input_dir, plot_input_list,
                  output_dir, file_tag, sample_tag, options, leg_labels, lumi_label):

    data = TFile.Open(datafile, "READ")
    mc = TFile.Open(mcfile, "READ")

    # create output dir(s) if necessary
    ensure_dir(output_dir+"/"+plot_input_dir+"/")

    # list of plots to make
    plot_names = [line.rstrip('\n') for line in open(plot_input_list) if (("#" not in line) or line == "")]

    # set basic plotting options
    options = options # pass string options to PlotOptions
    do_data = True # plot data as solid black markers
    do_data_mc = True # plot data_mc_ratio
    invert_data_mc_ratio = True # invert data mc ratio order
    do_ratio = True  # add a ratio plot under main plot
    do_norm = False # normalize all plots to unity
    do_norm_mc_to_data = False #
    do_auto_scale = False # scale y-axis appropriately
    draw_style = "hist" # pass args to h.Draw("args"), NB: not applicable for do_data
    marker_type = "data_mc" #sets the markers and colors
    if leg_labels == "":
        leg_labels = ["MC 2015", "MC15c"] #, "MC (Pythia JZ1W)"]
    atlas_label = "Internal"
    if lumi_label == "":
        lumi_label = "#sqrt{s} = 13 TeV, 1.6 nb^{-1}"

    for plot_name in plot_names:
        hdata = data.Get(plot_input_dir+"/"+plot_name)
        hmc = mc.Get(plot_input_dir+"/"+plot_name)

        plot_name_bg = plot_name
        plot_name_bg = re.sub("0_\d+_","",plot_name_bg)
        plot_name_bg = re.sub("(?<=_)[A-z0-9]+(?=_C)", "EM_BG", plot_name_bg)

        hdata_bg = data.Get(plot_input_dir+"/"+plot_name_bg)
        hmc_bg = mc.Get(plot_input_dir+"/"+plot_name_bg)

        if sample_tag != "":
            plot_tags = [sample_tag]
        else:
            plot_tags = []

        draw_avg_ratio = False
        if re.search("^eop.*_vs_", plot_name) or re.search("^p.*eop.*_vs_", plot_name):
            draw_avg_ratio = True

        if "TH2" in str(hdata):
            plot_opts = PlotOptions(plot_input_dir+"_"+plot_name,
                                    options,
                                    do_data,
                                    do_data_mc,
                                    invert_data_mc_ratio,
                                    do_ratio,
                                    do_norm,
                                    do_norm_mc_to_data,
                                    do_auto_scale,
                                    draw_style,
                                    marker_type,
                                    leg_labels,
                                    plot_tags,
                                    atlas_label,
                                    lumi_label,
                                    False,
                                    draw_avg_ratio)

            hdata.Rebin2D(plot_opts.rebinx, plot_opts.rebiny)
            hdata_p = tprofile2hist(hdata.ProfileX(hdata.GetName()+"_profile_data"))
            hmc.Rebin2D(plot_opts.rebinx, plot_opts.rebiny)
            hmc_p = tprofile2hist(hmc.ProfileX(hmc.GetName()+"_profile_mc"))

            hdata_bg.Rebin2D(plot_opts.rebinx, plot_opts.rebiny)
            hdata_bg_p = tprofile2hist(hdata_bg.ProfileX(hdata_bg.GetName()+"_profile_data_bg"))
            hmc_bg.Rebin2D(plot_opts.rebinx, plot_opts.rebiny)
            hmc_bg_p = tprofile2hist(hmc_bg.ProfileX(hmc_bg.GetName()+"_profile_mc_bg"))

            hdata_p.Add(hdata_bg_p,-4/3)
            hmc_p.Add(hmc_bg_p,-4/3)

            plot_output_path = output_dir+"/"+plot_input_dir+"/"+file_tag+"_"+plot_name

            save_hists_overlay([hdata_p, hmc_p], plot_opts, plot_output_path)

def plot_1d_hists_eop_eta_p_bins(datafile, mcfile, plot_input_dir, plot_input_list,
                                 output_dir, file_tag, sample_tag):

    data = TFile.Open(datafile, "READ")
    mc = TFile.Open(mcfile, "READ")

    # create output dir(s) if necessary
    ensure_dir(output_dir+"/"+plot_input_dir+"_eop_eta_p_bins")

    # list of plots to make
    plot_names = [line.rstrip('\n') for line in open(plot_input_list) if "#" not in line]

    # set basic plotting options
    options = "" # pass string options to PlotOptions
    do_data = True # plot data as solid black markers
    do_data_mc = True # plot data_mc_ratio
    invert_data_mc_ratio = True # invert data mc ratio order
    do_ratio = True  # add a ratio plot under main plot
    do_norm = False # normalize all plots to unity
    do_norm_mc_to_data = False #
    do_auto_scale = True # scale y-axis appropriately
    draw_style = "hist" # pass args to h.Draw("args"), NB: not applicable for do_data
    marker_type = "data_mc" #sets the markers and colors
    leg_labels = ["MC 2015", "MC15c"] #, "MC (Pythia JZ1W)"]
    atlas_label = "Internal"
    lumi_label = "#sqrt{s} = 13 TeV, 1.6 nb^{-1}"

    options = "eop_eta_p_bins"

    eta_bins = ["etaL06", "etaG06L11"]

    for eta_bin in eta_bins:
        for plot_name in plot_names:
            hdata = data.Get(plot_input_dir+"_"+eta_bin+"/"+plot_name)
            hmc = mc.Get(plot_input_dir+"_"+eta_bin+"/"+plot_name)
            print plot_input_dir+"_"+eta_bin+"/"+plot_name

            print hdata

            if "TH2" in str(hdata):
                for i in xrange(1,hdata.GetNbinsX()+1):
                    xlow = i-1
                    xhigh = i
                    xlow_val = hdata.GetXaxis().GetBinLowEdge(xlow+1)
                    xhigh_val = hdata.GetXaxis().GetBinLowEdge(xhigh)+hdata.GetXaxis().GetBinWidth(xhigh)
                    p_bin = "pG%dL%d" % (xlow_val*1000, xhigh_val*1000)

                    hdata_py = hdata.ProjectionY(hdata.GetName()+"_py_data",xlow,xhigh)
                    hmc_py = hmc.ProjectionY(hmc.GetName()+"_py_mc",xlow,xhigh)

                    plot_name_new = eta_bin+"_"+p_bin+"_"+re.search("^.*(?=_vs_trkP)", plot_name).group()
                    plot_output_path = output_dir+"/"+plot_input_dir+"_eop_eta_p_bins/"+file_tag+"_"+plot_name_new

                    if sample_tag != "":
                        plot_tags = [sample_tag]
                    else:
                        plot_tags = []

                    plot_opts = PlotOptions(plot_input_dir+"_"+plot_name_new,
                                            options,
                                            do_data,
                                            do_data_mc,
                                            invert_data_mc_ratio,
                                            do_ratio,
                                            do_norm,
                                            do_norm_mc_to_data,
                                            do_auto_scale,
                                            draw_style,
                                            marker_type,
                                            leg_labels,
                                            plot_tags,
                                            atlas_label,
                                            lumi_label)

                    save_hists_overlay([hdata_py, hmc_py], plot_opts, plot_output_path)

