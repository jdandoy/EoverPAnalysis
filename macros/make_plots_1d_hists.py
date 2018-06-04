#! /usr/bin/env python
import sys
import re
import os

def main(argv):

    # I/O
    # ---------------------------------
    # Modify based on local setup
    # run = "20170211_0"
    # run_pt_reweight = "20170211_0_pt_reweight"
    # run_pu = "20170210_0_pu_and_pt_reweight"

    # run_data = "20170527_0_tilepublic"
    # run_mc = "20170527_0_tilepublic"

    # filename_data = "condor_all_eop_lowmu_data_20170527_0_tilepublic.root"
    # filename_mc = "condor_all_eop_lowmu_mc_merged_20170527_0_tilepublic.root"

    # filename_data = "condor_all_eop_lowmu_data_20170528_0_tilepublic.root"
    # filename_mc = "condor_all_eop_lowmu_mc_merged_20170528_0_tilepublic_no_pt_reweight.root"

    # filename_data = "condor_all_eop_lowmu_data_20170528_0_tilepublic.root"
    # filename_mc = "condor_all_eop_lowmu_mc_merged_20170528_0_tilepublic.root"

    # filename_data = "condor_all_eop_pileup_data_20170529_0_tile_pileup.root"
    # filename_mc = "condor_all_eop_pileup_mc_merged_20170603_1_tile_pileup_no_pt_reweight.root"

    # filename_data = "condor_all_eop_pileup_data_20170529_0_tile_pileup.root"
    # filename_mc = "condor_all_eop_pileup_JZ1W_20170529_0_tile_pileup_no_pt_reweight.root"

    #filename_data = "condor_all_eop_pileup_data_20170529_0_tile_pileup.root"
    # filename_mc = "condor_all_eop_pileup_mc_merged_20170603_1_tile_pileup_no_pt_reweight.root"
    #filename_mc = "condor_all_eop_pileup_mc_merged_20170603_2_tile_pileup.root"

    filename_data = "condor_all_eop_lowmu_runII_general_data_20170529_0_deconvinputs.root"
    # filename_mc = "condor_all_eop_lowmu_runII_general_mc_merged_20170529_0_deconvinputs.root"
    filename_mc = "condor_all_eop_lowmu_runII_general_mc_merged_20170530_1_deconvinputs.root"

    eop_dir = "/Users/joakim/eoverp/"
    input_path = eop_dir+"results/histogramfiles/"
    output_path = eop_dir+"results/plots/"
    # ---------------------------------

    # Config
    # ---------------------------------
    do_pileup = False
    if (len(argv) > 0):
        if argv[1] == "pileup":
            do_pileup = True
    do_eop_logy = True
    do_eop_avg_g0 = True
    energy_calibs = ["ClusterEnergy"]
    # energy_calibs = ["ClusterEnergy", "ClusterEnergyLCW"]
    # energy_calibs = ["ClusterEnergy", "ClusterEnergyLCW", "CellEnergy"]
    # energy_calibs = ["ClusterEnergyLCW", "CellEnergy"]
    # energy_calibs = ["ClusterEnergyLCW"]
    # energy_calibs = ["CellEnergy"]
    # ---------------------------------

    # Default setup (low-mu samples)
    data = os.path.join(input_path, filename_data)
    mc = os.path.join(input_path, filename_mc)
    folder_tag = re.search("(?<=_)\d{1,8}_.+(?=\.)", data).group()
    plot_input_list = eop_dir+"/scripts/plotlists/plotlist_eop_1d_hists.txt"
    output_dir = output_path+"/eop_lowmu_data_mc_"+folder_tag+"/"
    file_tag = "overlay_lowmu_data_mc"
    # sample_tag = "13 TeV 2015 data from period B1 low-#mu run"
    sample_tag = ""
    leg_labels = ["Data 2015, Low-#LT#font[50]{#mu}#GT", "Pythia MinBias"]
    lumi_label = "#sqrt{s} = 13 TeV, 1.6 nb^{-1}"
    options = ""
    #options = "public"

    if (do_pileup):
        # mc = input_path+"/condor_all_eop_pileup_JZ1W_20170210_0_pu_reweight.root"
        # folder_tag = "eop_pileup_data_mc_"+runtag
        file_tag = "overlay_pileup_data_mc"
        output_dir = output_path+"/eop_pileup_data_mc_"+folder_tag+"/"
        # sample_tag = "13 TeV 2015 data from ZeroBias stream"
        lumi_label = "#sqrt{s} = 13 TeV, 0.42 pb^{-1}"
        leg_labels = ["Data 2015", "Pythia Multijet"]
        options += "pileup"

    if (do_eop_logy):
        file_tag_logy = file_tag+"_logy"
        plot_input_list_logy = eop_dir+"/scripts/plotlists/plotlist_eop_1d_hists_logy.txt"
        options_logy = "1d_eop_logy_"+options

    if (do_eop_avg_g0):
        file_tag_avg_g0 = file_tag+"_avg_g0"
        plot_input_list_avg_g0 = eop_dir+"/scripts/plotlists/plotlist_eop_1d_hists_avg_g0.txt"
        options_avg_g0 = "1d_eop_avg_g0_"+options


    for energy_calib in energy_calibs:

        selections = [# "EoverP_Run1paper_noSelection",
                      "EoverP_LoosePrimaryTrks_"+energy_calib+"_Run1paper",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Run1paper_pG1200L1800",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Run1paper_pG2800L3600",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Run1paper_pG4600L5600",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_notrkP_noLar_noTileEfrac",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_2GeVTrkP_noLar_noTileEfrac",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_2GeVTrkP_1GeVLar_noTileEfrac",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_TileEfrac000",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_TileEfrac010",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_TileEfrac030",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_TileEfrac050",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_TileEfrac060",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_TileEfrac070",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_TileEfrac075",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_TileEfrac080",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_pL4000",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_pG4000L8000",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_pG8000L12000",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_pG12000",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_etaL05",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_etaG05L07",
                      # "EoverP_LoosePrimaryTrks_"+energy_calib+"_Tile_defaultCuts_etaG07",
                      ]

        selection_tags = [#"run1paper_noselection",
                          "run1paper",
                          # "run1paper_pG1200L1800",
                          # "run1paper_pG2800L3600",
                          # "run1paper_pG4600L5600",
                          # "tile_notrkp_nolar_notileE",
                          # "tile_trkpG2_nolar_notileE",
                          # "tile_trkpG2_larL1_notileE",
                          # "tilecuts",
                          # "tilecuts_Efrac000",
                          # "tilecuts_Efrac010",
                          # "tilecuts_Efrac030",
                          # "tilecuts_Efrac050",
                          # "tilecuts_Efrac060",
                          # "tilecuts_Efrac070",
                          # "tilecuts_Efrac075",
                          # "tilecuts_Efrac080",
                          # "tilecuts_pL4000",
                          # "tilecuts_pG4000L8000",
                          # "tilecuts_pG8000L12000",
                          # "tilecuts_pG12000",
                          # "tilecuts_etaL05",
                          # "tilecuts_etaG05L07",
                          # "tilecuts_etaG07"
                          ]

        print "----> plot config: "
        print "selections:", selections
        print "data:", data
        print "mc:", mc
        print "folder_tag:", folder_tag
        print "plot_input_list:",plot_input_list
        print "output_dir:",output_dir
        print "file_tag:",file_tag
        print "sample_tag:",sample_tag
        print

        print "----> running script to save plots: "
        for i,selection in enumerate(selections):
            print "selections: "+selection
            file_tag_new = file_tag+"_"+selection_tags[i]
            from plot_maker import plot_1d_hists
            plot_1d_hists(data, mc, selection, plot_input_list,
                          output_dir, file_tag_new, sample_tag, options, leg_labels, lumi_label)
            if (do_eop_logy):
                file_tag_logy_new = file_tag_logy+"_"+selection_tags[i]
                plot_1d_hists(data, mc, selection, plot_input_list_logy,
                              output_dir, file_tag_logy_new, sample_tag, options_logy, leg_labels, lumi_label)
            if (do_eop_avg_g0):
                file_tag_logy_new = file_tag_avg_g0+"_"+selection_tags[i]
                plot_1d_hists(data, mc, selection, plot_input_list_avg_g0,
                              output_dir, file_tag_logy_new, sample_tag, options_avg_g0, leg_labels, lumi_label)

if __name__ == "__main__":
    main(sys.argv)
