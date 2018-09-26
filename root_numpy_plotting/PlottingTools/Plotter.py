from atlasplots import atlas_style as astyle
astyle.SetAtlasStyle()
import numpy as np
from variables.variables import calc_weight

from root_numpy import fill_hist

import ROOT

from DataPrep import FetchResults, FetchResultsSingleProcessor

def partitionDataset(N, partitionSize=2000000):
    '''given a dataset with length end, return tuples that partition N'''
    partitions=[(0, partitionSize)]
    counter = 1
    while partitions[-1][-1] < N:
        partitions.append((counter*partitionSize, (counter+1) * partitionSize))
        counter += 1

    partitions[-1] = ((counter-1)*partitionSize, N)
    return partitions

def DrawText(x, y, text, color=1, size=0.05):
    """Draw text.

    Parameters
    ----------
    x : float
        x position in NDC coordinates
    y : float
        y position in NDC coordinates
    text : string, optional
        The text
    color : int, optional
        Text colour (the default is 1, i.e. black).
        See https://root.cern.ch/doc/master/classTColor.html.
        If you know the hex code, rgb values, etc., use ``ROOT.TColor.GetColor()``
    size : float, optional
        Text size
        See https://root.cern.ch/doc/master/classTLatex.html
    """
    l = ROOT.TLatex()
    l.SetTextSize(size)
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x, y, text)

def handle_underflow_overflow(h):
    """
    add the underflow and overflow (underflow = 0, overflow=nbins+1)
    """

    # underflow
    underflow = h.GetBinContent(0)
    first_bin = h.GetBinContent(1)
    h.SetBinContent(1, first_bin+underflow)
    h.SetBinContent(0, 0) # remove underflow
    # underflow error
    underflow_error = h.GetBinError(0)
    first_bin_error = h.GetBinError(1)
    h.SetBinError(1, (underflow_error**2 + first_bin_error**2)**0.5)

    # overflow
    last_bin = h.GetNbinsX()
    last_bin_content = h.GetBinContent(last_bin)
    overflow = h.GetBinContent(last_bin+1)
    h.SetBinContent(last_bin, last_bin_content+overflow)
    h.SetBinContent(last_bin+1, 0) # remove overflow
    # overflow error
    last_bin_error = h.GetBinError(last_bin)
    overflow_error = h.GetBinError(last_bin+1)
    h.SetBinError(last_bin, (last_bin_error**2 + overflow_error**2)**0.5)

    return h

def batchSubmissions(list_for_batch, batch_size):
    '''take a list of submission configurations, and batch them into groups of batch-size. batch_size would typically be the number of processes available for the jop'''
    return [list_for_batch[x:x+batch_size] for x in xrange(0, len(list_for_batch), batch_size)]

#Keep things alive at the global scope
global_scope = []

def toGlobalScope(obj):
    global_scope.append(obj)

def DrawDataVsMC(histogram_dict, LegendLabels = {}, MCKey = "", DataKey = "", doLogy = True, ratio_min= 0.8, ratio_max = 1.2):
    ''' Draw the data vs MC ratio for the MC and data histograms'''

    MCHist = histogram_dict[MCKey]
    DataHist = histogram_dict[DataKey]

    legend = ROOT.TLegend(0.60, 0.65, 0.92, 0.89)
    legend.SetBorderSize(0)
    toGlobalScope(legend)

    canvas = ROOT.TCanvas(MCKey + DataKey + MCHist.GetTitle(), MCKey + DataKey + MCHist.GetTitle(), 1300, 800)

    top_pad = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    top_pad.Draw()
    top_pad.cd()
    top_pad.SetBottomMargin(0)

    filename = MCHist.GetTitle() + "histogram"

    maximum_bin = 0.0
    minimum_bin = 10000000000000000.0
    for bin in range(1, MCHist.GetNbinsX() + 1):
        content = MCHist.GetBinContent(bin)
        if content > maximum_bin:
            maximum_bin = content
        if content < minimum_bin and minimum_bin > 0.0:
            minimum_bin = content


    if doLogy:
        filename += "logy"

    if doLogy:
        if minimum_bin <= 0.0:
            minimum_bin = 0.011
        MCHist.SetMaximum(maximum_bin * 5000)
        MCHist.SetMinimum(minimum_bin * 1.0)

    else:
        MCHist.SetMaximum(maximum_bin * 1.8)
        MCHist.SetMinimum(0.0001)


    MCHist.Draw()
    DataHist.Draw("SAME")

    legend.AddEntry(MCHist, LegendLabels[MCKey])
    legend.AddEntry(DataHist, LegendLabels[DataKey])
    legend.SetTextSize(0.04)
    legend.Draw()

    if doLogy:
        top_pad.SetLogy()

    ROOT.gROOT.SetStyle("ATLAS")
    astyle.ATLASLabel(0.2, 0.87, "Internal")
    top_pad.Modified()
    top_pad.Update()
    toGlobalScope(top_pad)

    canvas.cd()
    bottom_pad = ROOT.TPad("pad2", "pad2", 0, 0.01, 1, 0.3)
    #bottom_pad.SetRightMargin(0.15)
    bottom_pad.Draw()
    bottom_pad.cd()
    toGlobalScope(bottom_pad)

    data_ratio = DataHist.Clone("data_histogram")
    data_ratio.Divide(MCHist)
    data_ratio.Draw()

    MCHist_label_size = MCHist.GetXaxis().GetLabelSize()
    variableLabel = MCHist.GetXaxis().GetTitle()

    data_ratio.GetYaxis().SetTitle("Data/MC")
    scale_ratio = (top_pad.GetWh()*top_pad.GetAbsHNDC())/(bottom_pad.GetWh() * bottom_pad.GetAbsHNDC())
    data_ratio.GetXaxis().SetLabelSize(MCHist_label_size*(scale_ratio))
    data_ratio.GetYaxis().SetLabelSize(MCHist_label_size*(scale_ratio))
    data_ratio.GetXaxis().SetTitle(variableLabel)
    data_ratio.GetXaxis().SetTitleOffset(1.1)
    data_ratio.SetMaximum(ratio_max - 0.0001)
    data_ratio.SetMinimum(ratio_min + 0.0001)
    data_ratio.GetXaxis().SetTitleSize(MCHist_label_size*scale_ratio)
    data_ratio.GetYaxis().SetTitleSize(MCHist_label_size*scale_ratio)
    data_ratio.GetYaxis().SetTitleOffset(0.25)
    data_ratio.GetYaxis().SetNdivisions(404);
    data_ratio.SetMaximum(ratio_max - 0.0001)
    data_ratio.SetMinimum(ratio_min + 0.0001)
    data_ratio.Draw()
    toGlobalScope(data_ratio)

    straight_line = ROOT.TF1("line", str(1.0) , -10e6, + 10e6)
    straight_line.SetLineWidth(2)
    straight_line.Draw("Same")
    toGlobalScope(straight_line)
    ##input("Does this look OK?")

    straight_line_up = ROOT.TF1("line",  str(1.0 + (1.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_up.SetLineWidth(1)
    straight_line_up.SetLineStyle(1)
    straight_line_up.Draw("Same")
    toGlobalScope(straight_line)

    straight_line_up2 = ROOT.TF1("line",  str(1.0 + (2.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_up2.SetLineWidth(1)
    straight_line_up2.SetLineStyle(3)
    straight_line_up2.Draw("Same")
    toGlobalScope(straight_line_up2)

    straight_line_up3 = ROOT.TF1("line",  str(1.0 + (3.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_up3.SetLineWidth(1)
    straight_line_up3.SetLineStyle(3)
    straight_line_up3.Draw("Same")
    toGlobalScope(straight_line_up3)

    straight_line_down3 = ROOT.TF1("line",  str(1.0 - (3.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_down3.SetLineWidth(1)
    straight_line_down3.SetLineStyle(3)
    straight_line_down3.Draw("Same")
    toGlobalScope(straight_line_down3)

    straight_line_down = ROOT.TF1("line",  str(1.0 - (2.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_down.SetLineWidth(1)
    straight_line_down.SetLineStyle(1)
    straight_line_down.Draw("Same")
    toGlobalScope(straight_line_down)

    straight_line_down2 = ROOT.TF1("line",  str(1.0 - (1.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_down2.SetLineWidth(1)
    straight_line_down2.SetLineStyle(3)
    straight_line_down2.Draw("Same")
    toGlobalScope(straight_line_down2)

    bottom_pad.Modified()
    bottom_pad.Update()

    bottom_pad.Draw()
    bottom_pad.cd()
    bottom_pad.SetTopMargin(0)
    bottom_pad.SetBottomMargin(0.45)

    return canvas


class Plotter:
    def __init__(self, inputs, treeName):
        self.channelFiles = {}
        self.channelLabels = {}
        self.AddInputDictionary(inputs)
        self.treeName = treeName
        self.partitionSize = 100000000
        self.Process = 1
        self.test = False

    def AddInputDictionary(self, dictionary):
        '''Take an input dictionary of channel names to: tuple of (channel descriptor for legend, list of [input root filename strings])'''
        self.channels = dictionary.keys()
        for channel in self.channels:
            self.channelFiles[channel] = dictionary[channel][1]
            self.channelLabels[channel] = dictionary[channel][0]

    def GetHistograms(self, variable, list_selections = [], nBins = 1, range_low = 0, range_high = 1, xlabel ="", ylabel = ""):
        '''given a variable, Draw the histogram for the given variable'''
        variableName = variable.name

        #First go and get all of the histograms that we need
        histogram_dictionary = {}
        LegendLabels = {}
        for channel in self.channels:
            histogram_dictionary[channel] = ROOT.TH1D(variableName + "_" + channel, variableName + "_" + channel, nBins, range_low, range_high)
            histogram_dictionary[channel].GetXaxis().SetTitle(xlabel)
            histogram_dictionary[channel].GetYaxis().SetTitle(ylabel)
            histogram_dictionary[channel].Sumw2()

        for channel in self.channels:
            for filename in self.channelFiles[channel]:
                print "Reading from file " + filename
                input_file = ROOT.TFile(filename, "READ")
                input_file.ls()

                print "Getting tree " + self.treeName
                input_tree = input_file.Get(self.treeName)
                entries = input_tree.GetEntries()
                input_file.Close()
                del input_tree
                del input_file

                if self.test:
                     entries = 10000

                partitions = partitionDataset(entries, self.partitionSize) #partition the reading of the tree to save disk space

                branches = []
                for branch in variable.branches:
                    if branch not in branches:
                        branches.append(branch)

                for selection in list_selections:
                    for branch in selection.branches:
                        if branch not in branches:
                            branches.append(branch)

                #The configuration for the fetch results function.
                configuration_dictionary = {}
                configuration_dictionary["selections"] = list_selections
                configuration_dictionary["branches"] = branches
                configuration_dictionary["variables"] = [variable]
                configuration_dictionary["filename"] = filename
                configuration_dictionary["treename"] = self.treeName
                configuration_dictionary["weightFunction"] = calc_weight

                submission_inputs = []
                for partition in partitions:
                    submission_inputs.append((partition, configuration_dictionary))

                submission_inputs = batchSubmissions(submission_inputs, self.Process)

                for submission_batch in submission_inputs:
                    if self.Process > 1:
                        #If you have configured the job to have more than one procress, call the job pool
                        results = FetchResults(submission_batch, self.Process)

                    else:
                        #Otherwise just use the map function. 
                        results = FetchResultsSingleProcessor(submission_batch)

                    for result in results:
                        selection_dict = result["selection_dict"]
                        variable_dict = result["variable_dict"]
                        weights = variable_dict["weights"]

                        #cutflows.Fill("no selection", len(weights))
                        #cutflowsWeighted.Fill("no selection", np.sum(weights))

                        #prepare to loop through the selections and apply all of them
                        print("\n\n\n\n\n===============\nPreselections are being applied")
                        print("Pre-selection event count" + str(len(weights)))
                        total_selection = np.ones(len(weights)) > 0.0
                        for selection in list_selections:
                            print("applying selection " + str(selection.name))
                            print(selection_dict[selection.name])
                            total_selection &= selection_dict[selection.name]
                            #cutflowsWeighted.Fill(selection_name, np.sum(weights[total_selection]))
                            #cutflows.Fill(selection_name, np.sum(1*total_selection))

                        print("post preselection event count" + str(np.sum(1*total_selection)))

                        #for the events that pass the selections, fill the events into the histograms
                        to_fill = variable_dict[variable.name][total_selection]
                        to_weight = weights[total_selection]
                        print(len(to_fill))
                        print(len(to_weight))
                        print to_fill
                        print to_weight
                        print("Filling Variable " + variable.name)
                        fill_hist(histogram_dictionary[channel], to_fill, to_weight)

        return histogram_dictionary








