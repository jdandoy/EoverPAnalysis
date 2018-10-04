from atlasplots import atlas_style as astyle
astyle.SetAtlasStyle()
import numpy as np
from variables.variables import calc_weight
import array
from root_numpy import fill_hist
import ROOT
from DataPrep import GetData
import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True)


from util import getWeightsFromBins

def WeightsToNormalizeToHistogram(variable_in_histogram, histogram):
    '''normalize the weights to the histogram'''
    low_edges = []
    high_edges = []
    normalizations = []
    #get the low bin edges of the histograms
    for bin in range(1, histogram.GetNbinsX() + 1):
        low_edges.append(histogram.GetBinLowEdge(bin))
        high_edges.append(histogram.GetBinLowEdge(bin + 1))
        normalizations.append(histogram.GetBinContent(bin))

    low_edges = np.array(low_edges, np.float)
    high_edges = np.array(high_edges, np.float)
    hist_weights = np.array(normalizations, np.float)
    ##There must be a very clever way to implement this in root numpy... I'm sure of it!
    weights = getWeightsFromBins(variable_in_histogram, low_edges, high_edges, hist_weights)

    return weights

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

#Keep things alive at the global scope
global_scope = []

def toGlobalScope(obj):
    global_scope.append(obj)

def cleanUpHistograms(h):
    h.SetLineStyle(1)
    h.SetLineWidth(3)
    h.SetFillStyle(0)
    h.SetFillColor(0)
    h.SetMarkerSize(0)
    return h

MCColor = ROOT.kRed
DataColor = ROOT.kBlack

def DrawDataVsMC(histogram_dict, LegendLabels = {}, MCKey = "", DataKey = "", doLogy = True, ratio_min= 0.0, ratio_max = 2.0, extra_description = None, extra_desx = None, extra_desy = None):
    ''' Draw the data vs MC ratio for the MC and data histograms'''

    MCHist = histogram_dict[MCKey]
    DataHist = histogram_dict[DataKey]

    MCHist = cleanUpHistograms(MCHist)

    MCHist.SetLineColor(MCColor)
    DataHist.SetLineColor(DataColor)

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
        MCHist.SetMaximum(maximum_bin * 500)
        MCHist.SetMinimum(minimum_bin * 1.0)

    else:
        MCHist.SetMaximum(maximum_bin * 1.8)
        MCHist.SetMinimum(0.0001)


    MCHist.Draw("HIST")
    DataHist.Draw("SAME")

    legend.AddEntry(MCHist, LegendLabels[MCKey])
    legend.AddEntry(DataHist, LegendLabels[DataKey])
    legend.SetTextSize(0.04)
    legend.Draw()

    if extra_description:
        pass

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
    data_ratio = cleanUpHistograms(data_ratio)
    data_ratio.SetLineColor(ROOT.kBlack)
    data_ratio.Draw("HIST E")

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
    data_ratio.GetYaxis().SetNdivisions(405);
    data_ratio.SetMaximum(ratio_max - 0.0001)
    data_ratio.SetMinimum(ratio_min + 0.0001)
    data_ratio.Draw("HIST E")
    toGlobalScope(data_ratio)

    straight_line = ROOT.TF1("line1", str(1.0) , -10e6, + 10e6)
    straight_line.SetLineWidth(2)
    straight_line.Draw("Same")
    toGlobalScope(straight_line)
    ##input("Does this look OK?")

    straight_line_up = ROOT.TF1("line2",  str(1.0 + (2.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_up.SetLineWidth(1)
    straight_line_up.SetLineStyle(1)
    straight_line_up.Draw("Same")
    toGlobalScope(straight_line_up)

    straight_line_up2 = ROOT.TF1("line3",  str(1.0 + (1.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_up2.SetLineWidth(1)
    straight_line_up2.SetLineStyle(3)
    straight_line_up2.Draw("Same")
    toGlobalScope(straight_line_up2)

    straight_line_up3 = ROOT.TF1("line4",  str(1.0 + (3.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_up3.SetLineWidth(1)
    straight_line_up3.SetLineStyle(3)
    straight_line_up3.Draw("Same")
    toGlobalScope(straight_line_up3)

    straight_line_down3 = ROOT.TF1("line5",  str(1.0 - (3.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_down3.SetLineWidth(1)
    straight_line_down3.SetLineStyle(3)
    straight_line_down3.Draw("Same")
    toGlobalScope(straight_line_down3)

    straight_line_down = ROOT.TF1("line6",  str(1.0 - (2.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
    straight_line_down.SetLineWidth(1)
    straight_line_down.SetLineStyle(1)
    straight_line_down.Draw("Same")
    toGlobalScope(straight_line_down)

    straight_line_down2 = ROOT.TF1("line7",  str(1.0 - (1.0 * (ratio_max - 1.0)/4)) , -10e6, + 10e6)
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

def GetListOfNeededBranches(variables, selections):
    branches = []

    for variable in variables:
        for branch in variable.branches:
            if branch not in branches:
                branches.append(branch)

    for selection in selections:
        for branch in selection.branches:
            if branch not in branches:
                branches.append(branch)

    return branches

def DivideHistograms(hist_dict1, hist_dict2):
    dict1_keys = hist_dict1.keys()
    dict2_keys = hist_dict2.keys()

    for key in dict1_keys:
        if key not in dict2_keys:
            raise ValueError("The two dictionaries do not have the same keys")
    for key in dict2_keys:
        if key not in dict1_keys:
            raise ValueError("The two dictionaries do not have the same keys")

    return_dict = {}
    for channel in dict1_keys:
        Hist_clone1 = dict1_keys[channel].Clone(channel + "divided")
        return_dict[channel] = Hist_clone1.Divide(dict2_keys[channel])

    return return_dict

class Plotter:
    def __init__(self, inputs, treeName, base_selections = ""):
        self.channelFiles = {}
        self.channelLabels = {}
        self.AddInputDictionary(inputs)
        self.treeName = treeName
        self.partitionSize = 10000000000
        self.Process = 1
        self.normalization_dictionary = {}
        self.test = False
        self.verbose = False
        self.base_selections = base_selections
        self.NormalizationWeightsDictionary = {} # A dictionary of channel to variable to the weights needed for the reweighting of this variable

    def AddInputDictionary(self, dictionary):
        '''Take an input dictionary of channel names to: tuple of (channel descriptor for legend, list of [input root filename strings])'''
        self.channels = dictionary.keys()
        for channel in self.channels:
            self.channelFiles[channel] = dictionary[channel][1]
            self.channelLabels[channel] = dictionary[channel][0]

    def CheckForNormalizationWeights(self, channel, filename):
        if channel not in self.NormalizationWeightsDictionary: 
            print "nor normalization weights found for channel " + filename
            return False
        if filename not in self.NormalizationWeightsDictionary[channel]: raise ValueError("There should have been a set of normalization weights for this file")

        return True

    def GetNumberOfTracks(self, channel):
        '''Get the total weighted number of Tracks. This is to normalize the distributions to each other'''

        list_selections = []
        branches = []
        variables = []

        weightedNumber = 0

        for filename in self.channelFiles[channel]:
            if self.verbose: print "Reading from file " + filename
            input_file = ROOT.TFile(filename, "READ")
            if self.verbose: input_file.ls()

            if self.verbose: print "Getting tree " + self.treeName
            input_tree = input_file.Get(self.treeName)
            total_entries = input_tree.GetEntries()
            selected_entries = input_tree.GetEntries(self.base_selections)
            input_file.Close()
            del input_tree
            del input_file

            branches = GetListOfNeededBranches(variables, list_selections)
            #The configuration for the fetch results function.
            result = GetData(partition = (0, total_entries), bare_branches = branches, filename = filename, treename = self.treeName, variables=variables, weightCalculator = calc_weight, selections = list_selections, selection_string = self.base_selections, verbose = self.verbose)

            variable_dict = result["variable_dict"]
            weights = variable_dict["weights"]
            if channel in self.NormalizationWeightsDictionary:
                weights = weights * self.NormalizationWeightsDictionary[channel][filename] #get the normalization weights to make distributions agree between data and MC
            weightedNumber += np.sum(weights)

        return weightedNumber

    def SetTotalTrackNumberNormalization(self, channel, normalization):
        if channel not in self.NormalizationWeightsDictionary:
            self.NormalizationWeightsDictionary[channel] = {}
            for filename in self.channelFiles[channel]:
                if filename not in self.NormalizationWeightsDictionary[channel]:
                   input_file = ROOT.TFile(filename, "READ")
                   if self.verbose: input_file.ls()
                   if self.verbose: print "Getting tree " + self.treeName
                   input_tree = input_file.Get(self.treeName)
                   total_entries = input_tree.GetEntries()
                   selected_entries = input_tree.GetEntries(self.base_selections)
                   input_file.Close()
                   del input_tree
                   del input_file
                   self.NormalizationWeightsDictionary[channel][filename] = np.ones(selected_entries) * normalization
        else:
            for filename in self.channelFiles[channel]:
                self.NormalizationWeightsDictionary[channel][filename] *= normalization


    def UseVariableAndHistogramToNormalize(self, variable, HistDict, ChannelToNormalize, ChannelToNormalizeTo):
         '''take a variable, and a histogram, and set up the plotter so that it always normalizes the ChannelToNormalize to the ChannelToNormalilzeTo'''
         #Take the ratio of the histograms that we want to use for normalization
         TargetHist = HistDict[ChannelToNormalizeTo]
         UnNormalizedHist = HistDict[ChannelToNormalize]

         ratio_hist = TargetHist.Clone("NormalizationHistogram" + variable.name + ChannelToNormalize)
         ratio_hist.Divide(UnNormalizedHist)

         if ChannelToNormalize not in self.NormalizationWeightsDictionary:
            self.NormalizationWeightsDictionary[ChannelToNormalize] = {}

         for filename in self.channelFiles[ChannelToNormalize]:
            print "Reading from file " + filename
            print "Corresponding to channel " + ChannelToNormalize
            input_file = ROOT.TFile(filename, "READ")

            if self.verbose: input_file.ls()

            if self.verbose: print "Getting tree " + self.treeName
            input_tree = input_file.Get(self.treeName)
            total_entries = input_tree.GetEntries()
            selected_entries = input_tree.GetEntries(self.base_selections)
            input_file.Close()
            del input_tree
            del input_file

            if filename not in self.NormalizationWeightsDictionary[ChannelToNormalize]:
                print "Filename not found in the weights dictionary"
                self.NormalizationWeightsDictionary[ChannelToNormalize][filename] = np.ones(selected_entries)

            variables = [variable]
            list_selections = []

            branches = GetListOfNeededBranches(variables, list_selections)
            #The configuration for the fetch results function.
            result = GetData(partition = (0, total_entries), bare_branches = branches, filename = filename, treename = self.treeName, variables=variables, weightCalculator = calc_weight, selections = list_selections, selection_string = self.base_selections,  verbose = self.verbose)
            variable_in_histogram = result["variable_dict"][variable.name]
            variable_in_histogram.shape
            normalization = WeightsToNormalizeToHistogram(variable_in_histogram, ratio_hist)
            print normalization.shape
            #store the normalization for this channel:
            self.NormalizationWeightsDictionary[ChannelToNormalize][filename] *= normalization


    def GetHistograms(self, variable, list_selections = [], bins = 1, range_low = 0, range_high=1,  xlabel ="", ylabel = ""):
        '''given a variable, Draw the histogram for the given variable'''
        variableNameToFill = variable.name
        variables = [variable]

        description_string = "Histogram"
        for variable in variables:
            description_string += variable.name
        for sel in list_selections:
            description_string += sel.name

        #First go and get all of the histograms that we need
        histogram_dictionary = {}
        LegendLabels = {}
        for channel in self.channels:
            if type(bins) == list:
                bins = array('d',bins)
                histogram_dictionary[channel] = ROOT.TH1D(description_string, description_string, len(bins)-1, bins)
            else:
                histogram_dictionary[channel] = ROOT.TH1D(description_string, description_string, bins, range_low, range_high)
            histogram_dictionary[channel].GetXaxis().SetTitle(xlabel)
            histogram_dictionary[channel].GetYaxis().SetTitle(ylabel)
            histogram_dictionary[channel].Sumw2()

        for channel in self.channels:
            for filename in self.channelFiles[channel]:
                print "Reading from file " + filename
                print "Corresponding to channel " + channel
                input_file = ROOT.TFile(filename, "READ")

                if self.verbose: input_file.ls()

                if self.verbose: print "Getting tree " + self.treeName
                input_tree = input_file.Get(self.treeName)
                total_entries = input_tree.GetEntries()
                selected_entries = input_tree.GetEntries(self.base_selections)
                input_file.Close()
                del input_tree
                del input_file

                branches = GetListOfNeededBranches(variables, list_selections)

                #The configuration for the fetch results function.
                result = GetData(partition = (0, total_entries), bare_branches = branches, filename = filename, treename = self.treeName, variables=variables, weightCalculator = calc_weight, selections = list_selections, selection_string = self.base_selections, verbose = self.verbose)

                selection_dict = result["selection_dict"]
                variable_dict = result["variable_dict"]
                weights = variable_dict["weights"]

                toNormalize = self.CheckForNormalizationWeights(channel, filename)
                if toNormalize:
                    weights = weights * self.NormalizationWeightsDictionary[channel][filename]

                #prepare to loop through the selections and apply all of them
                if self.verbose: print("\n\n\n\n\n===============\nPreselections are being applied")
                if self.verbose: print("Pre-selection event count" + str(len(weights)))
                total_selection = np.ones(len(weights)) > 0.0
                for selection in list_selections:
                    if self.verbose: print("applying selection " + str(selection.name))
                    if self.verbose: print(selection_dict[selection.name])
                    total_selection &= selection_dict[selection.name]
                    #cutflowsWeighted.Fill(selection_name, np.sum(weights[total_selection]))
                    #cutflows.Fill(selection_name, np.sum(1*total_selection))

                if self.verbose: print("post preselection event count" + str(np.sum(1*total_selection)))

                #for the events that pass the selections, fill the events into the histograms
                to_fill = variable_dict[variableNameToFill][total_selection]
                to_weight = weights[total_selection]
                if self.verbose: print(len(to_fill))
                if self.verbose: print(len(to_weight))
                if self.verbose: print to_fill
                if self.verbose: print to_weight
                if self.verbose: print("Filling Variable " + variable.name)
                fill_hist(histogram_dictionary[channel], to_fill, to_weight)

        return histogram_dictionary








