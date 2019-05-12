import numpy as np
from array import array
from calculation.calculation import calculation
import ROOT
import imp
import time
import os
import glob
import imp

def get_log_bins(minBin, maxBin, nBins):
    '''Get nBins logarithmically-evenly spaced bins ranging from minBin to maxBin'''
    bins = []
    base = (float(maxBin)/float(minBin)) ** (1.0/float(nBins))
    for i in range(0, nBins+1):
        bins.append(minBin * (base) ** i)
    return bins

def get_bins(minBin, maxBin, nBins):
    '''Get nBins evenly spaced bins ranging from minBin to maxBin'''
    bins = []
    step = float(maxBin - minBin)/float(nBins)
    for i in range(0, nBins+1):
        bins.append(minBin + (i*step))
    return bins

def get_p(Pt, eta):
    '''Given a value of pT and eta, return the momentum'''
    return Pt*np.cosh(eta)

def create_selection_function(template, branches, *args):
    if len(args) > 3:
        raise ValueError("Only up to three variables supported for the template function")
    if len(args) == 0:
        raise ValueError("You need to pass at least one argument")

    if len(args) == 1:
        function = lambda x, y=args[0]: template(x,y)
    if len(args) == 2:
        function = lambda x, y=args[0], z=args[1]: template(x,y,z)
    if len(args) == 3:
        function = lambda x, y=args[0], z=args[1], w=args[2]: template(x,y,z,w)
                
    function.__name__ = template.__name__ + "_".join([str(arg) for arg in args])
    selection_function = calculation(function, branches)
    return selection_function

class HistogramFiller:
    '''
    '''
    def __init__(self, inputs, treeName, weightCalculator, base_selections = "", partition_dictionary = None):
        self.channelFiles = {}
        self.channelLabels = {}
        self.binningHistograms = {}
        self.AddInputDictionary(inputs)
        self.treeName = treeName
        self.partition_dictionary = partition_dictionary
        self.normalization_dictionary = {}
        self.test = False
        self.verbose = False
        self.total_selections = []
        self.total_variables =[]
        self.HistogramCallDictionary = {}
        self.base_selections = base_selections
        self.NormalizationWeightsDictionary = {} # A dictionary of channel to variable to the weights needed for the reweighting of this variable
        self.object_counter = 0
        self.weightCalculator = weightCalculator

    def BookHistogramForBinning(self, histogram, histogramName):
        '''This saves a histogram that can be used later to determine bin sizes. This could be a track multiplicity distribution, for example'''
        histogram.SetDirectory(0)
        self.binningHistograms[histogramName] = histogram
        toGlobalScope(histogram) #keep the histogram alive at the global scope

    def AddInputDictionary(self, dictionary):
        '''Take an input dictionary of channel names to: tuple of (channel descriptor for legend, list of [input root filename strings])'''
        self.channels = dictionary.keys()
        for channel in self.channels:
            self.channelFiles[channel] = dictionary[channel][1]
            self.channelLabels[channel] = dictionary[channel][0]

    def GetVariablesAndWeights(self, channel, filename, variables, list_selections):
        '''
        Given a string channel, string filename, a list of calculation variables and a list of calculations list_selections, return a dictionary keys selection_dict, variable_dict and weights. selection_dict is a dictionary of key selection name to numpy array of bool. variable_dict is a dictionary of string variable name to numpy array variable. weights is a numpy array of floats
        '''
        branches = GetListOfNeededBranches(variables, list_selections)

        #get the parition of the ttree to be read
        partition = None
        if self.partition_dictionary == None:
            partition = (0, total_entries)
        else:
            partition = self.partition_dictionary[filename]
            if self.verbose: print("Found a partition")

        tree = self.tree_dict[filename]

        print("Getting data for partition " + str(partition))
        result = GetData(partition = partition, bare_branches = branches, channel = channel, filename = filename, tree = tree, treename = self.treeName, variables=variables, weightCalculator = self.weightCalculator, selections = list_selections, selection_string = self.base_selections, verbose = self.verbose)

        #Get the selections, variables and weights
        selection_dict = result["selection_dict"]
        variable_dict = result["variable_dict"]
        weights = result["weights"]

        if self.verbose: print "The following selections have been evaluated "
        for selection in selection_dict:
            if self.verbose: print selection
        if self.verbose: print "The following variables have be evaluated "
        for variable in variable_dict:
            if self.verbose: print variable

        return variable_dict, selection_dict, weights


    def GetHistograms(self, histogram_name, data_dictionary, variable, list_selections = [], bins = 1, range_low = 0.000001, range_high=1. - 0.00001,  xlabel ="", ylabel = "", HistogramPerFile=False, useWeights = True):
        '''
        Get the histogram for variable after list_selections is applied.
        '''

        variableNameToFill = variable.name
        variables = [variable]
        histogram_dictionary = {}
        if not HistogramPerFile:
            for channel in self.channels:
                if (type(bins) == list):
                    bins_array = array('d',bins)
                    histogram_dictionary[channel] = ROOT.TH1D(histogram_name + channel, histogram_name + channel, len(bins_array)-1, bins_array)
                else:
                    histogram_dictionary[channel] = ROOT.TH1D(histogram_name + channel, histogram_name + channel, bins, range_low + 0.0000001, range_high - 0.000001)
                histogram_dictionary[channel].GetXaxis().SetTitle(xlabel)
                histogram_dictionary[channel].GetYaxis().SetTitle(ylabel)
                histogram_dictionary[channel].Sumw2()

            for channel in self.channels:
                normalization_weight = 0.0
                for filename in self.channelFiles[channel]:
                    variable_dict, selection_dict, weights = data_dictionary[channel][filename]
                    total_selection = np.ones(len(weights)) > 0.0
                    for selection in list_selections:
                        total_selection &= selection_dict[selection.name]
                    to_fill = variable_dict[variableNameToFill][total_selection]
                    to_weight = weights[total_selection]
                    if self.verbose: print(len(to_fill))
                    if self.verbose: print(len(to_weight))
                    if self.verbose: print to_fill
                    if self.verbose: print to_weight
                    if self.verbose: print("Filling Variable " + variable.name)
                    if useWeights:
                        fill_hist(histogram_dictionary[channel], to_fill, to_weight)
                    else:
                        fill_hist(histogram_dictionary[channel], to_fill)

            return histogram_dictionary

    def Get2DHistograms(self, histogram_name, data_dictionary, variable_x, variable_y, list_selections = [], bins_x = 1, range_low_x = 0.000001, range_high_x=1. - 0.00001,  xlabel ="", bins_y=1, range_low_y=0.000001, range_high_y=1. - 0.00001, ylabel = "", zlabel="",):
        '''the 2-d histgram with variable_x and variable_y drawn'''
        variableNameToFill_x = variable_x.name
        variableNameToFill_y = variable_y.name
        variables = [variable_x, variable_y]
        histogram_dictionary = {}
        for channel in self.channels:
            if (type(bins_x) == list and type(bins_y) == list):
                bins_array_x = array('d',bins_x)
                bins_array_y = array('d',bins_y)
                histogram_dictionary[channel] = ROOT.TH2D(histogram_name + channel, histogram_name + channel, len(bins_array_x)-1, bins_array_x, len(bins_array_y)-1, bins_array_y)
            elif (type(bins_x) != list and type(bins_y) != list):
                histogram_dictionary[channel] = ROOT.TH2D(histogram_name + channel, histogram_name + channel, bins_x, range_low_x + 0.0000001, range_high_x - 0.000001, bins_y, range_low_y+0.0000001, range_high_y + 0.0000001)
            else:
                raise ValueError("both of the bins_x and bins_y variables need to be the same type. Both integers, or both lists")
            histogram_dictionary[channel].GetXaxis().SetTitle(xlabel)
            histogram_dictionary[channel].GetYaxis().SetTitle(ylabel)
            histogram_dictionary[channel].GetZaxis().SetTitle(zlabel)
            histogram_dictionary[channel].GetZaxis().SetTitleSize(0.035)
            histogram_dictionary[channel].GetZaxis().SetTitleOffset(1.35)
            histogram_dictionary[channel].Sumw2()

        for channel in self.channels:
            normalization_weight = 0.0
            for filename in self.channelFiles[channel]:
                variable_dict, selection_dict, weights = data_dictionary[channel][filename]
                total_selection = np.ones(len(weights)) > 0.0
                for selection in list_selections:
                    total_selection &= selection_dict[selection.name]
                to_weight = weights[total_selection]
                n_sel = len(to_weight)
                to_fill = np.zeros((n_sel,2))
                to_fill[:,0] = variable_dict[variableNameToFill_x][total_selection]
                to_fill[:,1] = variable_dict[variableNameToFill_y][total_selection]
                if self.verbose: print to_fill
                if self.verbose: print to_weight
                if self.verbose: print("Filling Variable " + variable.name)
                fill_hist(histogram_dictionary[channel], to_fill, to_weight)
        return histogram_dictionary


    def GetTProfileHistograms(self, histogram_name, data_dictionary, variable_x, variable_y, list_selections = [], bins = 1, range_low = 0.000001, range_high=1. - 0.00001,  xlabel ="", ylabel="",):
        '''Get a TProfile histogram with variable_y profiled against variable_x, after selections list_selections have been applied'''

        variableNameToFill_x = variable_x.name
        variableNameToFill_y = variable_y.name
        variables = [variable_x, variable_y]
        histogram_dictionary = {}
        for channel in self.channels:
            if (type(bins) == list):
                bins_array = array('d',bins)
                histogram_dictionary[channel] = ROOT.TProfile(histogram_name + channel, histogram_name + channel, len(bins_array)-1, bins_array)
            else:
                histogram_dictionary[channel] = ROOT.TProfile(histogram_name + channel, histogram_name + channel, bins, range_low + 0.0000001, range_high - 0.000001)
            histogram_dictionary[channel].Sumw2()

        for channel in self.channels:
            for filename in self.channelFiles[channel]:
                variable_dict, selection_dict, weights = data_dictionary[channel][filename]
                total_selection = np.ones(len(weights)) > 0.0
                for selection in list_selections:
                    total_selection &= selection_dict[selection.name]
                to_weight = weights[total_selection]
                n_sel = len(to_weight)
                to_fill = np.zeros((n_sel,2))
                to_fill[:,0] = variable_dict[variableNameToFill_x][total_selection]
                to_fill[:,1] = variable_dict[variableNameToFill_y][total_selection]
                if self.verbose: print to_fill
                if self.verbose: print to_weight
                if self.verbose: print("Filling Variable " + variable.name)
                if self.verbose: print("Filling Histogram")
                fill_profile(histogram_dictionary[channel], to_fill, to_weight)
                if self.verbose: print("Finished filling histogram")

            histogram_dictionary[channel].GetXaxis().SetTitle(xlabel)
            histogram_dictionary[channel].GetYaxis().SetTitle(ylabel)
        return histogram_dictionary

    def BookHistograms(self, histogram_name, variable, list_selections = [], bins = 1, range_low = 0.000001, range_high=1. - 0.00001,  xlabel ="", ylabel = "", HistogramPerFile=False, useWeights = True):
        if histogram_name not in self.HistogramCallDictionary:
            self.HistogramCallDictionary[histogram_name] = lambda data_dictionary : self.GetHistograms(histogram_name, data_dictionary, variable, list_selections = list_selections, bins = bins, range_low = range_low, range_high=range_high,  xlabel = xlabel, ylabel = ylabel, HistogramPerFile=HistogramPerFile, useWeights = useWeights)
        else:
            raise ValueError("histogram name already exists")
        for selection in list_selections:
            if selection.name not in [sel.name for sel in self.total_selections]:
                self.total_selections.append(selection)

        if variable.name not in [var.name for var in self.total_variables]:
            self.total_variables.append(variable)

    def Book2DHistograms(self, histogram_name, variable_x, variable_y, list_selections = [], bins_x = 1, range_low_x = 0.000001, range_high_x=1. - 0.00001,  xlabel ="", bins_y=1, range_low_y=0.000001, range_high_y=1. - 0.00001, ylabel = "", zlabel=""):
        if histogram_name not in self.HistogramCallDictionary:
            self.HistogramCallDictionary[histogram_name] = lambda data_dictionary : self.Get2DHistograms(histogram_name, data_dictionary, variable_x, variable_y, list_selections = list_selections, bins_x = bins_x, range_low_x =range_low_x, range_high_x=range_high_x,  xlabel =xlabel, bins_y=bins_y, range_low_y=range_low_y, range_high_y=range_high_y, ylabel = ylabel, zlabel=zlabel)
        else:
            raise ValueError("histogram name already exists")

        for selection in list_selections:
            if selection.name not in [sel.name for sel in self.total_selections]:
                self.total_selections.append(selection)

        if variable_x.name not in [var.name for var in self.total_variables]:
            self.total_variables.append(variable_x)

        if variable_y.name not in [var.name for var in self.total_variables]:
            self.total_variables.append(variable_y)

    def BookTProfileHistograms(self, histogram_name,  variable_x, variable_y, list_selections = [], bins = 1, range_low = 0.000001, range_high=1. - 0.00001,  xlabel ="", ylabel = ""):
        if histogram_name not in self.HistogramCallDictionary:
            self.HistogramCallDictionary[histogram_name] = lambda data_dictionary : self.GetTProfileHistograms(histogram_name, data_dictionary, variable_x, variable_y, list_selections = list_selections, bins = bins, range_low = range_low, range_high=range_high,  xlabel =xlabel, ylabel=ylabel)
        else:
            raise ValueError("histogram name already exists")

        for selection in list_selections:
            if selection.name not in [sel.name for sel in self.total_selections]:
                self.total_selections.append(selection)

        if variable_x.name not in [var.name for var in self.total_variables]:
            self.total_variables.append(variable_x)

        if variable_y.name not in [var.name for var in self.total_variables]:
            self.total_variables.append(variable_y)

    def DumpHistograms(self):
        data_dictionary = {}
        for channel in self.channels:
            data_dictionary[channel] = {}
            for filename in self.channelFiles[channel]:
                data_dictionary[channel][filename] = self.GetVariablesAndWeights(channel,filename, self.total_variables, self.total_selections)


        histogram_dictionary = {}
        for histogram_name in self.HistogramCallDictionary:
            histogram_dictionary[histogram_name] = self.HistogramCallDictionary[histogram_name](data_dictionary)

        return histogram_dictionary


#These are python JZW samples. I normalize to the number of generated events, the cross section and the filter efficiency
weight_dictionary = {\
"361020" : (1./(9999000.0)) * 78420000 * 0.9755,
"361021" : (1./(3995000.0)) * 78420000 * 0.00067143,
"361022" : (1./(1998000.0)) * 2433200 * 0.00033423,
"428001" : 1.0,
"428002" : 1.0,
}

try:
    imp.find_module('root_numpy')
    foundRootNumpy=True
    print "Found root_numpy"
except ImportError:
    foundRootNumpy=False
    print "Didn't find root_numpy module. Did NOT set atlas style. Continuing"

if foundRootNumpy:
    from root_numpy import fill_hist, fill_profile
    from root_numpy import tree2array

import os
import psutil
process = psutil.Process(os.getpid())

def branchDresser(branches):
    '''this is a function that dresses the branches to extend the branches'''
    return branches #No extension necessary

def getXSectionWeight(filename):
    weight = None
    for dsid in weight_dictionary:
        if dsid in filename:
            weight = weight_dictionary[dsid]

    if weight == None:
        raise ValueError("Can't find the dataset id in the weight file")
    return weight

def getIsData(filename):
    return "Data" in filename or "data" in filename

def GetData(partition = (0, 0), bare_branches = [], channel = "", filename = None, tree = None, treename = None, variables = [], weightCalculator = None, selections = [], selection_string = "",  verbose = False):
    '''a function used for multiprocessing. It makes all of the selections used by the analysis'''
    global process

    isData = getIsData(filename)
    for branch in weightCalculator.branches:
        if branch not in bare_branches:
            bare_branches.append(branch)
    branches = branchDresser(bare_branches)
    if verbose: print branches

    if verbose: print "Reading from file " + filename

    data = None
    for i in range(1, 50):
        try:
            data = tree2array(tree, branches, selection_string, start = partition[0], stop = partition[1])
        except Exception as e:
            print "Catching a failed attempt to retrieve data error. Trying agagin in 5 seconds"
            print e
            time.sleep(5) #try again in 5 seconds
        else:
            break

    if data == None:
        raise ValueError("Could not retrieve the data.")

    if verbose: print "Got the data for parition " + str(partition)

    # a dictionary of selections
    selection_dict = {}
    # a dictionary of variables = {}
    variable_dict = {}

    print("Evaluating weights")
    weights = weightCalculator.eval(data, isData, channel)

    if not isData:
        print("getting the xsection weight")
        xsec_weight = getXSectionWeight(filename)
        if verbose: print("X Section Weight Set To " + str(xsec_weight))
        weights = weights * xsec_weight

    ##calculate everything we need in one go!
    for variable in variables:
        if verbose: print("calculating variables for " + variable.name)
        variable_dict[variable.name] = variable.eval(data)

    #selection_dict is a dictionary of numpy arrays that have dimension # of events
    #each entry in the numpy array tells you if the event passed the selection
    for selection in selections:
        if verbose: print("calculating selection " + selection.name)
        if not selection.name in selection_dict:
            selection_dict[selection.name] = selection.eval(data)

    return_dict = {}
    return_dict["selection_dict"] = selection_dict
    return_dict["variable_dict"] = variable_dict
    return_dict["weights"] = weights
    f = tree.GetCurrentFile()
    tree.SetDirectory(0)
    f.Close()
    
    return return_dict

def GetListOfNeededBranches(variables, selections):
    '''given a list of variables and selections, get all of the branches that should be read from the tree'''
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
