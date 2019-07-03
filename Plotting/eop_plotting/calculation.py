from pyximport import install
import numpy as np
install(setup_args={"include_dirs":np.get_include()},reload_support=True)
from utils_cython import get_weights_from_bins

class CalculationDataMC:
    def __init__(self, function, list_of_branches):
        self.function = function
        self.name = function.__name__
        self.needsDataFlag = True
        self.branches = list_of_branches

    def eval(self, data, dataFlag):
        return self.function(data, dataFlag)

class Calculation:
    def __init__(self, function, list_of_branches):
        self.function = function
        self.name = function.__name__
        self.needsDataFlag = False
        self.branches = list_of_branches

    def eval(self, data):
        return self.function(data)

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

    weights = get_weights_from_bins(variable_in_histogram, low_edges, high_edges, hist_weights)

    weights[variable_in_histogram > high_edges[-1]] = 0.0 #set the weights of those events not in the reweighting histogram to 0

    return weights

class WeightCalculation:
    def __init__(self, function, list_of_branches):
        self.function = function
        self.name = function.__name__
        self.branches = list_of_branches
        self.reweightDictionary = {}

    def eval(self, data, isData, channel):
        weights = self.function(data, isData)
        if channel in self.reweightDictionary:
            for variable, histogram, selection in zip(self.reweightDictionary[channel]["variables"], self.reweightDictionary[channel]["histograms"], self.reweightDictionary[channel]["selections"]):
                print("Reweighting variable " + variable.name + " in channel " + channel)
                total_selection = np.ones(len(data)) > 0.5

                for s in selection:
                    total_selection &= s.eval(data)

                extra_weight = WeightsToNormalizeToHistogram(variable.eval(data), histogram)
                weights[total_selection] *= extra_weight[total_selection]
        return weights

    def add_reweight_histogram(self, channel, variable, histogram, selection=[]):
        histogram.SetDirectory(0)
        if channel not in self.reweightDictionary:
            self.reweightDictionary[channel] = {}
            self.reweightDictionary[channel]["variables"] = []
            self.reweightDictionary[channel]["histograms"] = []
            self.reweightDictionary[channel]["selections"] = []

            self.reweightDictionary[channel]["variables"].append(variable)
            self.reweightDictionary[channel]["histograms"].append(histogram)
            self.reweightDictionary[channel]["selections"].append(selection)
            ##make sure that we always read the variable that we need for the histogram reweighting
            for branch_name in variable.branches:
                if branch_name not in self.branches:
                    self.branches.append(branch_name)
            for s in selection:
                for branch_name in s.branches:
                    if branch_name not in self.branches:
                        self.branches.append(branch_name)
        else:
            self.reweightDictionary[channel]["variables"].append(variable)
            self.reweightDictionary[channel]["histograms"].append(histogram)
            self.reweightDictionary[channel]["selections"].append(selection)
            ##make sure that we always read the variable that we need for the histogram reweighting
            for branch_name in variable.branches:
                if branch_name not in self.branches:
                    self.branches.append(branch_name)
            for s in selection:
                for branch_name in s.branches:
                    if branch_name not in self.branches:
                        self.branches.append(branch_name)



