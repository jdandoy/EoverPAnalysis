from multiprocessing import Pool, cpu_count
from ROOT import TFile, TH1D, TTree
from numpy import loadtxt

#These are python JZW samples. I normalize to the number of generated events, the cross section and the filter efficiency
weight_dictionary = {\
"361020" : (1/(1000000.0)) * 78420000 * 0.9755,
"361021" : (1/(1000000.0)) * 78420000 * 0.00067143,
"361022" : (1/(999000.0)) *  2433200 * 0.00033423,
}

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

def GetData(partition = (0, 0), bare_branches = [], filename = None, treename = None, variables = [], weightCalculator = None, selections = [], selection_string = "",  verbose = False):
    '''a function used for multiprocessing. It makes all of the selections used by the analysis'''
    import root_numpy as rnp

    isData = getIsData(filename)
    for branch in weightCalculator.branches:
        if branch not in bare_branches:
            bare_branches.append(branch)
    branches = branchDresser(bare_branches)
    if verbose: print branches

    if verbose: print "Reading from file " + filename
    f = TFile(filename, "READ")
    t = f.Get(treename)
    data = rnp.tree2array(t, branches, selection_string, start = partition[0], stop = partition[1])
    if verbose: print "Got the data for parition " + str(partition)
    f.Close()
    del f
    del t

    # a dictionary of selections
    selection_dict = {}
    # a dictionary of variables = {}
    variable_dict = {}

    weights = weightCalculator.eval(data, isData)
    if not isData:
        xsec_weight = getXSectionWeight(filename)
        if verbose: print("X Section Weight Set To " + str(xsec_weight))
        weights = weights * xsec_weight
    variable_dict["weights"] = weights

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
    return return_dict
