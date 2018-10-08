
# a calculation that can be different in MC and in Data. E.G. reweighting the sample.
class calculationDataMC:
    def __init__(self, function, list_of_branches):
        self.function = function
        self.name = function.__name__
        self.needsDataFlag = True
        self.branches = list_of_branches

    def eval(self, data, dataFlag):
        return self.function(data, dataFlag)


class calculation:
    def __init__(self, function, list_of_branches):
        self.function = function
        self.name = function.__name__
        self.needsDataFlag = False
        self.branches = list_of_branches

    def eval(self, data):
        return self.function(data)

