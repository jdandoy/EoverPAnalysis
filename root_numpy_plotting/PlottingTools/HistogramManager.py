import ROOT
class HistogramManager:
    channels = []
    histograms = []
    filename = None

    def __init__(self, fileName):
       self.filename = fileName
       tFile = ROOT.TFile(self.filename, "READ")
       self.channels = [key.GetName() for key in tFile.GetListOfKeys()]
       dir = tFile.Get(self.channels[0])
       self.histograms = [key.GetName() for key in dir.GetListOfKeys()]
       for channel in self.channels:
           self.histograms = [key.replace(channel, "") for key in self.histograms]
       self.histograms = list(set(self.histograms))
       print "channels " + str(self.channels)
       print "histograms " + str(self.histograms)
       tFile.Close()

    def getHistograms(self, histogramName):
       tFile = ROOT.TFile(self.filename, "READ")
       histogram_dict = {}
       if not histogramName in self.histograms:
           raise ValueError("Couldn't find the histogram in the file")

       for channel in self.channels:
           tFile.cd(channel)
           histogram_dict[channel] = tFile.Get(channel + "/" + histogramName + channel)
           histogram_dict[channel].SetDirectory(0)
       tFile.Close()
       return histogram_dict
