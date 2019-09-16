import ROOT
import array
import uproot as ur
from histogram_manager import HistogramManager
from plotting_tools import DrawDataVsMC

def scotts_rule(histogram):
    N = histogram.Integral()
    sigma = histogram.GetRMS()
    width = (sigma)/(N**(1.0/3.0))
    return width

def generate_eop_var(low, high, sub_ranges = {}):
    eop = ROOT.RooRealVar("eop", "eop" ,low,high)
    eop.setRange("Full", low,high)
    for r in sub_ranges:
        eop.setRange(r, sub_ranges[r][0], sub_ranges[r][1])
    return eop

def prepare_for_fit():
    #minimize the output from roofit
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Fitting)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.InputArguments)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Eval)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.RooAbsReal.defaultIntegratorConfig().method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")  ## Better numerical integrator

def generate_dcb(x):
    eop_gaus_mean = ROOT.RooRealVar("eop_gaus_mean", "eop_gaus_mean", 1.0, 0.0, 2.0)
    eop_gaus_sigma = ROOT.RooRealVar("eop_gaus_sigma", "eop_gaus_sigma", 0.1, 0.0, 100.0)
    eop_gaus_alphaLo = ROOT.RooRealVar("eop_gaus_alphaLo", "eop_gaus_alphaLo", 1.0, 0.0, 100.0)
    eop_gaus_alphaHi = ROOT.RooRealVar("eop_gaus_alphaHi", "eop_gaus_alphaHi", 1.0, 0.0, 100.0)
    eop_gaus_nLo = ROOT.RooRealVar("eop_gaus_nLo", "eop_gaus_nLo", 10.0, 0.0, 100.0)
    eop_gaus_nHi = ROOT.RooRealVar("eop_gaus_nHi", "eop_gaus_nHi", 10.0, 0.0, 100.0)
    eop_gaus_model = ROOT.RooTwoSidedCBShape("signal_model", "signal_model", x, eop_gaus_mean, eop_gaus_sigma, eop_gaus_alphaLo, eop_gaus_nLo, eop_gaus_alphaHi, eop_gaus_nHi)
    var_list = [eop_gaus_mean, eop_gaus_sigma, eop_gaus_alphaLo, eop_gaus_alphaHi, eop_gaus_nLo, eop_gaus_nHi]
    return eop_gaus_model, var_list

keep_alive = []

def generate_landau_gaus(x):
    gaus, gaus_vars = generate_gaus(x)
    landau, landau_vars = generate_landau(x)
    lxg = ROOT.RooFFTConvPdf("lxg","landau (X) gauss",x,landau,gaus)
    keep_alive.append(gaus)
    keep_alive.append(landau)
    for var in gaus_vars:
        if "mean" in var.GetName():
            var.setVal(0.0)
            var.setConstant(True)
    return lxg, gaus_vars + landau_vars

def generate_gaus(x):
    eop_gaus_mean = ROOT.RooRealVar("eop_gaus_mean", "eop_gaus_mean", 0.0, 0.0, 2.0)
    eop_gaus_sigma = ROOT.RooRealVar("eop_gaus_sigma", "eop_gaus_sigma", 0.1, 0.0, 1.0)
    eop_gaus_model = ROOT.RooGaussian("gaus","gaus(x,mean,sigma)",x,eop_gaus_mean,eop_gaus_sigma)
    var_list = [eop_gaus_mean, eop_gaus_sigma]
    return eop_gaus_model, var_list

def generate_landau(x):
    eop_landau_mpv = ROOT.RooRealVar("eop_landau_mean", "eop_landau_mean", 0.0, 0.0, 1.3)
    eop_landau_sigma = ROOT.RooRealVar("eop_landau_sigma", "eop_landau_sigma", 0.1, 0.0, 1.0)
    eop_landau_model = ROOT.RooLandau('lx', 'lx', x, eop_landau_mpv, eop_landau_sigma)
    var_list = [eop_landau_mpv, eop_landau_sigma]
    return eop_landau_model, var_list

def do_fit(histograms, function="gaus"):
    mpvs = {}
    mpv_errs = {}
    fit_results = {}
    for channel in histograms:
        if channel != "LowMuData" and channel != "PythiaJetJet":
            continue
        print("Fitting the histogram in channel {}".format(channel))
        to_fit = histograms[channel]
        eop = generate_eop_var(-1.0, 3.0)
        eop_hist = ROOT.RooDataHist("eop_var", "eop_var", ROOT.RooArgList(eop), to_fit)
        mpv = -999
        mpv_entries = -999
        mpv_bin = -999
        sigma_up = -999
        sigma_down = -999

        #find the 10 bins with the highest average mpv
        bins = [i for i in range(1, 10)]
        for bin in range(5,to_fit.GetNbinsX()+1):
            print("for bin {} looking at bins {}".format(bin, bins))
            points = [to_fit.GetBinContent(b) for b in bins]
            average = sum(points)/float(len(points))
            if average > mpv_entries:
                mpv_bin = bins[4]
                mpv_entries = average
                mpv = to_fit.GetBinCenter(mpv_bin)
            bins = bins[1:] + [bin + 5]

        #find the integral for eop's below the mpv
        integral_left = to_fit.Integral(1, mpv_bin)

        #find the 68% quantile below the mpv
        lower_count = 0
        for lower_bin in range(mpv_bin-1, 0, -1):
            lower_count += to_fit.GetBinContent(lower_bin)
            if lower_count > (1.0 - 0.20) * integral_left:
                break
        lower_bin = lower_bin
        lower_content = to_fit.GetBinContent(lower_bin)

        #find the bin above the mpv that has the same number of entries as the lower bin at the 68% quantile
        up_bin = 0
        for upper_bin in range(mpv_bin +1, to_fit.GetNbinsX()+1):
             if to_fit.GetBinContent(upper_bin) <=  lower_content:
                 break
             up_bin = upper_bin
        upper_bin = up_bin

#       print("Integral {}".format(integral))
#       bins = [i for i in range(1, 10)]
#       for bin in range(5, to_fit.GetNbinsX()+1):
#           #skip those bins around the peak
#           print("for bin {} looking at bins {}".format(bin, bins))
#           points = [to_fit.GetBinContent(b) for b in bins]
#           average = sum(points)/float(len(points))
#           lower_count += to_fit.GetBinContent(bin)
#           if lower_count > 0.16 * integral and sigma_down < -100:
#               sigma_down = to_fit.GetBinCenter(bin)
#           #if abs(average - mpv_entries)/mpv_entries < 0.01:
#           #    bins = bins[1:] + [bin + 5]
#           #    continue
#           if sigma_down > -100:
#               upper_count += to_fit.GetBinContent(bin)
#           if upper_count > (1.0 - 0.16) * integral and sigma_up < -100:
#               sigma_up = to_fit.GetBinCenter(bin)
#           bins = bins[1:] + [bin + 5]

        #do the fit in this range.
        sigma_up = to_fit.GetBinCenter(upper_bin)
        sigma_down = to_fit.GetBinCenter(lower_bin)
        eop.setRange("Fit", sigma_down, sigma_up)
        print(integral_left)
        print("MPV: {}".format(mpv))
        print("Fitting in range [{},{}]".format(sigma_down, sigma_up))
        if function == "gaus":
            model, variables = generate_gaus(eop)
        elif function == "landau":
            model, variables = generate_landau(eop)
        elif function == "landauxgaus" or function == "gausxlandau":
            model_variables = generate_landau_gaus(eop)

        for var in variables:
            if "mean" in var.GetName():
                var.setVal(mpv)
        print(model)
        prepare_for_fit()
        result = model.fitTo(eop_hist, ROOT.RooFit.Range("Fit"))
        f=eop.frame()
        eop_hist.plotOn(f)
        model.plotOn(f)
        f.Draw()

        for var in variables:
            if "mean" in var.GetName():
                mpvs[channel]=var.getVal()
                mpv_errs[channel]=var.getError()
        fit_results[channel]=result

    return mpvs, mpv_errs, fit_results

def test_fit(f, histogram_base = "EOPDistribution", selection_name = "MIPSelectionHadFracAbove70", function = "gaus"):

    #get the binning vectors
    rf = ROOT.TFile(f, "READ")

    tree = rf.Get(selection_name + "BinningTree")
    for bins in tree:
        break

    eta_bins_low = getattr(bins, selection_name+"EtaBinsLow")
    eta_bins_high = getattr(bins, selection_name+"EtaBinsHigh")

    p_bins_low_for_eta_bin = []
    p_bins_high_for_eta_bin = []

    #get all of the binning information that we need
    for i in range(0, eta_bins_low.size()):
        p_bins_low_for_eta_bin.append(getattr(bins, selection_name+"PBinsLow_Eta"+str(i)))
        p_bins_high_for_eta_bin.append(getattr(bins, selection_name+"PBinsHigh_Eta"+str(i)))

    HM = HistogramManager(f)
    histograms_in_eta_bin = []
    for i, eta_low, eta_high in zip(range(0, len(eta_bins_low)),eta_bins_low, eta_bins_high):
        mpvs = {}
        mpv_errs = {}
        for j, p_low, p_high in zip(range(0, len(p_bins_low_for_eta_bin[i])),p_bins_low_for_eta_bin[i], p_bins_high_for_eta_bin[i]):
            histogram_name = histogram_base + "_" + selection_name + "_Eta_" + str(i) + "_Momentum_" + str(j)
            histograms = HM.getHistograms(histogram_name)
            mpv, mpv_err, fit_result = do_fit(histograms, function = function)
            for channel in mpv:
               if channel not in mpvs:
                   mpvs[channel]=[]
                   mpv_errs[channel]=[]
               mpvs[channel].append(mpv[channel])
               mpv_errs[channel].append(mpv_err[channel])
        hist_name = histogram_base + "{}_{}_".format(function, "fit") + selection_name + "_Eta_" + str(i)

        bins = [b for b in p_bins_low_for_eta_bin[i]] + [p_bins_high_for_eta_bin[i][-1]]
        bin_array = array.array('d', bins)

        histograms = {}
        for channel in mpvs:
            histograms[channel] = ROOT.TH1D(hist_name + channel, hist_name + channel, len(bins)-1, bin_array)
            for i in range(0, len(mpvs[channel])):
                histograms[channel].SetBinContent(i+1, mpvs[channel][i])
                histograms[channel].SetBinError(i+1, mpv_errs[channel][i])

        MCKeys = ["PythiaJetJet"]
        DataKey="LowMuData"
        base_description = ["P_{T} Reweighted"]
        description = base_description + ["MIP Selection", "{:.1f} < |#eta| < {:.1f}".format(eta_low, eta_high)]
        channelLabels = {"SinglePion": "Single Pion", "PythiaJetJet" : "#splitline{Pythia8}{MinBias and Dijet}", DataKey: "2017 Low-<#mu> Data", "PythiaJetJetPions    Reweighted":"Pythia8 MB+DJ Pions Only", "PythiaJetJetHardScatter":"Pythia8 MB+DJ Truth Matched", "PythiaJetJetTightIso": "#splitline{Pythia8}{MinBias and D    ijet}", "LowMuDataTightIso":"2017 Low-<#mu> Data"}
        ratio_min = 0.9
        ratio_max = 1.1
        #draw the fit results as a set of histograms
        canvas = DataVsMC1 = DrawDataVsMC(histograms,\
                       channelLabels,\
                       MCKeys = MCKeys,\
                       ratio_min = ratio_min,\
                       ratio_max = ratio_max,\
                       doLogy=False,\
                       doLogx=True,\
                       ylabel = "MPV(E/P)",\
                       DataKey=DataKey,\
                       extra_description = description)[0]

        canvas_name = hist_name + "_"  + selection_name + "_Eta_" + str(i) + ".png"
        canvas.Print(canvas_name)
        canvas.Close()

if __name__ == "__main__":
    f="pt_reweighted.root"
    for selection_name in ["MIPSelectionHadFracAbove70", "20TRTHitsNonZeroEnergy"]:
        test_fit(f, selection_name = selection_name)


def fitHistograms(histograms, fit_function, histogramName, channels=[], eta_low=-1,eta_high=-1,p_low=-1,p_high=-1, refit=False, rebin=False, rebin_rule = None):
    '''
    Fit function fit_function to the histograms
    fit_function can be gaus, landau or convolution
    '''

    low_value = 0.0
    max_x = {}
    max_x_err = {}
    chisq = {}

    low_fit = 1000.0
    high_fit = -1000.0
    low_rms = 0.0
    high_rms = 1.2
    means = {}
    sigmas = {}

    if rebin:
        histograms = rebin_histograms(histograms, binning_option = rebin_rule)

    if not refit:
       for channel in channels:
           histogram = histograms[channel]
           low_bin = histogram.FindBin(low_rms)
           high_bin = histogram.FindBin(high_rms)
           histogram.GetXaxis().SetRange(low_bin, high_bin)
           mean = histogram.GetMean()
           rms = histogram.GetRMS()
           sigmas[channel]=rms

           new_low_fit = mean - rms * 1.0
           new_high_fit = mean + rms * 1.0
           means[channel]=mean
           if new_low_fit < low_fit:
               low_fit = new_low_fit
           if new_high_fit > high_fit:
               high_fit = new_high_fit
    if refit:
       histogram = histograms["LowMuData"]
       histogram_name = histogram.GetName()
       fit_function_string = fit_function + histogram_name
       fit = histogram.GetFunction(fit_function_string)
       low_fit = fit.GetParameter(0) - fit.GetParameter(1) * 0.8
       high_fit = fit.GetParameter(0) + fit.GetParameter(1) * 0.8
       for channel in channels:
           means[channel] = fit.GetParameter(0)
           sigmas[channel] = fit.GetParameter(1)

    for channel in channels:
       histogram = histograms[channel]

       histogram_name = histogram.GetName()
       #this is the landau distribution that will be fit to the histogramsograms
       landau = ROOT.TF1("landau_" + channel +histogram_name, "[2]*TMath::Landau(x, [0], [1])", -1.0, 5.0)
       landau.SetParName(0, "mpv")
       landau.SetParameter(0, means[channel])
       landau.SetParLimits(0, 0.3, 1.1)
       landau.SetParName(1, "sigma")
       landau.SetName("landau" +  histogram_name)
       landau.SetParameter(1, sigmas[channel]/4.0)
       landau.SetParLimits(1, sigmas[channel]/100.0, sigmas[channel]*2.0)
       landau.SetParName(2, "Norm")
       landau.SetParameter(2, histogram.Integral())

       #this is the gaus distribution that will be fit to the histogramsograms
       gaus = ROOT.TF1("gaus_" + channel +histogram_name, "[2]*TMath::Gaus(x, [0], [1])", -1.0, 5.0)
       gaus.SetParName(0, "mu")
       gaus.SetParameter(0, means[channel])
       gaus.SetParLimits(0, 0.45, 0.95)
       gaus.SetParName(1, "sigma")
       gaus.SetName("gaus" + histogram_name)
       gaus.SetParameter(1, sigmas[channel])
       gaus.SetParLimits(1, sigmas[channel]/3.0, 1.2*sigmas[channel])
       gaus.SetParName(2, "Norm")
       gaus.SetParameter(2, histogram.Integral())

       #Create a gaussian convoluted with a landau histogramsogram
       gaus_forconvolution = ROOT.TF1("gaus_forconvolution_" + channel +histogram_name, "TMath::Gaus(x, 0.0, [0])", -10.0, +10.0)

       landau_forconvolution = ROOT.TF1("landau_forconvolution_" + channel +histogram_name, "[2]*TMath::Landau(x, [0], [1])", -1.0, 5.0)

       convolution = ROOT.TF1Convolution(gaus_forconvolution, landau_forconvolution,-1,6,True)
       convolution.SetRange(-1.,5.)
       convolution.SetNofPointsFFT(10000)
       convolution_tofit = ROOT.TF1("f",convolution, -1.0, 5., convolution.GetNpar())
       convolution_tofit.SetName("convolution" +histogram_name)

       convolution_tofit.SetParName(0, "SigmaSmear")
       convolution_tofit.SetParLimits(0, 0.0, 100.0)
       convolution_tofit.SetParameter(0, 0.5)

       convolution_tofit.SetParName(1, "mpv")
       convolution_tofit.SetParameter(1, means[channel])
       convolution_tofit.SetParLimits(1, 0.3, 1.1)
       convolution_tofit.SetParName(2, "sigma")
       convolution_tofit.SetParameter(2, sigmas[channel]/4.0)
       convolution_tofit.SetParLimits(2, sigmas[channel]/100.0, sigmas[channel]*2.0)
       convolution_tofit.SetParName(3, "Norm")
       convolution_tofit.SetParameter(3, histogram.Integral())
       print("Created convolution function")
       convolution_tofit.Print()

       #Good settings for a landau x gaus
       fit_function_string = fit_function + histogram_name
       print("convolution" +histogram_name)
       print(fit_function_string)
       histogram.GetXaxis().SetRange(histogram.FindBin(0.15), histogram.FindBin(1.3))
       histogram.GetXaxis().SetRange()
       histogram.Fit(fit_function_string, "", "", low_fit, high_fit)

       fit = histogram.GetFunction(fit_function_string)
       if  fit:
           fit.SetLineColor(histogram.GetLineColor())
           if fit_function == "convolution":
               max_x_tmp=0.0
               max_val_tmp=0.0
               max_x_tmp_err = 0.0
               chisq[channel]=fit.GetChisquare()
           elif fit_function == "gaus":
               max_x_tmp = fit.GetParameter(0)
               max_val_tmp = -1.0
               max_x_tmp_err = fit.GetParError(0)
               chisq[channel]=fit.GetChisquare()
       else:
           max_x_tmp = -1.0
           max_val_tmp = -1.0
           max_x_tmp_err = 0.0
           chisq[channel] = -1.0

       max_x[channel]=max_x_tmp
       max_x_err[channel]=max_x_tmp_err

       print("The maximum value in data was at eop " + str(max_x_tmp))

    return max_x, max_x_err, chisq

def FindMostProbableValue(fit_function, low, high, ndivisions=10000):
    x = np.linspace(low, high, ndivisions)
    vals = np.zeros(len(x))
    max_val = 0
    max_x = -999999.0
    for i in range(0, len(x)):
        vals[i] = fit_function.Eval(x[i])
        if vals[i] > max_val:
            max_val = vals[i]
            max_x = x[i]
    return max_x, max_val


