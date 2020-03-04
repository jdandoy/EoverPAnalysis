import ROOT
import array
import uproot as ur
from histogram_manager import HistogramManager
from plotting_tools import DrawDataVsMC, ProjectProfiles, SubtractHistograms
import numpy as np
from plotting_tools import DrawText
ROOT.gROOT.SetBatch(True)

ROOT.gSystem.Load("~/RooFitExtensions/build/libRooFitExtensions.dylib")

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

def generate_gaus(x, extra_str = ""):
    eop_gaus_mean = ROOT.RooRealVar("eop_gaus_mean{}".format(extra_str), "eop_gaus_mean{}".format(extra_str), 0.5, 0.0, 5.0)
    eop_gaus_sigma = ROOT.RooRealVar("eop_gaus_sigma{}".format(extra_str), "eop_gaus_sigma{}".format(extra_str), 0.1, 0.05, 5.0)
    eop_gaus_model = ROOT.RooGaussian("gaus{}".format(extra_str),"gaus(x,mean,sigma){}".format(extra_str),x,eop_gaus_mean,eop_gaus_sigma)
    var_list = [eop_gaus_mean, eop_gaus_sigma]
    return eop_gaus_model, var_list

def generate_landau(x):
    eop_landau_mpv = ROOT.RooRealVar("eop_landau_mean", "eop_landau_mean", 0.0, 0.0, 1.3)
    eop_landau_sigma = ROOT.RooRealVar("eop_landau_sigma", "eop_landau_sigma", 0.1, 0.05, 1.0)
    eop_landau_model = ROOT.RooLandau('landau', 'landau', x, eop_landau_mpv, eop_landau_sigma)
    var_list = [eop_landau_mpv, eop_landau_sigma]
    return eop_landau_model, var_list

def generate_two_gaus(x):
    gaus_one, vars_one = generate_gaus(x, extra_str="one")
    gaus_two, vars_two = generate_gaus(x, extra_str="two")
    for var in vars_one:
        if "mean" in var.GetName():
            var.setVal(0.5)
    for var in vars_two:
        if "mean" in var.GetName():
            var.setVal(0.7)
    coeff = ROOT.RooRealVar("frac", "frac", 0.5, 0.0, 1.0)
    keep_alive.append(coeff)
    keep_alive.append(gaus_one)
    keep_alive.append(gaus_two)
    pdf = ROOT.RooAddPdf("twogauss", "twogauss", gaus_one, gaus_two, coeff)
    return pdf,vars_one + vars_two + [coeff]

def generate_landau_plus_gaus(x):
    gaus, vars_one = generate_gaus(x, extra_str="one")
    landau, vars_two = generate_landau(x)
    for var in vars_one:
        if "mean" in var.GetName():
            var.setVal(0.5)
    for var in vars_two:
        if "mean" in var.GetName():
            var.setVal(0.7)
    coeff = ROOT.RooRealVar("frac", "frac", 0.5, 0.0, 1.0)
    keep_alive.append(coeff)
    keep_alive.append(landau)
    keep_alive.append(gaus)
    pdf = ROOT.RooAddPdf("landau+gauss", "landau+gauss", gaus, landau, coeff)
    return pdf,vars_one + vars_two + [coeff]

def montecarlo_uncertainties(model, variables, x_var, minimum, maximum):
    print("Using random sampling to estimate uncertainties")

    original_values = []
    original_errors = []
    for var in variables:
        original_values.append(var.getVal())
        original_errors.append(var.getError())
    original_values = np.array(original_values)
    original_errors = np.array(original_errors)

    f = model.asTF( ROOT.RooArgList(x_var) )
    original_xmax = f.GetMaximumX(minimum + 0.05, maximum - 0.05)

    montecarlo_maxima = []
    for montecarlo_round in range(0, 10000):
        if montecarlo_round % 1000 == 0:
            print("Bootstrap round {}".format(montecarlo_round))
        random_vals = np.random.normal(size=len(original_values))
        new_values = original_values + random_vals * original_errors #randomly sample the uncertainties on the parameters
        #set the values
        [variables[i].setVal(v) for i,v in enumerate(new_values)]
        #find the maximum
        f = model.asTF( ROOT.RooArgList(x_var) )
        xmax = f.GetMaximumX(minimum + 0.05, maximum - 0.05)
        montecarlo_maxima.append(xmax)

    montecarlo_maxima = np.array(montecarlo_maxima)
    print("The mean was {}".format(np.mean(montecarlo_maxima)))
    print("The original max was {}".format(original_xmax))
    print("The std dev was {}".format(np.std(montecarlo_maxima)))

    #reset the variables to their origin values
    [variables[i].setVal(v) for i,v in enumerate(original_values)]
    return original_xmax, np.std(montecarlo_maxima)


def do_fit(histograms, function="gaus", extra_str = "", montecarlo_errors = False, p_bin = 0, eta_bin=0,sel_type = ""):
    mpvs = {}
    mpv_errs = {}
    fit_results = {}
    for channel in histograms:
        if channel != "LowMuData" and channel != "PythiaJetJet" and channel != "LowMuDataTightIso" and channel != "PythiaJetJetTightIso":
            continue
        print("Fitting the histogram in channel {}".format(channel))
        to_fit = histograms[channel]
        if "PythiaJetJet" in channel:
            to_fit.Rebin(2)
        low = 0.0
        high = 2.0
        eop = generate_eop_var(0.0, 2.0)
        eop_hist = ROOT.RooDataHist("eop_var", "eop_var", ROOT.RooArgList(eop), to_fit)
        mpv = -999
        mpv_entries = -999
        mpv_bin = -999
        sigma_up = -999
        sigma_down = -999

        #find the 10 bins with the highest average mpv
        bins = [i for i in range(1, 10)]
        for bin in range(5,to_fit.GetNbinsX()-5):
            print("for bin {} looking at bins {}".format(bin, bins))
            points = [to_fit.GetBinContent(b) for b in bins]
            average = sum(points)/float(len(points))
            if average > mpv_entries:
                mpv_bin = bins[4]
                mpv_entries = average
                mpv = to_fit.GetBinCenter(mpv_bin)
            bins = bins[1:] + [bin + 5]

        low_limit = to_fit.FindBin(0.1)
        #if to_fit.Integral() > 2000:
        #find the integral for eop's below the mpv
        integral_left = to_fit.Integral(low_limit, mpv_bin)
        print("integral left {}".format(integral_left))

        #find the 68% quantile below the mpv
        lower_count = 0
        lower_bin = 0
        for lower_bin in range(mpv_bin-1, low_limit, -1):
            lower_count += to_fit.GetBinContent(lower_bin)
            if lower_count > (1.0 - 0.2) * integral_left:
                break
        lower_content = to_fit.GetBinContent(lower_bin)

        #find the bin above the mpv that has the same number of entries as the lower bin at the 68% quantile
        up_bin = 0
        max_bin = mpv_bin +1
        for upper_bin in range(mpv_bin +1, to_fit.GetNbinsX()+1):
             if to_fit.GetBinCenter(upper_bin) > high:
                 continue
             max_bin += 1
             print("{} {}".format(upper_bin, max_bin))
             if to_fit.GetBinContent(upper_bin) <=  lower_content:
                 break
             up_bin = upper_bin
             upper_bin -= 5
        upper_bin = up_bin
        if upper_bin == 0:
            upper_bin = max_bin

        if function == "gaus":
            old_min = to_fit.GetMinimum()
            old_max = to_fit.GetMaximum()
            to_fit.SetMinimum(0.15)
            to_fit.SetMaximum(1.1)
            mean = to_fit.GetMean()
            sig = to_fit.GetRMS()
            to_fit.SetMinimum(old_min)
            to_fit.SetMaximum(old_max)
            upper_bin = to_fit.FindBin(mean + sig)
            lower_bin = to_fit.FindBin(mean-sig)

            if p_bin > 7:
                upper_bin = to_fit.FindBin(mean + 1.2*sig)
                lower_bin = to_fit.FindBin(mean - 1.2*sig)
            #if p_bin > 6:
            #    to_fit.Rebin(2)
            #if p_bin > 10:
            #to_fit.Rebin(4)

            #hack in the fits:
            if "MIP" in sel_type:
               if eta_bin == 0:
                   if p_bin == 1:
                      if "Pythia" not in channel:
                          upper_bin = to_fit.FindBin(mean + (-0.225) *sig)
                          lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                      else:
                          upper_bin = to_fit.FindBin(mean + (-0.175) *sig)
                          lower_bin = to_fit.FindBin(mean - 0.95 *sig)

                   if p_bin == 2:
                      upper_bin = to_fit.FindBin(mean - 0.125 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.85 *sig)
                   if p_bin == 3:
                      upper_bin = to_fit.FindBin(mean + (-0.1) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                   if p_bin == 4:
                      upper_bin = to_fit.FindBin(mean + (-0.05) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)
                   if p_bin == 5:
                      upper_bin = to_fit.FindBin(mean + 0.1 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.1 *sig)
                   if p_bin == 6:
                      upper_bin = to_fit.FindBin(mean + 0.45 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.15 *sig)
                   if p_bin == 7:
                      upper_bin = to_fit.FindBin(mean + 0.475 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.1 *sig)

               if eta_bin == 1:
                   if p_bin == 1:
                      upper_bin = to_fit.FindBin(mean + (-0.15) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)
                   if p_bin == 2:
                      upper_bin = to_fit.FindBin(mean - 0.05 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.85 *sig)
                   if p_bin == 3:
                      upper_bin = to_fit.FindBin(mean + 0.05 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.90 *sig)
                   if p_bin == 4:
                      upper_bin = to_fit.FindBin(mean + 0.15 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)
                   if p_bin == 5:
                      upper_bin = to_fit.FindBin(mean + 0.275 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)
                   if p_bin == 6:
                      upper_bin = to_fit.FindBin(mean + 0.65 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)
                   if p_bin == 7:
                      upper_bin = to_fit.FindBin(mean + 0.7 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.975 *sig)
                   if p_bin == 8:
                      upper_bin = to_fit.FindBin(mean + 0.8 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                   if p_bin == 9:
                      upper_bin = to_fit.FindBin(mean + 0.85 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)
                   if p_bin == 10:
                      upper_bin = to_fit.FindBin(mean + 0.9 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)

               if eta_bin == 2:
                   if p_bin == 0:
                      upper_bin = to_fit.FindBin(mean + (-0.45) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)
                   if p_bin == 1:
                      upper_bin = to_fit.FindBin(mean + (-0.5) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)
                   if p_bin == 2:
                      upper_bin = to_fit.FindBin(mean - 0.4 * sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 * sig)
                   if p_bin == 3:
                      upper_bin = to_fit.FindBin(mean - 0.225 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                   if p_bin == 4:
                      upper_bin = to_fit.FindBin(mean - 0.1 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.85 *sig)
                   if p_bin == 5:
                      upper_bin = to_fit.FindBin(mean + 0.1 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.85 *sig)
                   if p_bin == 6:
                      upper_bin = to_fit.FindBin(mean + 0.45 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.15 *sig)
                   if p_bin == 7:
                      upper_bin = to_fit.FindBin(mean + 0.5 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)
                   if p_bin == 8:
                      upper_bin = to_fit.FindBin(mean + 0.65 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)
                   if p_bin == 9:
                      upper_bin = to_fit.FindBin(mean + 0.825 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)

               if eta_bin == 3:
                   if p_bin == 0:
                      upper_bin = to_fit.FindBin(mean + (-0.7) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.25 *sig)
                   if p_bin == 1:
                      upper_bin = to_fit.FindBin(mean + (-0.4) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)
                   if p_bin == 2:
                      upper_bin = to_fit.FindBin(mean - 0.2 * sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 * sig)
                   if p_bin == 3:
                      upper_bin = to_fit.FindBin(mean - 0.0 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.875 *sig)
                   if p_bin == 4:
                      upper_bin = to_fit.FindBin(mean + 0.05 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.925 *sig)
                   if p_bin == 5:
                      upper_bin = to_fit.FindBin(mean + 0.15 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)
                   if p_bin == 6:
                      upper_bin = to_fit.FindBin(mean + 0.5 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)
                   if p_bin == 7:
                      upper_bin = to_fit.FindBin(mean + 0.7 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.8 *sig)
                   if p_bin == 8:
                      upper_bin = to_fit.FindBin(mean + 0.85 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                   if p_bin == 9:
                      upper_bin = to_fit.FindBin(mean + 1.0 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)

               if eta_bin == 4:
                   if p_bin == 0:
                      upper_bin = to_fit.FindBin(mean + (0.7) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.7 *sig)
                   if p_bin == 1:
                      upper_bin = to_fit.FindBin(mean + (0.7) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.7 *sig)
                   if p_bin == 2:
                      upper_bin = to_fit.FindBin(mean + 0.6 * sig)
                      lower_bin = to_fit.FindBin(mean - 0.8 * sig)
                   if p_bin == 3:
                      upper_bin = to_fit.FindBin(mean + 0.72 * sig)
                      lower_bin = to_fit.FindBin(mean - 0.85 *sig)
                   if p_bin == 4:
                      upper_bin = to_fit.FindBin(mean + 0.75 * sig)
                      lower_bin = to_fit.FindBin(mean - 0.85 * sig)
                   if p_bin == 5:
                      upper_bin = to_fit.FindBin(mean + 0.7 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                   if p_bin == 6:
                      upper_bin = to_fit.FindBin(mean + 0.7 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                   if p_bin == 7:
                      upper_bin = to_fit.FindBin(mean + 0.7 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                   if p_bin == 8:
                      upper_bin = to_fit.FindBin(mean + 0.7 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.8 *sig)
                   if p_bin == 9:
                      upper_bin = to_fit.FindBin(mean + 0.8 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)
                   if p_bin == 10:
                      upper_bin = to_fit.FindBin(mean + 0.9 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)

            if "20TRT" in sel_type:
               if eta_bin == 0:
                   if p_bin == 1:
                      upper_bin = to_fit.FindBin(mean + (0.25) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.75 *sig)
                   if p_bin == 2:
                      upper_bin = to_fit.FindBin(mean + 0.15 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.8 *sig)
                   if p_bin == 3:
                      upper_bin = to_fit.FindBin(mean + 0.1 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.8 *sig)
                   if p_bin == 4:
                      upper_bin = to_fit.FindBin(mean + 0.075 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                   if p_bin == 5:
                      upper_bin = to_fit.FindBin(mean + 0.35 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)
                   if p_bin == 6:
                      upper_bin = to_fit.FindBin(mean + 0.4 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.10 *sig)
                   if p_bin == 7:
                      upper_bin = to_fit.FindBin(mean + 0.60 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.1 *sig)
                   if p_bin == 8:
                      upper_bin = to_fit.FindBin(mean + 0.65 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.1 *sig)
                   if p_bin == 9:
                      upper_bin = to_fit.FindBin(mean + 0.75 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)

               if eta_bin == 1:
                   if p_bin == 1:
                      upper_bin = to_fit.FindBin(mean + (0.1) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                   if p_bin == 2:
                      upper_bin = to_fit.FindBin(mean + 0.075 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.75 *sig)
                   if p_bin == 3:
                      upper_bin = to_fit.FindBin(mean + 0.15 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.775 *sig)
                   if p_bin == 4:
                      upper_bin = to_fit.FindBin(mean + 0.15 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.85 *sig)
                   if p_bin == 5:
                      upper_bin = to_fit.FindBin(mean + 0.2 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)
                   if p_bin == 6:
                      upper_bin = to_fit.FindBin(mean + 0.25 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)
                   if p_bin == 7:
                      upper_bin = to_fit.FindBin(mean + 0.4 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.975 *sig)
                   if p_bin == 8:
                      upper_bin = to_fit.FindBin(mean + 0.7 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.025 *sig)
                   if p_bin == 9:
                      upper_bin = to_fit.FindBin(mean + 0.5 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)
                   if p_bin == 10:
                      upper_bin = to_fit.FindBin(mean + 0.9 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)

               if eta_bin == 2:
                   if p_bin == 0:
                      upper_bin = to_fit.FindBin(mean + (-0.15) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.6 *sig)
                   if p_bin == 1:
                      upper_bin = to_fit.FindBin(mean + (-0.3) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 *sig)
                   if p_bin == 2:
                      upper_bin = to_fit.FindBin(mean + 0.1 * sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 * sig)
                   if p_bin == 3:
                      upper_bin = to_fit.FindBin(mean - 0.1 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)
                   if p_bin == 4:
                      upper_bin = to_fit.FindBin(mean + 0.05 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.75 *sig)
                   if p_bin == 5:
                      upper_bin = to_fit.FindBin(mean + 0.2 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.76 *sig)
                   if p_bin == 6:
                      upper_bin = to_fit.FindBin(mean + 0.25 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)
                   if p_bin == 7:
                      upper_bin = to_fit.FindBin(mean + 0.35 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.0 *sig)
                   if p_bin == 8:
                      upper_bin = to_fit.FindBin(mean + 0.4 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)
                   if p_bin == 9:
                      upper_bin = to_fit.FindBin(mean + 0.4 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)
                   if p_bin == 10:
                      upper_bin = to_fit.FindBin(mean + 0.6 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.95 *sig)

               if eta_bin == 3:
                   if p_bin == 0:
                      upper_bin = to_fit.FindBin(mean + (0.1) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.7 *sig)
                   if p_bin == 1:
                      upper_bin = to_fit.FindBin(mean + (0.1) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.8 *sig)
                   if p_bin == 2:
                      upper_bin = to_fit.FindBin(mean + 0.15 * sig)
                      lower_bin = to_fit.FindBin(mean - 0.9 * sig)
                   if p_bin == 3:
                      upper_bin = to_fit.FindBin(mean + 0.125 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.725 *sig)
                   if p_bin == 4:
                      upper_bin = to_fit.FindBin(mean + 0.125 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.775 *sig)
                   if p_bin == 5:
                      upper_bin = to_fit.FindBin(mean + 0.075 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.875 *sig)
                   if p_bin == 6:
                      upper_bin = to_fit.FindBin(mean + 0.125 *sig)
                      lower_bin = to_fit.FindBin(mean - 0.975 *sig)
                   if p_bin == 7:
                      upper_bin = to_fit.FindBin(mean + 0.15 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)
                   if p_bin == 8:
                      upper_bin = to_fit.FindBin(mean + 0.15 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.10 *sig)
                   if p_bin == 9:
                      upper_bin = to_fit.FindBin(mean + 0.55 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)
                   if p_bin == 10:
                      upper_bin = to_fit.FindBin(mean + 0.7 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)
                   if p_bin == 11:
                      upper_bin = to_fit.FindBin(mean + 0.8 *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)

               if eta_bin == 4:
                   if p_bin == 0:
                      upper_bin = to_fit.FindBin(mean + (-0.25) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.8 *sig)
                   if p_bin == 1:
                      upper_bin = to_fit.FindBin(mean + (-0.2) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.8 *sig)
                   if p_bin == 2:
                      upper_bin = to_fit.FindBin(mean + (-0.2) * sig)
                      lower_bin = to_fit.FindBin(mean - 0.8 * sig)
                   if p_bin == 3:
                      upper_bin = to_fit.FindBin(mean + (-0.15) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.8 *sig)
                   if p_bin == 4:
                      upper_bin = to_fit.FindBin(mean + (-0.25) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.75 *sig)
                   if p_bin == 5:
                      upper_bin = to_fit.FindBin(mean + (-0.05) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.85 *sig)
                   if p_bin == 6:
                      upper_bin = to_fit.FindBin(mean + (-0.05) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.90 *sig)
                   if p_bin == 7:
                      upper_bin = to_fit.FindBin(mean + (0.0) *sig)
                      lower_bin = to_fit.FindBin(mean - 0.90 *sig)
                   if p_bin == 8:
                      upper_bin = to_fit.FindBin(mean + (0.0) * sig)
                      lower_bin = to_fit.FindBin(mean - 0.90 * sig)
                   if p_bin == 9:
                      upper_bin = to_fit.FindBin(mean + (0.15) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)
                   if p_bin == 10:
                      upper_bin = to_fit.FindBin(mean + (0.2) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.05 *sig)
                   if p_bin == 11:
                      upper_bin = to_fit.FindBin(mean + (0.3) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.1 *sig)
                   if p_bin == 12:
                      upper_bin = to_fit.FindBin(mean + (0.4) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.15 *sig)
                   if p_bin == 13:
                      upper_bin = to_fit.FindBin(mean + (0.5) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.15 *sig)
                   if p_bin == 14:
                      upper_bin = to_fit.FindBin(mean + (0.65) *sig)
                      lower_bin = to_fit.FindBin(mean - 1.2 *sig)

        #do the fit in this range.
        sigma_up = to_fit.GetBinCenter(upper_bin)
        sigma_down = to_fit.GetBinCenter(lower_bin)
        eop.setRange("Fit", sigma_down, sigma_up)
        print("MPV: {}".format(mpv))
        print("Fitting in range [{},{}]".format(sigma_down, sigma_up))
        if function == "gaus":
            model, variables = generate_gaus(eop)
        if function == "two_gaus":
            model, variables = generate_two_gaus(eop)
        elif function == "landau":
            model, variables = generate_landau(eop)
        elif function == "dcb":
            model, variables = generate_dcb(eop)
        elif function == "landauxgaus" or function == "gausxlandau":
            model, variables = generate_landau_gaus(eop)
        elif function == "landau+gaus" or function == "gaus+landau":
             model, variables = generate_landau_plus_gaus(eop)

        for var in variables:
            if "mean" in var.GetName():
                var.setVal(mpv)
                break

        print(model)
        prepare_for_fit()
        result = model.fitTo(eop_hist, ROOT.RooFit.Range("Fit"), ROOT.RooFit.Save(True))
        fit_results[channel]=result

        if not montecarlo_errors:
            for var in variables:
                if "mean" in var.GetName():
                    mpvs[channel]=var.getVal()
                    mpv_errs[channel]=var.getError()
        else:
            mpvs[channel], mpv_errs[channel]=montecarlo_uncertainties(model, variables,eop, to_fit.GetBinCenter(lower_bin), to_fit.GetBinCenter(upper_bin))

        #draw the fit
        c = ROOT.TCanvas("canv", "canv")
        c.Draw()

        top = ROOT.TPad("top", "top", 0.0, 0.3, 1.0, 1.0)
        top.SetLeftMargin(0.15)
        top.SetBottomMargin(0.0)
        top.Draw()
        top.cd()

        f=eop.frame()
        eop_hist.plotOn(f)
        model.plotOn(f)
        pull = f.pullHist()
        model.plotOn(f, ROOT.RooFit.Components("gausone"), ROOT.RooFit.LineStyle(ROOT.kDashed))
        model.plotOn(f, ROOT.RooFit.Components("landau"), ROOT.RooFit.LineStyle(ROOT.kDashed))
        model.plotOn(f, ROOT.RooFit.Components("gaustwo"), ROOT.RooFit.LineStyle(ROOT.kDashed))

        n_param = result.floatParsFinal().getSize()
        reduced_chi_square = f.chiSquare(n_param)

        f.Draw()
        legend = ROOT.TLegend(0.5, 0.7, 0.89, 0.89)
        legend.SetBorderSize(0)  # no border
        legend.SetFillStyle(0)  # make transparent
        legend.Draw()
        legend.AddEntry(None,\
                '#chi^{2}' + ' / {} = {:.3f}'.format( "nDOF",\
                reduced_chi_square), '')

        f.GetXaxis().SetTitleFont(43)
        f.GetXaxis().SetTitleSize(25)
        f.GetXaxis().SetTitleOffset(1.0)
        f.GetXaxis().SetLabelSize(20)
        f.GetXaxis().SetLabelFont(43)
        f.GetYaxis().SetLabelSize(13)
        f.GetYaxis().SetLabelFont(43)
        f.GetYaxis().SetTitleSize(20)
        f.GetYaxis().SetTitleFont(43)

        DrawText(0.60, 0.5,"Bin: {}".format(p_bin) ,  size = 0.18)

        pull.SetMarkerSize(0.3)
        bottom=ROOT.TPad("bottom", "bottom", 0.0, 0.0, 1.0, 0.3)
        bottom.SetBottomMargin(0.4)
        bottom.SetLeftMargin(0.15)
        bottom.SetTopMargin(0.0)
        c.cd()
        bottom.Draw()
        bottom.cd()
        frame2 = eop.frame()
        frame2.addPlotable(pull, "P")
        frame2.Draw()
        frame2.GetXaxis().SetTitle("E/P")
        frame2.GetYaxis().SetTitle("#frac{fit - data}/{#sigma_{data}}")
        frame2.GetXaxis().SetTitleFont(43)
        frame2.GetXaxis().SetTitleSize(25)
        frame2.GetXaxis().SetTitleOffset(2.7)
        frame2.GetXaxis().SetLabelSize(20)
        frame2.GetXaxis().SetLabelFont(43)
        frame2.GetYaxis().SetLabelSize(13)
        frame2.GetYaxis().SetLabelFont(43)
        frame2.GetYaxis().SetTitleSize(20)
        frame2.GetYaxis().SetTitleFont(43)

        bottom.Update()
        bottom.Modified()

        c.Update()
        c.Modified()
        if extra_str != "":
           c.Print("{}_{}_fit_plot.png".format(channel,extra_str))
        top.Close()
        bottom.Close()
        c.Close()

    return mpvs, mpv_errs, fit_results

def test_fit(f, histogram_base = "EOPDistribution", selection_name = "MIPSelectionHadFracAbove70", function = "gaus", montecarlo_errors = True, sel_type=""):

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
        bins = []
        bin_numbers = []
        for j, p_low, p_high in zip(range(0, len(p_bins_low_for_eta_bin[i])),p_bins_low_for_eta_bin[i], p_bins_high_for_eta_bin[i]):
            if j==0:
                continue
            if "MIP" in selection_name:
                if i == 0 and j == 10:
                    break
                if i == 1 and j == 10:
                    break
                if i == 2 and j == 11:
                    break
                if i == 3 and j  == 12:
                    break
                if i == 4 and j == 20:
                    break
            if "20TRT" in selection_name:
                if i == 0 and j == 10:
                    break
                if i == 1 and j == 11:
                    break
                if i == 2 and j == 12:
                    break
                if i == 3 and j == 12:
                    break
                if i == 4 and j == 14:
                    pass
>>>>>>> Update plotting code
            histogram_name = histogram_base + "_" + selection_name + "_Eta_" + str(i) + "_Momentum_" + str(j)
            histograms = HM.getHistograms(histogram_name)
            to_fit_function = function
            mpv, mpv_err, fit_result = do_fit(histograms, function = to_fit_function, extra_str = "Eta_{}_P_{}_{}".format(i,j, selection_name), montecarlo_errors = montecarlo_errors, p_bin = j, eta_bin=i, sel_type = sel_type)
            for channel in mpv:
               if channel not in mpvs:
                   mpvs[channel]=[]
                   mpv_errs[channel]=[]
               mpvs[channel].append(mpv[channel])
               mpv_errs[channel].append(mpv_err[channel])
            bins.append(p_low)
            bin_numbers.append(j)
        bins.append(p_high)

        hist_name = histogram_base + "{}_{}_".format(function, "fit") + selection_name + "_Eta_" + str(i)
        bin_array = array.array('d', bins)
        histograms = {}
        for channel in mpvs:
            histograms[channel] = ROOT.TH1D(hist_name + channel, hist_name + channel, len(bins)-1, bin_array)
            for n in range(0, len(mpvs[channel])):
                histograms[channel].SetBinContent(n+1, mpvs[channel][n])
                histograms[channel].SetBinError(n+1, mpv_errs[channel][n])

        MCKeys = ["PythiaJetJet"]
        DataKey="LowMuData"
        if "20TRT" in selection_name:
            MCKeys = ["PythiaJetJetTightIso"]
            DataKey = "LowMuDataTightIso"

        base_description = ["P_{T} Reweighted"]
        if "MIP" in selection_name:
            description = base_description + ["MIP Selection", "{:.1f} < |#eta| < {:.1f}".format(eta_low, eta_high)]
        elif "20TRT" in selection_name:
            description = base_description + ["N_{TRT} > 20","Tight Isolation", "E_{Total} != 0", "{:.1f} < |#eta| < {:.1f}".format(eta_low, eta_high)]
        channelLabels = {"SinglePion": "Single Pion", "PythiaJetJet" : "#splitline{Pythia8}{MinBias and Dijet}", DataKey: "2017 Low-<#mu> Data", "PythiaJetJetPionsReweighted":"Pythia8 MB+DJ Pions Only", "PythiaJetJetHardScatter":"Pythia8 MB+DJ Truth Matched", "PythiaJetJetTightIso": "#splitline{Pythia8}{MinBias and Dijet}", "LowMuDataTightIso":"2017 Low-<#mu> Data"}
        ratio_min = 0.9
        ratio_max = 1.1
        if "20TRT" in selection_name:
            ratio_min = 0.9
            ratio_max = 1.1
        if "MIP" in selection_name and (i == 0 or i == 1 or i == 3):
            ratio_min = 0.95
            ratio_max = 1.05
        clones = {}
        for c in histograms:
            clones[c] = histograms[c].Clone(histograms[c].GetName() + "Clone")
        #draw the fit results as a set of histograms
        canvas, top, bottom = DrawDataVsMC(histograms,\
                       channelLabels,\
                       MCKeys = MCKeys,\
                       ratio_min = ratio_min,\
                       ratio_max = ratio_max,\
                       doLogy=False,\
                       doLogx=True,\
                       ylabel = "MPV(E/P)",\
                       xlabel = "P [GeV]",\
                       DataKey=DataKey,\
                       extra_description = description)
        canvas_name = hist_name + "_"  + selection_name + "_Eta_" + str(i) + ".png"
        canvas.Print(canvas_name)
        top.Close()
        bottom.Close()
        canvas.Close()

        if "MIP" in selection_name:
            bkg_eop = HM.getHistograms("EnergyBkgProfileVsMomentum" +"__" + selection_name + "_Eta_"+str(i))
        elif "20TRT" in selection_name:
            bkg_eop = HM.getHistograms("EnergyBigBkgProfileVsMomentum" +"__" + selection_name + "_Eta_"+str(i))
        bkg_eop_clone = {}
        for channel in mpvs:
            bkg_eop_clone[channel] = bkg_eop[channel]
        bkg_eop = ProjectProfiles(bkg_eop_clone)
        bkg_eop_subrange = {}
        for channel in bkg_eop:
            bkg_eop_subrange[channel] = ROOT.TH1D(hist_name + channel + "bkg", hist_name + channel, len(bins)-1, bin_array)
            for a,b in enumerate(bin_numbers):
                bkg_eop_subrange[channel].SetBinContent(a+1,bkg_eop[channel].GetBinContent(b+1))
                bkg_eop_subrange[channel].SetBinError(a+1,bkg_eop[channel].GetBinError(b+1))
        corr_eop = SubtractHistograms(clones, bkg_eop_subrange)

        stuff = DrawDataVsMC(corr_eop,\
                                channelLabels,\
                                MCKeys = MCKeys,\
                                DataKey=DataKey,\
                                doLogy=False,\
                                doLogx=True,\
                                ratio_min=ratio_min,\
                                ratio_max=ratio_max,\
                                ylabel="MPV(E/P)_{CORR}",\
                                xlabel = "P [GeV]",\
                                extra_description = description)


        for s in stuff:
            s.Update()
            s.Modified()
        canvas = stuff[0]
        canvas.Update()
        canvas.Modified()
        canvas.Print(hist_name + "corrected_"  + selection_name + ".png")
        canvas.Draw()
        histogram_name += ("_".join(MCKeys) + "_{}".format(DataKey))
        stuff[1].Close()
        stuff[2].Close()
        canvas.Close()

        ratio_min = 0.9
        ratio_max = 1.1
        if "20TRT" in selection_name:
            ratio_min = 0.9
            ratio_max = 1.1
        if "MIP" in selection_name and (i == 0 or i == 1 or i == 2):
            ratio_min = 0.95
            ratio_max = 1.05

        average_eop = HM.getHistograms("EOPProfileVsMomentum" +"__" + selection_name + "_Eta_"+str(i))
        average_eop_clone = {}
        for channel in mpvs:
            average_eop_clone[channel] =  ROOT.TH1D(hist_name + channel + "bkg", hist_name + channel, len(bins)-1, bin_array)
            for a,b in enumerate(bin_numbers):
                average_eop_clone[channel].SetBinContent(a+1,average_eop[channel].GetBinContent(b+1))
                average_eop_clone[channel].SetBinError(a+1,average_eop[channel].GetBinError(b+1))
        corr_average_eop = SubtractHistograms(average_eop_clone, bkg_eop_subrange)

        stuff = DrawDataVsMC(corr_average_eop,\
                                channelLabels,\
                                MCKeys = MCKeys,\
                                DataKey=DataKey,\
                                doLogy=False,\
                                doLogx=True,\
                                ratio_min=ratio_min,\
                                ratio_max=ratio_max,\
                                ylabel="<E/P>_{CORR}",\
                                xlabel = "P [GeV]",\
                                extra_description = description)


        for s in stuff:
            s.Update()
            s.Modified()
        canvas = stuff[0]
        canvas.Update()
        canvas.Modified()
        canvas.Print("EOPProfileVsMomentum_" + "corrected_"  + selection_name +  "_Eta_"+str(i) + ".png")
        canvas.Draw()
        histogram_name += ("_".join(MCKeys) + "_{}".format(DataKey))
        stuff[1].Close()
        stuff[2].Close()
        canvas.Close()

        stuff = DrawDataVsMC(bkg_eop_subrange,\
                                channelLabels,\
                                MCKeys = MCKeys,\
                                DataKey=DataKey,\
                                doLogy=False,\
                                doLogx=True,\
                                ratio_min=0.5,\
                                ratio_max=1.5,\
                                ylabel="<E/P>_{BKG}",\
                                xlabel = "P [GeV]",\
                                extra_description = description)


        for s in stuff:
            s.Update()
            s.Modified()
        canvas = stuff[0]
        canvas.Update()
        canvas.Modified()
        canvas.Print("EOPBkgEstimateProfileVsMomentum_" + selection_name +  "_Eta_"+str(i) + ".png")
        canvas.Draw()
        histogram_name += ("_".join(MCKeys) + "_{}".format(DataKey))
        stuff[1].Close()
        stuff[2].Close()
        canvas.Close()


if __name__ == "__main__":
    f="pt_reweighted.root"
    #for selection_name in ["MIPSelectionHadFracAbove70", "20TRTHitsNonZeroEnergy"][::-1]:
    #for selection_name in ["MIPSelectionHadFracAbove70"]:
    for selection_name in ["20TRTHitsNonZeroEnergy"]:
        test_fit(f, selection_name = selection_name, function="gaus", sel_type = selection_name, montecarlo_errors = False)

