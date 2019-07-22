import ROOT
import array

def scotts_rule(histogram):
    N = histogram.Integral()
    sigma = histogram.GetRMS()
    width = (sigma)/(N**(1.0/3.0))
    return width


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


