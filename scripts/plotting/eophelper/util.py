import os
import numpy as np
from array import array
from ROOT import *

def ensure_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

def draw_atlas_label(x, y, text, color = kBlack):
    l = TLatex()
    l.SetNDC()
    l.SetTextAlign(13)
    l.SetTextFont(72)
    l.SetTextColor(color)

    delx = 0.115*696*gPad.GetWh()/(472*gPad.GetWw())

    l.DrawLatex(x,y,"ATLAS")
    if text:
        p = TLatex()
        p.SetNDC()
        p.SetTextAlign(13)
        p.SetTextFont(42)
        p.SetTextColor(color)
        p.DrawLatex(x+delx,y,text)

def draw_latex(x, y, text, size, color = kBlack):
    l = TLatex()
    l.SetNDC()
    l.SetTextFont(42)
    l.SetTextAlign(13)
    l.SetTextSize(size)
    l.SetTextColor(color)
    l.DrawLatex(x,y,text)


def tprofile2hist(p):
    p.Sumw2()
    nbins = p.GetNbinsX()
    xlabel = p.GetXaxis().GetTitle()
    xbins = np.zeros(nbins+1)
    for i in xrange(nbins):
        xbins[i] = p.GetXaxis().GetBinLowEdge(i+1)
    xbins[-1] = p.GetXaxis().GetBinLowEdge(nbins)+p.GetXaxis().GetBinWidth(nbins)
    # print "---------------> ", xbins
    # print len(xbins)
    # print nbins
    h = TH1F(p.GetName()+"_hist", xlabel, nbins, array("d", xbins))
    h.Sumw2()
    for i in xrange(1,nbins+1):
        h.SetBinContent(i, p.GetBinContent(i))
        h.SetBinError(i, p.GetBinError(i))
    return h

def get_bg_cor_hist(h, h_bg):
    mean_y = get_mean_y(h)
    bg_mean_y = get_mean_y(h_bg)
    cor_y = mean_y - 4/3*bg_mean_y
    nbins = h.GetNbinsX()
    xlabel = h.GetXaxis().GetTitle()
    xmin = h.GetXaxis().GetBinCenter(1)
    xmax = h.GetXaxis().GetBinCenter(nbins)
    h_cor = TH1F(h.GetName()+"_cor", h.GetXaxis().GetTitle(), nbins, xmin, xmax)
    h_cor.Sumw2()
    return h_cor

def get_mean_x(h2):
    nbinsx = h2.GetNbinsX()
    nbinsy = h2.GetNbinsY()
    array = np.zeros(nbinsx)
    for i in xrange(1,nbinsy):
        mean = 0.
        n = 0.
        for j in xrange(1,nbinsx):
            mean = mean + h2.GetBinContent(i,j)*h2.GetXaxis().GetBinCenter(j)
            n = n + h2.GetBinContent(i,j)
        if n != 0:
            mean = mean/n
        else:
            mean = 0

        array[i] = mean
    return array

def get_mean_y(h2):
    nbinsx = h2.GetNbinsX()
    nbinsy = h2.GetNbinsY()
    array = np.zeros(nbinsx)
    for i in xrange(1,nbinsx):
        mean = 0.
        n = 0.
        for j in xrange(1,nbinsy):
            mean = mean + h2.GetBinContent(i,j)*h2.GetYaxis().GetBinCenter(j)
            n = n + h2.GetBinContent(i,j)
        if n != 0:
            mean = mean/n
        else:
            mean = 0

        array[i] = mean
    return array

def get_1dhist_y(h2, name, xmin, xmax):
    nbinsy = h2.GetNbinsY()
    ylabel = h2.GetYaxis().GetTitle()
    ymin = h2.GetYaxis().GetBinLowEdge(1)
    ymax = h2.GetYaxis().GetBinLowEdge(nbinsy)+h2.GetYaxis().GetBinWidth(nbinsy)
    h = TH1F(name, ylabel, nbinsy, ymin, ymax)
    print nbinsy
    print ylabel
    print ymin
    print ymax
    for i in xrange(1,nbinsy+1):
        sumx = 0
        errx = 0
        for xbin in xrange(xmin,xmax+1):
            sumx += h2.GetBinContent(xbin,i)
            errx += h2.GetBinError(xbin,i)**2
        errx = np.sqrt(errx)
        h.SetBinContent(i, sumx)
        h.SetBinError(i, errx)
    return h
