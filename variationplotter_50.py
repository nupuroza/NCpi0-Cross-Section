## ########################################################
## Written by: Cricket Bergner
## Date Created: 06/28/24
## Purpose: to plot first 50 numerator
##          reco signal universes with
##          the cv numerator reco signal
## ########################################################

## ############################
## imports
## ############################

import ROOT
import os
from customHistAndPlotMethods import Rebin as r                 ## r(rebin, reference, name)
from customHistAndPlotMethods import DrawWithOverflow as d      ## d(hist, canvas, draw options)
ROOT.gStyle.SetOptStat(0)

## ############################
## variables
## ############################

cvs_filepath = "/app/users/crbergner/data/variation_spectra/SBNfit_variation_spectra_inclusive_2gXp.root"
mvs_filepath = "/mnt/morrigan/NCPi0_XS_data/SBNfit_variation_spectra_inclusive_2gXp.root"

cvs_file = ROOT.TFile(cvs_filepath)
mvs_file = ROOT.TFile(mvs_filepath)

## ############################
## creating the histograms
## ############################

## cv histograms
cvs_cv = cvs_file.Get("inclusive_2gXp_CV_Dir/Sys2gXp_numerator_reco_Signal")
mvs_cv = mvs_file.Get("inclusive_2gXp_CV_Dir/Sys2gXp_numerator_reco_Signal")

cvs_cv = ROOT.TH1D(cvs_cv)
mvs_cv = ROOT.TH1D(mvs_cv)

## ############################
## canvases
## ############################

cvs_canvas = ROOT.TCanvas("cvs_canvas", "cvs_cv", 1000, 750)
mvs_canvas = ROOT.TCanvas("mvs_canvas", "mvs_cv", 1000, 750)

# overflow_cvs_cv = d(cvs_cv, cvs_canvas, "hist")
# overflow_mvs_cv = d(mvs_cv, mvs_canvas, "hist")

hists = []

for i in range(1,51):
    cvs_ubgenie = cvs_file.Get("inclusive_2gXp_All_UBGenie_Dir/Sys2gXp_numerator_reco_universe_" + str(i) + "_Signal")
    cvs_ubgenie = ROOT.TH1D(cvs_ubgenie)
    overflow_cvs_ubgenie = d(cvs_ubgenie, cvs_canvas, "same hist")
    hists.append(cvs_ubgenie)
    hists.append(overflow_cvs_ubgenie)

    mvs_ubgenie = mvs_file.Get("inclusive_2gXp_All_UBGenie_Dir/Sys2gXp_numerator_reco_universe_" + str(i) + "_Signal")
    mvs_ubgenie = ROOT.TH1D(mvs_ubgenie)
    mvs_ubgenie_rebin = r(mvs_ubgenie, cvs_ubgenie, ("mvs_rebin" + str(i)))
    mvs_ubgenie_rebin.GetYaxis().SetRangeUser(0, 200)
    overflow_mvs_ubgenie = d(mvs_ubgenie_rebin, mvs_canvas, "same hist")
    hists.append(mvs_ubgenie_rebin)
    hists.append(overflow_mvs_ubgenie)

rebin_cvs_cv = r(cvs_cv, cvs_ubgenie, "cvs_rebin")
rebin_mvs_cv = r(mvs_cv, mvs_ubgenie, "mvs_rebin")

rebin_cvs_cv.SetLineColor(ROOT.kBlack)
rebin_mvs_cv.SetLineColor(ROOT.kBlack)

rebin_cvs_cv.SetLineWidth(3)
rebin_mvs_cv.SetLineWidth(3)

hists.append(rebin_cvs_cv)
hists.append(rebin_mvs_cv)

overflow_cvs_cv = d(rebin_cvs_cv, cvs_canvas, "same hist")
overflow_mvs_cv = d(rebin_mvs_cv, mvs_canvas, "same hist")

## ############################
## printing the histograms
## ############################

fout_cvs = "/app/users/crbergner/work/NCpi0-Cross-Section/pictures/cvs.png"
fout_mvs = "/app/users/crbergner/work/NCpi0-Cross-Section/pictures/mvs.png"

cvs_canvas.Print(fout_cvs)
mvs_canvas.Print(fout_mvs)

## ########################################################


