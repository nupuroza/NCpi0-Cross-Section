## ########################################################
## Written by: Cricket Bergner
## Date Created: 06/26/24
## Purpose: to plot the difference 
##          between the response matrix 
##          and variation spectra histograms
## ########################################################

## ############################
## imports
## ############################

import ROOT
import array
from customHistAndPlotMethods import Rebin as r                 ## r(rebin, reference, name)
from customHistAndPlotMethods import DrawWithOverflow as d      ## d(hist, canvas, draw options)
ROOT.gStyle.SetOptStat(0)

## ############################
## variables
## ############################

## 2g1p
filepath_rm_2g1p = "/app/users/crbergner/data/response_matrices/response_matrices_exclusive.root"
filepath_vs_2g1p = "/app/users/crbergner/data/variation_spectra/SBNfit_variation_spectra_exclusive_2g1p.root"
file_rm_2g1p = ROOT.TFile(filepath_rm_2g1p)
file_vs_2g1p = ROOT.TFile(filepath_vs_2g1p)

## 2g0p
filepath_rm_2g0p = "/app/users/crbergner/data/response_matrices/response_matrices_exclusive.root"
filepath_vs_2g0p = "/app/users/crbergner/data/variation_spectra/SBNfit_variation_spectra_exclusive_2g0p.root"
file_rm_2g0p = ROOT.TFile(filepath_rm_2g0p)
file_vs_2g0p = ROOT.TFile(filepath_vs_2g0p)

## 2gNp
filepath_rm_2gXp = "/app/users/crbergner/data/response_matrices/response_matrices_exclusive.root"
filepath_vs_2gXp = "/app/users/crbergner/data/variation_spectra/SBNfit_variation_spectra_inclusive_2gXp.root"
file_rm_2gXp = ROOT.TFile(filepath_rm_2gXp)
file_vs_2gXp = ROOT.TFile(filepath_vs_2gXp)

## ############################
## 2g1p histograms
## ############################

true_rm_2g1p = file_rm_2g1p.Get("htrue_2g1p")
true_rm_2g1p = ROOT.TH1D(true_rm_2g1p)
true_vs_2g1p = file_vs_2g1p.Get("exclusive_2g1p_CV_Dir/Sys2g1p_denominator_truth_Signal")
true_vs_2g1p = ROOT.TH1D(true_vs_2g1p)
true_vs_rebin_2g1p = r(true_vs_2g1p, true_rm_2g1p, "true_vs_rebin_2g1p")

reco_rm_2g1p = file_rm_2g1p.Get("hreco_2g1p")
reco_rm_2g1p = ROOT.TH1D(reco_rm_2g1p)
reco_vs_2g1p = file_vs_2g1p.Get("exclusive_2g1p_CV_Dir/Sys2g1p_numerator_reco_Signal")
reco_vs_2g1p = ROOT.TH1D(reco_vs_2g1p)
reco_vs_rebin_2g1p = r(reco_vs_2g1p, reco_rm_2g1p, "true_vs_rebin_2g1p")

## subtract hists
diff_true_2g1p = true_rm_2g1p - true_vs_rebin_2g1p
diff_reco_2g1p = reco_rm_2g1p - reco_vs_rebin_2g1p

## ############################
## 2g0p histograms
## ############################

true_rm_2g0p = file_rm_2g0p.Get("htrue_2g0p")
true_rm_2g0p = ROOT.TH1D(true_rm_2g0p)
true_vs_2g0p = file_vs_2g0p.Get("exclusive_2g0p_CV_Dir/Sys2g0p_denominator_truth_Signal")
true_vs_2g0p = ROOT.TH1D(true_vs_2g0p)
true_vs_rebin_2g0p = r(true_vs_2g0p, true_rm_2g0p, "true_vs_rebin_2g0p")

reco_rm_2g0p = file_rm_2g0p.Get("hreco_2g0p")
reco_rm_2g0p = ROOT.TH1D(reco_rm_2g0p)
reco_vs_2g0p = file_vs_2g0p.Get("exclusive_2g0p_CV_Dir/Sys2g0p_numerator_reco_Signal")
reco_vs_2g0p = ROOT.TH1D(reco_vs_2g0p)
reco_vs_rebin_2g0p = r(reco_vs_2g0p, reco_rm_2g0p, "reco_vs_rebin_2g0p")

## subtract hists
diff_true_2g0p = true_rm_2g0p - true_vs_rebin_2g0p
diff_reco_2g0p = reco_rm_2g0p - reco_vs_rebin_2g0p

## ############################
## 2gXp histograms
## ############################

true_rm_2gXp = file_rm_2gXp.Get("htrue_2gXp")
true_rm_2gXp = ROOT.TH1D(true_rm_2gXp)
true_vs_2gXp = file_vs_2gXp.Get("inclusive_2gXp_CV_Dir/Sys2gXp_denominator_truth_Signal")
true_vs_2gXp = ROOT.TH1D(true_vs_2gXp)
true_vs_rebin_2gXp = r(true_vs_2gXp, true_rm_2gXp, "true_vs_rebin_2gXp")

reco_rm_2gXp = file_rm_2gXp.Get("hreco_2gXp")
reco_rm_2gXp = ROOT.TH1D(reco_rm_2gXp)
reco_vs_2gXp = file_vs_2gXp.Get("inclusive_2gXp_CV_Dir/Sys2gXp_numerator_reco_Signal")
reco_vs_2gXp = ROOT.TH1D(reco_vs_2gXp)
reco_vs_rebin_2gXp = r(reco_vs_2gXp, reco_rm_2gXp, "reco_vs_rebin_2gXp")


## subtract hists
diff_true_2gXp = true_rm_2gXp - true_vs_rebin_2gXp
# diff_reco_2gXp = reco_rm_2gXp - reco_vs_rebin_2gXp

## ############################
## canvases
## ############################

diff_true_2g1p_canvas = ROOT.TCanvas("diff_true_2g1p_canvas", "diff_true_2g1p", 1000, 750)
diff_true_2g0p_canvas = ROOT.TCanvas("diff_true_2g0p_canvas", "diff_true_2g0p", 1000, 750)
diff_true_2gXp_canvas = ROOT.TCanvas("diff_true_2gXp_canvas", "diff_true_2gXp", 1000, 750)
diff_reco_2g1p_canvas = ROOT.TCanvas("diff_reco_2g1p_canvas", "diff_reco_2g1p", 1000, 750)
diff_reco_2g0p_canvas = ROOT.TCanvas("diff_reco_2g0p_canvas", "diff_reco_2g0p", 1000, 750)
# diff_reco_2gXp_canvas = ROOT.TCanvas("diff_reco_2gXp_canvas", "diff_reco_2gXp", 1000, 750)
rm_reco_2g1p_canvas = ROOT.TCanvas("rm_reco_2g1p_canvas", "rm_reco_2g1p", 1000, 750)
vs_reco_2g1p_canvas = ROOT.TCanvas("vs_reco_2g1p_canvas", "vs_reco_2g1p", 1000, 750)
rm_reco_2g0p_canvas = ROOT.TCanvas("rm_reco_2g0p_canvas", "rm_reco_2g0p", 1000, 750)
vs_reco_2g0p_canvas = ROOT.TCanvas("vs_reco_2g0p_canvas", "vs_reco_2g0p", 1000, 750)
rm_reco_2gXp_canvas = ROOT.TCanvas("rm_reco_2gXp_canvas", "rm_reco_2gXp", 1000, 750)
vs_reco_2gXp_canvas = ROOT.TCanvas("vs_reco_2gXp_canvas", "vs_reco_2gXp", 1000, 750)

## ############################
## d function
## ############################

overflow_diff_true_2g1p = d(diff_true_2g1p, diff_true_2g1p_canvas, "hist")
overflow_diff_reco_2g1p = d(diff_reco_2g1p, diff_reco_2g1p_canvas, "hist")
overflow_diff_true_2g0p = d(diff_true_2g0p, diff_true_2g0p_canvas, "hist")
overflow_diff_reco_2g0p = d(diff_reco_2g0p, diff_reco_2g0p_canvas, "hist")
overflow_diff_true_2gXp = d(diff_true_2gXp, diff_true_2gXp_canvas, "hist")
# overflow_diff_reco_2gXp = d(diff_reco_2gXp, diff_reco_2gXp_canvas, "hist")
overflow_reco_2g1p_rm = d(reco_rm_2g1p, rm_reco_2g1p_canvas, "hist")
overflow_reco_2g1p_vs = d(reco_vs_rebin_2g1p, vs_reco_2g1p_canvas, "hist")
overflow_reco_2g0p_rm = d(reco_rm_2g0p, rm_reco_2g0p_canvas, "hist")
overflow_reco_2g0p_vs = d(reco_vs_rebin_2g0p, vs_reco_2g0p_canvas, "hist")
overflow_reco_2gXp_rm = d(reco_rm_2gXp, rm_reco_2gXp_canvas, "hist")
overflow_reco_2gXp_vs = d(reco_vs_rebin_2gXp, vs_reco_2gXp_canvas, "hist")

## ############################
## printing images
## ############################

fout_true_2g1p = "/app/users/crbergner/work/NCpi0-Cross-Section/diff_true_2g1p.png"
fout_reco_2g1p = "/app/users/crbergner/work/NCpi0-Cross-Section/diff_reco_2g1p.png"
fout_true_2g0p = "/app/users/crbergner/work/NCpi0-Cross-Section/diff_true_2g0p.png"
fout_reco_2g0p = "/app/users/crbergner/work/NCpi0-Cross-Section/diff_reco_2g0p.png"
fout_true_2gXp = "/app/users/crbergner/work/NCpi0-Cross-Section/diff_true_2gXp.png"
# fout_reco_2gXp = "/app/users/crbergner/work/NCpi0-Cross-Section/diff_reco_2gXp.png"
fout_rm_reco_2g1p = "/app/users/crbergner/work/NCpi0-Cross-Section/rm_reco_2g1p.png"
fout_vs_reco_2g1p = "/app/users/crbergner/work/NCpi0-Cross-Section/vs_reco_2g1p.png"
fout_rm_reco_2g0p = "/app/users/crbergner/work/NCpi0-Cross-Section/rm_reco_2g0p.png"
fout_vs_reco_2g0p = "/app/users/crbergner/work/NCpi0-Cross-Section/vs_reco_2g0p.png"
fout_rm_reco_2gXp = "/app/users/crbergner/work/NCpi0-Cross-Section/rm_reco_2gXp.png"
fout_vs_reco_2gXp = "/app/users/crbergner/work/NCpi0-Cross-Section/vs_reco_2gXp.png"


diff_true_2g1p_canvas.Print(fout_true_2g1p)
diff_reco_2g1p_canvas.Print(fout_reco_2g1p)
diff_true_2g0p_canvas.Print(fout_true_2g0p)
diff_reco_2g0p_canvas.Print(fout_reco_2g0p)
diff_true_2gXp_canvas.Print(fout_true_2gXp)
# diff_reco_2gXp_canvas.Print(fout_reco_2gXp)
rm_reco_2g1p_canvas.Print(fout_rm_reco_2g1p)
vs_reco_2g1p_canvas.Print(fout_vs_reco_2g1p)
rm_reco_2g0p_canvas.Print(fout_rm_reco_2g0p)
vs_reco_2g0p_canvas.Print(fout_vs_reco_2g0p)
rm_reco_2gXp_canvas.Print(fout_rm_reco_2gXp)
vs_reco_2gXp_canvas.Print(fout_vs_reco_2gXp)

## ########################################################