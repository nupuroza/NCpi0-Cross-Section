import ROOT
import argparse
import os
from customHistAndPlotMethods import makeEnv_TCanvas
from calculateChi2 import *

## Plots GENIE and NuWro signal, background, and combined 2g1p and 2g0p reco distributions with associated errors.

## Default GENIE input file and output directory for Leon. Feel free to change it locally. in_file_path must point to output from translatehists.py for proper errors.
in_file_path = "/uboone/data/users/ltong/gLEE/NCPi0/2024-03-06_fixed-minmax-nuwro/2024-03-06_out.root"
out_dir = "/uboone/data/users/ltong/gLEE/NCPi0/GENIEvsNuWro"

## Use command line arguments for input and output if provided.
parser = argparse.ArgumentParser(description='Script to plot GENIE and NuWro signal, background, and combined 2g1p and 2g0p reco distributions with associated errors')
parser.add_argument('in_file_path', help='Path to input directory', type=str,nargs='?')
parser.add_argument('out_dir', help='Path to ouput directory', type=str,nargs='?')
p = parser.parse_args()
if p.in_file_path > 0:
  in_file_path = p.in_file_path
if p.out_dir > 0:
  out_dir = p.out_dir

## Create output directory if it doesn't exist
if not os.path.isdir(out_dir):
  print "Making output directory {0}".format(out_dir)
  os.system( "mkdir %s" % out_dir )

## Set batch mode
ROOT.gROOT.SetBatch()

## This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

## GENIE and NuWro paths and files
in_file = ROOT.TFile(in_file_path, "read")

## Retrieve GENIE and NuWro reco signal, background, and fake data histograms. GENIE "fakedata" is simply the sum of its reco signal and background histograms
mHist_GENIE_2g1p_exclusive_signal = in_file.Get("effNum_reco_2g1p_exclusive")
mHist_GENIE_2g0p_exclusive_signal = in_file.Get("effNum_reco_2g0p_exclusive")
mHist_GENIE_2g1p_exclusive_background = in_file.Get("background_2g1p_exclusive")
mHist_GENIE_2g0p_exclusive_background = in_file.Get("background_2g0p_exclusive")
mHist_GENIE_2g1p_exclusive_fakedata = mHist_GENIE_2g1p_exclusive_signal.Clone("fakedata_2g1p_exclusive")
mHist_GENIE_2g1p_exclusive_fakedata.Add(mHist_GENIE_2g1p_exclusive_background)
mHist_GENIE_2g0p_exclusive_fakedata = mHist_GENIE_2g0p_exclusive_signal.Clone("fakedata_2g0p_exclusive")
mHist_GENIE_2g0p_exclusive_fakedata.Add(mHist_GENIE_2g0p_exclusive_background)
mHist_GENIE_2g1p_exclusive_folded = mHist_GENIE_2g1p_exclusive_signal.Clone("folded_2g1p_exclusive")
mHist_GENIE_2g0p_exclusive_folded = mHist_GENIE_2g0p_exclusive_signal.Clone("folded_2g0p_exclusive")
GENIE_2g1p_exclusive_signal = mHist_GENIE_2g1p_exclusive_signal.GetCVHistoWithError()
GENIE_2g0p_exclusive_signal = mHist_GENIE_2g0p_exclusive_signal.GetCVHistoWithError()
GENIE_2g1p_exclusive_background = mHist_GENIE_2g1p_exclusive_background.GetCVHistoWithError()
GENIE_2g0p_exclusive_background = mHist_GENIE_2g0p_exclusive_background.GetCVHistoWithError()
GENIE_2g1p_exclusive_fakedata = mHist_GENIE_2g1p_exclusive_fakedata.GetCVHistoWithError()
GENIE_2g0p_exclusive_fakedata = mHist_GENIE_2g0p_exclusive_fakedata.GetCVHistoWithError()
GENIE_2g1p_exclusive_folded = mHist_GENIE_2g1p_exclusive_folded.GetCVHistoWithError()
GENIE_2g0p_exclusive_folded = mHist_GENIE_2g0p_exclusive_folded.GetCVHistoWithError()
NuWro_2g1p_exclusive_signal = in_file.Get("nu_uBooNE_breakdown_2g1psig")
NuWro_2g0p_exclusive_signal = in_file.Get("nu_uBooNE_breakdown_2g0psig")
NuWro_2g1p_exclusive_background = in_file.Get("nu_uBooNE_breakdown_2g1pbkg")
NuWro_2g0p_exclusive_background = in_file.Get("nu_uBooNE_breakdown_2g0pbkg")
NuWro_2g1p_exclusive_fakedata = in_file.Get("nu_uBooNE_fakedata_2g1p")
NuWro_2g0p_exclusive_fakedata = in_file.Get("nu_uBooNE_fakedata_2g0p")

## Scale NuWro histograms to GENIE POT
#for np in ["2g1p", "2g0p"]:
#    for datatype in ["signal", "background", "fakedata"]:
#        exec("NuWro_{0}_exclusive_{1}.Scale(6.868/3.0041393)".format(np, datatype))

mHist_NuWro_2g1p_exclusive_folded = in_file.Get("evtRate_2g1p_exclusive")
mHist_NuWro_2g0p_exclusive_folded = in_file.Get("evtRate_2g0p_exclusive")
NuWro_2g1p_exclusive_folded = mHist_NuWro_2g1p_exclusive_folded.GetCVHistoWithError()
NuWro_2g0p_exclusive_folded = mHist_NuWro_2g0p_exclusive_folded.GetCVHistoWithError()

## Loop over all plots to be made
for np in ["2g1p", "2g0p"]:
    for datatype in ["signal", "background", "fakedata", "folded"]:
        with makeEnv_TCanvas(out_dir + "/{0}_exclusive_{1}.png".format(np, datatype)) as canvas:

            ## Clone GENIE histograms with different variable names to the same variable within the for loop for conciseness of code
            exec("local_GENIE_exclusive = GENIE_{0}_exclusive_{1}.Clone(\"local_Sys{0}_numerator_reco_{1}\")".format(np, datatype))

            ## Convert datatype to words used in plot titles 
            if datatype == "signal":
                DT = "Signal"
            elif datatype == "background":
                DT = "Background"
            elif datatype == "fakedata":
                DT = "Selected Data"
            elif datatype == "folded":
                DT = "NuWro Selected Data - GENIE Background"

            ## Set GENIE histogram Draw properties
            local_GENIE_exclusive.SetLineColor(ROOT.kRed)
            local_GENIE_exclusive.SetMarkerColor(ROOT.kRed)
            local_GENIE_exclusive.SetFillColor(ROOT.kRed)
            local_GENIE_exclusive.SetFillStyle(0)
            local_GENIE_exclusive.SetLineWidth(2)
            local_GENIE_exclusive.SetTitle("{0} Exclusive {1}".format(np, DT))
            local_GENIE_exclusive.GetXaxis().SetTitle("Reco #pi^{0} momentum (GeV)")
            local_GENIE_exclusive.GetYaxis().SetTitle("# Reco {0} Events".format(DT))
            local_GENIE_exclusive.SetFillStyle(0)
            ROOT.gStyle.SetOptStat(0)

            ## Get content and error of GENIE bin with max y value
            GENIE_maxybin = local_GENIE_exclusive.GetMaximumBin()
            GENIE_maxy = local_GENIE_exclusive.GetBinContent(GENIE_maxybin)
            GENIE_maxyerror = local_GENIE_exclusive.GetBinError(GENIE_maxybin)

            ## Clone NuWro histograms with different variable names to the same variable within the for loop for conciseness of code
            exec("local_NuWro_exclusive = NuWro_{0}_exclusive_{1}.Clone(\"local_nu_uBooNE_reco_{1}_{0}\")".format(np, datatype))

            ## Set NuWro histogram Draw properties. Plotted on same canvas as GENIE histogram, so doesn't need own titles.
            local_NuWro_exclusive.SetLineWidth(2)
            local_NuWro_exclusive.SetMarkerSize(2)
            local_NuWro_exclusive.SetMarkerStyle(20)

            ## Get content and error of NuWro bin with max y value
            NuWro_maxybin = local_NuWro_exclusive.GetMaximumBin()
            NuWro_maxy = local_NuWro_exclusive.GetBinContent(NuWro_maxybin)
            NuWro_maxyerror = local_NuWro_exclusive.GetBinError(NuWro_maxybin)

            ## Choose maximum of GENIE or NuWro maxy + maxyerror and multiply by 1.1 for maximum y axis value of combined histogram. Set minimum to 0.
            actual_maxy = max(GENIE_maxy + GENIE_maxyerror, NuWro_maxy + NuWro_maxyerror)*1.1
            local_GENIE_exclusive.SetMaximum(actual_maxy)
            local_GENIE_exclusive.SetMinimum(0)

            ## Draw GENIE histogram, GENIE error, and NuWro histogram with error. GENIE histogram drawn first drawn without errors; then error is drawn with hashed fill style.
            ## DrawCopy() prevents first drawn local_GENIE_exclusive histogram from changing its fill style once drawn.
            local_GENIE_exclusive.DrawCopy("HIST")
            local_GENIE_exclusive.SetFillStyle(3545)
            local_GENIE_exclusive.Draw("E2 SAME")
            local_NuWro_exclusive.Draw("E1P SAME")

            ## Draw borderless legend.
            legend = ROOT.TLegend(0.555, 0.5, 0.90, 0.7)
            legend.SetBorderSize(0)
            legend.AddEntry(local_GENIE_exclusive, "GENIE")
            legend.AddEntry(local_NuWro_exclusive, "NuWro Fake Data")
            legend.Draw()

            ## Calculate and draw Chi^2
            ## Statistical uncertainties probably incorrect. (!!!)
            #exec("local_covhist_exclusive = in_file.Get(\"cov_evtRate_{0}_exclusive\")".format(np))
            exec("local_covmat_exclusive = mHist_GENIE_{0}_exclusive_{1}.GetTotalErrorMatrix()".format(np, datatype))
            if datatype == "folded":
                exec("local_covmat_exclusive += mHist_NuWro_{0}_exclusive_folded.GetTotalErrorMatrix(False)".format(np))
            nbins = local_GENIE_exclusive.GetNbinsX()
            local_covmat_exclusive = local_covmat_exclusive.GetSub(1, nbins, 1, nbins)
            #local_covmat_exclusive = ROOT.TMatrixD(nbins, nbins)
            #for row in range(0, nbins):
            #    for col in range(0, nbins):
            #        local_covmat_exclusive[row][col] = local_covhist_exclusive.GetBinContent(row + 1, col + 1)
            chi2 = calculateChi2(local_GENIE_exclusive, local_NuWro_exclusive, local_covmat_exclusive, False)
            #pt = ROOT.TPaveText()
            pt = ROOT.TPaveText(0.555, 0.75, 0.90, 0.90, "NDC")
            pt.SetBorderSize(0)
            pt.SetFillColorAlpha(0, 0)
            pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2, local_GENIE_exclusive.GetNbinsX()))
            pt.Draw()
