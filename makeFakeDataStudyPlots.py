### Script to take make plots, etc.
### using the MINERvA Analysis Toolkit

import ROOT
import datetime as dt
import argparse
import os,math
from customHistAndPlotMethods import makeEnv_TCanvas,localDrawErrorSummary

## Set ROOT to batch mode
ROOT.gROOT.SetBatch()

## This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

## Setup MnvPlotter, which has all of the plotting utilities
plotter = ROOT.PlotUtils.MnvPlotter()

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################
parser = argparse.ArgumentParser(description='Script to make cross-section plots using the MINERvA Analysis Toolkit')
parser.add_argument('in_file', help='Path to input file', type=str,nargs='?')
parser.add_argument('out_dir', help='Path to output directory. Defaults to input directory', type=str,nargs='?')
p = parser.parse_args()

## If in_file is not provided, exit
if p.in_file < 0:
  print "ERROR: Input directory argument not provided"
  parser.print_help()
  exit(1)

histFilePath = p.in_file
histFile = ROOT.TFile(histFilePath)

## If out_dir is not provided, default to using in_dir
if p.out_dir < 0:
  outFileDir = os.path.dirname(histFilePath)
  plotDir = outFileDir+"/{0}_xsec-plots".format(dt.date.today())
else:
  plotDir = p.out_dir

## Create output directory if it doesn't exist
if not os.path.isdir(plotDir):
  print "Making plot directory {0}".format(plotDir)
  os.system( "mkdir %s" % plotDir )

#############################################################################################################
### Pull out relevant objects from input file ###############################################################
#############################################################################################################

sigDef = "2g1p_exclusive"

### Get THnDs
################
for histCat in ["cov_evtRate","unfolded_cov_evtRate","unfolded_evtRate","unfolded_xSection","add_smear_matrix"]:

  exec("tHist_{0}_{1} = histFile.Get(\"{0}_{1}\")".format(histCat,sigDef))

### Get MnvHnDs
################
for histCat in ["xSection_mc"]:

  exec("mHist_{0}_{1} = histFile.Get(\"{0}_{1}\")".format(histCat,sigDef))
  exec("tHist_{0}_{1} = mHist_{0}_{1}.GetCVHistoWithStatError()".format(histCat,sigDef))

### Construct correlation matrices
##################################
tHist_corr_evtRate_2g1p_exclusive = tHist_cov_evtRate_2g1p_exclusive.Clone("tHist_corr_evtRate_2g1p_exclusive")
tHist_unfolded_corr_evtRate_2g1p_exclusive = tHist_unfolded_cov_evtRate_2g1p_exclusive.Clone("tHist_unfolded_corr_evtRate_2g1p_exclusive")

nBins = tHist_unfolded_evtRate_2g1p_exclusive.GetNbinsX()

for histCat in ["","unfolded_"]:

  # Loop over each bin in the covariance matrix
  for i in range(nBins):
      for j in range(nBins):
          exec("cov_ij = tHist_{0}cov_evtRate_{1}.GetBinContent(i+1, j+1)".format(histCat,sigDef))
          exec("cov_ii = tHist_{0}cov_evtRate_{1}.GetBinContent(i+1, i+1)".format(histCat,sigDef))
          exec("cov_jj = tHist_{0}cov_evtRate_{1}.GetBinContent(j+1, j+1)".format(histCat,sigDef))
          
          # Calculate the correlation coefficient
          if cov_ii > 0 and cov_jj > 0:
              corr_ij = cov_ij / (cov_ii**0.5 * cov_jj**0.5)
          else:
              corr_ij = 0.0
          
          # Set the correlation coefficient in the correlation matrix
          exec("tHist_{0}corr_evtRate_{1}.SetBinContent(i+1, j+1, corr_ij)".format(histCat,sigDef))

#############################################################################################################
### Plot response and event rate covariance and correlation #################################################
#############################################################################################################
plotter = ROOT.PlotUtils.MnvPlotter()
plotter.SetROOT6Palette(54)
ROOT.gStyle.SetNumberContours(200)

with makeEnv_TCanvas('{0}/response_2g1p_exclusive.png'.format(plotDir)) as canvas:
  
  local_tHist_Response = histFile.Get('Response')
  local_tHist_Response.Draw("colz")
  local_tHist_Response.GetXaxis().SetTitle("True #pi^{0} momentum")
  local_tHist_Response.GetYaxis().SetTitle("Reconstructed #pi^{0} momentum")
  canvas.canvas.SetLogz()

with makeEnv_TCanvas('{0}/add_smear_matrix_2g1p_exclusive.png'.format(plotDir)) as canvas:
  
  local_tHist_add_smear_matrix = histFile.Get('add_smear_matrix_2g1p_exclusive')
  local_tHist_add_smear_matrix.Draw("colz")
  local_tHist_add_smear_matrix.GetXaxis().SetTitle("True #pi^{0} momentum")
  local_tHist_add_smear_matrix.GetYaxis().SetTitle("Reconstructed #pi^{0} momentum")
  #canvas.canvas.SetLogz()

with makeEnv_TCanvas('{0}/cov_evtRate_2g1p_exclusive.png'.format(plotDir)) as canvas:

  local_tHist_cov_evtRate_2g1p_exclusive = tHist_cov_evtRate_2g1p_exclusive.Clone("local_tHist_cov_evtRate_2g1p_exclusive")
  local_tHist_cov_evtRate_2g1p_exclusive.Draw("colz")
  local_tHist_cov_evtRate_2g1p_exclusive.GetXaxis().SetTitle("Reconstructed #pi^{0} momentum")
  local_tHist_cov_evtRate_2g1p_exclusive.GetYaxis().SetTitle("Reconstructed #pi^{0} momentum")
  canvas.canvas.SetLogz()

with makeEnv_TCanvas('{0}/unfolded_cov_evtRate_2g1p_exclusive.png'.format(plotDir)) as canvas:

  local_tHist_unfolded_cov_evtRate_2g1p_exclusive = tHist_unfolded_cov_evtRate_2g1p_exclusive.Clone("local_tHist_unfolded_cov_evtRate_2g1p_exclusive")
  local_tHist_unfolded_cov_evtRate_2g1p_exclusive.Draw("colz")
  local_tHist_unfolded_cov_evtRate_2g1p_exclusive.GetXaxis().SetTitle("True #pi^{0} momentum")
  local_tHist_unfolded_cov_evtRate_2g1p_exclusive.GetYaxis().SetTitle("True #pi^{0} momentum")
  canvas.canvas.SetLogz()

## Change to different color palette for correlation matrix
plotter.SetROOT6Palette(87)
ROOT.gStyle.SetNumberContours(200)

with makeEnv_TCanvas('{0}/corr_evtRate_2g1p_exclusive.png'.format(plotDir)) as canvas:

  local_tHist_corr_evtRate_2g1p_exclusive = tHist_corr_evtRate_2g1p_exclusive.Clone("local_tHist_corr_evtRate_2g1p_exclusive")
  local_tHist_corr_evtRate_2g1p_exclusive.GetZaxis().SetRangeUser(-1,1)
  local_tHist_corr_evtRate_2g1p_exclusive.Draw("colz")
  local_tHist_corr_evtRate_2g1p_exclusive.GetXaxis().SetTitle("Reconstructed #pi^{0} momentum")
  local_tHist_corr_evtRate_2g1p_exclusive.GetYaxis().SetTitle("Reconstructed #pi^{0} momentum")

with makeEnv_TCanvas('{0}/unfolded_corr_evtRate_2g1p_exclusive.png'.format(plotDir)) as canvas:

  local_tHist_unfolded_corr_evtRate_2g1p_exclusive = tHist_unfolded_corr_evtRate_2g1p_exclusive.Clone("local_tHist_unfolded_corr_evtRate_2g1p_exclusive")
  local_tHist_unfolded_corr_evtRate_2g1p_exclusive.GetZaxis().SetRangeUser(-1,1)
  local_tHist_unfolded_corr_evtRate_2g1p_exclusive.Draw("colz")
  local_tHist_unfolded_corr_evtRate_2g1p_exclusive.GetXaxis().SetTitle("True #pi^{0} momentum")
  local_tHist_unfolded_corr_evtRate_2g1p_exclusive.GetYaxis().SetTitle("True #pi^{0} momentum")

#############################################################################################################
### Make fake data study comparison plots ###################################################################
#############################################################################################################
with makeEnv_TCanvas('{0}/fakedatavsgenie_xSection.png'.format(plotDir)):

  local_tHist_unfolded_xSection_2g1p_exclusive = tHist_unfolded_xSection_2g1p_exclusive.Clone("local_tHist_unfolded_xSection_2g1p_exclusive")
  local_tHist_xSection_mc_2g1p_exclusive = tHist_xSection_mc_2g1p_exclusive.Clone("local_tHist_xSection_mc_2g1p_exclusive")

  local_tHist_unfolded_xSection_2g1p_exclusive.SetMarkerSize(0.5)
  local_tHist_unfolded_xSection_2g1p_exclusive.Draw()

  local_tHist_xSection_mc_2g1p_exclusive.GetYaxis().SetTitleSize(0.05)
  local_tHist_xSection_mc_2g1p_exclusive.GetXaxis().SetTitleSize(0.05)
  local_tHist_xSection_mc_2g1p_exclusive.GetYaxis().SetTitle("#sigma_{NC 1 #pi^{0}}[10^{-38} cm^{2}/Atom/GeV]")
  local_tHist_xSection_mc_2g1p_exclusive.GetXaxis().SetTitle("#pi^{0} momentum")

  local_tHist_xSection_mc_2g1p_exclusive.SetLineColor(ROOT.kRed)
  local_tHist_xSection_mc_2g1p_exclusive.SetMarkerColor(ROOT.kRed)
  local_tHist_xSection_mc_2g1p_exclusive.SetMarkerSize(0.5)
  local_tHist_xSection_mc_2g1p_exclusive.Draw("SAME")

  legend = ROOT.TLegend(0.55,0.7,0.8,0.85, "")
  legend.SetBorderSize(0);
  legend.AddEntry(local_tHist_unfolded_xSection_2g1p_exclusive,"NuWro Fake Data","lep")
  legend.AddEntry(local_tHist_xSection_mc_2g1p_exclusive,"GENIE Prediction","lep")
  legend.Draw()

## Import NuWro Truth
NuWroTruthFile = ROOT.TFile("/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/MajorMerge_GGE_mark/working_dir/ToTH1D/NuWro_FakeData_Generation_NEW/NuWro_Apr2023_v2_CV.SBNspec.root")

with makeEnv_TCanvas('{0}/fakedatavsgenie_evtRate.png'.format(plotDir)):

  local_tHist_unfolded_evtRate_2g1p_exclusive = tHist_unfolded_evtRate_2g1p_exclusive.Clone("local_tHist_unfolded_evtRate_2g1p_exclusive")
  local_tHist_genie_evtRate_smeared_2g1p_exclusive = histFile.Get("smeared_true_signal_2g1p_exclusive")
  local_tHist_genie_evtRate_2g1p_exclusive = histFile.Get("prior_true_signal_2g1p_exclusive")
  local_tHist_nuwro_truth_2g1p_exclusive = NuWroTruthFile.Get("nu_uBooNE_denom_2g1p")
  local_tHist_nuwro_truth_2g1p_exclusive.Scale(6.7873/3.0041393)

  for i in range(1,nBins+1):## Loop over bins
    cov_binContent = tHist_unfolded_cov_evtRate_2g1p_exclusive.GetBinContent(i,i)
    local_tHist_unfolded_evtRate_2g1p_exclusive.SetBinError(i,math.sqrt(cov_binContent))

  local_tHist_unfolded_evtRate_2g1p_exclusive.SetMarkerSize(0.5)
  local_tHist_unfolded_evtRate_2g1p_exclusive.GetYaxis().SetRangeUser(0,1500)
  local_tHist_unfolded_evtRate_2g1p_exclusive.Draw()

  local_tHist_genie_evtRate_2g1p_exclusive.SetMarkerColor(ROOT.kGreen+2)
  local_tHist_genie_evtRate_2g1p_exclusive.SetLineColor(ROOT.kGreen+2)
  local_tHist_genie_evtRate_2g1p_exclusive.Draw("SAME")
  
  local_tHist_nuwro_truth_2g1p_exclusive.SetMarkerColor(ROOT.kViolet)
  local_tHist_nuwro_truth_2g1p_exclusive.SetLineColor(ROOT.kViolet)
  local_tHist_nuwro_truth_2g1p_exclusive.Draw("SAME")

  local_tHist_genie_evtRate_smeared_2g1p_exclusive.GetYaxis().SetTitleSize(0.05)
  local_tHist_genie_evtRate_smeared_2g1p_exclusive.GetXaxis().SetTitleSize(0.05)
  local_tHist_genie_evtRate_smeared_2g1p_exclusive.GetYaxis().SetTitle("\# Events")
  local_tHist_genie_evtRate_smeared_2g1p_exclusive.GetXaxis().SetTitle("#pi^{0} momentum")

  local_tHist_genie_evtRate_smeared_2g1p_exclusive.SetLineColor(ROOT.kRed)
  local_tHist_genie_evtRate_smeared_2g1p_exclusive.SetMarkerColor(ROOT.kRed)
  local_tHist_genie_evtRate_smeared_2g1p_exclusive.SetMarkerSize(0.5)
  local_tHist_genie_evtRate_smeared_2g1p_exclusive.Draw("SAME")

  legend = ROOT.TLegend(0.5,0.7,0.845,0.9, "")
  legend.SetBorderSize(0);
  legend.AddEntry(local_tHist_unfolded_evtRate_2g1p_exclusive,"NuWro Fake Data","lep")
  legend.AddEntry(local_tHist_nuwro_truth_2g1p_exclusive,"NuWro Truth","lep")
  legend.AddEntry(local_tHist_genie_evtRate_2g1p_exclusive,"GENIE Prediction","lep")
  legend.AddEntry(local_tHist_genie_evtRate_smeared_2g1p_exclusive,"GENIE Prediction, smeared","lep")
  legend.Draw()
