### Script to take make plots, etc.
### using the MINERvA Analysis Toolkit

import ROOT
import datetime as dt
import argparse
import os,math
from customHistAndPlotMethods import *
from calculateChi2 import *
from errorMaps import *

## Set ROOT to batch mode
ROOT.gROOT.SetBatch()

## This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

## Setup MnvPlotter, which has all of the plotting utilities
plotter = ROOT.PlotUtils.MnvPlotter()
plotter.error_summary_group_map.clear()
for group in error_bands:
  for error in error_bands[group]:
    plotter.error_summary_group_map[group].push_back(error)
plotter.stat_error_name = "Data Statistical"
ROOT.gStyle.SetNumberContours(200)

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################
parser = argparse.ArgumentParser(description='Script to make cross-section plots using the MINERvA Analysis Toolkit')
parser.add_argument('in_file', help='Path to input file', type=str,nargs='?')
parser.add_argument('out_dir', help='Path to output directory. Defaults to input directory', type=str,nargs='?')
parser.add_argument('--closureTest',help='Input file corresponds to closure test',action='store_true')
parser.add_argument('--test', help = 'Produce plots in test mode; useful for debugging', action = 'store_true')
p = parser.parse_args()

## If in_file is not provided, exit
if not p.in_file:
  print("ERROR: Input directory argument not provided")
  parser.print_help()
  exit(1)

histFilePath = p.in_file
histFile = ROOT.TFile(histFilePath)

## If out_dir is not provided, default to using in_dir
if not p.out_dir:
  outFileDir = os.path.dirname(histFilePath)
  plotDir = outFileDir+"/{0}_xsec-plots".format(dt.date.today())
else:
  plotDir = p.out_dir

## Create output directory if it doesn't exist
if not os.path.isdir(plotDir):
  print("Making plot directory {0}".format(plotDir))
  os.system( "mkdir %s" % plotDir )

#############################################################################################################
### Set Configurations ######################################################################################
#############################################################################################################

## Prescription is slightly different for closure test input 
is_closure_test = True if p.closureTest>0 else False

## Running in the test configuration adds unsmeared distributions to the plots and chi2 values to closure tests.
is_test = True if p.test > 0 else False

## Check whether this run is for fake data.
is_fake_data = histFile.Get("is_fake_data")
data_string = "NuWro Fake Data" if is_fake_data else "Data"

#############################################################################################################
### Pull out relevant objects from input file ###############################################################
#############################################################################################################

#for sigDef in ["2g1p","2g0p","2gXp"]:
for sigDefnp in ["2g1p","2g0p","2gXp"]:
  #for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["exclusive","inclusive"]:

    sigDef = sigDefnp + "_" + sigDefexcl

    if sigDefnp == "2gXp" and sigDefexcl == "exclusive":
      continue

    if sigDefnp == "2g1p" and sigDefexcl == "inclusive":
      continue

    if sigDefnp == "2g0p" and sigDefexcl == "inclusive":
      continue

    ### Get THnDs
    ################
    for histCat in ["cov_evtRate","unfolded_cov_evtRate","unfolded_xSection","add_smear_matrix", "smeared_xSection_mc", "smeared_nuwro_signal"]:
    
      exec("tHist_{0}_{1} = histFile.Get(\"{0}_{1}\")".format(histCat,sigDef))
     
    ### Get MnvHnDs
    ################
    for histCat in ["unfolded_evtRate","xSection_mc", "background"]:
    
      exec("mHist_{0}_{1} = histFile.Get(\"{0}_{1}\")".format(histCat,sigDef))
      exec("tHist_{0}_{1} = mHist_{0}_{1}.GetCVHistoWithStatError()".format(histCat,sigDef))

    ### Construct correlation matrices
    ##################################
    exec("tHist_corr_evtRate_{0} = tHist_cov_evtRate_{0}.Clone(\"tHist_corr_evtRate_{0}\")".format(sigDef))
    exec("tHist_unfolded_corr_evtRate_{0} = tHist_unfolded_cov_evtRate_{0}.Clone(\"tHist_unfolded_corr_evtRate_{0}\")".format(sigDef))
    
    exec("nBins = tHist_unfolded_evtRate_{0}.GetNbinsX()".format(sigDef))
    
    for histCat in ["","unfolded_"]:
    
      # Loop over each bin in the covariance matrix
      for i in range(nBins+2):
          for j in range(nBins+2):
              exec("cov_ij = tHist_{0}cov_evtRate_{1}.GetBinContent(i, j)".format(histCat,sigDef))
              exec("cov_ii = tHist_{0}cov_evtRate_{1}.GetBinContent(i, i)".format(histCat,sigDef))
              exec("cov_jj = tHist_{0}cov_evtRate_{1}.GetBinContent(j, j)".format(histCat,sigDef))
              
              # Calculate the correlation coefficient
              if cov_ii > 0 and cov_jj > 0:
                  corr_ij = cov_ij / (cov_ii**0.5 * cov_jj**0.5)
              else:
                  corr_ij = 0.0
              
              # Set the correlation coefficient in the correlation matrix
              exec("tHist_{0}corr_evtRate_{1}.SetBinContent(i, j, corr_ij)".format(histCat,sigDef))

### Get x-axis titles
######################
reco_xlabel = mHist_background_2g1p_exclusive.GetXaxis().GetTitle()
true_xlabel = mHist_xSection_mc_2g1p_exclusive.GetXaxis().GetTitle()

mHist_POT_scaling = histFile.Get("POT_scaling")
tHist_POT_scaling = mHist_POT_scaling.GetCVHistoWithStatError()
  
#############################################################################################################
### Plot response and event rate covariance and correlation #################################################
#############################################################################################################

# Text to include on all plots
ptall = ROOT.TPaveText(0.605, 0.15, 0.85, 0.25, "NDC")
ptall.SetBorderSize(0)
ptall.SetFillColorAlpha(0, 0)
ptall.AddText("MicroBooNE Simulation In Progress")

#for sigDef in ["2g1p","2g0p","2gnp"]:
for sigDefnp in ["2g1p","2g0p","2gXp"]:
  #for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["exclusive","inclusive"]:

    sigDef = sigDefnp + "_" + sigDefexcl

    if sigDefnp == "2gXp" and sigDefexcl == "exclusive":
      continue

    if sigDefnp == "2g1p" and sigDefexcl == "inclusive":
      continue

    if sigDefnp == "2g0p" and sigDefexcl == "inclusive":
      continue

    #overflow and underflow for unfolded_evtRate
    threshold = 10e-6
    exec("local_tHist_unfolded_xSection = tHist_unfolded_xSection_{0}".format(sigDef))
    nBins = local_tHist_unfolded_xSection.GetNbinsX()

    binVal_overflow_true = abs(local_tHist_unfolded_xSection.GetBinContent(nBins + 1))
    include_overflow_true = 1 if binVal_overflow_true > threshold else 0

    binVal_underflow_true = abs(local_tHist_unfolded_xSection.GetBinContent(0))
    include_underflow_true = 1 if binVal_underflow_true > threshold else 0

    #upper and lower bounds for later use
    lowerBound = 0 if binVal_underflow_true else 1
    upperBound = nBins + 1 if binVal_overflow_true else nBins

    plotter.SetROOT6Palette(54)
    
    with makeEnv_TCanvas('{0}/response_{1}.png'.format(plotDir,sigDef)) as canvas:
      
      local_tHist_Response = histFile.Get("Response_{0}".format(sigDefnp))
      ROOT.gStyle.SetOptTitle(1)
      local_tHist_Response.SetTitle(sigDefnp + " Exclusive Response Matrix")
      local_tHist_Response.GetXaxis().SetTitle(true_xlabel)
      local_tHist_Response.GetYaxis().SetTitle(reco_xlabel)
      #canvas.canvas.SetLogz()
      local_tHist_Response.Draw("colz")
      ptall.Draw()

    with makeEnv_TCanvas('{0}/efficiency_{1}.png'.format(plotDir, sigDef)) as canvas:
      local_mHist_eff = histFile.Get("eff_{0}".format(sigDef))
      local_covMat_eff = local_mHist_eff.GetTotalErrorMatrix(True)
      for i in range(nBins + 2):
        local_mHist_eff.SetBinError(i, math.sqrt(local_covMat_eff(i, i)))
      local_mHist_eff.SetTitle(sigDefnp + " Exclusive Efficiency")
      local_mHist_eff.GetXaxis().SetTitle(true_xlabel)
      local_mHist_eff.GetYaxis().SetTitle("Efficiency")
      local_mHist_eff.SetFillStyle(3545)
      local_mHist_eff.SetFillColor(local_mHist_eff.GetLineColor())
      local_mHist_eff.DrawCopy("E2")
      local_mHist_eff.SetFillColorAlpha(ROOT.kWhite, 0.00)
      local_mHist_eff.Draw("HIST SAME")
    
    with makeEnv_TCanvas('{0}/add_smear_matrix_{1}.png'.format(plotDir, sigDef)) as canvas:
      local_tHist_add_smear_matrix = histFile.Get("add_smear_matrix_{0}".format(sigDef))
      local_tHist_add_smear_matrix.Draw("colz")
      local_tHist_add_smear_matrix.SetTitle(sigDefnp + " Exclusive Additional Smearing Matrix")
      local_tHist_add_smear_matrix.GetXaxis().SetTitle(true_xlabel)
      local_tHist_add_smear_matrix.GetYaxis().SetTitle("Smeared " + true_xlabel)
      ptall.SetTextColor(ROOT.kWhite)
      ptall.Draw()
      #canvas.canvas.SetLogz()
    
    with makeEnv_TCanvas('{0}/cov_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_cov_evtRate_{0} = tHist_cov_evtRate_{0}.Clone(\"local_tHist_cov_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_cov_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_cov_evtRate_{0}.SetTitle(\"{1} Exclusive Folded Covariance Matrix\")".format(sigDef, sigDefnp))
      exec("local_tHist_cov_evtRate_{0}.GetXaxis().SetTitle(reco_xlabel)".format(sigDef))
      exec("local_tHist_cov_evtRate_{0}.GetYaxis().SetTitle(reco_xlabel)".format(sigDef))
      ptall.Draw()
      canvas.canvas.SetLogz()
    
    with makeEnv_TCanvas('{0}/unfolded_cov_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_unfolded_cov_evtRate_{0} = tHist_unfolded_cov_evtRate_{0}.Clone(\"local_tHist_unfolded_cov_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_unfolded_cov_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_unfolded_cov_evtRate_{0}.SetTitle(\"{1} Exclusive Unfolded Covariance Matrix\")".format(sigDef, sigDefnp))
      exec("local_tHist_unfolded_cov_evtRate_{0}.GetXaxis().SetTitle(true_xlabel)".format(sigDef))
      exec("local_tHist_unfolded_cov_evtRate_{0}.GetYaxis().SetTitle(true_xlabel)".format(sigDef))
      ptall.Draw()
      #canvas.canvas.SetLogz()
    
    ## Change to different color palette for correlation matrix. Must set to 0 before setting to 87 or it will stay at 54 after the first iteration for some reason.
    plotter.SetROOT6Palette(0)
    plotter.SetROOT6Palette(87)

    with makeEnv_TCanvas('{0}/corr_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_corr_evtRate_{0} = tHist_corr_evtRate_{0}.Clone(\"local_tHist_corr_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_corr_evtRate_{0}.GetZaxis().SetRangeUser(-1,1)".format(sigDef))
      exec("local_tHist_corr_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_corr_evtRate_{0}.SetTitle(\"{1} Exclusive Folded Correlation Matrix\")".format(sigDef, sigDefnp))
      exec("local_tHist_corr_evtRate_{0}.GetXaxis().SetTitle(reco_xlabel)".format(sigDef))
      exec("local_tHist_corr_evtRate_{0}.GetYaxis().SetTitle(reco_xlabel)".format(sigDef))
      ptall.SetTextColor(ROOT.kBlack)
      ptall.Draw()
    
    with makeEnv_TCanvas('{0}/unfolded_corr_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_unfolded_corr_evtRate_{0} = tHist_unfolded_corr_evtRate_{0}.Clone(\"local_tHist_unfolded_corr_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.GetZaxis().SetRangeUser(-1,1)".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.SetTitle(\"{1} Unfolded Correlation Matrix\")".format(sigDef, sigDefnp))
      exec("local_tHist_unfolded_corr_evtRate_{0}.GetXaxis().SetTitle(true_xlabel)".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.GetYaxis().SetTitle(true_xlabel)".format(sigDef))
      ptall.Draw()
    
#############################################################################################################
### Make fake data study comparison plots ###################################################################
#############################################################################################################

# Move TPaveText for TH1D plots.
ptall.SetX1NDC(0.65)
ptall.SetY1NDC(0.4)
ptall.SetX2NDC(0.925)
ptall.SetY2NDC(0.5)

#for sigDef in ["2g1p","2g0p","2gnp"]:
for sigDefnp in ["2g1p","2g0p","2gXp"]:
  #for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["exclusive","inclusive"]:

    sigDef = sigDefnp + "_" + sigDefexcl

    if sigDefnp == "2gXp" and sigDefexcl == "exclusive":
      continue

    if sigDefnp == "2g1p" and sigDefexcl == "inclusive":
      continue

    if sigDefnp == "2g0p" and sigDefexcl == "inclusive":
      continue

    #overflow and underflow for unfolded_evtRate
    threshold = 10e-6
    exec("local_tHist_unfolded_xSection = tHist_unfolded_xSection_{0}".format(sigDef))
    nBins = local_tHist_unfolded_xSection.GetNbinsX()

    binVal_overflow_true = abs(local_tHist_unfolded_xSection.GetBinContent(nBins + 1))
    include_overflow_true = 1 if binVal_overflow_true > threshold else 0
    include_overflow_reco = include_overflow_true # Change if this changes, although it seems unlikely.

    binVal_underflow_true = abs(local_tHist_unfolded_xSection.GetBinContent(0))
    include_underflow_true = 1 if binVal_underflow_true > threshold else 0
    include_underflow_reco = include_underflow_true

    #upper and lower bounds for later use
    lowerBound = 0 if binVal_underflow_true else 1
    upperBound = nBins + 1 if binVal_overflow_true else nBins

    ### Plot of NuWro fake data unfolded events vs. GENIE generator prediction
    ##########################################################################
    with makeEnv_TCanvas('{0}/fakedatavsgenie_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      ## Create local copies of all needed hists
      exec("local_tHist_unfolded_evtRate = tHist_unfolded_evtRate_{0}.Clone(\"local_tHist_unfolded_evtRate_{0}\")".format(sigDef))
      local_tHist_genie_evtRate_smeared = histFile.Get("smeared_true_signal_{0}".format(sigDef))
      local_tHist_genie_evtRate = histFile.Get("prior_true_signal_{0}".format(sigDef))
      local_mHist_effDenom = histFile.Get("effDenom_{0}".format(sigDef))
      local_tHist_effDenom =local_mHist_effDenom.GetCVHistoWithError()
      exec("local_tHist_nuwro_truth = histFile.Get(\"nu_uBooNE_denom_{0}\")".format(sigDefnp))
      exec("local_tHist_smeared_nuwro_truth = tHist_smeared_nuwro_signal_{0}.Clone(\"smeared_nuwro_truth\")".format(sigDef))
    
      for i in range(nBins+2):## Loop over bins
        exec("cov_binContent = tHist_unfolded_cov_evtRate_{0}.GetBinContent(i,i)".format(sigDef))
        local_tHist_unfolded_evtRate.SetBinError(i,math.sqrt(cov_binContent))
        if is_closure_test:
          local_tHist_unfolded_evtRate.SetBinError(i, 0.00001)
      
      
      ## Calculate chi2 before scaling.
      local_tMat_unfolded_cov_evtRate = ROOT.TMatrixD(upperBound + include_underflow_true, upperBound + include_underflow_true)
      for i in range(lowerBound, upperBound+1):
        for j in range(lowerBound, upperBound+1):
          exec("local_tMat_unfolded_cov_evtRate[i-1][j-1] = tHist_unfolded_cov_evtRate_{0}.GetBinContent(j, i)".format(sigDef))
      chi2_NuWro = calculateChi2(local_tHist_nuwro_truth, local_tHist_unfolded_evtRate, local_tMat_unfolded_cov_evtRate, True)[0]
      chi2_NuWro_smeared = calculateChi2(local_tHist_smeared_nuwro_truth, local_tHist_unfolded_evtRate, local_tMat_unfolded_cov_evtRate, True)[0]
      chi2_GENIE = calculateChi2(local_tHist_genie_evtRate, local_tHist_unfolded_evtRate, local_tMat_unfolded_cov_evtRate, True)[0]
      chi2_GENIE_smeared = calculateChi2(local_tHist_genie_evtRate_smeared, local_tHist_unfolded_evtRate, local_tMat_unfolded_cov_evtRate, True)[0]
    
      ## Scale and bin-width-normalize all hists (and POT normalize NuWro truth)
      #local_tHist_unfolded_evtRate.Scale(1e-3)
      #local_tHist_genie_evtRate_smeared.Scale(1e-3)
      #local_tHist_genie_evtRate.Scale(1e-3)
      #local_tHist_effDenom.Scale(1e-3)
      #local_tHist_nuwro_truth.Scale(1e-3)
      #local_tHist_smeared_nuwro_truth.Scale(1e-3)
           
      ## Calculate yrange
      minybin_local_tHist_unfolded_evtRate = 0
      minycontent_local_tHist_unfolded_evtRate = 100000
      for i in range(1, nBins + 2):
        if local_tHist_unfolded_evtRate.GetBinContent(i) < minycontent_local_tHist_unfolded_evtRate:
          minybin_local_tHist_unfolded_evtRate = i
          minycontent_local_tHist_unfolded_evtRate = local_tHist_unfolded_evtRate.GetBinContent(i)
      miny_local_tHist_unfolded_evtRate = minycontent_local_tHist_unfolded_evtRate - (local_tHist_unfolded_evtRate.GetBinError(minybin_local_tHist_unfolded_evtRate) if minycontent_local_tHist_unfolded_evtRate < 0 else 0)
      maxy_local_tHist_unfolded_evtRate = local_tHist_unfolded_evtRate.GetMaximum() + local_tHist_unfolded_evtRate.GetBinError(local_tHist_unfolded_evtRate.GetMaximumBin())
      miny_local_tHist_genie_evtRate = local_tHist_genie_evtRate.GetMinimum()
      maxy_local_tHist_genie_evtRate = local_tHist_genie_evtRate.GetMaximum()
      miny_local_tHist_genie_evtRate_smeared = local_tHist_genie_evtRate_smeared.GetMinimum()
      maxy_local_tHist_genie_evtRate_smeared = local_tHist_genie_evtRate_smeared.GetMaximum()
      miny = min(0, miny_local_tHist_unfolded_evtRate, (miny_local_tHist_genie_evtRate if not is_closure_test else 0), miny_local_tHist_genie_evtRate_smeared)*1.1
      maxy = max(0, maxy_local_tHist_unfolded_evtRate, (maxy_local_tHist_genie_evtRate if not is_closure_test else 0), maxy_local_tHist_genie_evtRate_smeared)*1.1

      ## Set plot formatting 
      local_tHist_unfolded_evtRate.SetMarkerSize(0.5)
      local_tHist_unfolded_evtRate.SetLineColor(ROOT.kBlue+2)
      local_tHist_unfolded_evtRate.SetMarkerColor(ROOT.kBlack)
      local_tHist_unfolded_evtRate.SetMarkerStyle(ROOT.kFullCircle)
      local_tHist_unfolded_evtRate.SetFillStyle(0)
      local_tHist_unfolded_evtRate.SetLineWidth(1)

      local_tHist_effDenom.SetLineColor(ROOT.kGreen+2)
      local_tHist_effDenom.SetMarkerColor(ROOT.kGreen+2)
      local_tHist_effDenom.SetMarkerSize(0.5)
      #local_tHist_effDenom.Draw("SAME")

      local_tHist_genie_evtRate_smeared.SetLineColor(ROOT.kRed)
      local_tHist_genie_evtRate_smeared.SetMarkerColor(ROOT.kRed)
      local_tHist_genie_evtRate_smeared.SetFillColor(ROOT.kRed)
      local_tHist_genie_evtRate_smeared.SetFillStyle(3554)
      local_tHist_genie_evtRate_smeared.SetTitle(sigDefnp + " Exclusive Unfolded Events") 
      local_tHist_genie_evtRate_smeared.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_genie_evtRate_smeared.GetYaxis().SetTitleSize(0.05)
      local_tHist_genie_evtRate_smeared.GetXaxis().SetTitleSize(0.05)
      local_tHist_genie_evtRate_smeared.GetYaxis().SetTitle("Events")
      local_tHist_genie_evtRate_smeared.GetXaxis().SetTitle(true_xlabel)

      overflow1 = DrawWithOverflow(local_tHist_genie_evtRate_smeared, canvas.canvas, "HIST")
      if not is_test:
        canvas.canvas.cd(1) # Place chi^2 pavetext and legend relative to first pad.
        ptall.Draw()
        canvas.canvas.cd(0)

      # If it's a closure test, change the historgram title to reflect that.
      if is_closure_test:
        local_tHist_genie_evtRate_smeared.SetTitle(sigDefnp + " Exclusive Closure Test")
      
      # Print chi2 value of closure test using unfolded (fake) data covariance if in test mode.
      # If not in test mode, print only the smeared chi2 value regardless of whether or not it's a closure test.
      if (is_closure_test and is_test) or (not is_closure_test and not is_test):
        canvas.canvas.cd(1)
        pt = ROOT.TPaveText(0.53, 0.6, 0.775, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        if(is_closure_test):
          pt.AddText("#chi^{{2}}/DOF: {0:.5f}/{1}".format(chi2_GENIE_smeared, nBins + 1))
        else:
          pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2_GENIE_smeared, nBins + 1))
        pt.Draw()
        canvas.canvas.cd(0)
    
      # If configured to test mode, also draw and print chi2 of the unsmeared distribution for non-closure tests.
      elif is_test:
        local_tHist_genie_evtRate.SetMarkerColor(ROOT.kCyan - 3)
        local_tHist_genie_evtRate.SetLineColor(ROOT.kCyan - 3)
        local_tHist_genie_evtRate.SetFillColor(ROOT.kCyan - 3)
        local_tHist_genie_evtRate.SetFillStyle(3545)
        overflow2 = DrawWithOverflow(local_tHist_genie_evtRate, canvas.canvas, "SAME HIST")
        
        canvas.canvas.cd(1)
        pt = ROOT.TPaveText(0.53, 0.5, 0.875, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}_{{unsmeared}}/DOF: {0:.0f}/{1}".format(chi2_GENIE, nBins + 1))
        pt.AddText("#chi^{{2}}_{{smeared}}/DOF: {0:.2f}/{1}".format(chi2_GENIE_smeared, nBins + 1))
        pt.Draw()
        canvas.canvas.cd(0)

      overflow3 = DrawWithOverflow(local_tHist_unfolded_evtRate, canvas.canvas, "SAME")

      canvas.canvas.cd(1)
      # Draw legend with appropriate labels.
      legend = ROOT.TLegend(0.52,0.7,0.945,0.9, "")
      legend.SetBorderSize(0)
      if is_closure_test:
        legend.AddEntry(local_tHist_unfolded_evtRate,"GENIE Reconstruction Unfolded","lep")
        legend.AddEntry(local_tHist_genie_evtRate_smeared,"GENIE Truth, Smeared","f")
      else:
        legend.AddEntry(local_tHist_unfolded_evtRate,"Unfolded " + data_string,"lep")
        #legend.AddEntry(local_tHist_nuwro_truth,"NuWro Truth","f")
        #legend.AddEntry(local_tHist_effDenom,"effDenom","lep")
        if is_test:
          legend.AddEntry(local_tHist_genie_evtRate,"GENIE Truth","f")
        legend.AddEntry(local_tHist_genie_evtRate_smeared,"GENIE Truth, Smeared","f")
      legend.Draw()
      canvas.canvas.cd(0)

    ### Plot of NuWro fake data cross section vs. GENIE generator prediction
    ########################################################################
    with makeEnv_TCanvas('{0}/fakedatavsgenie_xSection_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      # Cross sections are on the order of 1e-38. Scale by 1e38 and move the exponent to the units label.
      exec("local_tHist_unfolded_xSection_scaled = tHist_unfolded_xSection_{0}.Clone(\"local_tHist_unfolded_xSection_scaled_{0}\")".format(sigDef))
      local_tHist_unfolded_xSection_scaled.Scale(1e40) 
      exec("local_tHist_xSection_mc_scaled = tHist_xSection_mc_{0}.Clone(\"local_tHist_xSection_mc_scaled_{0}\")".format(sigDef))
      local_tHist_xSection_mc_scaled.Scale(1e40)
      exec("local_tHist_smeared_xSection_mc_scaled = tHist_smeared_xSection_mc_{0}.Clone(\"local_tHist_smeared_xSection_mc_scaled_{0}\")".format(sigDef))
      local_tHist_smeared_xSection_mc_scaled.Scale(1e40)

      # Prescription for determining minimum and maximum y values from individual histograms.
      minycontent_local_tHist_unfolded_xSection_scaled = local_tHist_unfolded_xSection_scaled.GetMinimum()
      miny_local_tHist_unfolded_xSection_scaled = minycontent_local_tHist_unfolded_xSection_scaled - (local_tHist_unfolded_xSection_scaled.GetBinError(local_tHist_unfolded_xSection_scaled.GetMinimumBin()) if minycontent_local_tHist_unfolded_xSection_scaled < 0 else 0)
      maxy_local_tHist_unfolded_xSection_scaled = local_tHist_unfolded_xSection_scaled.GetMaximum() + local_tHist_unfolded_xSection_scaled.GetBinError(local_tHist_unfolded_xSection_scaled.GetMaximumBin())
      miny_local_tHist_xSection_mc_scaled = local_tHist_xSection_mc_scaled.GetMinimum()
      maxy_local_tHist_xSection_mc_scaled = local_tHist_xSection_mc_scaled.GetMaximum()
      miny_local_tHist_smeared_xSection_mc_scaled = local_tHist_smeared_xSection_mc_scaled.GetMinimum()
      maxy_local_tHist_smeared_xSection_mc_scaled = local_tHist_smeared_xSection_mc_scaled.GetMaximum()
      miny = min(0, miny_local_tHist_unfolded_xSection_scaled, miny_local_tHist_xSection_mc_scaled, miny_local_tHist_smeared_xSection_mc_scaled)*1.1
      maxy = max(0, maxy_local_tHist_unfolded_xSection_scaled, maxy_local_tHist_xSection_mc_scaled, maxy_local_tHist_smeared_xSection_mc_scaled)*1.1
    
      # Set plot styles.
      local_tHist_unfolded_xSection_scaled.SetMarkerSize(0.5)
      local_tHist_unfolded_xSection_scaled.SetMarkerStyle(ROOT.kFullCircle)
      local_tHist_unfolded_xSection_scaled.SetMarkerColor(ROOT.kBlack)
    
      local_tHist_smeared_xSection_mc_scaled.SetTitle(sigDefnp + " Exclusive Cross Section")
      local_tHist_smeared_xSection_mc_scaled.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_smeared_xSection_mc_scaled.GetYaxis().SetTitleSize(0.05)
      local_tHist_smeared_xSection_mc_scaled.GetXaxis().SetTitleSize(0.05)
      local_tHist_smeared_xSection_mc_scaled.GetYaxis().SetTitle("#sigma_{NC 1 #pi^{0}}[10^{-40} cm^{2}/Atom]")
      local_tHist_smeared_xSection_mc_scaled.GetXaxis().SetTitle(true_xlabel)
     
      local_tHist_xSection_mc_scaled.SetLineColor(ROOT.kCyan - 3)
      local_tHist_xSection_mc_scaled.SetMarkerColor(ROOT.kCyan - 3)
      local_tHist_xSection_mc_scaled.SetFillColor(ROOT.kCyan - 3)
      local_tHist_xSection_mc_scaled.SetFillStyle(3545)

      local_tHist_smeared_xSection_mc_scaled.SetLineColor(ROOT.kRed)
      local_tHist_smeared_xSection_mc_scaled.SetMarkerColor(ROOT.kRed)
      local_tHist_smeared_xSection_mc_scaled.SetFillColor(ROOT.kRed)
      local_tHist_smeared_xSection_mc_scaled.SetFillStyle(3554)

      # Draw only the smeared GENIE prediction if not in test configuration.
      overflow1 = DrawWithOverflow(local_tHist_smeared_xSection_mc_scaled, canvas.canvas, "HIST")
      if is_test:
        overflow2 = DrawWithOverflow(local_tHist_xSection_mc_scaled, canvas.canvas, "SAME HIST")
      overflow3 = DrawWithOverflow(local_tHist_unfolded_xSection_scaled, canvas.canvas, "SAME")

      # Draw the legend
      canvas.canvas.cd(1)
      legend = ROOT.TLegend(0.52,0.7,0.945,0.9, "")
      legend.SetBorderSize(0)
      legend.AddEntry(local_tHist_unfolded_xSection_scaled,data_string,"lep")
      if is_test:
        legend.AddEntry(local_tHist_xSection_mc_scaled,"GENIE Prediction","f")
      legend.AddEntry(local_tHist_smeared_xSection_mc_scaled, "GENIE Prediction, Smeared", "f")
      legend.Draw()

      if is_test:
        pt = ROOT.TPaveText(0.53, 0.5, 0.875, 0.7, "NDC")
        pt.AddText("#chi^{{2}}_{{unsmeared}}/DOF: {0:.0f}/{1}".format(chi2_GENIE, nBins + 1))
        pt.AddText("#chi^{{2}}_{{smeared}}/DOF: {0:.2f}/{1}".format(chi2_GENIE_smeared, nBins + 1))
      else:
        pt = ROOT.TPaveText(0.53, 0.6, 0.775, 0.7, "NDC")
        pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2_GENIE_smeared, nBins + 1))
        ptall.Draw()
      pt.SetBorderSize(0)
      pt.SetFillColorAlpha(0, 0)
      pt.Draw()
      canvas.canvas.cd(0)

    ### Plot of NuWro fake data unfolded events vs. NuWro generator prediction
    ##########################################################################
    with makeEnv_TCanvas('{0}/nuwro_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
      # Prescription for determining minimum and maximum y values from individual histograms.
      truth_max = local_tHist_nuwro_truth.GetMaximum()
      truth_smeared_max = local_tHist_smeared_nuwro_truth.GetMaximum()
      data_maxbin = local_tHist_unfolded_evtRate.GetMaximumBin()
      local_tHist_unfolded_evtRate.GetBinContent(data_maxbin)
      data_max = local_tHist_unfolded_evtRate.GetBinContent(data_maxbin)
      data_max += local_tHist_unfolded_evtRate.GetBinError(data_maxbin)
      truth_smeared_min = local_tHist_smeared_nuwro_truth.GetMinimum()
      data_minbin = local_tHist_unfolded_evtRate.GetMinimumBin()
      if local_tHist_unfolded_evtRate.GetBinContent(upperBound) < local_tHist_unfolded_evtRate.GetMinimum():
        data_minbin = upperBound;
      data_min = local_tHist_unfolded_evtRate.GetBinContent(data_minbin)
      data_min -= local_tHist_unfolded_evtRate.GetBinError(data_minbin) if data_min < 0 else 0
      maxy = max(truth_max, truth_smeared_max, data_max)*1.1
      miny = min(0, truth_smeared_min, data_min)*1.1

      # Set plot styles.
      local_tHist_smeared_nuwro_truth.SetLineColor(ROOT.kViolet)
      local_tHist_smeared_nuwro_truth.SetFillColor(ROOT.kViolet)
      local_tHist_smeared_nuwro_truth.SetFillStyle(3554)
      local_tHist_smeared_nuwro_truth.SetLineWidth(1)
      local_tHist_smeared_nuwro_truth.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_smeared_nuwro_truth.GetYaxis().SetTitleSize(0.05)
      local_tHist_smeared_nuwro_truth.GetXaxis().SetTitleSize(0.05)
      local_tHist_smeared_nuwro_truth.SetTitle(sigDefnp + " Exclusive Unfolded Events")
      local_tHist_smeared_nuwro_truth.GetYaxis().SetTitle("Events")
      local_tHist_smeared_nuwro_truth.GetXaxis().SetTitle(true_xlabel)
      overflow1 = DrawWithOverflow(local_tHist_smeared_nuwro_truth, canvas.canvas, "HIST")

      local_tHist_nuwro_truth.SetMarkerColor(ROOT.kGreen -3)
      local_tHist_nuwro_truth.SetLineColor(ROOT.kGreen - 3)
      local_tHist_nuwro_truth.SetFillColor(ROOT.kGreen -3)
      local_tHist_nuwro_truth.SetFillStyle(3545)
      local_tHist_nuwro_truth.SetLineWidth(1)

      # Draw unsmeared distribution and corresponding labels only in test mode.
      if is_test:
        overflow2 = DrawWithOverflow(local_tHist_nuwro_truth, canvas.canvas, "HIST SAME")
      overflow3 = DrawWithOverflow(local_tHist_unfolded_evtRate, canvas.canvas, "SAME")
      canvas.canvas.cd(1)
      legend = ROOT.TLegend(0.52,0.7,0.945,0.9, "")
      legend.SetBorderSize(0)
      legend.AddEntry(local_tHist_unfolded_evtRate,"Unfolded " + data_string,"lep")
      if is_test:
        legend.AddEntry(local_tHist_nuwro_truth,"NuWro Truth","f")
      legend.AddEntry(local_tHist_smeared_nuwro_truth, "NuWro Truth, Smeared", "f")
      legend.Draw()
      if is_test:
        pt = ROOT.TPaveText(0.53, 0.5, 0.875, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}_{{unsmeared}}/DOF: {0:.0f}/{1}".format(chi2_NuWro, nBins + 1))
        pt.AddText("#chi^{{2}}_{{smeared}}/DOF: {0:.2f}/{1}".format(chi2_NuWro_smeared, nBins + 1))
        pt.Draw()
      else:
        pt = ROOT.TPaveText(0.53, 0.6, 0.775, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2_NuWro_smeared, nBins))
        pt.Draw()
        ptall.Draw()
      canvas.canvas.cd(0)
      #nbins = local_tHist_nuwro_truth.GetNbinsX()
      #exec("local_tHist_unfolded_cov_evtRate = local_tHist_unfolded_cov_evtRate_{0}".format(sigDef))
      #local_tMat_unfolded_cov_evtRate = 

    ### Plot of NuWro fake data folded events vs. GENIE generator prediction
    ########################################################################
    with makeEnv_TCanvas('{0}/fakedatavsgenie_evtRate_folded_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      ## Create local copies of all needed hists
      exec("local_mHist_evtRate_reco = histFile.Get(\"evtRate_{0}\")".format(sigDef))
      local_tHist_evtRate_reco = local_mHist_evtRate_reco.GetCVHistoWithError()
    
      exec("local_tHist_genie_evtRate_reco = histFile.Get(\"folded_true_signal_{0}\")".format(sigDef))
     
      exec("local_mHist_hreco = histFile.Get(\"hreco_{0}\")".format(sigDefnp))
      local_tHist_hreco = local_mHist_hreco.GetCVHistoWithStatError()
      local_tHist_hreco.Multiply(tHist_POT_scaling)
     
      exec("local_mHist_effNum_reco = histFile.Get(\"effNum_reco_{0}\")".format(sigDef))
      local_tHist_effNum_reco = local_mHist_effNum_reco.GetCVHistoWithError()
      local_tMat_cov_effNum_reco = local_mHist_effNum_reco.GetTotalErrorMatrix(True)

      # Get folded covariance matrix to compare folded distributions. Remove underflow and overflow bins for now.
      # Does this make sense? See ref C's comments on Ben Bogart's numu CC paper. (???)
      # On the one hand, the fake data should not have signal systematic uncertainties prior to unfolding since it's "data".
      # On the other hand, do we want the covariance matrix we're unfolding to be different from the one we use for our reco distributions?
      # Currently subtracting signal systematic component of the folded covariance matrix.
      exec("local_tMat_cov_evtRate = ROOT.TMatrixD(nBins + 2, nBins + 2, tHist_cov_evtRate_{0}.GetArray())".format(sigDef))
      local_tMat_cov_evtRate = local_tMat_cov_evtRate.GetSub(lowerBound, upperBound, lowerBound, upperBound)
      local_tMat_cov_effNum_reco = local_tMat_cov_effNum_reco.GetSub(lowerBound, upperBound, lowerBound, upperBound)
      local_tMat_cov_reco_evtRate = ROOT.TMatrixD(upperBound + include_underflow_reco, upperBound + include_underflow_reco)
      local_tMat_cov_reco_evtRate.Minus(local_tMat_cov_evtRate, local_tMat_cov_effNum_reco)

      exec("local_mHist_background = mHist_background_{0}.Clone(\"local_mHist_background\")".format(sigDef))
      # Transfer all uncertainties onto the folded fake data (local_mHist_evtRate_reco) so that fractional uncertainties are calculated relative to its bin content.
      # Folded fake data is fake data minus MC background. The only systematic uncertainties should come from MC background.
      local_mHist_background.TransferErrorBands(local_mHist_evtRate_reco, False)
      # Calculate data statistical uncertainty using prescription in calculateChi2.py
      local_mHist_data_selected = histFile.Get("data_selected_" + sigDefnp)
      local_mHist_fakedata_mc = histFile.Get("fakedata_mc_" + sigDef)
      exec("local_tMat_folded_covariance = ROOT.TMatrixD(nBins + 2, nBins + 2, tHist_cov_evtRate_{0}.GetArray())".format(sigDef))
      local_tMat_folded_covariance_sub = local_tMat_folded_covariance.GetSub(1, nBins + 1, 1, nBins + 1)
      local_tMat_folded_covariance_with_stat_sub = calculateChi2(local_mHist_data_selected, local_mHist_fakedata_mc, local_tMat_folded_covariance_sub, False)[1] 
      local_tMat_folded_covariance_with_stat = local_tMat_folded_covariance.Clone()
      for i in range(nBins + 1):
        for j in range(nBins + 1):
          local_tMat_folded_covariance_with_stat[i + 1][j + 1] = local_tMat_folded_covariance_with_stat_sub[i][j]
        local_mHist_evtRate_reco.SetBinError(i, local_mHist_background.GetBinError(i))
      local_tMat_stat_error = ROOT.TMatrixD(nBins+2, nBins + 2)
      local_tMat_stat_error.Minus(local_tMat_folded_covariance_with_stat, local_tMat_folded_covariance)
      local_mHist_evtRate_reco.FillSysErrorMatrix("data_statistical", local_tMat_stat_error)
      local_tMat_cov_reco_evtRate = local_mHist_evtRate_reco.GetTotalErrorMatrix(False)
      local_tMat_cov_reco_evtRate = local_tMat_cov_reco_evtRate.GetSub(1, nBins + 1, 1, nBins + 1)
      
      # Set errors to square root of diagonal of folded covariance matrix.
      for i in range(upperBound + include_underflow_reco):
        local_tHist_effNum_reco.SetBinError(i + lowerBound, math.sqrt(local_tMat_cov_reco_evtRate[i][i]))
        local_tHist_evtRate_reco.SetBinError(i + lowerBound, math.sqrt(local_tMat_cov_reco_evtRate[i][i]))

      # Calculate chi2 before scaling.
      chi2_effNum_hreco = calculateChi2(local_tHist_hreco, local_tHist_effNum_reco, local_tMat_cov_reco_evtRate, True)[0]
      chi2_effNum_genie_evtRate = calculateChi2(local_tHist_genie_evtRate_reco, local_tHist_effNum_reco, local_tMat_cov_reco_evtRate, True)[0]
      chi2_evtRate_reco_effNum = calculateChi2(local_tHist_effNum_reco, local_tHist_evtRate_reco, local_tMat_cov_reco_evtRate, True)[0]

      #local_tHist_evtRate_reco.Scale(1e-2,"width")
      #local_tHist_genie_evtRate_reco.Scale(1e-2,"width")
      #local_tHist_hreco.Scale(1e-2,"width")
      #local_tHist_effNum_reco.Scale(1e-2,"width")

      # Prescription for determining minimum and maximum y values from individual histograms.
      miny_hreco = local_tHist_hreco.GetMinimum()
      miny_genie_evtRate_reco = local_tHist_genie_evtRate_reco.GetMinimum()
      minycont_effNum_reco = local_tHist_effNum_reco.GetMinimum()
      miny_effNum_reco = minycont_effNum_reco - (local_tHist_effNum_reco.GetBinError(local_tHist_effNum_reco.GetMinimumBin()) if minycont_effNum_reco < 0 and is_closure_test else 0)
      minycont_evtRate_reco = local_tHist_evtRate_reco.GetMinimum()
      miny_evtRate_reco = minycont_evtRate_reco - (local_tHist_evtRate_reco.GetBinError(local_tHist_evtRate_reco.GetMinimumBin()) if minycont_evtRate_reco < 0 else 0)
      maxy_hreco = local_tHist_hreco.GetMaximum()
      maxy_genie_evtRate_reco = local_tHist_genie_evtRate_reco.GetMaximum()
      maxybin_effNum_reco = local_tHist_effNum_reco.GetMaximumBin()
      maxy_effNum_reco = local_tHist_effNum_reco.GetBinContent(maxybin_effNum_reco) + (local_tHist_effNum_reco.GetBinError(maxybin_effNum_reco) if is_closure_test else 0)
      maxybin_evtRate_reco = local_tHist_evtRate_reco.GetMaximumBin()
      maxy_evtRate_reco = local_tHist_evtRate_reco.GetBinContent(maxybin_evtRate_reco) + local_tHist_evtRate_reco.GetBinError(maxybin_evtRate_reco)
      miny = min(0, (miny_hreco if is_closure_test else 0), (miny_genie_evtRate_reco if is_closure_test else 0), miny_effNum_reco, (miny_evtRate_reco if not is_closure_test else 0))*1.1
      maxy = max(0, (maxy_hreco if is_closure_test else 0), (maxy_genie_evtRate_reco if is_closure_test else 0), maxy_effNum_reco, (maxy_evtRate_reco if not is_closure_test else 0))*1.1

      ## Set plot formatting 
      local_tHist_evtRate_reco.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_evtRate_reco.SetMarkerSize(0.5)
      local_tHist_evtRate_reco.SetMarkerStyle(ROOT.kFullCircle)
      local_tHist_evtRate_reco.SetMarkerColor(ROOT.kBlack)
      local_tHist_evtRate_reco.GetYaxis().SetTitleSize(0.05)
      local_tHist_evtRate_reco.GetXaxis().SetTitleSize(0.05)
      local_tHist_evtRate_reco.SetTitle(sigDefnp + " Exclusive Folded Events") 
      local_tHist_evtRate_reco.GetYaxis().SetTitle("Events")
      local_tHist_evtRate_reco.GetXaxis().SetTitle(reco_xlabel)
      if is_closure_test:
        local_tHist_evtRate_reco.SetLineColor(ROOT.kCyan-3)
      else:
        local_tHist_evtRate_reco.SetLineColor(ROOT.kBlue + 2)
    
      local_tHist_hreco.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_hreco.GetYaxis().SetTitleSize(0.05)
      local_tHist_hreco.GetXaxis().SetTitleSize(0.05)
      local_tHist_hreco.SetTitle(sigDefnp + " Exclusive Folded Closure Test") 
      local_tHist_hreco.GetYaxis().SetTitle("Events")
      local_tHist_hreco.GetXaxis().SetTitle(reco_xlabel)

      if is_closure_test:
        local_tHist_effNum_reco.SetLineColor(ROOT.kGreen+2)
        local_tHist_effNum_reco.SetMarkerColor(ROOT.kGreen+2)
        local_tHist_effNum_reco.SetFillColor(ROOT.kGreen+2)
      else:
        local_tHist_effNum_reco.SetLineColor(ROOT.kRed)
        local_tHist_effNum_reco.SetMarkerColor(ROOT.kRed)
        local_tHist_effNum_reco.SetFillColor(ROOT.kRed)
      local_tHist_effNum_reco.SetMarkerSize(0.5)
      local_tHist_effNum_reco.SetFillStyle(3545)
    
      local_tHist_genie_evtRate_reco.SetLineColor(ROOT.kRed)
      local_tHist_genie_evtRate_reco.SetFillColor(ROOT.kRed)
      local_tHist_genie_evtRate_reco.SetMarkerColor(ROOT.kRed)
      local_tHist_genie_evtRate_reco.SetMarkerSize(0.5)
      local_tHist_genie_evtRate_reco.SetFillStyle(3554)
    
      local_tHist_hreco.SetLineColor(ROOT.kCyan-3)
      local_tHist_hreco.SetFillColor(ROOT.kCyan-3)
      local_tHist_hreco.SetMarkerColor(ROOT.kCyan-3)
      local_tHist_hreco.SetMarkerSize(0.5)
      local_tHist_hreco.SetFillStyle(3545)
      
      legend = ROOT.TLegend(0.57,0.7,0.945,0.9, "")
      legend.SetBorderSize(0)
      # If it's a closure test, draw the same GENIE reco distribution calculated using three different methods.
      if is_closure_test:
        overflow1 = DrawWithOverflow(local_tHist_hreco, canvas.canvas, "HIST")
        overflow2 = DrawWithOverflow(local_tHist_genie_evtRate_reco, canvas.canvas, "SAME HIST")
        overflow3 = DrawWithOverflow(local_tHist_effNum_reco, canvas.canvas, "SAME")
        legend.AddEntry(local_tHist_effNum_reco,"GENIE Reconstructed (from Sbnfit)","lep")
        legend.AddEntry(local_tHist_hreco,"GENIE Reconstructed (from ResponseMaker)","f")
        legend.AddEntry(local_tHist_genie_evtRate_reco,"GENIE Truth (from Sbnfit) Folded","f")
        canvas.canvas.cd(1)
        pt = ROOT.TPaveText(0.58, 0.5, 0.875, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}_{{ResponseMaker}}/DOF: {0:.2f}/{1}".format(chi2_effNum_hreco, nBins + 1))
        pt.AddText("#chi^{{2}}_{{Folded}}/DOF: {0:.2f}/{1}".format(chi2_effNum_genie_evtRate, nBins + 1))
        pt.Draw()
      else:
        overflow4 = DrawWithOverflow(local_tHist_evtRate_reco, canvas.canvas, "")
        legend.AddEntry(local_tHist_evtRate_reco,data_string,"lep")
        legend.AddEntry(local_tHist_effNum_reco,"GENIE Reco Signal","f")
        overflow5 = DrawWithOverflow(local_tHist_effNum_reco, canvas.canvas, "HIST SAME") 
        canvas.canvas.cd(1)
        pt = ROOT.TPaveText(0.58, 0.6, 0.775, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2_evtRate_reco_effNum, nBins + 1))
        pt.Draw()
      legend.Draw()
      ptall.Draw()
      canvas.canvas.cd(0)
   
    ### Plot of NuWro fake data folded events vs. NuWro generator prediction
    ########################################################################
    with makeEnv_TCanvas('{0}/nuwro_evtRate_folded_{1}.png'.format(plotDir, sigDef)) as canvas:
      # The only difference between NuWro fake data and NuWro reco signal is the background, so we extract and only use the background uncertainty.
      local_tHist_nuwro_signal = histFile.Get("nu_uBooNE_breakdown_" + sigDefnp + "sig")
      local_tHist_nuwro_background = histFile.Get("nu_uBooNE_breakdown_" + sigDefnp + "bkg")
      local_tHist_evtRate_reco = local_mHist_evtRate_reco.GetCVHistoWithError()
      local_tHist_background = local_mHist_background.GetCVHistoWithError()
      local_tMat_background_covMat = local_mHist_background.GetTotalErrorMatrix(True) 
      local_tMat_background_covMat = local_tMat_background_covMat.GetSub(lowerBound, upperBound, lowerBound, upperBound)

      # Update covariance with data statistical uncertainties and set error bars to the diagonal.
      chi2_evtRate_reco_nuwro_signal, local_tMat_background_covMat = calculateChi2(local_tHist_nuwro_background, local_tHist_background, local_tMat_background_covMat, False)
      for i in range(upperBound):
        local_tHist_evtRate_reco.SetBinError(i + 1, math.sqrt(local_tMat_background_covMat[i][i]))
      
      # Scale and set max y.
      #local_tHist_nuwro_signal.Scale(1e-2, "width")
      #local_tHist_evtRate_reco.Scale(1e-2, "width")

      maxy_nuwro_signal = local_tHist_nuwro_signal.GetMaximum()
      maxy = max(maxy_nuwro_signal, maxy_evtRate_reco)*1.1

      # Set plot formatting.
      local_tHist_nuwro_signal.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_nuwro_signal.GetYaxis().SetTitleSize(0.05)
      local_tHist_nuwro_signal.GetXaxis().SetTitleSize(0.05)
      local_tHist_nuwro_signal.SetTitle(sigDefnp + " Exclusive Folded Events")
      local_tHist_nuwro_signal.GetYaxis().SetTitle("Events")
      local_tHist_nuwro_signal.GetXaxis().SetTitle(reco_xlabel)
      local_tHist_nuwro_signal.SetMarkerColor(ROOT.kViolet)
      local_tHist_nuwro_signal.SetLineColor(ROOT.kViolet)
      local_tHist_nuwro_signal.SetFillColor(ROOT.kViolet)
      local_tHist_nuwro_signal.SetFillStyle(3545)
      local_tHist_nuwro_signal.SetLineWidth(1)

      local_tHist_evtRate_reco.SetMarkerSize(0.5)
      local_tHist_evtRate_reco.SetLineColor(ROOT.kBlue+2)
      local_tHist_evtRate_reco.SetMarkerStyle(ROOT.kFullCircle)
      local_tHist_evtRate_reco.SetMarkerColor(ROOT.kBlack)
      local_tHist_evtRate_reco.SetFillStyle(0)
      local_tHist_evtRate_reco.SetLineWidth(1)

      # Draw plot.
      overflow1 = DrawWithOverflow(local_tHist_nuwro_signal, canvas.canvas, "HIST")
      overflow2 = DrawWithOverflow(local_tHist_evtRate_reco, canvas.canvas, "SAME E0")

      canvas.canvas.cd(1)
      legend = ROOT.TLegend(0.57,0.7,0.965,0.9, "")
      legend.SetBorderSize(0)
      legend.AddEntry(local_tHist_evtRate_reco, data_string + " - GENIE Bkgd.", "lep")
      legend.AddEntry(local_tHist_nuwro_signal, "NuWro Reco Signal", "f")
      legend.Draw()

      pt = ROOT.TPaveText(0.58, 0.6, 0.825, 0.7, "NDC")
      pt.SetBorderSize(0)
      pt.SetFillColorAlpha(0, 0)
      pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2_evtRate_reco_nuwro_signal, nBins + 1))
      pt.Draw()
      ptall.Draw()
      canvas.canvas.cd(0)
  

    ### Plot of error breakdown for Monte Carlo Background
    ######################################################
    exec("localDrawErrorSummary(plotter, local_mHist_background, \"{0} Exclusive Background Error Summary\", reco_xlabel, \"{2}/errorSummary_background_{1}.png\")".format(sigDefnp, sigDef, plotDir))

    ### Plot of error breakdown for Folded Fake Data
    ######################################################
    exec("localDrawErrorSummary(plotter, local_mHist_evtRate_reco, \"{0} Exclusive Folded Events Error Summary\", reco_xlabel, \"{2}/errorSummary_evtRate_folded_{1}.png\")".format(sigDefnp, sigDef, plotDir))

    ### Plot of error breakdown for Unfolded Fake Data
    ######################################################
    exec("localDrawErrorSummary(plotter, mHist_unfolded_evtRate_{1}, \"{0} Exclusive Unfolded Events Error Summary\", true_xlabel, \"{2}/errorSummary_evtRate_{1}.png\")".format(sigDefnp, sigDef, plotDir))

