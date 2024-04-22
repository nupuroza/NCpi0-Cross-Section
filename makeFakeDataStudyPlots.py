### Script to take make plots, etc.
### using the MINERvA Analysis Toolkit

import ROOT
import datetime as dt
import argparse
import os,math
from customHistAndPlotMethods import makeEnv_TCanvas,localDrawErrorSummary
from calculateChi2 import *

## Set ROOT to batch mode
ROOT.gROOT.SetBatch()

## This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

## Setup MnvPlotter, which has all of the plotting utilities
plotter = ROOT.PlotUtils.MnvPlotter()
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
### Set Configurations ######################################################################################
#############################################################################################################

## Prescription is slightly different for closure test input 
is_closure_test = True if p.closureTest>0 else False

## Running in the test configuration adds unsmeared distributions to the plots and chi2 values to closure tests.
is_test = True if p.test > 0 else False

#############################################################################################################
### Pull out relevant objects from input file ###############################################################
#############################################################################################################

#for sigDef in ["2g1p","2g0p","2gnp"]:
for sigDefnp in ["2g1p","2g0p"]:
  #for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["exclusive"]:

    if sigDefnp == "2gnp" and sigDefexcl == "exclusive":
      continue
  
    sigDef = sigDefnp+"_"+sigDefexcl

    ### Get THnDs
    ################
    for histCat in ["cov_evtRate","unfolded_cov_evtRate","unfolded_evtRate","unfolded_xSection","add_smear_matrix", "smeared_xSection_mc", "smeared_nuwro_signal"]:
    
      exec("tHist_{0}_{1} = histFile.Get(\"{0}_{1}\")".format(histCat,sigDef))
    
    ### Get MnvHnDs
    ################
    for histCat in ["xSection_mc", "background"]:
    
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
for sigDefnp in ["2g1p", "2g0p"]:
  #for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["exclusive"]:

    if sigDefnp == "2gnp" and sigDefexcl == "exclusive":
      continue
  
    sigDef = sigDefnp+"_"+sigDefexcl

    plotter.SetROOT6Palette(54)
    
    with makeEnv_TCanvas('{0}/response_{1}.png'.format(plotDir,sigDef)) as canvas:
      
      local_tHist_Response = histFile.Get("Response_{0}".format(sigDefnp))
      ROOT.gStyle.SetOptTitle(1)
      local_tHist_Response.SetTitle(sigDefnp + " Exclusive Response Matrix")
      local_tHist_Response.GetXaxis().SetTitle("True #pi^{0} Momentum")
      local_tHist_Response.GetYaxis().SetTitle("Reco #pi^{0} Momentum")
      #canvas.canvas.SetLogz()
      local_tHist_Response.Draw("colz")
      ptall.Draw()

    # Probably not accurate (see translateHists.py). (!!!)
    with makeEnv_TCanvas('{0}/efficiency_{1}.png'.format(plotDir, sigDef)) as canvas:
      local_tHist_eff = histFile.Get("eff_{0}".format(sigDef))
      local_tHist_eff.SetTitle(sigDefnp + " Exclusive Efficiency")
      local_tHist_eff.GetXaxis().SetTitle("True #pi^{0} Momentum")
      local_tHist_eff.GetYaxis().SetTitle("Efficiency")
      local_tHist_eff.Draw("HIST")
      local_tHist_eff.SetFillStyle(3545)
      local_tHist_eff.Draw("E2 SAME")
    
    with makeEnv_TCanvas('{0}/add_smear_matrix_{1}.png'.format(plotDir, sigDef)) as canvas:
      local_tHist_add_smear_matrix = histFile.Get("add_smear_matrix_{0}".format(sigDef))
      local_tHist_add_smear_matrix.Draw("colz")
      local_tHist_add_smear_matrix.SetTitle(sigDefnp + " Exclusive Additional Smearing Matrix")
      local_tHist_add_smear_matrix.GetXaxis().SetTitle("True #pi^{0} Momentum")
      local_tHist_add_smear_matrix.GetYaxis().SetTitle("Reco #pi^{0} Momentum")
      ptall.SetTextColor(ROOT.kWhite)
      ptall.Draw()
      #canvas.canvas.SetLogz()
    
    with makeEnv_TCanvas('{0}/cov_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_cov_evtRate_{0} = tHist_cov_evtRate_{0}.Clone(\"local_tHist_cov_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_cov_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_cov_evtRate_{0}.SetTitle(\"{1} Exclusive Folded Covariance Matrix\")".format(sigDef, sigDefnp))
      exec("local_tHist_cov_evtRate_{0}.GetXaxis().SetTitle(\"Reco #pi^{{0}} Momentum\")".format(sigDef))
      exec("local_tHist_cov_evtRate_{0}.GetYaxis().SetTitle(\"Reco #pi^{{0}} Momentum\")".format(sigDef))
      ptall.Draw()
      #canvas.canvas.SetLogz()
    
    with makeEnv_TCanvas('{0}/unfolded_cov_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_unfolded_cov_evtRate_{0} = tHist_unfolded_cov_evtRate_{0}.Clone(\"local_tHist_unfolded_cov_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_unfolded_cov_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_unfolded_cov_evtRate_{0}.SetTitle(\"{1} Exclusive Unfolded Covariance Matrix\")".format(sigDef, sigDefnp))
      exec("local_tHist_unfolded_cov_evtRate_{0}.GetXaxis().SetTitle(\"True #pi^{{0}} Momentum\")".format(sigDef))
      exec("local_tHist_unfolded_cov_evtRate_{0}.GetYaxis().SetTitle(\"True #pi^{{0}} Momentum\")".format(sigDef))
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
      exec("local_tHist_corr_evtRate_{0}.GetXaxis().SetTitle(\"Reco #pi^{{0}} Momentum\")".format(sigDef))
      exec("local_tHist_corr_evtRate_{0}.GetYaxis().SetTitle(\"Reco #pi^{{0}} Momentum\")".format(sigDef))
      ptall.SetTextColor(ROOT.kBlack)
      ptall.Draw()
    
    with makeEnv_TCanvas('{0}/unfolded_corr_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_unfolded_corr_evtRate_{0} = tHist_unfolded_corr_evtRate_{0}.Clone(\"local_tHist_unfolded_corr_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.GetZaxis().SetRangeUser(-1,1)".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.SetTitle(\"{1} Unfolded Correlation Matrix\")".format(sigDef, sigDefnp))
      exec("local_tHist_unfolded_corr_evtRate_{0}.GetXaxis().SetTitle(\"True #pi^{{0}} Momentum\")".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.GetYaxis().SetTitle(\"True #pi^{{0}} Momentum\")".format(sigDef))
      ptall.Draw()
    
#############################################################################################################
### Make fake data study comparison plots ###################################################################
#############################################################################################################

# Move TPaveText for TH1D plots.
ptall.SetX1NDC(0.55)
ptall.SetY1NDC(0.5)
ptall.SetX2NDC(0.825)
ptall.SetY2NDC(0.6)

#for sigDef in ["2g1p","2g0p","2gnp"]:
for sigDefnp in ["2g1p","2g0p"]:
  #for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["exclusive"]:

    if sigDefnp == "2gnp" and sigDefexcl == "exclusive":
      continue

    sigDef = sigDefnp+"_"+sigDefexcl

        
    ### Plot of NuWro fake data unfolded events vs. GENIE generator prediction
    ##########################################################################
    with makeEnv_TCanvas('{0}/fakedatavsgenie_evtRate_{1}.png'.format(plotDir,sigDef)):
    
      ## Create local copies of all needed hists
      exec("local_tHist_unfolded_evtRate_scaled = tHist_unfolded_evtRate_{0}.Clone(\"local_tHist_unfolded_evtRate_{0}\")".format(sigDef))
      local_tHist_genie_evtRate_smeared_scaled = histFile.Get("smeared_true_signal_{0}".format(sigDef))
      local_tHist_genie_evtRate_scaled = histFile.Get("prior_true_signal_{0}".format(sigDef))
      local_mHist_effDenom_scaled = histFile.Get("effDenom_{0}".format(sigDef))
      local_tHist_effDenom_scaled =local_mHist_effDenom_scaled.GetCVHistoWithError()
      exec("local_tHist_nuwro_truth_scaled = histFile.Get(\"nu_uBooNE_denom_{0}\")".format(sigDefnp))
      exec("local_tHist_smeared_nuwro_truth_scaled = tHist_smeared_nuwro_signal_{0}.Clone(\"smeared_nuwro_truth_scaled\")".format(sigDef))
    
      for i in range(nBins+2):## Loop over bins
        exec("cov_binContent = tHist_unfolded_cov_evtRate_{0}.GetBinContent(i,i)".format(sigDef))
        local_tHist_unfolded_evtRate_scaled.SetBinError(i,math.sqrt(cov_binContent))
        if is_closure_test:
          local_tHist_unfolded_evtRate_scaled.SetBinError(i, 0.00001)
      
      
      ## Calculate chi2 before scaling.
      local_tMat_unfolded_cov_evtRate = ROOT.TMatrixD(nBins, nBins)
      for i in range(1, nBins+1):
        for j in range(1, nBins+1):
          exec("local_tMat_unfolded_cov_evtRate[i-1][j-1] = tHist_unfolded_cov_evtRate_{0}.GetBinContent(i, j)".format(sigDef))
      chi2_NuWro = calculateChi2(local_tHist_nuwro_truth_scaled, local_tHist_unfolded_evtRate_scaled, local_tMat_unfolded_cov_evtRate, True)
      chi2_NuWro_smeared = calculateChi2(local_tHist_smeared_nuwro_truth_scaled, local_tHist_unfolded_evtRate_scaled, local_tMat_unfolded_cov_evtRate, True)
      chi2_GENIE = calculateChi2(local_tHist_genie_evtRate_scaled, local_tHist_unfolded_evtRate_scaled, local_tMat_unfolded_cov_evtRate, True)
      chi2_GENIE_smeared = calculateChi2(local_tHist_genie_evtRate_smeared_scaled, local_tHist_unfolded_evtRate_scaled, local_tMat_unfolded_cov_evtRate, True)
    
      ## Scale and bin-width-normalize all hists (and POT normalize NuWro truth)
      local_tHist_unfolded_evtRate_scaled.Scale(1e-3,"width")
      local_tHist_genie_evtRate_smeared_scaled.Scale(1e-3,"width")
      local_tHist_genie_evtRate_scaled.Scale(1e-3,"width")
      local_tHist_effDenom_scaled.Scale(1e-3,"width")
      local_tHist_nuwro_truth_scaled.Scale(1e-3,"width")
      local_tHist_smeared_nuwro_truth_scaled.Scale(1e-3, "width")
           
      ## Calculate yrange
      minycontent_local_tHist_unfolded_evtRate_scaled = local_tHist_unfolded_evtRate_scaled.GetMinimum()
      miny_local_tHist_unfolded_evtRate_scaled = minycontent_local_tHist_unfolded_evtRate_scaled - (local_tHist_unfolded_evtRate_scaled.GetBinError(local_tHist_unfolded_evtRate_scaled.GetMinimumBin()) if minycontent_local_tHist_unfolded_evtRate_scaled < 0 else 0)
      maxy_local_tHist_unfolded_evtRate_scaled = local_tHist_unfolded_evtRate_scaled.GetMaximum() + local_tHist_unfolded_evtRate_scaled.GetBinError(local_tHist_unfolded_evtRate_scaled.GetMaximumBin())
      miny_local_tHist_genie_evtRate_scaled = local_tHist_genie_evtRate_scaled.GetMinimum()
      maxy_local_tHist_genie_evtRate_scaled = local_tHist_genie_evtRate_scaled.GetMaximum()
      miny_local_tHist_genie_evtRate_smeared_scaled = local_tHist_genie_evtRate_smeared_scaled.GetMinimum()
      maxy_local_tHist_genie_evtRate_smeared_scaled = local_tHist_genie_evtRate_smeared_scaled.GetMaximum()
      miny = min(0, miny_local_tHist_unfolded_evtRate_scaled, (miny_local_tHist_genie_evtRate_scaled if not is_closure_test else 0), miny_local_tHist_genie_evtRate_smeared_scaled)*1.1
      maxy = max(0, maxy_local_tHist_unfolded_evtRate_scaled, (maxy_local_tHist_genie_evtRate_scaled if not is_closure_test else 0), maxy_local_tHist_genie_evtRate_smeared_scaled)*1.1

      ## Set plot formatting 
      local_tHist_unfolded_evtRate_scaled.SetMarkerSize(0.5)
      local_tHist_unfolded_evtRate_scaled.SetLineColor(ROOT.kBlue+2)
      local_tHist_unfolded_evtRate_scaled.SetMarkerColor(ROOT.kBlack)
      local_tHist_unfolded_evtRate_scaled.SetMarkerStyle(ROOT.kFullCircle)
      local_tHist_unfolded_evtRate_scaled.SetFillStyle(0)
      local_tHist_unfolded_evtRate_scaled.SetLineWidth(1)

      local_tHist_effDenom_scaled.SetLineColor(ROOT.kGreen+2)
      local_tHist_effDenom_scaled.SetMarkerColor(ROOT.kGreen+2)
      local_tHist_effDenom_scaled.SetMarkerSize(0.5)
      #local_tHist_effDenom_scaled.Draw("SAME")

      local_tHist_genie_evtRate_smeared_scaled.SetLineColor(ROOT.kRed)
      local_tHist_genie_evtRate_smeared_scaled.SetMarkerColor(ROOT.kRed)
      local_tHist_genie_evtRate_smeared_scaled.SetFillColor(ROOT.kRed)
      local_tHist_genie_evtRate_smeared_scaled.SetFillStyle(3554)
      local_tHist_genie_evtRate_smeared_scaled.SetTitle(sigDefnp + " Exclusive Unfolded Events") 
      local_tHist_genie_evtRate_smeared_scaled.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_genie_evtRate_smeared_scaled.GetYaxis().SetTitleSize(0.05)
      local_tHist_genie_evtRate_smeared_scaled.GetXaxis().SetTitleSize(0.05)
      local_tHist_genie_evtRate_smeared_scaled.GetYaxis().SetTitle("Events [10^{3}/GeV]")
      local_tHist_genie_evtRate_smeared_scaled.GetXaxis().SetTitle("True #pi^{0} Momentum [GeV]")

      local_tHist_genie_evtRate_smeared_scaled.Draw("HIST")

      ptall.Draw()

      # If it's a closure test, change the historgram title to reflect that.
      if is_closure_test:
        local_tHist_genie_evtRate_smeared_scaled.SetTitle(sigDefnp + " Exclusive Closure Test")
      
      # Print chi2 value of closure test using unfolded (fake) data covariance if in test mode.
      # If not in test mode, print only the smeared chi2 value regardless of whether or not it's a closure test.
      if (is_closure_test and is_test) or (not is_closure_test and not is_test):
        pt = ROOT.TPaveText(0.55, 0.6, 0.795, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2_GENIE_smeared, nBins))
        pt.Draw()
    
      # If configured to test mode, also draw and print chi2 of the unsmeared distribution for non-closure tests.
      elif is_test:
        local_tHist_genie_evtRate_scaled.SetMarkerColor(ROOT.kCyan - 3)
        local_tHist_genie_evtRate_scaled.SetLineColor(ROOT.kCyan - 3)
        local_tHist_genie_evtRate_scaled.SetFillColor(ROOT.kCyan - 3)
        local_tHist_genie_evtRate_scaled.SetFillStyle(3545)
        local_tHist_genie_evtRate_scaled.Draw("SAME HIST")
        
        pt = ROOT.TPaveText(0.5, 0.5, 0.845, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}_{{unsmeared}}/DOF: {0:.2f}/{1}".format(chi2_GENIE, nBins))
        pt.AddText("#chi^{{2}}_{{smeared}}/DOF: {0:.2f}/{1}".format(chi2_GENIE_smeared, nBins))
        pt.Draw()

      local_tHist_unfolded_evtRate_scaled.Draw("SAME")

      # Draw legend with appropriate labels.
      legend = ROOT.TLegend(0.42,0.7,0.845,0.9, "")
      legend.SetBorderSize(0)
      if is_closure_test:
        legend.AddEntry(local_tHist_unfolded_evtRate_scaled,"GENIE Reconstruction Unfolded","lep")
        legend.AddEntry(local_tHist_genie_evtRate_smeared_scaled,"GENIE Truth, Smeared","f")
      else:
        legend.AddEntry(local_tHist_unfolded_evtRate_scaled,"NuWro Fake Data Unfolded","lep")
        #legend.AddEntry(local_tHist_nuwro_truth_scaled,"NuWro Truth","f")
        #legend.AddEntry(local_tHist_effDenom_scaled,"effDenom","lep")
        if is_test:
          legend.AddEntry(local_tHist_genie_evtRate_scaled,"GENIE Truth","f")
        legend.AddEntry(local_tHist_genie_evtRate_smeared_scaled,"GENIE Truth, Smeared","f")
      legend.Draw()

    ### Plot of NuWro fake data cross section vs. GENIE generator prediction
    ########################################################################
    with makeEnv_TCanvas('{0}/fakedatavsgenie_xSection_{1}.png'.format(plotDir,sigDef)):
    
      # Cross sections are on the order of 1e-38. Scale by 1e38 and move the exponent to the units label.
      exec("local_tHist_unfolded_xSection_scaled = tHist_unfolded_xSection_{0}.Clone(\"local_tHist_unfolded_xSection_{0}\")".format(sigDef))
      local_tHist_unfolded_xSection_scaled.Scale(1e38,"width") 
      exec("local_tHist_xSection_mc_scaled = tHist_xSection_mc_{0}.Clone(\"local_tHist_xSection_mc_{0}\")".format(sigDef))
      local_tHist_xSection_mc_scaled.Scale(1e38,"width")
      exec("local_tHist_smeared_xSection_mc_scaled = tHist_smeared_xSection_mc_{0}.Clone(\"local_tHist_smeared_xSection_mc_{0}\")".format(sigDef))
      local_tHist_smeared_xSection_mc_scaled.Scale(1e38, "width")

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
      local_tHist_smeared_xSection_mc_scaled.GetYaxis().SetTitle("#sigma_{NC 1 #pi^{0}}[10^{-38} cm^{2}/Atom/GeV]")
      local_tHist_smeared_xSection_mc_scaled.GetXaxis().SetTitle("True #pi^{0} momentum [GeV]")
     
      local_tHist_xSection_mc_scaled.SetLineColor(ROOT.kCyan - 3)
      local_tHist_xSection_mc_scaled.SetMarkerColor(ROOT.kCyan - 3)
      local_tHist_xSection_mc_scaled.SetFillColor(ROOT.kCyan - 3)
      local_tHist_xSection_mc_scaled.SetFillStyle(3545)

      local_tHist_smeared_xSection_mc_scaled.SetLineColor(ROOT.kRed)
      local_tHist_smeared_xSection_mc_scaled.SetMarkerColor(ROOT.kRed)
      local_tHist_smeared_xSection_mc_scaled.SetFillColor(ROOT.kRed)
      local_tHist_smeared_xSection_mc_scaled.SetFillStyle(3554)

      # Draw only the smeared GENIE prediction if not in test configuration.
      local_tHist_smeared_xSection_mc_scaled.Draw("HIST")
      if is_test:
        local_tHist_xSection_mc_scaled.Draw("SAME HIST")
      local_tHist_unfolded_xSection_scaled.Draw("SAME")
      ptall.Draw()

      # Draw the legend
      legend = ROOT.TLegend(0.55,0.7,0.8,0.85, "")
      legend.SetBorderSize(0);
      legend.AddEntry(local_tHist_unfolded_xSection_scaled,"NuWro Fake Data","lep")
      if is_test:
        legend.AddEntry(local_tHist_xSection_mc_scaled,"GENIE Prediction","f")
      legend.AddEntry(local_tHist_smeared_xSection_mc_scaled, "GENIE Prediction, Smeared", "f")
      legend.Draw()

      if is_test:
        pt = ROOT.TPaveText(0.5, 0.5, 0.845, 0.7, "NDC")
      else:
        pt = ROOT.TPaveText(0.55, 0.6, 0.795, 0.7, "NDC")
      pt.SetBorderSize(0)
      pt.SetFillColorAlpha(0, 0)
      pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2_GENIE_smeared, nBins))
      pt.Draw()

    ### Plot of NuWro fake data unfolded events vs. NuWro generator prediction
    ##########################################################################
    with makeEnv_TCanvas('{0}/nuwro_evtRate_{1}.png'.format(plotDir,sigDef)):
      # Prescription for determining minimum and maximum y values from individual histograms.
      truth_max = local_tHist_nuwro_truth_scaled.GetMaximum()
      truth_smeared_max = local_tHist_smeared_nuwro_truth_scaled.GetMaximum()
      reco_maxbin = local_tHist_unfolded_evtRate_scaled.GetMaximumBin() local_tHist_unfolded_evtRate_scaled.GetBinContent(reco_maxbin)
      reco_max = local_tHist_unfolded_evtRate_scaled.GetBinContent(reco_maxbin)
      reco_max += local_tHist_unfolded_evtRate_scaled.GetBinError(reco_maxbin)
      truth_smeared_min = local_tHist_smeared_nuwro_truth_scaled.GetMinimum()
      reco_minbin = local_tHist_unfolded_evtRate_scaled.GetMinimumBin()
      reco_min = local_tHist_unfolded_evtRate_scaled.GetBinContent(reco_minbin)
      reco_min -= local_tHist_unfolded_evtRate_scaled.GetBinError(reco_minbin) if reco_min < 0 else 0
      maxy = max(truth_max, truth_smeared_max, reco_max)*1.1
      miny = min(0, truth_smeared_min, reco_min)*1.1

      # Set plot styles.
      local_tHist_smeared_nuwro_truth_scaled.SetLineColor(ROOT.kViolet)
      local_tHist_smeared_nuwro_truth_scaled.SetFillColor(ROOT.kViolet)
      local_tHist_smeared_nuwro_truth_scaled.SetFillStyle(3554)
      local_tHist_smeared_nuwro_truth_scaled.SetLineWidth(1)
      local_tHist_smeared_nuwro_truth_scaled.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_smeared_nuwro_truth_scaled.GetYaxis().SetTitleSize(0.05)
      local_tHist_smeared_nuwro_truth_scaled.GetXaxis().SetTitleSize(0.05)
      local_tHist_smeared_nuwro_truth_scaled.SetTitle(sigDefnp + " Exclusive Unfolded Events")
      local_tHist_smeared_nuwro_truth_scaled.GetYaxis().SetTitle("Events [10^{3}/GeV]")
      local_tHist_smeared_nuwro_truth_scaled.GetXaxis().SetTitle("True #pi^{0} momentum [GeV]")
      local_tHist_smeared_nuwro_truth_scaled.Draw("HIST")

      local_tHist_nuwro_truth_scaled.SetMarkerColor(ROOT.kGreen -3)
      local_tHist_nuwro_truth_scaled.SetLineColor(ROOT.kGreen - 3)
      local_tHist_nuwro_truth_scaled.SetFillColor(ROOT.kGreen -3)
      local_tHist_nuwro_truth_scaled.SetFillStyle(3545)
      local_tHist_nuwro_truth_scaled.SetLineWidth(1)

      # Draw unsmeared distribution and corresponding labels only in test mode.
      if is_test:
        local_tHist_nuwro_truth_scaled.Draw("HIST SAME")
      local_tHist_unfolded_evtRate_scaled.Draw("SAME")
      legend = ROOT.TLegend(0.42,0.7,0.845,0.9, "")
      legend.SetBorderSize(0);
      legend.AddEntry(local_tHist_unfolded_evtRate_scaled,"NuWro Fake Data Unfolded","lep")
      if is_test:
        legend.AddEntry(local_tHist_nuwro_truth_scaled,"NuWro Truth","f")
      legend.AddEntry(local_tHist_smeared_nuwro_truth_scaled, "NuWro Truth, Smeared", "f")
      legend.Draw()
      if is_test:
        pt = ROOT.TPaveText(0.5, 0.6, 0.845, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}_{{unsmeared}}/DOF: {0:.2f}/{1}".format(chi2_NuWro, nBins))
        pt.AddText("#chi^{{2}}_{{smeared}}/DOF: {0:.2f}/{1}".format(chi2_NuWro_smeared, nBins))
        pt.Draw()
      else:
        pt = ROOT.TPaveText(0.55, 0.6, 0.795, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2_NuWro_smeared, nBins))
        pt.Draw()
      ptall.Draw()
      #nbins = local_tHist_nuwro_truth_scaled.GetNbinsX()
      #exec("local_tHist_unfolded_cov_evtRate = local_tHist_unfolded_cov_evtRate_{0}".format(sigDef))
      #local_tMat_unfolded_cov_evtRate = 

    ### Plot of NuWro fake data folded events vs. GENIE generator prediction
    ########################################################################
    with makeEnv_TCanvas('{0}/fakedatavsgenie_evtRate_folded_{1}.png'.format(plotDir,sigDef)):
    
      ## Create local copies of all needed hists
      exec("local_mHist_evtRate_reco = histFile.Get(\"evtRate_{0}\")".format(sigDef))
      local_tHist_evtRate_reco_scaled = local_mHist_evtRate_reco.GetCVHistoWithError()
    
      exec("local_tHist_genie_evtRate_reco_scaled = histFile.Get(\"folded_true_signal_{0}\")".format(sigDef))
     
      exec("local_mHist_hreco_scaled = histFile.Get(\"hreco_{0}\")".format(sigDefnp))
      local_tHist_hreco_scaled = local_mHist_hreco_scaled.GetCVHistoWithStatError()
      local_tHist_hreco_scaled.Multiply(tHist_POT_scaling)
     
      exec("local_mHist_effNum_reco = histFile.Get(\"effNum_reco_{0}\")".format(sigDef))
      local_tHist_effNum_reco_scaled = local_mHist_effNum_reco.GetCVHistoWithError()

      # Get folded covariance matrix to compare folded distributions. Remove underflow and overflow bins for now.
      # Does this make sense? See ref C's comments on Ben Bogart's numu CC paper. (???)
      # On the one hand, the fake data should not have signal systematic uncertainties prior to unfolding since it's "data".
      # On the other hand, do we want the covariance matrix we're unfolding to be different from the one we use for our reco distributions?
      exec("local_tMat_cov_evtRate = ROOT.TMatrixD(nBins + 2, nBins + 2, tHist_cov_evtRate_{0}.GetArray())".format(sigDef))
      local_tMat_cov_evtRate = local_tMat_cov_evtRate.GetSub(1, nBins, 1, nBins)
      
      # Set errors to square root of diagonal of folded covariance matrix.
      for i in range(nBins):
        local_tHist_effNum_reco_scaled.SetBinError(i + 1, math.sqrt(local_tMat_cov_evtRate[i][i]))
        local_tHist_evtRate_reco_scaled.SetBinError(i + 1, math.sqrt(local_tMat_cov_evtRate[i][i]))

      # Calculate chi2 before scaling.
      chi2_effNum_hreco = calculateChi2(local_tHist_hreco_scaled, local_tHist_effNum_reco_scaled, local_tMat_cov_evtRate, True)
      chi2_effNum_genie_evtRate = calculateChi2(local_tHist_genie_evtRate_reco_scaled, local_tHist_effNum_reco_scaled, local_tMat_cov_evtRate, True)
      chi2_evtRate_reco_effNum = calculateChi2(local_tHist_effNum_reco_scaled, local_tHist_evtRate_reco_scaled, local_tMat_cov_evtRate, True)

      local_tHist_evtRate_reco_scaled.Scale(1e-2,"width")
      local_tHist_genie_evtRate_reco_scaled.Scale(1e-2,"width")
      local_tHist_hreco_scaled.Scale(1e-2,"width")
      local_tHist_effNum_reco_scaled.Scale(1e-2,"width")

      # Prescription for determining minimum and maximum y values from individual histograms.
      miny_hreco = local_tHist_hreco_scaled.GetMinimum()
      miny_genie_evtRate_reco = local_tHist_genie_evtRate_reco_scaled.GetMinimum()
      minycont_effNum_reco = local_tHist_effNum_reco_scaled.GetMinimum()
      miny_effNum_reco = minycont_effNum_reco - (local_tHist_effNum_reco_scaled.GetBinError(local_tHist_effNum_reco_scaled.GetMinimumBin()) if minycont_effNum_reco < 0 and is_closure_test else 0)
      minycont_evtRate_reco = local_tHist_evtRate_reco_scaled.GetMinimum()
      miny_evtRate_reco = minycont_evtRate_reco - (local_tHist_evtRate_reco_scaled.GetBinError(local_tHist_evtRate_reco_scaled.GetMinimumBin()) if minycont_evtRate_reco < 0 else 0)
      maxy_hreco = local_tHist_hreco_scaled.GetMaximum()
      maxy_genie_evtRate_reco = local_tHist_genie_evtRate_reco_scaled.GetMaximum()
      maxybin_effNum_reco = local_tHist_effNum_reco_scaled.GetMaximumBin()
      maxy_effNum_reco = local_tHist_effNum_reco_scaled.GetBinContent(maxybin_effNum_reco) + (local_tHist_effNum_reco_scaled.GetBinError(maxybin_effNum_reco) if is_closure_test else 0)
      maxybin_evtRate_reco = local_tHist_evtRate_reco_scaled.GetMaximumBin()
      maxy_evtRate_reco = local_tHist_evtRate_reco_scaled.GetBinContent(maxybin_evtRate_reco) + local_tHist_evtRate_reco_scaled.GetBinError(maxybin_evtRate_reco)
      miny = min(0, (miny_hreco if is_closure_test else 0), (miny_genie_evtRate_reco if is_closure_test else 0), miny_effNum_reco, (miny_evtRate_reco if not is_closure_test else 0))*1.1
      maxy = max(0, (maxy_hreco if is_closure_test else 0), (maxy_genie_evtRate_reco if is_closure_test else 0), maxy_effNum_reco, (maxy_evtRate_reco if not is_closure_test else 0))*1.1

      ## Set plot formatting 
      local_tHist_evtRate_reco_scaled.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_evtRate_reco_scaled.SetMarkerSize(0.5)
      local_tHist_evtRate_reco_scaled.SetMarkerStyle(ROOT.kFullCircle)
      local_tHist_evtRate_reco_scaled.SetMarkerColor(ROOT.kBlack)
      local_tHist_evtRate_reco_scaled.GetYaxis().SetTitleSize(0.05)
      local_tHist_evtRate_reco_scaled.GetXaxis().SetTitleSize(0.05)
      local_tHist_evtRate_reco_scaled.SetTitle(sigDefnp + " Exclusive Folded Events") 
      local_tHist_evtRate_reco_scaled.GetYaxis().SetTitle("Events [10^{2}/GeV]")
      local_tHist_evtRate_reco_scaled.GetXaxis().SetTitle("Reco #pi^{{0}} momentum [GeV]".format(sigDef))
      if is_closure_test:
        local_tHist_evtRate_reco_scaled.SetLineColor(ROOT.kCyan-3)
      else:
        local_tHist_evtRate_reco_scaled.SetLineColor(ROOT.kBlue + 2)
    
      local_tHist_hreco_scaled.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_hreco_scaled.GetYaxis().SetTitleSize(0.05)
      local_tHist_hreco_scaled.GetXaxis().SetTitleSize(0.05)
      local_tHist_hreco_scaled.SetTitle(sigDefnp + " Exclusive Folded Closure Test") 
      local_tHist_hreco_scaled.GetYaxis().SetTitle("Events [10^{2}/GeV]")
      local_tHist_hreco_scaled.GetXaxis().SetTitle("Reco #pi^{{0}} momentum [GeV]".format(sigDef))

      if is_closure_test:
        local_tHist_effNum_reco_scaled.SetLineColor(ROOT.kGreen+2)
        local_tHist_effNum_reco_scaled.SetMarkerColor(ROOT.kGreen+2)
        local_tHist_effNum_reco_scaled.SetFillColor(ROOT.kGreen+2)
      else:
        local_tHist_effNum_reco_scaled.SetLineColor(ROOT.kRed)
        local_tHist_effNum_reco_scaled.SetMarkerColor(ROOT.kRed)
        local_tHist_effNum_reco_scaled.SetFillColor(ROOT.kRed)
      local_tHist_effNum_reco_scaled.SetMarkerSize(0.5)
      local_tHist_effNum_reco_scaled.SetFillStyle(3545)
    
      local_tHist_genie_evtRate_reco_scaled.SetLineColor(ROOT.kRed)
      local_tHist_genie_evtRate_reco_scaled.SetFillColor(ROOT.kRed)
      local_tHist_genie_evtRate_reco_scaled.SetMarkerColor(ROOT.kRed)
      local_tHist_genie_evtRate_reco_scaled.SetMarkerSize(0.5)
      local_tHist_genie_evtRate_reco_scaled.SetFillStyle(3554)
    
      local_tHist_hreco_scaled.SetLineColor(ROOT.kCyan-3)
      local_tHist_hreco_scaled.SetFillColor(ROOT.kCyan-3)
      local_tHist_hreco_scaled.SetMarkerColor(ROOT.kCyan-3)
      local_tHist_hreco_scaled.SetMarkerSize(0.5)
      
      legend = ROOT.TLegend(0.5,0.7,0.845,0.9, "")
      legend.SetBorderSize(0)
      # If it's a closure test, draw the same GENIE reco distribution calculated using three different methods.
      if is_closure_test:
        local_tHist_hreco_scaled.Draw("HIST")
        local_tHist_genie_evtRate_reco_scaled.Draw("SAME HIST")
        local_tHist_effNum_reco_scaled.Draw("SAME")
        legend.AddEntry(local_tHist_effNum_reco_scaled,"GENIE Reconstructed (from Sbnfit)","lep")
        legend.AddEntry(local_tHist_hreco_scaled,"GENIE Reconstructed (from ResponseMaker)","f")
        legend.AddEntry(local_tHist_genie_evtRate_reco_scaled,"GENIE Truth (from Sbnfit) Folded","f")
        local_tHist_hreco_scaled.SetFillStyle(3545)
        pt = ROOT.TPaveText(0.5, 0.5, 0.845, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}_{{ResponseMaker}}/DOF: {0:.2f}/{1}".format(chi2_effNum_hreco, nBins))
        pt.AddText("#chi^{{2}}_{{Folded}}/DOF: {0:.2f}/{1}".format(chi2_effNum_genie_evtRate, nBins))
        pt.Draw()
      else:
        local_tHist_evtRate_reco_scaled.Draw()
        legend.AddEntry(local_tHist_evtRate_reco_scaled,"NuWro Fake Data","lep")
        legend.AddEntry(local_tHist_effNum_reco_scaled,"GENIE Reco Signal","f")
        local_tHist_effNum_reco_scaled.Draw("HIST SAME") 
        pt = ROOT.TPaveText(0.5, 0.5, 0.845, 0.7, "NDC")
        pt.SetBorderSize(0)
        pt.SetFillColorAlpha(0, 0)
        pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2_evtRate_reco_effNum, nBins))
        pt.Draw()
      legend.Draw()
      ptall.Draw()
   
    ### Plot of NuWro fake data folded events vs. NuWro generator prediction
    ########################################################################
    with makeEnv_TCanvas('{0}/nuwro_evtRate_folded_{1}.png'.format(plotDir, sigDef)):
      # The only difference between NuWro fake data and NuWro reco signal is the background, so we extract and only use the background uncertainty.
      local_tHist_nuwro_signal_scaled = histFile.Get("nu_uBooNE_breakdown_" + sigDefnp + "sig")
      exec("local_mHist_background = mHist_background_{0}.Clone(\"local_mHist_background\")".format(sigDef))
      local_tHist_evtRate_reco_scaled = local_mHist_evtRate_reco.GetCVHistoWithError()
      local_tMat_background_covMat = local_mHist_background.GetTotalErrorMatrix(True)
      
      # Bug: The data statistical uncertainty should be on the background, not the remaining signal. (!!!)
      for i in range(1, nBins + 1):
        local_tHist_evtRate_reco_scaled.SetBinError(i, math.sqrt(local_tMat_background_covMat[i][i]))
      local_tMat_background_covMat = local_tMat_background_covMat.GetSub(1, nBins, 1, nBins)
      chi2_evtRate_reco_nuwro_signal = calculateChi2(local_tHist_nuwro_signal_scaled, local_tHist_evtRate_reco_scaled, local_tMat_background_covMat, False)
      
      # Scale and set max y.
      local_tHist_nuwro_signal_scaled.Scale(1e-2, "width")
      local_tHist_evtRate_reco_scaled.Scale(1e-2, "width")

      maxy_nuwro_signal_scaled = local_tHist_nuwro_signal_scaled.GetMaximum()
      maxy = max(maxy_nuwro_signal_scaled, maxy_evtRate_reco)*1.1

      # Set plot formatting.
      local_tHist_nuwro_signal_scaled.GetYaxis().SetRangeUser(miny, maxy)
      local_tHist_nuwro_signal_scaled.GetYaxis().SetTitleSize(0.05)
      local_tHist_nuwro_signal_scaled.GetXaxis().SetTitleSize(0.05)
      local_tHist_nuwro_signal_scaled.SetTitle(sigDefnp + " Exclusive Folded Events")
      local_tHist_nuwro_signal_scaled.GetYaxis().SetTitle("Events [10^{2}/GeV]")
      local_tHist_nuwro_signal_scaled.GetXaxis().SetTitle("Reco #pi^{0} momentum [GeV]")
      local_tHist_nuwro_signal_scaled.SetMarkerColor(ROOT.kViolet)
      local_tHist_nuwro_signal_scaled.SetLineColor(ROOT.kViolet)
      local_tHist_nuwro_signal_scaled.SetFillColor(ROOT.kViolet)
      local_tHist_nuwro_signal_scaled.SetFillStyle(3545)
      local_tHist_nuwro_signal_scaled.SetLineWidth(1)

      local_tHist_evtRate_reco_scaled.SetMarkerSize(0.5)
      local_tHist_evtRate_reco_scaled.SetLineColor(ROOT.kBlue+2)
      local_tHist_evtRate_reco_scaled.SetMarkerStyle(ROOT.kFullCircle)
      local_tHist_evtRate_reco_scaled.SetMarkerColor(ROOT.kBlack)
      local_tHist_evtRate_reco_scaled.SetFillStyle(0)
      local_tHist_evtRate_reco_scaled.SetLineWidth(1)

      # Draw plot.
      local_tHist_nuwro_signal_scaled.Draw("HIST")
      local_tHist_evtRate_reco_scaled.Draw("SAME E0")

      legend = ROOT.TLegend(0.5,0.7,0.845,0.9, "")
      legend.SetBorderSize(0)
      legend.AddEntry(local_tHist_evtRate_reco_scaled, "NuWro Fake Data - GENIE Bkgd.", "lep")
      legend.AddEntry(local_tHist_nuwro_signal_scaled, "NuWro Reco Signal", "f")
      legend.Draw()

      pt = ROOT.TPaveText(0.55, 0.5, 0.795, 0.7, "NDC")
      pt.SetBorderSize(0)
      pt.SetFillColorAlpha(0, 0)
      pt.AddText("#chi^{{2}}/DOF: {0:.2f}/{1}".format(chi2_evtRate_reco_nuwro_signal, nBins))
      pt.Draw()
      ptall.Draw()
