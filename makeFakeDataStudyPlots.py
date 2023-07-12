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
parser.add_argument('--closureTest',help='Input file corresponds to closure test',action='store_true')
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
### Is this a closure test input? ###########################################################################
#############################################################################################################

## Prescription is slightly different for closure test input 
is_closure_test = True if p.closureTest>0 else False

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
    for histCat in ["cov_evtRate","unfolded_cov_evtRate","unfolded_evtRate","unfolded_xSection","add_smear_matrix"]:
    
      exec("tHist_{0}_{1} = histFile.Get(\"{0}_{1}\")".format(histCat,sigDef))
    
    ### Get MnvHnDs
    ################
    for histCat in ["xSection_mc"]:
    
      exec("mHist_{0}_{1} = histFile.Get(\"{0}_{1}\")".format(histCat,sigDef))
      exec("tHist_{0}_{1} = mHist_{0}_{1}.GetCVHistoWithStatError()".format(histCat,sigDef))
    
    ### Construct correlation matrices
    ##################################
    exec("tHist_corr_evtRate_{0} = tHist_cov_evtRate_{0}.Clone(\"tHist_corr_evtRate_{0}\")".format(sigDef))
    exec("tHist_unfolded_corr_evtRate_{0} = tHist_unfolded_cov_evtRate_{0}.Clone(\"tHist_unfolded_corr_evtRate_{0}\")".format(sigDef))
    
    exec("nBins = tHist_unfolded_evtRate_{0}.GetNbinsX()".format(sigDef))
    
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

#for sigDef in ["2g1p","2g0p","2gnp"]:
for sigDefnp in ["2g1p","2g0p"]:
  #for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["exclusive"]:

    if sigDefnp == "2gnp" and sigDefexcl == "exclusive":
      continue
  
    sigDef = sigDefnp+"_"+sigDefexcl

    plotter = ROOT.PlotUtils.MnvPlotter()
    plotter.SetROOT6Palette(54)
    ROOT.gStyle.SetNumberContours(200)
    
    with makeEnv_TCanvas('{0}/response_{1}.png'.format(plotDir,sigDef)) as canvas:
      
      local_tHist_Response = histFile.Get("Response_{0}".format(sigDefnp))
      local_tHist_Response.Draw("colz")
      local_tHist_Response.GetXaxis().SetTitle("True #pi^{0} momentum")
      local_tHist_Response.GetYaxis().SetTitle("Reconstructed #pi^{0} momentum")
      canvas.canvas.SetLogz()
    
    with makeEnv_TCanvas('{0}/add_smear_matrix_{1}.png'.format(plotDir,sigDef)) as canvas:
      
      local_tHist_add_smear_matrix = histFile.Get("add_smear_matrix_{0}".format(sigDef))
      local_tHist_add_smear_matrix.Draw("colz")
      local_tHist_add_smear_matrix.GetXaxis().SetTitle("True #pi^{0} momentum")
      local_tHist_add_smear_matrix.GetYaxis().SetTitle("Reconstructed #pi^{0} momentum")
      #canvas.canvas.SetLogz()
    
    with makeEnv_TCanvas('{0}/cov_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_cov_evtRate_{0} = tHist_cov_evtRate_{0}.Clone(\"local_tHist_cov_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_cov_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_cov_evtRate_{0}.GetXaxis().SetTitle(\"Reconstructed #pi^{0} momentum\")".format(sigDef))
      exec("local_tHist_cov_evtRate_{0}.GetYaxis().SetTitle(\"Reconstructed #pi^{0} momentum\")".format(sigDef))
      canvas.canvas.SetLogz()
    
    with makeEnv_TCanvas('{0}/unfolded_cov_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_unfolded_cov_evtRate_{0} = tHist_unfolded_cov_evtRate_{0}.Clone(\"local_tHist_unfolded_cov_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_unfolded_cov_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_unfolded_cov_evtRate_{0}.GetXaxis().SetTitle(\"True #pi^{0} momentum\")".format(sigDef))
      exec("local_tHist_unfolded_cov_evtRate_{0}.GetYaxis().SetTitle(\"True #pi^{0} momentum\")".format(sigDef))
      canvas.canvas.SetLogz()
    
    ## Change to different color palette for correlation matrix
    plotter.SetROOT6Palette(87)
    ROOT.gStyle.SetNumberContours(200)
    
    with makeEnv_TCanvas('{0}/corr_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_corr_evtRate_{0} = tHist_corr_evtRate_{0}.Clone(\"local_tHist_corr_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_corr_evtRate_{0}.GetZaxis().SetRangeUser(-1,1)".format(sigDef))
      exec("local_tHist_corr_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_corr_evtRate_{0}.GetXaxis().SetTitle(\"Reconstructed #pi^{0} momentum\")".format(sigDef))
      exec("local_tHist_corr_evtRate_{0}.GetYaxis().SetTitle(\"Reconstructed #pi^{0} momentum\")".format(sigDef))
    
    with makeEnv_TCanvas('{0}/unfolded_corr_evtRate_{1}.png'.format(plotDir,sigDef)) as canvas:
    
      exec("local_tHist_unfolded_corr_evtRate_{0} = tHist_unfolded_corr_evtRate_{0}.Clone(\"local_tHist_unfolded_corr_evtRate_{0}\")".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.GetZaxis().SetRangeUser(-1,1)".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.Draw(\"colz\")".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.GetXaxis().SetTitle(\"True #pi^{0} momentum\")".format(sigDef))
      exec("local_tHist_unfolded_corr_evtRate_{0}.GetYaxis().SetTitle(\"True #pi^{0} momentum\")".format(sigDef))
    
#############################################################################################################
### Make fake data study comparison plots ###################################################################
#############################################################################################################

#for sigDef in ["2g1p","2g0p","2gnp"]:
for sigDefnp in ["2g1p","2g0p"]:
  #for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["exclusive"]:

    if sigDefnp == "2gnp" and sigDefexcl == "exclusive":
      continue

    sigDef = sigDefnp+"_"+sigDefexcl

    with makeEnv_TCanvas('{0}/fakedatavsgenie_xSection_{1}.png'.format(plotDir,sigDef)):
    
      exec("local_tHist_unfolded_xSection_scaled = tHist_unfolded_xSection_{0}.Clone(\"local_tHist_unfolded_xSection_{0}\")".format(sigDef))
      local_tHist_unfolded_xSection_scaled.Scale(1e38,"width") 
      exec("local_tHist_xSection_mc_scaled = tHist_xSection_mc_{0}.Clone(\"local_tHist_xSection_mc_{0}\")".format(sigDef))
      local_tHist_xSection_mc_scaled.Scale(1e38,"width") 
    
      local_tHist_unfolded_xSection_scaled.SetMarkerSize(0.5)
    
      local_tHist_unfolded_xSection_scaled.GetYaxis().SetRangeUser(0,2.8)
      local_tHist_unfolded_xSection_scaled.GetYaxis().SetTitleSize(0.05)
      local_tHist_unfolded_xSection_scaled.GetXaxis().SetTitleSize(0.05)
      local_tHist_unfolded_xSection_scaled.GetYaxis().SetTitle("#sigma_{NC 1 #pi^{0}}[10^{-38} cm^{2}/Atom/GeV]")
      local_tHist_unfolded_xSection_scaled.GetXaxis().SetTitle("True #pi^{{0}} momentum (GeV) [{0}]".format(sigDef))
    
      local_tHist_unfolded_xSection_scaled.Draw()
    
      local_tHist_xSection_mc_scaled.SetLineColor(ROOT.kRed)
      local_tHist_xSection_mc_scaled.SetMarkerColor(ROOT.kRed)
      local_tHist_xSection_mc_scaled.SetMarkerSize(0.5)
      local_tHist_xSection_mc_scaled.Draw("SAME")
    
      legend = ROOT.TLegend(0.55,0.7,0.8,0.85, "")
      legend.SetBorderSize(0);
      legend.AddEntry(local_tHist_unfolded_xSection_scaled,"NuWro Fake Data","lep")
      legend.AddEntry(local_tHist_xSection_mc_scaled,"GENIE Prediction","lep")
      legend.Draw()
    
    ## Import NuWro Truth
    NuWroTruthFile = ROOT.TFile("/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/MajorMerge_GGE_mark/working_dir/ToTH1D/NuWro_FakeData_Generation_NEW/NuWro_Apr2023_v2_CV.SBNspec.root")
    
    with makeEnv_TCanvas('{0}/fakedatavsgenie_evtRate_{1}.png'.format(plotDir,sigDef)):
    
      ## Create local copies of all needed hists
      exec("local_tHist_unfolded_evtRate_scaled = tHist_unfolded_evtRate_{0}.Clone(\"local_tHist_unfolded_evtRate_{0}\")".format(sigDef))
      local_tHist_genie_evtRate_smeared_scaled = histFile.Get("smeared_true_signal_{0}".format(sigDef))
      local_tHist_genie_evtRate_scaled = histFile.Get("prior_true_signal_{0}".format(sigDef))
      local_mHist_effDenom_scaled = histFile.Get("effDenom_{0}".format(sigDef))
      local_tHist_effDenom_scaled =local_mHist_effDenom_scaled.GetCVHistoWithError()
      exec("local_tHist_nuwro_truth_scaled = NuWroTruthFile.Get(\"nu_uBooNE_denom_{0}\")".format(sigDefnp))
    
      for i in range(1,nBins+1):## Loop over bins
        exec("cov_binContent = tHist_unfolded_cov_evtRate_{0}.GetBinContent(i,i)".format(sigDef))
        local_tHist_unfolded_evtRate_scaled.SetBinError(i,math.sqrt(cov_binContent))
    
      ## Scale and bin-width-normalize all hists (and POT normalize NuWro truth)
      local_tHist_unfolded_evtRate_scaled.Scale(1e-3,"width")
      local_tHist_genie_evtRate_smeared_scaled.Scale(1e-3,"width")
      local_tHist_genie_evtRate_scaled.Scale(1e-3,"width")
      local_tHist_effDenom_scaled.Scale(1e-3,"width")
      local_tHist_nuwro_truth_scaled.Scale(1e-3,"width")
      local_tHist_nuwro_truth_scaled.Scale(6.7873/3.0041393) ## POT normalization
     
      ## Set plot formatting 
      local_tHist_unfolded_evtRate_scaled.SetMarkerSize(0.5)
      if is_closure_test:
        local_tHist_unfolded_evtRate_scaled.SetMarkerColor(ROOT.kCyan-3)
        local_tHist_unfolded_evtRate_scaled.SetLineColor(ROOT.kCyan-3)
    
      local_tHist_unfolded_evtRate_scaled.GetYaxis().SetRangeUser(0,18)
      local_tHist_unfolded_evtRate_scaled.GetYaxis().SetTitleSize(0.05)
      local_tHist_unfolded_evtRate_scaled.GetXaxis().SetTitleSize(0.05)
      local_tHist_unfolded_evtRate_scaled.GetYaxis().SetTitle("# Events [10^{3}/GeV]")
      local_tHist_unfolded_evtRate_scaled.GetXaxis().SetTitle("True #pi^{{0}} momentum (GeV) [{0}]".format(sigDef))
    
      local_tHist_unfolded_evtRate_scaled.Draw()
    
      local_tHist_effDenom_scaled.SetLineColor(ROOT.kGreen+2)
      local_tHist_effDenom_scaled.SetMarkerColor(ROOT.kGreen+2)
      local_tHist_effDenom_scaled.SetMarkerSize(0.5)
      local_tHist_effDenom_scaled.Draw("SAME")
    
      if not is_closure_test:
        local_tHist_genie_evtRate_scaled.SetMarkerColor(ROOT.kRed)
        local_tHist_genie_evtRate_scaled.SetLineColor(ROOT.kRed)
        local_tHist_genie_evtRate_scaled.Draw("SAME")
        
        local_tHist_nuwro_truth_scaled.SetMarkerColor(ROOT.kViolet)
        local_tHist_nuwro_truth_scaled.SetLineColor(ROOT.kViolet)
        local_tHist_nuwro_truth_scaled.Draw("SAME")
    
      local_tHist_genie_evtRate_smeared_scaled.SetLineColor(ROOT.kCyan-3)
      local_tHist_genie_evtRate_smeared_scaled.SetMarkerColor(ROOT.kCyan-3)
      local_tHist_genie_evtRate_smeared_scaled.SetMarkerSize(0.5)
      local_tHist_genie_evtRate_smeared_scaled.Draw("SAME")
    
      legend = ROOT.TLegend(0.5,0.7,0.845,0.9, "")
      legend.SetBorderSize(0);
      if is_closure_test:
        legend.AddEntry(local_tHist_unfolded_evtRate_scaled,"GENIE Fake Data (closure test)","lep")
        legend.AddEntry(local_tHist_genie_evtRate_smeared_scaled,"GENIE Prediction","lep")
      else:
        legend.AddEntry(local_tHist_unfolded_evtRate_scaled,"NuWro Fake Data","lep")
        legend.AddEntry(local_tHist_nuwro_truth_scaled,"NuWro Truth","lep")
        legend.AddEntry(local_tHist_effDenom_scaled,"effDenom","lep")
        legend.AddEntry(local_tHist_genie_evtRate_scaled,"GENIE Prediction","lep")
        legend.AddEntry(local_tHist_genie_evtRate_smeared_scaled,"GENIE Prediction, smeared","lep")
      legend.Draw()
    
    with makeEnv_TCanvas('{0}/fakedatavsgenie_evtRate_folded_{1}.png'.format(plotDir,sigDef)):
    
      ## Create local copies of all needed hists
      exec("local_mHist_evtRate_reco = histFile.Get(\"evtRate_{0}\")".format(sigDef))
      local_tHist_evtRate_reco_scaled = local_mHist_evtRate_reco.GetCVHistoWithError()
      local_tHist_evtRate_reco_scaled.Scale(1e-3,"width")
    
      exec("local_tHist_genie_evtRate_reco_scaled = histFile.Get(\"prior_reco_signal_{0}\")".format(sigDef))
      local_tHist_genie_evtRate_reco_scaled.Scale(1e-3,"width")
     
      exec("local_tHist_hreco_scaled = histFile.Get(\"hreco_{0}\")".format(sigDefnp))
      local_tHist_hreco_scaled.Scale(1e-3,"width")
     
      exec("local_mHist_effNum_reco = histFile.Get(\"effNum_reco_{0}\")".format(sigDef))
      local_tHist_effNum_reco_scaled = local_mHist_effNum_reco.GetCVHistoWithError() 
      local_tHist_effNum_reco_scaled.Scale(1e-3,"width")

      ## Set plot formatting 
      local_tHist_evtRate_reco_scaled.SetMarkerSize(0.5)
      if is_closure_test:
        local_tHist_evtRate_reco_scaled.SetMarkerColor(ROOT.kCyan-3)
        local_tHist_evtRate_reco_scaled.SetLineColor(ROOT.kCyan-3)
    
      local_tHist_evtRate_reco_scaled.GetYaxis().SetRangeUser(0,4)
      local_tHist_evtRate_reco_scaled.GetYaxis().SetTitleSize(0.05)
      local_tHist_evtRate_reco_scaled.GetXaxis().SetTitleSize(0.05)
      local_tHist_evtRate_reco_scaled.GetYaxis().SetTitle("# Events [10^{3}/GeV]")
      local_tHist_evtRate_reco_scaled.GetXaxis().SetTitle("Reco #pi^{{0}} momentum (GeV) [{0}]".format(sigDef))
    
      local_tHist_evtRate_reco_scaled.Draw()
    
      local_tHist_effNum_reco_scaled.SetLineColor(ROOT.kGreen+2)
      local_tHist_effNum_reco_scaled.SetMarkerColor(ROOT.kGreen+2)
      local_tHist_effNum_reco_scaled.SetMarkerSize(0.5)
      local_tHist_effNum_reco_scaled.Draw("SAME")
    
      local_tHist_genie_evtRate_reco_scaled.SetLineColor(ROOT.kRed)
      local_tHist_genie_evtRate_reco_scaled.SetMarkerColor(ROOT.kRed)
      local_tHist_genie_evtRate_reco_scaled.SetMarkerSize(0.5)
      local_tHist_genie_evtRate_reco_scaled.Draw("SAME")
    
      local_tHist_hreco_scaled.SetLineColor(ROOT.kCyan-3)
      local_tHist_hreco_scaled.SetMarkerColor(ROOT.kCyan-3)
      local_tHist_hreco_scaled.SetMarkerSize(0.5)
      local_tHist_hreco_scaled.Draw("SAME")
    
      legend = ROOT.TLegend(0.5,0.7,0.845,0.9, "")
      legend.SetBorderSize(0);
      if is_closure_test:
        legend.AddEntry(local_tHist_evtRate_reco_scaled,"GENIE Fake Data (closure test)","lep")
      else:
        legend.AddEntry(local_tHist_evtRate_reco_scaled,"NuWro Fake Data","lep")
      legend.AddEntry(local_tHist_genie_evtRate_reco_scaled,"GENIE Prediction (derived)","lep")
      legend.AddEntry(local_tHist_hreco_scaled,"GENIE Prediction (raw)","lep")
      legend.AddEntry(local_tHist_effNum_reco_scaled,"effNum_reco","lep")
      legend.Draw()
    
