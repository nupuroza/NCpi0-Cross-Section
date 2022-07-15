### Script to take make plots, etc.
### using the MINERvA Analysis Toolkit

import ROOT
import datetime as dt
import argparse
from plottingClasses import *
from errorMaps import *
#from errorMaps_gLEE_GENIE_breakout import *
#from errorMaps_gLEE_Detector_breakout import *

## Set ROOT to batch mode
ROOT.gROOT.SetBatch()

## #Load and implement Phil's plot style header file
#ROOT.gROOT.ProcessLine(".L ../style/myPlotStyle.h")
#ROOT.myPlotStyle()

## This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

## Setup MnvPlotter, which has all of the plotting utilities
plotter = ROOT.PlotUtils.MnvPlotter()

## Manually override default error summary groups
plotter.error_summary_group_map.clear()
for group in error_bands:
  for error in error_bands[group]:
    plotter.error_summary_group_map[group].push_back(error)
plotter.SetLegendNColumns(2)

#############################################################################################################
### Systematic Universes ####################################################################################
#############################################################################################################

## List of flux systematics
FLUX_SYSTS = [
  ## [syst_name,universe_prefix,n_universes]
  "expskin_FluxUnisim",
  "horncurrent_FluxUnisim",
  "kminus_PrimaryHadronNormalization",
  "kplus_PrimaryHadronFeynmanScaling",
  "kzero_PrimaryHadronSanfordWang",
  "nucleoninexsec_FluxUnisim",
  "nucleonqexsec_FluxUnisim",
  "nucleontotxsec_FluxUnisim",
  "piminus_PrimaryHadronSWCentralSplineVariation",
  "pioninexsec_FluxUnisim",
  "pionqexsec_FluxUnisim",
  "piontotxsec_FluxUnisim",
  "piplus_PrimaryHadronSWCentralSplineVariation"
]

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################
parser = argparse.ArgumentParser(description='Script to make cross-section plots using the MINERvA Analysis Toolkit')
parser.add_argument('in_dir', help='Path to input directory', type=str,nargs='?')
parser.add_argument('in_date', help='Creation date of input file (yyyy-mm-dd). Defaults to a file dated today if it exists', type=str,nargs='?')
parser.add_argument('out_dir', help='Path to output directory. Defaults to input directory', type=str,nargs='?')
p = parser.parse_args()

## If in_dir is not provided, exit
if p.in_dir < 0:
  print "ERROR: Input directory argument not provided"
  parser.print_help()
  exit(1)

## If in_date is not provided, search file created today
if p.in_date < 0:
  histFileLocation = p.in_dir+"/{0}_out_final.root".format(dt.date.today())
  if not os.path.exists(histFileLocation):
    print "ERROR: An input ROOT file created today does not exist. Specify input date argument"
    parser.print_help()
    exit(1)
  else:
    print "This is the input file I'm opening: {0}".format(histFileLocation)
    histFile = ROOT.TFile(histFileLocation)
else: 
  histFileLocation = p.in_dir+"/"+p.in_date+"_out_final.root"
  if not os.path.exists(histFileLocation):
    print "ERROR: An input ROOT file created on "+p.in_date+" does not exist"
    parser.print_help()
    exit(1)
  else:
    print "This is the input file I'm opening: {0}".format(histFileLocation)
    histFile = ROOT.TFile(histFileLocation)  

## If out_dir is not provided, default to using in_dir
if p.out_dir < 0:
  plotDir = p.in_dir+"/{0}_xsec-plots".format(dt.date.today())
else:
  ## Create p.out_dir if it doesn't exist
  if not os.path.isdir(p.out_dir):
    print "Making plot directory {0}".format(p.out_dir)
    os.system( "mkdir %s" % p.out_dir )
  plotDir = p.out_dir+"/{0}_xsec-plots".format(dt.date.today())

## Create output directory if it doesn't exist
if not os.path.isdir(plotDir):
  print "Making plot directory {0}".format(plotDir)
  os.system( "mkdir %s" % plotDir )

#############################################################################################################
### Make simple xsec, xsec component plots ##################################################################
#############################################################################################################

#for sigDef in ["2g0p","2g1p","2gnp"]:
for sigDef in ["2g1p"]:  
  ## for sigDefexcl in ["exclusive", "inclusive"]:
  for sigDefexcl in ["inclusive"]:
  
    if sigDef == "2gnp" and sigDefexcl == "exclusive":
      continue

    else:
      for histCat in ["background","evtRate","effNum","effDenom", "eff", "unfolded_xSection", "xSection_mc","data_selected", "unfolded_evtRate"]:
        exec("{0}_{1}_{2} = histFile.Get(\"{0}_{1}_{2}\")".format(histCat,sigDef,sigDefexcl))

        with makeEnv_TCanvas("{0}/{1}_{2}_{3}.png".format(plotDir,histCat,sigDef,sigDefexcl)):
          exec("tHist_{0}_{1}_{2} = {0}_{1}_{2}.GetCVHistoWithError()".format(histCat,sigDef,sigDefexcl))
          ## Set axis label sizes
          exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitleSize(0.05)".format(histCat,sigDef,sigDefexcl))
          exec("tHist_{0}_{1}_{2}.GetYaxis().SetTitleSize(0.05)".format(histCat,sigDef,sigDefexcl))
          ## Set vertical axis label
          if histCat == "eff":
            exec("tHist_{0}_{1}_{2}.GetYaxis().SetTitle(\"Efficiency/GeV\")".format(histCat,sigDef,sigDefexcl))
          elif histCat == "unfolded_xSection" or histCat == "xSection_mc":
            exec("tHist_{0}_{1}_{2}.GetYaxis().SetTitle(\"#sigma_{{NC 1 #pi^{{0}}}}[10^{{-38}} cm^{{2}}/Atom]/GeV\")".format(histCat,sigDef,sigDefexcl))
            exec("tHist_{0}_{1}_{2}.Scale(10**38)".format(histCat,sigDef,sigDefexcl))
          else:
            exec("tHist_{0}_{1}_{2}.GetYaxis().SetTitle(\"Number of Events/GeV\")".format(histCat,sigDef,sigDefexcl))
          ## Set horizontal axis label
          if histCat == "background" or histCat == "evtRate" or histCat == "data_selected" or histCat == "BNB_ext":
            exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitle(\"Reconstructed #pi^{{0}} momentum\")".format(histCat,sigDef,sigDefexcl))
          else: 
            exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitle(\"#pi^{{0}} momentum\")".format(histCat,sigDef,sigDefexcl))
          exec("tHist_{0}_{1}_{2}.Scale(1.0,\"width\")".format(histCat,sigDef,sigDefexcl))
          exec("tHist_{0}_{1}_{2}.Draw()".format(histCat,sigDef,sigDefexcl))

        if histCat == "background" or histCat == "evtRate" or histCat == "data_selected" or histCat == "BNB_ext":
          with makeEnv_TCanvas("{0}/errorSummary_{1}_{2}_{3}.png".format(plotDir,histCat,sigDef,sigDefexcl)):
            exec("localDrawErrorSummary(plotter,{0}_{1}_{2},\"Reconstructed #pi^{{0}} momentum\")".format(histCat,sigDef,sigDefexcl))
        else:
          with makeEnv_TCanvas("{0}/errorSummary_{1}_{2}_{3}.png".format(plotDir,histCat,sigDef,sigDefexcl)):
            exec("localDrawErrorSummary(plotter,{0}_{1}_{2},\"#pi^{{0}} momentum\")".format(histCat,sigDef,sigDefexcl))

        ## Print out value and error for Mark to package into table
        exec("tHist = {0}_{1}_{2}.GetCVHistoWithError()".format(histCat,sigDef,sigDefexcl))
        cv_val = tHist.GetBinContent(1)
        err_val = tHist.GetBinError(1)
        exec("tHist_statError = {0}_{1}_{2}.GetCVHistoWithStatError()".format(histCat,sigDef,sigDefexcl))
        exec("tHist_systError = {0}_{1}_{2}.GetCVHistoWithError(False)".format(histCat,sigDef,sigDefexcl))
        err_val_stat = tHist_statError.GetBinError(1)
        err_val_syst = tHist_systError.GetBinError(1)
        print "sigDef: {0}\thistCat: {1}_{2}\tcv_val: {3}\terr_val_stat: {4}\terr_val_syst: {5}\terr_val_tot: {6}".format(histCat,sigDef,sigDefexcl,cv_val,err_val_stat,err_val_syst,err_val)
		
  for histCat in ["POT"]:
    print "I'm producing plots for histCat: {0}".format(histCat)
    exec("{0}_{1} = histFile.Get(\"{0}_{1}\")".format(histCat,sigDef))

    with makeEnv_TCanvas("{0}/{1}_{2}.png".format(plotDir,histCat,sigDef)):
      exec("tHist_{0}_{1} = {0}_{1}.GetCVHistoWithError()".format(histCat,sigDef))
      ## Set axis label sizes
      exec("tHist_{0}_{1}.GetXaxis().SetTitleSize(0.05)".format(histCat,sigDef))
      exec("tHist_{0}_{1}.GetYaxis().SetTitleSize(0.05)".format(histCat,sigDef)) 
      ## Set vertical axis label
      exec("tHist_{0}_{1}.GetYaxis().SetTitle(\"Number of Events/GeV\")".format(histCat,sigDef))
      exec("tHist_{0}_{1}.Scale(1.0,\"width\")".format(histCat,sigDef))
      exec("tHist_{0}_{1}.Draw()".format(histCat,sigDef))
      ## Set horizontal axis label
      exec("tHist_{0}_{1}.GetXaxis().SetTitle(\"#pi^{{0}} momentum\")".format(histCat,sigDef))

    with makeEnv_TCanvas("{0}/errorSummary_{1}_{2}.png".format(plotDir,histCat,sigDef)):
      exec("localDrawErrorSummary(plotter,{0}_{1},\"#pi^{{0}} momentum\")".format(histCat,sigDef))

    ## Print out value and error for Mark to package into table
    exec("tHist = {0}_{1}.GetCVHistoWithError()".format(histCat,sigDef))
    cv_val = tHist.GetBinContent(1)
    err_val = tHist.GetBinError(1)
    exec("tHist_statError = {0}_{1}.GetCVHistoWithStatError()".format(histCat,sigDef))
    exec("tHist_systError = {0}_{1}.GetCVHistoWithError(False)".format(histCat,sigDef))
    err_val_stat = tHist_statError.GetBinError(1)
    err_val_syst = tHist_systError.GetBinError(1)
    print "sigDef: {0}\thistCat: {1}\tcv_val: {2}\terr_val_stat: {3}\terr_val_syst: {4}\terr_val_tot: {5}".format(histCat,sigDef,cv_val,err_val_stat,err_val_syst,err_val)
    
#############################################################################################################
### Make simple flux plots ##################################################################################
#############################################################################################################

flux = histFile.Get("flux")
flux.GetXaxis().SetRangeUser(0,5)

integratedFlux = histFile.Get("integratedFlux")
nTargets = histFile.Get("nTargets")

for histCat in ["flux","integratedFlux","nTargets"]:

  with makeEnv_TCanvas("{0}/{1}.png".format(plotDir,histCat)) as can:
    exec("tHist_{0} = {0}.GetCVHistoWithError()".format(histCat))
    ## Set axes labels
    exec("tHist_{0}.GetXaxis().SetTitleSize(0.05)".format(histCat))
    exec("tHist_{0}.GetYaxis().SetTitleSize(0.05)".format(histCat))
    if histCat == "flux":
      exec("tHist_{0}.GetXaxis().SetTitle(\"E_{{#nu}}(GeV)\")".format(histCat))
      exec("tHist_{0}.GetYaxis().SetTitle(\"#nu/10^{{-6}} POT/cm^{{2}}/GeV\")".format(histCat))
      exec("tHist_{0}.Scale(10**-6)".format(histCat))
    elif histCat == "integratedFlux":
      exec("tHist_{0}.GetXaxis().SetTitle(\"#pi^{{0}} momentum\")".format(histCat))
      exec("tHist_{0}.GetYaxis().SetTitle(\"#nu/10^{{9}} POT/cm^{{2}}/GeV\")".format(histCat))
      exec("tHist_{0}.Scale(10**9)".format(histCat))
    elif histCat == "nTargets":
      exec("tHist_{0}.GetXaxis().SetTitle(\"#pi^{{0}} momentum\")".format(histCat))
      exec("tHist_{0}.GetYaxis().SetTitle(\"Number of Targets[10^{{28}} atoms]/GeV\")".format(histCat))
      exec("tHist_{0}.Scale(10**-28)".format(histCat)) 
    exec("tHist_{0}.Scale(1.0,\"width\")".format(histCat))
    exec("tHist_{0}.Draw()".format(histCat))

  with makeEnv_TCanvas("{0}/errorSummary_{1}.png".format(plotDir,histCat)):
    if histCat == "flux":
      exec("localDrawErrorSummary(plotter,{0},\"#nu/10^{{-6}} POT/cm^{{2}}/GeV\")".format(histCat))
    else:  
      exec("localDrawErrorSummary(plotter,{0},\"#pi^{{0}} momentum\")".format(histCat))

#############################################################################################################
### Plot response and event rate covariance #################################################################
#############################################################################################################
plotter = ROOT.PlotUtils.MnvPlotter()
plotter.SetROOT6Palette(54)
ROOT.gStyle.SetNumberContours(200)

with makeEnv_TCanvas('{0}/cov_evtRate_2g1p_inclusive.png'.format(plotDir)) as canvas:
  tmp = histFile.Get('cov_evtRate_2g1p_inclusive')
  tmp.Draw("colz")
  tmp.GetXaxis().SetTitle("Reconstructed bins")
  tmp.GetYaxis().SetTitle("Reconstructed bins")
  canvas.canvas.SetLogz()

with makeEnv_TCanvas('{0}/response_2g1p_inclusive.png'.format(plotDir)) as canvas:
  tmp = histFile.Get('response_2g1p_inclusive')
  tmp.Draw("colz")
  tmp.GetXaxis().SetTitle("True #pi^{0} momentum")
  tmp.GetYaxis().SetTitle("Reconstructed #pi^{0} momentum")
  canvas.canvas.SetLogz()

with makeEnv_TCanvas('{0}/unfolded_cov_evtRate_2g1p_inclusive.png'.format(plotDir)) as canvas:
  tmp = histFile.Get('unfolded_cov_evtRate_2g1p_inclusive')
  tmp.Draw("colz")
  tmp.GetXaxis().SetTitle("True bins")
  tmp.GetYaxis().SetTitle("True bins")
  canvas.canvas.SetLogz()
