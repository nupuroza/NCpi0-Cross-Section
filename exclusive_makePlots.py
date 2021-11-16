### Script to take make plots, etc.
### using the MINERvA Analysis Toolkit

import ROOT
from plottingClasses import *
from errorMaps_gLEE import *
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

histFileLocation = "2021-11-15_out_exclusive_20MeV.root"
#histFileLocation = "2021-11-15_out_exclusive_50MeV.root"
histFile = ROOT.TFile(histFileLocation)

plotDir = "/uboone/data/users/finer/gLEE/NCPi0/2021-11-15_exclusive_20MeV_xsec-plots"
#plotDir = "/uboone/data/users/finer/gLEE/NCPi0/2021-11-15_exclusive_50MeV_xsec-plots"
if not os.path.isdir(plotDir):
  print "Making plot directory {0}".format(plotDir)
  os.system( "mkdir %s" % plotDir )

#############################################################################################################
### Make simple xsec, xsec component plots ##################################################################
#############################################################################################################

for sigDef in ["2g0p","2g1p"]:
  for histCat in ["effDenom","effNum","eff","data_selected","BNB_ext","background","evtRate","xSection","xSection_mc","POT"]:

    exec("{0}_{1} = histFile.Get(\"{0}_{1}\")".format(histCat,sigDef))

    with makeEnv_TCanvas("{0}/{1}_{2}.png".format(plotDir,histCat,sigDef)):
      exec("local_{0}_{1} = {0}_{1}.GetCVHistoWithError()".format(histCat,sigDef))
      ## Set horizontal axis label
      exec("local_{0}_{1}.GetXaxis().SetTitle(\"cos(#theta_{{#pi^{{0}}}})\")".format(histCat,sigDef))
      exec("local_{0}_{1}.GetXaxis().SetTitleSize(0.05)".format(histCat,sigDef))
      exec("local_{0}_{1}.GetYaxis().SetTitleSize(0.05)".format(histCat,sigDef))
      ## Set vertical axis label
      if histCat == "eff":
        exec("local_{0}_{1}.GetYaxis().SetTitle(\"Efficiency\")".format(histCat,sigDef))
      elif histCat == "xSection" or histCat == "xSection_mc":
        exec("local_{0}_{1}.GetYaxis().SetTitle(\"#sigma_{{NC 1 #pi^{{0}}}} [10^{{-38}} cm^{{2}}/Atom]\")".format(histCat,sigDef))
      else:
        exec("local_{0}_{1}.GetYaxis().SetTitle(\"Number of Events\")".format(histCat,sigDef))
      exec("local_{0}_{1}.Draw()".format(histCat,sigDef))

    if not histCat == "xSection_mc":
      for flux_uni in FLUX_SYSTS:
        exec("local_{0}_{1} = {0}_{1}.GetVertErrorBand(\"{2}\")".format(histCat,sigDef,flux_uni))
        with makeEnv_TCanvas("{0}/breakout_{1}_{2}_{3}.png".format(plotDir,flux_uni,histCat,sigDef)):
          exec("local_{0}_{1}.DrawAll(\"\",True)".format(histCat,sigDef))

    with makeEnv_TCanvas("{0}/errorSummary_{1}_{2}.png".format(plotDir,histCat,sigDef)):
      exec("localDrawErrorSummary(plotter,{0}_{1})".format(histCat,sigDef))

    ## Print out value and error for Mark to package into table
    exec("local_hist = {0}_{1}.GetCVHistoWithError()".format(histCat,sigDef))
    cv_val = local_hist.GetBinContent(1)
    err_val = local_hist.GetBinError(1)
    exec("local_hist_statError = {0}_{1}.GetCVHistoWithStatError()".format(histCat,sigDef))
    exec("local_hist_systError = {0}_{1}.GetCVHistoWithError(False)".format(histCat,sigDef))
    err_val_stat = local_hist_statError.GetBinError(1)
    err_val_syst = local_hist_systError.GetBinError(1)
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
    exec("local_{0} = {0}.GetCVHistoWithError()".format(histCat))
    ## Set horizontal axis label
    if histCat == "flux":
      exec("local_{0}.GetXaxis().SetTitle(\"E_{{#nu}} (GeV))\")".format(histCat))
    else:
      exec("local_{0}.GetXaxis().SetTitle(\"cos(#theta_{{#pi^{{0}}}})\")".format(histCat))
    exec("local_{0}.GetXaxis().SetTitleSize(0.05)".format(histCat))
    exec("local_{0}.GetYaxis().SetTitleSize(0.05)".format(histCat))
    ## Set vertical axis label
    exec("local_{0}.GetYaxis().SetTitle(\"#nu/POT/cm^{{2}}\")".format(histCat))
    exec("local_{0}.Draw()".format(histCat))

  for flux_uni in FLUX_SYSTS:
    exec("local_{0} = {0}.GetVertErrorBand(\"{1}\")".format(histCat,flux_uni))
    with makeEnv_TCanvas("{0}/breakout_{1}_{2}.png".format(plotDir,flux_uni,histCat)):
      exec("local_{0}.DrawAll(\"\",True)".format(histCat))

  with makeEnv_TCanvas("{0}/errorSummary_{1}.png".format(plotDir,histCat)):
    exec("localDrawErrorSummary(plotter,{0})".format(histCat))

#############################################################################################################
### Make downstream breakout plots ##########################################################################
#############################################################################################################

for syst_uni in ["piplus_PrimaryHadronSWCentralSplineVariation","All_UBGenie"]:

  exec("local_xSection_2g1p_breakout = xSection_2g1p.GetVertErrorBand(\"{0}\")".format(syst_uni))
  with makeEnv_TCanvas("{0}/breakout_{1}_xSection.png".format(plotDir,syst_uni)):
    exec("local_xSection_2g1p_breakout.DrawAll(\"\",True)".format(sigDef))

