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

#histFileLocation = "2021-08-26_out.root"
histFileLocation = "2021-09-06_out.root"
histFile = ROOT.TFile(histFileLocation)

plotDir = "/uboone/data/users/noza/gLEE/xsection/2021-10-18_xsec-plots_test"
#plotDir = "/uboone/data/users/finer/gLEE/NCPi0/2021-09-06_GENIE-breakout_xsec-plots"
#plotDir = "/uboone/data/users/finer/gLEE/NCPi0/2021-09-06_Detector-breakout_xsec-plots"
#plotDir = "/uboone/data/users/finer/gLEE/NCPi0/2021-08-05_xsec-plots_Zarko_MCC9Flux"
#plotDir = "2021-07-07_xsec-plots"
if not os.path.isdir(plotDir):
  print "Making plot directory {0}".format(plotDir)
  os.system( "mkdir %s" % plotDir )

#############################################################################################################
### Make simple xsec, xsec component plots ##################################################################
#############################################################################################################

for sigDef in ["2g0p","2g1p","2gnp"]:
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

#############################################################################################################
### Make pretty cross section plot ##########################################################################
#############################################################################################################
with makeEnv_TCanvas("{0}/xSection_pretty.png".format(plotDir)):

  local_xSection_2g0p.Scale(10**38)
  local_xSection_2g0p.GetYaxis().SetTitle("#sigma_{NC 1 #pi^{0}} [10^{-38} cm^{2}/Atom]")
  local_xSection_2g0p.GetYaxis().SetRangeUser(0,2)
  local_xSection_2g0p.Draw()
  
  local_xSection_2g1p.Scale(10**38)
  local_xSection_2g1p.SetMarkerColor(ROOT.kRed)
  local_xSection_2g1p.SetLineColor(ROOT.kRed)
  local_xSection_2g1p.Draw("same")

  local_xSection_2gnp.Scale(10**38)
  local_xSection_2gnp.SetMarkerColor(ROOT.kBlue)
  local_xSection_2gnp.SetLineColor(ROOT.kBlue)
  local_xSection_2gnp.Draw("same")

  ############################################################### Legend
  ######################################################################

  leg = ROOT.TLegend(0.55,0.7,0.95,0.9)
  setPlotSpecs_legend(leg)
  leg.AddEntry(local_xSection_2g0p,"2g0p","lep")
  leg.AddEntry(local_xSection_2g1p,"2g1p","lep")
  leg.AddEntry(local_xSection_2gnp,"2g(0+1)p","lep")
  leg.Draw()

#############################################################################################################
### Mark's pretty cross section plot code ###################################################################
#############################################################################################################

with makeEnv_TCanvas("{0}/xSection_pretty_markFormat.png".format(plotDir)):

  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetEndErrorSize(7)
  ROOT.gStyle.SetErrorX(0.00001)

  local_xSection_2g0p.SetLineColor(ROOT.kGreen-3)
  local_xSection_2g0p.SetMarkerColor(ROOT.kGreen-3)
  local_xSection_2g0p.SetLineWidth(2)
  local_xSection_2g0p.SetMarkerStyle(20)

  local_xSection_2g0p.Draw("E1p")
  
  local_xSection_2g0p.SetMaximum(3.0)
  local_xSection_2g0p.SetMinimum(1.0)
  local_xSection_2g0p.GetYaxis().SetTitle("#sigma_{NC 1 #pi^{0} +0 p} [10^{-38} cm^{2}/Atom]")
  local_xSection_2g0p.GetYaxis().SetTitleOffset(0.8)
  local_xSection_2g0p.GetYaxis().SetTitleSize(0.05)
  local_xSection_2g0p.GetXaxis().SetLabelOffset(999)
  local_xSection_2g0p.GetXaxis().SetLabelSize(0)
  local_xSection_2g0p.GetXaxis().SetTickLength(0.)

  local_xSection_mc_2g0p.Scale(10**38)


#    int n = 4;
#    double mcerr = E_genie_scale;
#    double mcval = xs_genie_scale;
#    std::vector<double> ymax(n,mcval+mcerr);
#    std::vector<double> ymin(n,mcval-mcerr);
#    std::vector<double> y(n,mcval);
#    std::vector<double> x = {-1,0,1,2};
#
#    TGraph *gr    = new TGraph(n,&x[0],&y[0]);
#    TGraph *grshade = new TGraph(2*n);
#
#    for (int i=0;i<n;i++) {
#        grshade->SetPoint(i,x[i],ymax[i]);
#        grshade->SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
#    }
#
#    grshade->SetFillStyle(3013);
#    grshade->SetFillColor(kRed-7);
#    grshade->SetLineColor(kRed-7);
#    grshade->Draw("f");
#    gr->SetLineWidth(3);
#    gr->SetLineColor(kRed-7);
#    gr->Draw("l");
#
#    Dat_2g0p->DrawCopy("E1p same");


  local_xSection_mc_2g0p.SetFillStyle(3013)
  local_xSection_mc_2g0p.SetFillColor(ROOT.kRed-7)
  local_xSection_mc_2g0p.SetLineColor(ROOT.kRed-7)
  #local_xSection_mc_2g0p.SetMarkerColor(ROOT.kRed)
  #local_xSection_mc_2g0p.SetMarkerStyle(1)
  local_xSection_mc_2g0p.Draw("f same")

  #lineAtOne = ROOT.TLine(horizontalAxis_lowerBound,1,horizontalAxis_upperBound,1)
  #lineAtOne.SetLineColor(ROOT.kRed)
  #lineAtOne.Draw("same")

  #local_xSection_2g1p.Scale(10**42)
  #local_xSection_2g1p.SetMarkerColor(ROOT.kRed)
  #local_xSection_2g1p.SetLineColor(ROOT.kRed)
  #local_xSection_2g1p.Draw("same")

  #local_xSection_2gnp.Scale(10**42)
  #local_xSection_2gnp.SetMarkerColor(ROOT.kBlue)
  #local_xSection_2gnp.SetLineColor(ROOT.kBlue)
  #local_xSection_2gnp.Draw("same")

  ############################################################### Legend
  ######################################################################

  #leg = ROOT.TLegend(0.11,0.7,0.89,0.89)
  leg = ROOT.TLegend(0.15,0.7,0.95,0.9)
  setPlotSpecs_legend(leg)
  leg.AddEntry(local_xSection_2g0p,"Run 1-3 Data (2#gamma0p Extracted)","l")
  leg.AddEntry(local_xSection_mc_2g0p,"Genie (G18_10a_02_11)","lf")
  #leg.SetLineColor(ROOT.kWhite)
  #leg.AddEntry(local_xSection_2g1p,"2g1p","lep")
  #leg.AddEntry(local_xSection_2gnp,"2g(0+1)p","lep")
  leg.Draw()

