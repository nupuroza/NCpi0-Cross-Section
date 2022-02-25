### Script to take make plots, etc.
### using the MINERvA Analysis Toolkit

import ROOT
import datetime as dt
import argparse
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

## Define localDrawErrorSummary with x-axis arg
## TO-DO: Fix in plottingClasses.py and remove
def localDrawErrorSummary( plotter , hist , xaxis_label ):

  # Not working at the moment; meant to block out the 1-2 GeV Enu bin  
  box = ROOT.TBox(1,0,2,0.16)
  box.SetFillColor(ROOT.kGray)
  box.SetFillStyle(3001)
  box.Draw()

  hist.GetXaxis().SetTitle(xaxis_label)
  hist.GetXaxis().SetTitleSize(0.05)

  plotter.DrawErrorSummary(hist,"TL",True,True,0.00001,False,"",True,"",True)

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
parser.add_argument('tag', help='String to append to end of output directory name. e.g. "test"', type=str,nargs='?')
p = parser.parse_args()

## If in_dir is not provided, exit
if p.in_dir < 0:
  print "ERROR: Input directory argument not provided"
  parser.print_help()
  exit(1)

## If in_date is not provided, search file created today
if p.in_date < 0:
  histFileLocation = p.in_dir+"/{0}_out.root".format(dt.date.today())
  if not os.path.exists(histFileLocation):
    print "ERROR: An input ROOT file created today does not exist. Specify input date argument"
    parser.print_help()
    exit(1)
  else:
    print "This is the input file I'm opening: {0}".format(histFileLocation)
    histFile = ROOT.TFile(histFileLocation)
else: 
  histFileLocation = p.in_dir+"/"+p.in_date+"_out.root"
  if not os.path.exists(histFileLocation):
    print "ERROR: An input ROOT file created on "+p.in_date+" does not exist"
    parser.print_help()
    exit(1)
  else:
    print "This is the input file I'm opening: {0}".format(histFileLocation)
    histFile = ROOT.TFile(histFileLocation)  

## If tag is provided, prepend it with an underscore
if p.tag >= 0:
  tag = "_{0}".format(p.tag)
else:
  tag = ''

## If out_dir is not provided, default to using in_dir
if p.out_dir < 0:
  plotDir = p.in_dir+"/{0}_xsec-plots{1}".format(dt.date.today(),tag)
else:
  ## Create p.out_dir if it doesn't exist
  if not os.path.isdir(p.out_dir):
    print "Making plot directory {0}".format(p.out_dir)
    os.system( "mkdir %s" % p.out_dir )
  plotDir = p.out_dir+"/{0}_xsec-plots{1}".format(dt.date.today(),tag)

## Create output directory if it doesn't exist
if not os.path.isdir(plotDir):
  print "Making plot directory {0}".format(plotDir)
  os.system( "mkdir %s" % plotDir )

#############################################################################################################
### Make simple xsec, xsec component plots ##################################################################
#############################################################################################################

for sigDef in ["2g0p","2g1p","2gnp"]:
  
  ## OMIT EXCL:
  ## for sigDefexcl in ["exclusive", "inclusive"]:
  for sigDefexcl in ["inclusive"]:
  
    if sigDef == "2gnp" and sigDefexcl == "exclusive":
      continue

    else:
      for histCat in ["background","evtRate","effNum","effDenom", "eff", "xSection", "xSection_mc"]:
        exec("{0}_{1}_{2} = histFile.Get(\"{0}_{1}_{2}\")".format(histCat,sigDef,sigDefexcl))
      
        with makeEnv_TCanvas("{0}/{1}_{2}_{3}.png".format(plotDir,histCat,sigDef,sigDefexcl)):
          exec("tHist_{0}_{1}_{2} = {0}_{1}_{2}.GetCVHistoWithError()".format(histCat,sigDef,sigDefexcl))
          ## Set horizontal axis label
          exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitle(\"Reconstructed #pi^{{0}} momentum\")".format(histCat,sigDef,sigDefexcl))
          exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitleSize(0.05)".format(histCat,sigDef,sigDefexcl))
          exec("tHist_{0}_{1}_{2}.GetYaxis().SetTitleSize(0.05)".format(histCat,sigDef,sigDefexcl))
          ## Set vertical axis label
          if histCat == "eff":
            exec("tHist_{0}_{1}_{2}.GetYaxis().SetTitle(\"Efficiency/GeV\")".format(histCat,sigDef,sigDefexcl))
            exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitle(\"Reconstructed #pi^{{0}} momentum\")".format(histCat,sigDef,sigDefexcl))
            exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitleSize(0.05)".format(histCat,sigDef,sigDefexcl))
          elif histCat == "xSection" or histCat == "xSection_mc":
            exec("tHist_{0}_{1}_{2}.GetYaxis().SetTitle(\"#sigma_{{NC 1 #pi^{{0}}}}[10^{{-38}} cm^{{2}}/Atom]/GeV\")".format(histCat,sigDef,sigDefexcl))
            exec("tHist_{0}_{1}_{2}.Scale(10**38)".format(histCat,sigDef,sigDefexcl))
            exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitle(\"Reconstructed #pi^{{0}} momentum\")".format(histCat,sigDef,sigDefexcl))
            exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitleSize(0.05)".format(histCat,sigDef,sigDefexcl))
          else:
            exec("tHist_{0}_{1}_{2}.GetYaxis().SetTitle(\"Number of Events/GeV\")".format(histCat,sigDef,sigDefexcl))
            exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitle(\"Reconstructed #pi^{{0}} momentum\")".format(histCat,sigDef,sigDefexcl))
            exec("tHist_{0}_{1}_{2}.GetXaxis().SetTitleSize(0.05)".format(histCat,sigDef,sigDefexcl))
          exec("tHist_{0}_{1}_{2}.Scale(1.0,\"width\")".format(histCat,sigDef,sigDefexcl))
          exec("tHist_{0}_{1}_{2}.Draw()".format(histCat,sigDef,sigDefexcl))

        if not histCat == "xSection_mc":
          for flux_uni in FLUX_SYSTS:
            exec("veb_{0}_{1}_{2} = {0}_{1}_{2}.GetVertErrorBand(\"{3}\")".format(histCat,sigDef,sigDefexcl,flux_uni))
            with makeEnv_TCanvas("{0}/breakout_{1}_{2}_{3}_{4}.png".format(plotDir,flux_uni,histCat,sigDef,sigDefexcl)):
              exec("veb_{0}_{1}_{2}.DrawAll(\"\",True)".format(histCat,sigDef,sigDefexcl))

        with makeEnv_TCanvas("{0}/errorSummary_{1}_{2}_{3}.png".format(plotDir,histCat,sigDef,sigDefexcl)):
          exec("localDrawErrorSummary(plotter,{0}_{1}_{2},\"Reconstructed #pi^{{0}} momentum\")".format(histCat,sigDef,sigDefexcl))

        ## Print out value and error for Mark to package into table
        exec("tHist = {0}_{1}_{2}.GetCVHistoWithError()".format(histCat,sigDef,sigDefexcl))
        cv_val = tHist.GetBinContent(1)
        err_val = tHist.GetBinError(1)
        exec("tHist_statError = {0}_{1}_{2}.GetCVHistoWithStatError()".format(histCat,sigDef,sigDefexcl))
        exec("tHist_systError = {0}_{1}_{2}.GetCVHistoWithError(False)".format(histCat,sigDef,sigDefexcl))
        err_val_stat = tHist_statError.GetBinError(1)
        err_val_syst = tHist_systError.GetBinError(1)
        print "sigDef: {0}\thistCat: {1}_{2}\tcv_val: {3}\terr_val_stat: {4}\terr_val_syst: {5}\terr_val_tot: {6}".format(histCat,sigDef,sigDefexcl,cv_val,err_val_stat,err_val_syst,err_val)
		
  for histCat in ["data_selected","BNB_ext","POT"]:
    print "I'm producing plots for histCat: {0}".format(histCat)
    exec("{0}_{1} = histFile.Get(\"{0}_{1}\")".format(histCat,sigDef))

    with makeEnv_TCanvas("{0}/{1}_{2}.png".format(plotDir,histCat,sigDef)):
      exec("tHist_{0}_{1} = {0}_{1}.GetCVHistoWithError()".format(histCat,sigDef))
      ## Set horizontal axis label
      exec("tHist_{0}_{1}.GetXaxis().SetTitle(\"Reconstructed #pi^{{0}} momentum\")".format(histCat,sigDef,sigDefexcl))
      exec("tHist_{0}_{1}.GetXaxis().SetTitleSize(0.05)".format(histCat,sigDef,sigDefexcl))
      exec("tHist_{0}_{1}.GetYaxis().SetTitleSize(0.05)".format(histCat,sigDef))
      ## Set vertical axis label
      exec("tHist_{0}_{1}.GetYaxis().SetTitle(\"Number of Events/GeV\")".format(histCat,sigDef))
      exec("tHist_{0}_{1}.Scale(1.0,\"width\")".format(histCat,sigDef))
      exec("tHist_{0}_{1}.Draw()".format(histCat,sigDef))
      
    for flux_uni in FLUX_SYSTS:
      exec("veb_{0}_{1} = {0}_{1}.GetVertErrorBand(\"{2}\")".format(histCat,sigDef,flux_uni))
      with makeEnv_TCanvas("{0}/breakout_{1}_{2}_{3}.png".format(plotDir,flux_uni,histCat,sigDef)):
        exec("veb_{0}_{1}.DrawAll(\"\",True)".format(histCat,sigDef))

    with makeEnv_TCanvas("{0}/errorSummary_{1}_{2}.png".format(plotDir,histCat,sigDef)):
      exec("localDrawErrorSummary(plotter,{0}_{1},\"Reconstructed #pi^{{0}} momentum\")".format(histCat,sigDef))

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
    ## Set horizontal axis label
    if histCat == "flux":
      exec("tHist_{0}.GetXaxis().SetTitle(\"E_{{#nu}}/GeV\")".format(histCat))
    else:
      exec("tHist_{0}.GetXaxis().SetTitle(\"Reconstructed #pi^{{0}} momentum)\")".format(histCat))
    exec("tHist_{0}.GetXaxis().SetTitleSize(0.05)".format(histCat))
    ## Set vertical axis label
    if histCat == "nTargets":
      exec("tHist_{0}.GetYaxis().SetTitle(\"Number of Targets[10^{{28}} atoms]/GeV\")".format(histCat))
      exec("tHist_{0}.Scale(10**-28)".format(histCat))
    else:
      exec("tHist_{0}.GetYaxis().SetTitle(\"#Phi[10^{{-9}} x #nu/POT/cm^{{2}}]/GeV\")".format(histCat))
      exec("tHist_{0}.Scale(10**9)".format(histCat))
    exec("tHist_{0}.GetYaxis().SetTitleSize(0.05)".format(histCat))
    exec("tHist_{0}.Scale(1.0,\"width\")".format(histCat))
    exec("tHist_{0}.Draw()".format(histCat))

  for flux_uni in FLUX_SYSTS:
    exec("veb_{0} = {0}.GetVertErrorBand(\"{1}\")".format(histCat,flux_uni))
    with makeEnv_TCanvas("{0}/breakout_{1}_{2}.png".format(plotDir,flux_uni,histCat)):
      exec("veb_{0}.DrawAll(\"\",True)".format(histCat))

  with makeEnv_TCanvas("{0}/errorSummary_{1}.png".format(plotDir,histCat)):
    exec("localDrawErrorSummary(plotter,{0},\"Reconstructed #pi^{{0}} momentum\")".format(histCat))

#############################################################################################################
### Make downstream breakout plots ##########################################################################
#############################################################################################################

for syst_uni in ["piplus_PrimaryHadronSWCentralSplineVariation","All_UBGenie"]:

  ## OMIT EXCL:
  ## for sigDefexcl in ["exclusive", "inclusive"]:
  for sigDefexcl in ["inclusive"]:

    exec("veb_xSection_2g1p_{1}_breakout = xSection_2g1p_{1}.GetVertErrorBand(\"{0}\")".format(syst_uni,sigDefexcl))
    with makeEnv_TCanvas("{0}/breakout_{1}_xSection.png".format(plotDir,syst_uni)):
      exec("veb_xSection_2g1p_{0}_breakout.DrawAll(\"\",True)".format(sigDefexcl))

