### Script to take make plots, etc.
### using the MINERvA Analysis Toolkit

import ROOT
import argparse,math
from plottingClasses import *
from errorMaps_gLEE import *

## Set ROOT to batch mode
ROOT.gROOT.SetBatch()

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
### File Management #########################################################################################
#############################################################################################################

parser = argparse.ArgumentParser(description='Script to make cross-section plots using the MINERvA Analysis Toolkit')
parser.add_argument('in_file', help='Path to input file', type=str,nargs='?')
p = parser.parse_args()

## If in_dir is not provided, exit
if p.in_file < 0:
  print "ERROR: Input file argument not provided"
  parser.print_help()
  exit(1)

histFileLocation = p.in_file
print "This is the input file I'm opening: {0}".format(histFileLocation)
histFile = ROOT.TFile(histFileLocation)  

#############################################################################################################
### Pull out systematic uncertainties for xsecs #############################################################
#############################################################################################################

for sigDef in ["2g0p","2g1p","2gnp"]:
  
  for sigDefexcl in ["exclusive", "inclusive"]:
    
    if sigDef == "2gnp" and sigDefexcl == "exclusive":
      continue

    else:
  
      print "sigDef: {0}\tsigDefexcl: {1}".format(sigDef,sigDefexcl)
      print "-----------------------------------------"

      # Grab relevant histogram
      exec("xsec_{0}_{1} = histFile.Get(\"xSection_{0}_{1}\")".format(sigDef,sigDefexcl))
      # Total and stat errors first
      exec("tHist_totalError = xsec_{0}_{1}.GetCVHistoWithError()".format(sigDef,sigDefexcl))
      exec("tHist_statError = xsec_{0}_{1}.GetCVHistoWithStatError()".format(sigDef,sigDefexcl))
      cv_val = tHist_totalError.GetBinContent(1)
      err_val_total = tHist_totalError.GetBinError(1)/cv_val
      err_val_stat = tHist_statError.GetBinError(1)/cv_val
      err_val_syst = 0.0
      print "cv_val: {0}".format(cv_val)
      print "err_val_total: {0}".format(err_val_total)
      print "err_val_stat: {0}".format(err_val_stat)
      for err_to_keep in error_bands:
        # Make a copy of the hist that we'll pop error bands out of
        exec("local_xsec_only_{0} = xsec_{1}_{2}.Clone(\"local_xsec_only_{0}\")".format(err_to_keep,sigDef,sigDefexcl))
        # These are the error bands we'll pop
        local_err_list = list(error_bands)
        # We want to pop all but a single error band
        local_err_list.remove(err_to_keep)
        # Loop over error bands to pop
        for err_to_remove in local_err_list:
          for error_band in error_bands[err_to_remove]:
            exec("local_xsec_only_{0}.PopVertErrorBand(\"{1}\")".format(err_to_keep,error_band))
        # Now there should only be one error band left. Pull out CV and read off error 
        exec("tHist = local_xsec_only_{0}.GetCVHistoWithError(False)".format(err_to_keep))
        err_val = tHist.GetBinError(1)/cv_val
        print "err_to_keep: {1}\t{2}".format(cv_val,err_to_keep,err_val)
        err_val_syst = math.sqrt(err_val_syst*err_val_syst+err_val*err_val)

      print "err_val_syst: {0}".format(err_val_syst)
      print "-----------------------------------------"

#############################################################################################################
### Pull out systematic uncertainties for xsec ratio ########################################################
#############################################################################################################

print "xsec_ratio"
print "-----------------------------------------"

# Grab relevant histogram
xsec_ratio = histFile.Get("xSecRatio")
# Total and stat errors first
tHist_totalError = xsec_ratio.GetCVHistoWithError()
tHist_statError = xsec_ratio.GetCVHistoWithStatError()
cv_val = tHist_totalError.GetBinContent(1)
err_val_total = tHist_totalError.GetBinError(1)/cv_val
err_val_stat = tHist_statError.GetBinError(1)/cv_val
err_val_syst = 0.0
print "cv_val: {0}".format(cv_val)
print "err_val_total: {0}".format(err_val_total)
print "err_val_stat: {0}".format(err_val_stat)
for err_to_keep in error_bands:
  # Make a copy of the hist that we'll pop error bands out of
  exec("local_xsec_only_{0} = xsec_ratio.Clone(\"local_xsec_only_{0}\")".format(err_to_keep))
  # These are the error bands we'll pop
  local_err_list = list(error_bands)
  # We want to pop all but a single error band
  local_err_list.remove(err_to_keep)
  # Loop over error bands to pop
  for err_to_remove in local_err_list:
    for error_band in error_bands[err_to_remove]:
      exec("local_xsec_only_{0}.PopVertErrorBand(\"{1}\")".format(err_to_keep,error_band))
  # Now there should only be one error band left. Pull out CV and read off error 
  exec("tHist = local_xsec_only_{0}.GetCVHistoWithError(False)".format(err_to_keep))
  err_val = tHist.GetBinError(1)/cv_val
  print "err_to_keep: {1}\t{2}".format(cv_val,err_to_keep,err_val)
  err_val_syst = math.sqrt(err_val_syst*err_val_syst+err_val*err_val)

print "err_val_syst: {0}".format(err_val_syst)
print "-----------------------------------------"

