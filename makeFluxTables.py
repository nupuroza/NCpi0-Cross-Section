### Script to take make plots, etc.
### using the MINERvA Analysis Toolkit

import ROOT
import datetime as dt
import argparse
import os
from customHistAndPlotMethods import makeEnv_TCanvas,localDrawErrorSummary

## Set ROOT to batch mode
ROOT.gROOT.SetBatch()

## #Load and implement Phil's plot style header file
#ROOT.gROOT.ProcessLine(".L ../style/myPlotStyle.h")
#ROOT.myPlotStyle()

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
  plotDir = outFileDir+"/{0}_flux-tables".format(dt.date.today())
else:
  ## Create p.out_dir if it doesn't exist
  if not os.path.isdir(p.out_dir):
    print "Making plot directory {0}".format(p.out_dir)
    os.system( "mkdir %s" % p.out_dir )
  plotDir = p.out_dir

#############################################################################################################
### Make simple flux plots ##################################################################################
#############################################################################################################

for histCat in ["flux","flux_numu","flux_numubar","flux_nue","flux_nuebar"]:

  # Open a file for writing the output
  with open("{0}/values_{1}.txt".format(plotDir,histCat), "w") as f:

    f.write("{0}:\n".format(histCat))
    f.write("-----------------------------\n")
    f.write("BIN\tVAL\tABS ERR\tFRAC ERROR\n")
  
    local_mHist = histFile.Get(histCat)
    local_tHist = local_mHist.GetCVHistoWithError()
    nBins_flux = local_tHist.GetNbinsX()
  
    # Create a dictionary to store the bin values and errors
    local_binDict = {}
  
    # Loop over bins and store the bin values and errors in the dictionary
    for i in range(nBins_flux):
        if i == 0: continue
        if i > 60: continue
        local_binDict["binVal_{0}_{1}".format(histCat, i)] = local_tHist.GetBinContent(i)*1./(4997.*5*10**8)/(256.35*233.) # Units of POT*cm^2; from fluxFileDir/readme.txt
        local_binDict["binErr_{0}_{1}".format(histCat, i)] = local_tHist.GetBinError(i)*1./(4997.*5*10**8)/(256.35*233.) # Units of POT*cm^2; from fluxFileDir/readme.txt
        local_binDict["binErr_frac_{0}_{1}".format(histCat, i)] = local_tHist.GetBinError(i)/local_tHist.GetBinContent(i)
  
        # Construct the output string using the values from the dictionary
        f.write("{0}\t{1}\t{2}\t{3}\n".format(i, local_binDict["binVal_{0}_{1}".format(histCat, i)], local_binDict["binErr_{0}_{1}".format(histCat, i)], local_binDict["binErr_frac_{0}_{1}".format(histCat, i)]))

  exec("{0} = histFile.Get(\"{0}\")".format(histCat))
  exec("{0}.GetXaxis().SetRangeUser(0,3)".format(histCat))

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

#for histCat in ["integratedFlux","integratedFlux_numu","integratedFlux_numubar","integratedFlux_nue","integratedFlux_nuebar"]:
#
#  exec("{0} = histFile.Get(\"{0}\")".format(histCat))
#
#  with makeEnv_TCanvas("{0}/{1}.png".format(plotDir,histCat)) as can:
#    exec("tHist_{0} = {0}.GetCVHistoWithError()".format(histCat))
#    ## Set axes labels
#    exec("tHist_{0}.GetXaxis().SetTitleSize(0.05)".format(histCat))
#    exec("tHist_{0}.GetYaxis().SetTitleSize(0.05)".format(histCat))
#    if histCat == "flux":
#      exec("tHist_{0}.GetXaxis().SetTitle(\"E_{{#nu}}(GeV)\")".format(histCat))
#      exec("tHist_{0}.GetYaxis().SetTitle(\"#nu/10^{{-6}} POT/cm^{{2}}/GeV\")".format(histCat))
#      exec("tHist_{0}.Scale(10**-6)".format(histCat))
#    elif histCat == "integratedFlux":
#      exec("tHist_{0}.GetXaxis().SetTitle(\"#pi^{{0}} momentum\")".format(histCat))
#      exec("tHist_{0}.GetYaxis().SetTitle(\"#nu/10^{{9}} POT/cm^{{2}}/GeV\")".format(histCat))
#      exec("tHist_{0}.Scale(10**9)".format(histCat))
#    elif histCat == "nTargets":
#      exec("tHist_{0}.GetXaxis().SetTitle(\"#pi^{{0}} momentum\")".format(histCat))
#      exec("tHist_{0}.GetYaxis().SetTitle(\"Number of Targets[10^{{28}} atoms]/GeV\")".format(histCat))
#      exec("tHist_{0}.Scale(10**-28)".format(histCat)) 
#    exec("tHist_{0}.Scale(1.0,\"width\")".format(histCat))
#    exec("tHist_{0}.Draw()".format(histCat))

