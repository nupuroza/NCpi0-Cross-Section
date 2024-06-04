## Script to take unfolded signal and output from translateHists.py perform the cross-section calculation

import ROOT
import argparse
import os
import math
from customHistAndPlotMethods import writeHist

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################

parser = argparse.ArgumentParser(description='Script to take unfolded signal and output from translateHists.py perform the cross-section calculation')
parser.add_argument('in_file', help='Path to input file', type=str,nargs='?')
p = parser.parse_args()

## If in_file is not provided, exit
if p.in_file < 0:
  print "ERROR: Input file argument not provided"
  parser.print_help()
  exit(1)

inFilePath = p.in_file
inFile = ROOT.TFile(inFilePath)

outFileDir = os.path.dirname(inFilePath)
fileBaseName,fileExtension = os.path.splitext(os.path.basename(inFilePath))
outFilePath = "{0}/{1}_xsec-extracted{2}".format(outFileDir,fileBaseName,fileExtension)

inFile.Cp(outFilePath)
inFile.Close()

outFile = ROOT.TFile(outFilePath, "UPDATE")

#############################################################################################################
### Cross section calculation ###############################################################################
#############################################################################################################

## Get xsec ingredients that are common to all calculations
mHist_flux_integral = outFile.Get("integratedFlux")
mHist_nTargets = outFile.Get("nTargets")
mHist_POT_data = outFile.Get("POT_data")

#for sigDef in ["2g1p","2g0p","2gnp"]:
for sigDefnp in ["2g1p","2g0p"]:
  #for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["exclusive"]:

    sigDef = sigDefnp + "_" + sigDefexcl

    if sigDef == "2gnp" and sigDefexcl == "exclusive":
      continue
    
    else:
      
      ## Get xsec ingredients that are unique to this sigDef/sigDefexcl
      tHist_unfolded_cov_evtRate = outFile.Get("unfolded_cov_evtRate_" + sigDef)
      tHist_evtRate_unfolded = outFile.Get("unfolded_evtRate_" + sigDef)
      mHist_xSection_mc = outFile.Get("xSection_mc_" + sigDef)
      tHist2D_add_smear_matrix = outFile.Get("add_smear_matrix_" + sigDef)
      nBins = tHist_evtRate_unfolded.GetNbinsX()
      for i in range(0, nBins+2):## Loop over bins
        cov_binContent = tHist_unfolded_cov_evtRate.GetBinContent(i,i)
        tHist_evtRate_unfolded.SetBinError(i,math.sqrt(cov_binContent))

      ## Calculate cross section
      tHist_xSection_data = tHist_evtRate_unfolded.Clone("unfolded_xSection_" + sigDef)
      ## The efficiency should *not* be divided out here because it is already corrected for in the unfolding through the smearcept matrix
      tHist_xSection_data.Divide(tHist_xSection_data,mHist_flux_integral)
      tHist_xSection_data.Divide(tHist_xSection_data,mHist_POT_data) # Remove units of per POT
      tHist_xSection_data.Divide(tHist_xSection_data,mHist_nTargets)

      ## Smear mc cross section
      tMat_xSection_mc = ROOT.TMatrixD(nBins + 2, 1, mHist_xSection_mc.GetArray())
      tMat_add_smear_matrix = ROOT.TMatrixD(nBins + 2, nBins + 2, tHist2D_add_smear_matrix.GetArray())
      tMat_smeared_xSection_mc = ROOT.TMatrixD(tMat_add_smear_matrix, ROOT.TMatrixD.kMult, tMat_xSection_mc)
      tHist_smeared_xSection_mc = mHist_xSection_mc.Clone("smeared_xSection_mc_" + sigDef)
      for i in range(nBins + 2):
        tHist_smeared_xSection_mc.SetBinContent(i, tMat_smeared_xSection_mc[i][0])

      ## Write xsec to output file
      writeHist(tHist_xSection_data,outFile)
      writeHist(tHist_smeared_xSection_mc, outFile)

outFile.Close()
