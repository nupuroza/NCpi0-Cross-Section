## Script to take unfolded signal and output from translateHists.py perform the cross-section calculation

import ROOT
import argparse
import os
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

mHist_evtRate_2g1p_inclusive_unfolded = outFile.Get("unfolded_evtRate_2g1p_inclusive")
mHist_eff_2g1p_inclusive = outFile.Get("eff_2g1p_inclusive")
mHist_flux_integral = outFile.Get("integratedFlux")
mHist_POT_2g1p = outFile.Get("POT_2g1p")
mHist_nTargets = outFile.Get("nTargets")

## Cross section calculation
mHist_xSection_2g1p_inclusive_unfolded = mHist_evtRate_2g1p_inclusive_unfolded.Clone("unfolded_xSection_2g1p_inclusive")
mHist_xSection_2g1p_inclusive_unfolded.Divide(mHist_xSection_2g1p_inclusive_unfolded,mHist_eff_2g1p_inclusive)
mHist_xSection_2g1p_inclusive_unfolded.Divide(mHist_xSection_2g1p_inclusive_unfolded,mHist_flux_integral)
mHist_xSection_2g1p_inclusive_unfolded.Divide(mHist_xSection_2g1p_inclusive_unfolded,mHist_POT_2g1p)# Remove units of per POT
mHist_xSection_2g1p_inclusive_unfolded.Divide(mHist_xSection_2g1p_inclusive_unfolded,mHist_nTargets)

writeHist(mHist_xSection_2g1p_inclusive_unfolded,outFile)

outFile.Close()
