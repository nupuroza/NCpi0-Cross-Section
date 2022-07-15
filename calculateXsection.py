## Script to take unfolded signal and output from translateHists.py perform the cross-section calculation

import ROOT
import datetime as dt
import argparse
import os
from array import *

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

def writeHist(hist,outFile):
  outFile.cd()
  print 'Writing 2g1p to output file'.format(hist.GetName())
  hist.Write()

## Set arguments
parser = argparse.ArgumentParser(description='Script to take unfolded signal and output from translateHists.py perform the cross-section calculation')
parser.add_argument('inout_dir', help='Path to input/ouput directory', type=str)
parser.add_argument('in_date', help='Creation date of input file (yyyy-mm-dd)', type=str)
p = parser.parse_args()

histFileLocation = p.inout_dir+"/"+p.in_date+"_out_unfolded.root"
histFile = ROOT.TFile(histFileLocation, "UPDATE")

mHist_eff_2g1p_inclusive = histFile.Get("eff_2g1p_inclusive")
mHist_flux_integral = histFile.Get("integratedFlux")
mHist_POT_2g1p = histFile.Get("POT_2g1p")
mHist_nTargets = histFile.Get("nTargets")
mHist_evtRate_2g1p_inclusive_unfolded = histFile.Get("unfolded_evtRate_2g1p_inclusive")

## Cross section calculation
mHist_xSection_2g1p_inclusive_unfolded = mHist_evtRate_2g1p_inclusive_unfolded.Clone("unfolded_xSection_2g1p_inclusive")
mHist_xSection_2g1p_inclusive_unfolded.Divide(mHist_xSection_2g1p_inclusive_unfolded,mHist_eff_2g1p_inclusive)
mHist_xSection_2g1p_inclusive_unfolded.Divide(mHist_xSection_2g1p_inclusive_unfolded,mHist_flux_integral)
mHist_xSection_2g1p_inclusive_unfolded.Divide(mHist_xSection_2g1p_inclusive_unfolded,mHist_POT_2g1p)# Remove units of per POT
mHist_xSection_2g1p_inclusive_unfolded.Divide(mHist_xSection_2g1p_inclusive_unfolded,mHist_nTargets)

writeHist(mHist_xSection_2g1p_inclusive_unfolded,histFile)

## MC signal prediction
mHist_xSection_mc_2g1p_inclusive = mHist_effNum_2g1p_inclusive.Clone("xSection_mc_2g1p_inclusive")
## We don't want to subtract off the background prediction because the MC signal doesn't include backgrounds
mHist_xSection_mc_2g1p_inclusive.Divide(mHist_xSection_mc_2g1p_inclusive,mHist_eff_2g1p_inclusive)
mHist_xSection_mc_2g1p_inclusive.Divide(mHist_xSection_mc_2g1p_inclusive,mHist_flux_integral)
mHist_xSection_mc_2g1p_inclusive.Divide(mHist_xSection_mc_2g1p_inclusive,mHist_POT_2g1p) # Remove units of per POT
mHist_xSection_mc_2g1p_inclusive.Divide(mHist_xSection_mc_2g1p_inclusive,mHist_nTargets)

## Pop out all error bands except GENIE from the MC xsection
for systName,universePrefix,nUniverses in FLUX_SYSTS + DETECTOR_SYSTS + G4_SYSTS + OTHER_SYSTS:
  exec("mHist_xSection_mc_2g1p_inclusive.PopVertErrorBand(\"{2}\")".format(systName))
  writeHist(mHist_xSection_mc_2g1p_inclusive,outFile)

histFile.Close()
