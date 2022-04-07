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
  print 'Writing {0} to output file'.format(hist.GetName())
  hist.Write()

## Set arguments
parser = argparse.ArgumentParser(description='Script to take unfolded signal and output from translateHists.py perform the cross-section calculation')
parser.add_argument('inout_dir', help='Path to input/ouput directory', type=str)
parser.add_argument('in_date', help='Creation date of input file (yyyy-mm-dd)', type=str)
p = parser.parse_args()

histFileLocation = p.inout_dir+"/"+p.in_date+"_out.root"
histFile = ROOT.TFile(histFileLocation, "UPDATE")

unfoldedFileLocation = p.inout_dir+"/"+p.in_date+"_unfolded.root"
unfoldedFile = ROOT.TFile(unfoldedFileLocation, "READ") 

mHist_eff_2gnp_inclusive = histFile.Get("eff_2gnp_inclusive")
mHist_flux_integral = histFile.Get("flux_integral")
mHist_POT_2gnp = histFile.Get("POT_2gnp")
mHist_nTargets = histFile.Get("nTargets")

ref_mHist = mHist_eff_2gnp_inclusive.Clone("ref_mHist")
tHist_evtRate_2gnp_inclusive_unfolded = unfoldedFile.Get("unf_signal")
mHist_evtRate_2gnp_inclusive_unfolded = ROOT.PlotUtils.MnvH1D(tHist_evtRate_2gnp_inclusive_unfolded)
mHist_evtRate_2gnp_inclusive_unfolded.AddMissingErrorBandsAndFillWithCV(ref_mHist)
mHist_evtRate_2gnp_inclusive_unfolded.SetName("evtRate_2gnp_inclusive_unfolded")
writeHist(mHist_evtRate_2gnp_inclusive_unfolded,histFile)

## Cross section calculation
mHist_xSection_2gnp_inclusive_unfolded = mHist_evtRate_2gnp_inclusive_unfolded.Clone("xSection_2gnp_inclusive_unfolded")
mHist_xSection_2gnp_inclusive_unfolded.Divide(mHist_xSection_2gnp_inclusive_unfolded,mHist_eff_2gnp_inclusive)
mHist_xSection_2gnp_inclusive_unfolded.Divide(mHist_xSection_2gnp_inclusive_unfolded,mHist_flux_integral)
mHist_xSection_2gnp_inclusive_unfolded.Divide(mHist_xSection_2gnp_inclusive_unfolded,mHist_POT_2gnp)# Remove units of per POT
mHist_xSection_2gnp_inclusive_unfolded.Divide(mHist_xSection_2gnp_inclusive_unfolded,mHist_nTargets)

writeHist(mHist_xSection_2gnp_inclusive_unfolded,histFile)

histFile.Close()
