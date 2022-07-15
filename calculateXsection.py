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

### File Management #########################################################################################
#############################################################################################################
parser = argparse.ArgumentParser(description='Script to take unfolded signal and output from translateHists.py perform the cross-section calculation')
parser.add_argument('in_dir', help='Path to input directory', type=str,nargs='?')
parser.add_argument('in_date', help='Creation date of input file (yyyy-mm-dd). Defaults to a file dated today if it exists', type=str,nargs='?')
p = parser.parse_args()

## If in_dir is not provided, exit
if p.in_dir < 0:
  print "ERROR: Input directory argument not provided"
  parser.print_help()
  exit(1)

## If in_date is not provided, search file created today
if p.in_date < 0:
  inFile_Location = p.in_dir+"/{0}_out_unfolded.root".format(dt.date.today())
  if not os.path.exists(inFile_Location):
    print "ERROR: An input ROOT file created today does not exist. Specify input date argument"
    parser.print_help()
    exit(1)
  else:
    print "This is the input file I'm opening: {0}".format(inFile_Location)
    inFile = ROOT.TFile(inFile_Location)
else: 
  inFile_Location = p.in_dir+"/"+p.in_date+"_out_unfolded.root"
  if not os.path.exists(inFile_Location):
    print "ERROR: An input ROOT file created on "+p.in_date+" does not exist"
    parser.print_help()
    exit(1)
  else:
    print "This is the input file I'm opening: {0}".format(inFile_Location)
    inFile = ROOT.TFile(inFile_Location)

## Set arguments
outFile_Location = p.in_dir+"/{0}_out_final.root".format(dt.date.today())
inFile.Cp(outFile_Location)
inFile.Close()
outFile = ROOT.TFile(outFile_Location, "UPDATE")

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
