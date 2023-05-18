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

## Get xsec ingredients that are common to all calculations
mHist_flux_integral = outFile.Get("integratedFlux")
mHist_nTargets = outFile.Get("nTargets")

#for sigDef in ["2g1p","2g0p","2gnp"]:
for sigDef in ["2g1p","2g0p"]:
  #for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["exclusive"]:

    if sigDef == "2gnp" and sigDefexcl == "exclusive":
      continue
    
    else:
      
      ## Get xsec ingredients that are unique to this sigDef/sigDefexcl
      exec("mHist_POT_{0} = outFile.Get(\"POT_{0}\")".format(sigDef))
      exec("tHist_evtRate_unfolded_{0}_{1} = outFile.Get(\"unfolded_evtRate_{0}_{1}\")".format(sigDef,sigDefexcl))
      exec("mHist_evtRate_unfolded_{0}_{1} = ROOT.MnvH1D(tHist_evtRate_unfolded_{0}_{1})".format(sigDef,sigDefexcl))
      exec("mHist_eff_{0}_{1} = outFile.Get(\"eff_{0}_{1}\")".format(sigDef,sigDefexcl))

      ## Calculate cross section
      exec("mHist_xSection_{0}_{1} = mHist_evtRate_unfolded_{0}_{1}.Clone(\"unfolded_xSection_{0}_{1}\")".format(sigDef,sigDefexcl))
      ## The efficiency should *not* be divided out here because it is already corrected for in the unfolding through the smearcept matrix
      exec("mHist_xSection_{0}_{1}.Divide(mHist_xSection_{0}_{1},mHist_flux_integral)".format(sigDef,sigDefexcl))
      exec("mHist_xSection_{0}_{1}.Divide(mHist_xSection_{0}_{1},mHist_POT_{0})".format(sigDef,sigDefexcl)) # Remove units of per POT
      exec("mHist_xSection_{0}_{1}.Divide(mHist_xSection_{0}_{1},mHist_nTargets)".format(sigDef,sigDefexcl))

      ## Write xsec to output file
      exec("writeHist(mHist_xSection_{0}_{1},outFile)".format(sigDef,sigDefexcl))

outFile.Close()
