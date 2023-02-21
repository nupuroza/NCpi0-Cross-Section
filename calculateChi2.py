### Script to calculate chi2 using functionality of the MINERvA Analysis Toolkit
### using the MINERvA Analysis Toolkit

import ROOT
import datetime as dt
import argparse
import numpy
from array import *

from plottingClasses import *

## Set ROOT to batch mode
ROOT.gROOT.SetBatch()

## This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

## Setup MnvPlotter, which has all of the plotting utilities
plotter = ROOT.PlotUtils.MnvPlotter()

plotter.SetROOT6Palette(54) # I think this needs to be done again after instantiating MnvPlotter, because I guess MnvPlotter sets the palette as well
ROOT.gStyle.SetNumberContours(200)

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################

parser = argparse.ArgumentParser(description='Script to calculate chi2 using the MINERvA Analysis Toolkit')
parser.add_argument('in_file', help='Path to input file', type=str,nargs='?')
parser.add_argument('out_dir', help='Path to ouput directory', type=str)
p = parser.parse_args()

## If in_file is not provided, exit
if p.in_file < 0:
  print "ERROR: Input file argument not provided"
  parser.print_help()
  exit(1)

histFile = ROOT.TFile(p.in_file)

## If output_dir is not provided, exit
if p.out_dir < 0:
  print "ERROR: Output directory argument not provided"
  parser.print_help()
  exit(1)

## Create output directory if it doesn't exist
if not os.path.isdir(p.out_dir):
  print "Making output directory {0}".format(p.out_dir)
  os.system( "mkdir %s" % p.out_dir )

plotDir = p.out_dir

#############################################################################################################
### Calculate Data/GENIE Chi2 ###############################################################################
#############################################################################################################

mHist_xsec_data = histFile.Get("xSection_2gnp_inclusive")
mHist_xsec_mc   = histFile.Get("xSection_mc_2gnp_inclusive")

chi2 = plotter.Chi2DataMC(mHist_xsec_data,mHist_xsec_mc,1.0,True,False,False)

print "chi2: {0}".format(chi2)

tHist_data = mHist_xsec_data.GetCVHistoWithError()
data_cv = tHist_data.GetBinContent(1)
data_err = tHist_data.GetBinError(1)

print "data CV: {0}".format(data_cv)
print "data err: {0}".format(data_err)

tHist_mc = mHist_xsec_mc.GetCVHistoWithError()
mc_cv = tHist_mc.GetBinContent(1)

print "mc CV: {0}".format(mc_cv)

#### Manual chi2 calculation
chi2_manual = ((mc_cv-data_cv)*(mc_cv-data_cv))/(data_err*data_err)

print "chi2 manual: {0}".format(chi2_manual)

#############################################################################################################
### Repackage exclusive xsecs into combined 2-bin MnvH1D ####################################################
#############################################################################################################

## Pull MnvH1Ds out of input file
mHist_xsec_2g1p_data = histFile.Get("xSection_2g1p_exclusive")
mHist_xsec_2g0p_data = histFile.Get("xSection_2g0p_exclusive")

mHist_xsec_2g1p_mc   = histFile.Get("xSection_mc_2g1p_exclusive")
mHist_xsec_2g0p_mc   = histFile.Get("xSection_mc_2g0p_exclusive")

## Define reference hist with complete set of error bands (using mHist_xsec_2g1p_data, but any of them would work for this)
referenceHist_data = mHist_xsec_2g1p_data
referenceHist_mc   = mHist_xsec_2g1p_mc

## Merge 2g1p and 2g0p hists into single two-bin "exclusive_xsec" hist for each of data and mc
for datamc in ["data","mc"]:

  exec("referenceHist = referenceHist_{0}".format(datamc))

  ## Create empty 2-bin TH1D to be CV in MnvH1D
  exec("tHist_exclusive_xsecs_{0} = ROOT.TH1D(\"h_exclusive_xsecs_{0}\",\"h_exclusive_xsecs_{0}\",2,array('d',[0,1,2]))".format(datamc))
  
  ## Fill 2 bins of CV TH1D 
  exec("tHist_exclusive_xsecs_{0}.SetBinContent(1,mHist_xsec_2g1p_{0}.GetCVHistoWithStatError().GetBinContent(1))".format(datamc))
  exec("tHist_exclusive_xsecs_{0}.SetBinContent(2,mHist_xsec_2g0p_{0}.GetCVHistoWithStatError().GetBinContent(1))".format(datamc))
  exec("tHist_exclusive_xsecs_{0}.SetBinError(1,mHist_xsec_2g1p_{0}.GetCVHistoWithStatError().GetBinError(1))".format(datamc))
  exec("tHist_exclusive_xsecs_{0}.SetBinError(2,mHist_xsec_2g0p_{0}.GetCVHistoWithStatError().GetBinError(1))".format(datamc))

  ## Create MnvH1D from TH1D
  exec("mHist_exclusive_xsecs_{0} = ROOT.PlotUtils.MnvH1D(tHist_exclusive_xsecs_{0})".format(datamc))
  exec("mHist_exclusive_xsecs_{0}.SetName(\"exclusive_xsecs_{0}\")".format(datamc))

  ## Loop over error bands
  for eb_name in referenceHist.GetErrorBandNames():

    exec("nUniverses = referenceHist.GetVertErrorBand(\"{1}\").GetNHists()".format(datamc,eb_name))

    ## Add blank error band
    exec("mHist_exclusive_xsecs_{0}.AddVertErrorBand(\"{1}\",{2})".format(datamc,eb_name,nUniverses))

    ## Loop over universes in error band and fill hists
    for i in range(nUniverses):

      ## Grab values from separate universe hists
      exec("universeVal_2g1p = mHist_xsec_2g1p_{0}.GetVertErrorBand(\"{1}\").GetHist({2}).GetBinContent(1)".format(datamc,eb_name,i))
      exec("universeVal_2g0p = mHist_xsec_2g0p_{0}.GetVertErrorBand(\"{1}\").GetHist({2}).GetBinContent(1)".format(datamc,eb_name,i))

      ## Fill 2 bins of universe hist
      exec("mHist_exclusive_xsecs_{0}.GetVertErrorBand(\"{1}\").GetHist({2}).SetBinContent(1,universeVal_2g1p)".format(datamc,eb_name,i))
      exec("mHist_exclusive_xsecs_{0}.GetVertErrorBand(\"{1}\").GetHist({2}).SetBinContent(2,universeVal_2g0p)".format(datamc,eb_name,i))

#############################################################################################################
### DEBUG ###################################################################################################
#############################################################################################################

bins_exclusive_combo_data = [mHist_exclusive_xsecs_data.GetBinContent(1),mHist_exclusive_xsecs_data.GetBinContent(2)]
bins_exclusive_combo_mc   = [mHist_exclusive_xsecs_mc.GetBinContent(1),mHist_exclusive_xsecs_mc.GetBinContent(2)]
print "bins_exclusive_combo_data: {0}".format(bins_exclusive_combo_data)
print "bins_exclusive_combo_mc: {0}".format(bins_exclusive_combo_mc)

bins_exclusive_combo_data_err = [mHist_exclusive_xsecs_data.GetBinError(1),mHist_exclusive_xsecs_data.GetBinError(2)]
bins_exclusive_combo_mc_err   = [mHist_exclusive_xsecs_mc.GetBinError(1),mHist_exclusive_xsecs_mc.GetBinError(2)]
print "bins_exclusive_combo_data_err: {0}".format(bins_exclusive_combo_data_err)
print "bins_exclusive_combo_mc_err: {0}".format(bins_exclusive_combo_mc_err)

tMatrix_cov_data = mHist_exclusive_xsecs_data.GetTotalErrorMatrix(True,False,False)
list_cov_data = numpy.array([[tMatrix_cov_data[1][1],tMatrix_cov_data[1][2]],[tMatrix_cov_data[2][1],tMatrix_cov_data[2][2]]])
#list_cov_data = numpy.array([[tMatrix_cov_data[1][1],0.],[0.,tMatrix_cov_data[2][2]]]) ## DEBUG: remove correlations
print "list_cov_data:\n {0}".format(list_cov_data)
list_errMat_data = numpy.linalg.inv(list_cov_data)
print "list_errMat_data:\n {0}".format(list_errMat_data)
unity_test = numpy.dot(list_cov_data,list_errMat_data)
print "unity_test:\n {0}".format(unity_test)


chi2 = [[0,0],[0,0]]
chi2_tot = 0.0
for i in range(2):
  for j in range(2):
    tmp = (bins_exclusive_combo_data[i]-bins_exclusive_combo_mc[i])*list_errMat_data[i][j]*(bins_exclusive_combo_data[j]-bins_exclusive_combo_mc[j])
    chi2[i][j] = tmp
    chi2_tot += tmp
    print "chi2[{0}][{1}]: {2}\tchi2_tot: {3}".format(i,j,tmp,chi2_tot)

#############################################################################################################
### Stuff other model predictions into 2-bin MnvH1Ds (the universes don't matter for the calculation) #######
#############################################################################################################

## GENIE v2
mHist_exclusive_xsecs_genieV2 = mHist_exclusive_xsecs_mc.Clone("mHist_exclusive_xsecs_genieV2")
mHist_exclusive_xsecs_genieV2.SetBinContent(1,0.72318e-38) ## 2g1p
mHist_exclusive_xsecs_genieV2.SetBinContent(2,1.0033e-38)  ## 2g0p

## NEUT 
mHist_exclusive_xsecs_neut = mHist_exclusive_xsecs_mc.Clone("mHist_exclusive_xsecs_neut")
mHist_exclusive_xsecs_neut.SetBinContent(1,0.61092e-38) ## 2g1p
mHist_exclusive_xsecs_neut.SetBinContent(2,0.76273e-38) ## 2g0p

## NuWro 
mHist_exclusive_xsecs_nuwro = mHist_exclusive_xsecs_mc.Clone("mHist_exclusive_xsecs_nuwro")
mHist_exclusive_xsecs_nuwro.SetBinContent(1,0.58938e-38) ## 2g1p
mHist_exclusive_xsecs_nuwro.SetBinContent(2,0.96381e-38) ## 2g0p

## GiBUU 
mHist_exclusive_xsecs_gibuu = mHist_exclusive_xsecs_mc.Clone("mHist_exclusive_xsecs_gibuu")
mHist_exclusive_xsecs_gibuu.SetBinContent(1,0.73226e-38) ## 2g1p
mHist_exclusive_xsecs_gibuu.SetBinContent(2,0.77836e-38) ## 2g0p

#############################################################################################################
### Calculate Data/GENIE Chi2 again, but now for 2g1p and 2g0p exclusive xsecs simultaneously ###############
#############################################################################################################

chi2_combo_genieV3 = plotter.Chi2DataMC(mHist_exclusive_xsecs_data,mHist_exclusive_xsecs_mc,1.0,True,False,False)
print "chi2 combo data vs geniev3: {0}".format(chi2_combo_genieV3)

for generator in ["genieV2","neut","nuwro","gibuu"]:
  exec("chi2_combo_tmp = plotter.Chi2DataMC(mHist_exclusive_xsecs_data,mHist_exclusive_xsecs_{0},1.0,True,False,False)".format(generator))
  print "chi2 combo data vs {0}: {1}".format(generator,chi2_combo_tmp)

#############################################################################################################
### Draw Covariance Matrices ################################################################################
#############################################################################################################

with makeEnv_TCanvas("{0}/covarianceMatrix_exclusive_xsecs_data.png".format(plotDir)) as canvas:
  mHist_exclusive_xsecs_data_scaled = mHist_exclusive_xsecs_data.Clone("mHist_exclusive_xsecs_data_scaled")
  mHist_exclusive_xsecs_data_scaled.Scale(1e38)
  tmp_matrix = mHist_exclusive_xsecs_data_scaled.GetTotalErrorMatrix(True,False,False)
  tmp_hist = ROOT.TH2D(tmp_matrix)
  tmp_hist.GetZaxis().SetRangeUser(0,0.05)
  tmp_hist.Draw("colzTEXT")

with makeEnv_TCanvas("{0}/covarianceMatrix_exclusive_xsecs_mc.png".format(plotDir)) as canvas:
  mHist_exclusive_xsecs_mc_scaled = mHist_exclusive_xsecs_mc.Clone("mHist_exclusive_xsecs_mc_scaled")
  mHist_exclusive_xsecs_mc_scaled.Scale(1e38)
  tmp_matrix = mHist_exclusive_xsecs_mc_scaled.GetTotalErrorMatrix(True,False,False)
  tmp_hist = ROOT.TH2D(tmp_matrix)
  tmp_hist.GetZaxis().SetRangeUser(0,0.05)
  tmp_hist.Draw("colzTEXT")

## Correlation matrices should use a different color palette
plotter.SetROOT6Palette(87)

with makeEnv_TCanvas("{0}/correlationMatrix_exclusive_xsecs_data.png".format(plotDir)) as canvas:
  tmp_matrix = mHist_exclusive_xsecs_data.GetTotalCorrelationMatrix(True,True)
  tmp_hist = ROOT.TH2D(tmp_matrix)
  tmp_hist.GetZaxis().SetRangeUser(-1.0,1.0)
  tmp_hist.Draw("colzTEXT")

with makeEnv_TCanvas("{0}/correlationMatrix_exclusive_xsecs_mc.png".format(plotDir)) as canvas:
  tmp_matrix = mHist_exclusive_xsecs_mc.GetTotalCorrelationMatrix(True,True)
  tmp_hist = ROOT.TH2D(tmp_matrix)
  tmp_hist.GetZaxis().SetRangeUser(-1.0,1.0)
  tmp_hist.Draw("colzTEXT")

