### Script to take TH1Ds evaluated in various systematic universes (provided by Guanqun circa 2021 Jun 8)
### and package them into MnvH1Ds using the MINERvA Analysis Toolkit

import ROOT
import datetime as dt
import argparse
import os
from array import *

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

#############################################################################################################
### Custom classes ##########################################################################################
#############################################################################################################

def writeHist(hist,outFile):
  outFile.cd()
  print 'Writing {0} to output file'.format(hist.GetName())
  hist.Write()

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################

## CV
cvFilePath_2g1p = "/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/XS_calculation_July2021/MomemtumVersions/XS2g1p_Momentum_FinalSelection_CV.SBNspec.root"
cvFile_2g1p = ROOT.TFile(cvFilePath_2g1p)
cvFilePath_2g0p = "/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/XS_calculation_July2021/MomemtumVersions/XS2g0p_Momentum_FinalSelection_CV.SBNspec.root"
cvFile_2g0p = ROOT.TFile(cvFilePath_2g0p)

## Efficiency Denominators
## TO-DO: Change to "true" momemtum binning
effDenomFilePath_2g1p_inclusive = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/earlier_stage/Inclusive/MultipleBin/2g1p/variation_spectra/SBNfit_variation_spectra_Flux_XS.root"
effDenomFile_2g1p_inclusive = ROOT.TFile(effDenomFilePath_2g1p_inclusive)

effDenomFilePath_2g0p_inclusive = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/earlier_stage/Inclusive/MultipleBin/2g0p/variation_spectra/SBNfit_variation_spectra_Flux_XS.root"
effDenomFile_2g0p_inclusive = ROOT.TFile(effDenomFilePath_2g0p_inclusive)

## OMIT EXCL:
## effDenomFilePath_2g1p_exclusive = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/earlier_stage/Exclusive/2g1p/KE_50MeV/variation_spectra/SBNfit_variation_spectra_Flux_XS.root"
## effDenomFile_2g1p_exclusive = ROOT.TFile(effDenomFilePath_2g1p_exclusive)

## effDenomFilePath_2g0p_exclusive = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/earlier_stage/Exclusive/2g0p/KE_50MeV/variation_spectra/SBNfit_variation_spectra_Flux_XS.root"
## effDenomFile_2g0p_exclusive = ROOT.TFile(effDenomFilePath_2g0p_exclusive)

## Efficiency numerators; backgrounds
## TO-DO: Change to "true" momemtum binning
## Final stage; flux, XS, Det systematics included
effNumFilePath_inclusive = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/final_stage/Inclusive/MultipleBin/variation_spectra/Merged_SBNfit_variation_spectra_FluxXSDet.root"
effNumFile_inclusive = ROOT.TFile(effNumFilePath_inclusive)
# 2g0p and 2g1p effNum hists come from the same file for the inclusive analysis
effNumFile_2g0p_inclusive = effNumFile_inclusive 
effNumFile_2g1p_inclusive = effNumFile_inclusive 

## OMIT EXCL:
## effNumFilePath_2g0p_exclusive = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/final_stage/Exclusive/NCPi0NoVisProton/KE_50MeV/variation_spectra/Merged_SBNfit_variation_spectra_FluxXSDet.root"
## effNumFile_2g0p_exclusive = ROOT.TFile(effNumFilePath_2g0p_exclusive)
## effNumFilePath_2g1p_exclusive = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/final_stage/Exclusive/NCPi0OneProton/KE_50MeV/variation_spectra/Merged_SBNfit_variation_spectra_FluxXSDet.root"
## effNumFile_2g1p_exclusive = ROOT.TFile(effNumFilePath_2g1p_exclusive)

## g4 variations
g4FilePath_inclusive = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/final_stage/Inclusive/MultipleBin/variation_spectra/SBNfit_variation_spectra_GEANT4.root"
g4File_inclusive = ROOT.TFile(g4FilePath_inclusive)
# 2g0p and 2g1p g4 hists come from the same file for the inclusive analysis
g4File_2g0p_inclusive = g4File_inclusive 
g4File_2g1p_inclusive = g4File_inclusive 

## OMIT EXCL:
## g4FilePath_2g0p_exclusive = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/final_stage/Exclusive/NCPi0NoVisProton/KE_50MeV/variation_spectra/SBNfit_variation_spectra_GEANT4.root"
## g4File_2g0p_exclusive = ROOT.TFile(g4FilePath_2g0p_exclusive)
## g4FilePath_2g1p_exclusive = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/final_stage/Exclusive/NCPi0OneProton/KE_50MeV/variation_spectra/SBNfit_variation_spectra_GEANT4.root"
## g4File_2g1p_exclusive = ROOT.TFile(g4FilePath_2g1p_exclusive)

## Output file
parser = argparse.ArgumentParser(description='Script to take TH1Ds evaluated in various systematic universes and package them into MnvH1Ds using the MINERvA Analysis Toolkit')
parser.add_argument('output_dir', help='Path to ouput directory', type=str,nargs='?')
parser.add_argument('test',help='Run in test mode using smaller number of syst universes (faster)',type=str,nargs='?')
p = parser.parse_args()

## If output_dir is not provided, exit
if p.output_dir < 0:
  print "ERROR: Output directory argument not provided"
  parser.print_help()
  exit(1)

## Create output directory if it doesn't exist
if not os.path.isdir(p.output_dir):
  print "Making output directory {0}".format(p.output_dir)
  os.system( "mkdir %s" % p.output_dir )

outputFilePath = p.output_dir+"/{0}_out.root".format(dt.date.today())
outFile = ROOT.TFile(outputFilePath,"recreate")

#############################################################################################################
### Create Reference Hist ###################################################################################
#############################################################################################################

## Create reference Hist that will be a template for whatever input binning is being used
histToBeCloned = cvFile_2g1p.Get("nu_uBooNE_2g1p_Data")
referenceHist = histToBeCloned.Clone("referenceHist")
nBins_analysis = referenceHist.GetNbinsX()
bin_lowEdges = []
for i in range(1,nBins_analysis+1):
  referenceHist.SetBinContent(i,-999.)
  referenceHist.SetBinError(i,0.)
  bin_lowEdges.append(referenceHist.GetBinLowEdge(i))
  binWidth = referenceHist.GetBinWidth(i)
bin_lowEdges.append(bin_lowEdges[-1] + binWidth)

#############################################################################################################
### Systematic Universes ####################################################################################
#############################################################################################################

## Number of Universes
## If running in test mode use a smaller number of flux universes
if p.test>0:
  nMultiverses = 10
else:
  nMultiverses = 1000

nMinMaxUniverses = 2

## List of cross section systematics
XS_SYSTS = [
  ## [syst_name,universe_prefix,n_universes]
  ["All_UBGenie","universe",nMultiverses],
  ["AxFFCCQEshape_UBGenie","minmax", nMinMaxUniverses],
  ["DecayAngMEC_UBGenie","minmax",nMinMaxUniverses],
  ["NormCCCOH_UBGenie","minmax",nMinMaxUniverses],
  ["NormNCCOH_UBGenie","minmax",nMinMaxUniverses],
  ["RPA_CCQE_UBGenie","universe",nMinMaxUniverses],
  ["Theta_Delta2Npi_UBGenie","minmax",nMinMaxUniverses],
  ["VecFFCCQEshape_UBGenie","minmax",nMinMaxUniverses],
]

## List of flux systematics
FLUX_SYSTS = [
  ## [syst_name,universe_prefix,n_universes]
  ["expskin_FluxUnisim","universe",nMultiverses],
  ["horncurrent_FluxUnisim","universe",nMultiverses],
  ["kminus_PrimaryHadronNormalization","universe",nMultiverses],
  ["kplus_PrimaryHadronFeynmanScaling","universe",nMultiverses],
  ["kzero_PrimaryHadronSanfordWang","universe",nMultiverses],
  ["nucleoninexsec_FluxUnisim","universe",nMultiverses],
  ["nucleonqexsec_FluxUnisim","universe",nMultiverses],
  ["nucleontotxsec_FluxUnisim","universe",nMultiverses],
  ["piminus_PrimaryHadronSWCentralSplineVariation","universe",nMultiverses],
  ["pioninexsec_FluxUnisim","universe",nMultiverses],
  ["pionqexsec_FluxUnisim","universe",nMultiverses],
  ["piontotxsec_FluxUnisim","universe",nMultiverses],
  ["piplus_PrimaryHadronSWCentralSplineVariation","universe",nMultiverses]
]

## List of detector systematics
DETECTOR_SYSTS = [
  ## [syst_name,universe_prefix,n_universes]
  ["AngleXZ","minmax", nMinMaxUniverses],
  ["AngleYZ","minmax", nMinMaxUniverses],
  ["LYAtt","minmax", nMinMaxUniverses],
  ["LYRay","minmax", nMinMaxUniverses],
  ["LY","minmax", nMinMaxUniverses],
  ["Recom2","minmax", nMinMaxUniverses],
  ["SCE","minmax", nMinMaxUniverses],
  ["WireX","minmax", nMinMaxUniverses],
  ["WireYZ","minmax", nMinMaxUniverses]
]

G4_SYSTS = [
  ## [syst_name,universe_prefix,n_universes]
  ["piminus","universe",nMultiverses],
  ["piplus","universe",nMultiverses],
  ["proton","universe",nMultiverses]
]
 
OTHER_SYSTS = [
  ## [syst_name,universe_prefix,n_universes]
  ["POT_variation","minmax", nMinMaxUniverses],
  ["target_variation","minmax", nMinMaxUniverses]
]
 
#############################################################################################################
### Calculate N_targets #####################################################################################
#############################################################################################################

## From Mark's script

fid_vol = 56408336.14 #cm^3 5cm
density = 1.3954 #from MC https://files.slack.com/files-pri/T0LFJ3CE5-F01ER1SK6G4/screen_shot_2020-11-11_at_4.18.34_pm.png
#At the temperature 89.2K  (https://lar.bnl.gov/properties/) is 1.3837 [g/cm3]
avo = 6.02214e23
molar_mass = 39.948
n_targets = density*fid_vol*avo/molar_mass

print "Number of targets: {0}".format(n_targets)

## Put scalar into TH1D
tHist_nTargets = referenceHist.Clone("tHist_nTargets")
for i in range(1,nBins_analysis+1):
  tHist_nTargets.SetBinContent(i,n_targets)

## Create MnvH1D from TH1D
mHist_nTargets = ROOT.PlotUtils.MnvH1D(tHist_nTargets)
mHist_nTargets.SetName("nTargets")
## Populate error bands
for systName,universePrefix,nUniverses in XS_SYSTS + FLUX_SYSTS + DETECTOR_SYSTS + G4_SYSTS + OTHER_SYSTS:
  mHist_nTargets.AddVertErrorBandAndFillWithCV(systName,nUniverses)


## Steven Gardiner told us to use a 1% variation
for i in range(1,nBins_analysis+1):
  mHist_nTargets.GetVertErrorBand("target_variation").GetHist(0).SetBinContent(i,n_targets*0.99)
  mHist_nTargets.GetVertErrorBand("target_variation").GetHist(1).SetBinContent(i,n_targets*1.01)

writeHist(mHist_nTargets,outFile)

#############################################################################################################
### Calculate POT ###########################################################################################
#############################################################################################################

POT_2g1p = 5.8447*10**20
POT_2g0p = 5.8930*10**20
POT_scaling = POT_2g1p/POT_2g0p

## For 2gnp, scale 2g0p POT to 2g1p POT before combining samples
POT_2gnp = 5.8447*10**20

for sigDef in ["2g0p","2g1p","2gnp"]:
  ## Put scalar into TH1D
  exec("tHist_POT_{0} = referenceHist.Clone(\"tHist_POT_{0}\")".format(sigDef))
  for i in range(1,nBins_analysis+1):
    exec("tHist_POT_{0}.SetBinContent(i,POT_{0})".format(sigDef))
  ## Create MnvH1D from TH1D
  exec("mHist_POT_{0} = ROOT.PlotUtils.MnvH1D(tHist_POT_{0})".format(sigDef))
  exec("mHist_POT_{0}.SetName(\"POT_{0}\")".format(sigDef))
  ## Populate error bands
  for systName,universePrefix,nUniverses in XS_SYSTS + FLUX_SYSTS + DETECTOR_SYSTS + G4_SYSTS + OTHER_SYSTS:
    exec("mHist_POT_{0}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(sigDef))

  ## Steven Gardiner told us to use a 2% variation
  for i in range(1,nBins_analysis+1):
    exec("mHist_POT_{0}.GetVertErrorBand(\"POT_variation\").GetHist(0).SetBinContent(i,POT_{0}*0.98)".format(sigDef))
    exec("mHist_POT_{0}.GetVertErrorBand(\"POT_variation\").GetHist(1).SetBinContent(i,POT_{0}*1.02)".format(sigDef))

  ## This mHist_POT_{0} object is only used in the xsec calculation, where we need the correct units
  #exec("mHist_POT_{0}.Scale(10**20)".format(sigDef))

  exec("writeHist(mHist_POT_{0},outFile)".format(sigDef))

#############################################################################################################
### Assemble Flux MnvH1D ####################################################################################
#############################################################################################################

## MCC9 (official files)
fluxFileDir = "/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_fluxes_01.09.2019_463_hist"
fluxFilePath = "{0}/MCC9_FluxHist_volTPCActive.root".format(fluxFileDir)
fluxFile = ROOT.TFile(fluxFilePath)

## Loop over neutrino species
for nuSpec in ["numu","numubar","nue","nuebar"]:
  
  ### Start with flux spectrum #####################
  ##################################################
  
  ## CV
  # Pull out the TH1D
  exec("tHist_flux_{0}_CV = fluxFile.Get(\"hE{0}_cv\")".format(nuSpec))
  # Copy this into an MnvH1D (no systs yet)
  exec("mHist_flux_{0} = ROOT.PlotUtils.MnvH1D(tHist_flux_{0}_CV)".format(nuSpec))
  exec("mHist_flux_{0}.SetName(\"flux_{0}\")".format(nuSpec))

  # We'll need this below
  exec("nBins_flux = tHist_flux_{0}_CV.GetNbinsX()".format(nuSpec))

  ## Loop over cross section, detector, GEANT4, and other systematics
  for systName,universePrefix,nUniverses in XS_SYSTS + DETECTOR_SYSTS + G4_SYSTS + OTHER_SYSTS:
  
    # Create the appropriate error band in the MnvH1D
    exec("mHist_flux_{0}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(nuSpec))
  
    ## Cross section variations don't affect the flux, so leave it at that
  
  ## Loop over flux systematics
  for systName,universePrefix,nUniverses in FLUX_SYSTS:

    print "systName: {0}".format(systName)   
    # Create the appropriate error band in the MnvH1D
    exec("mHist_flux_{0}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(nuSpec))
  
    # Loop over universes in this category of systematic
    for i in range(nUniverses):
      # Pull out the relevant TH1D
      exec("tHist_flux_{0}_{1}_{2} = fluxFile.Get(\"{0}_ms_{1}/hE{0}_{1}_ms_{2}\")".format(nuSpec,systName,i))
      # Loop over bins of flux
      for j in range(nBins_flux):
        # Pull out content of the relevant bin
        exec("binVal_{1}_{2} = tHist_flux_{0}_{1}_{2}.GetBinContent(j)".format(nuSpec,systName,i))
        # Populate corresponding hist in error band with bin content
        exec("mHist_flux_{0}.GetVertErrorBand(systName).GetHist(i).SetBinContent(j,binVal_{1}_{2})".format(nuSpec,systName,i))
    
  exec("writeHist(mHist_flux_{0},outFile)".format(nuSpec))
  
  ### Now integrate flux ###########################
  ##################################################
  
  ## Calculate CV flux integral
  exec("flux_integral_CV = tHist_flux_{0}_CV.Integral()".format(nuSpec))
  scale_factor = 1./(4997.*5*10**8)/(256.35*233.) # Units of POT*cm^2; from fluxFileDir/readme.txt

  ## Put scalar integrated flux into TH1D
  # Note that this hist has the binninng of the analysis, not the binning of the flux
  exec("tHist_flux_{0}_integral = referenceHist.Clone(\"tHist_flux_{0}_integral\")".format(nuSpec))
  for i in range(1,nBins_analysis+1):
    exec("tHist_flux_{0}_integral.SetBinContent(i,flux_integral_CV)".format(nuSpec))
  
  ## Create MnvH1D from TH1D
  exec("mHist_flux_{0}_integral = ROOT.PlotUtils.MnvH1D(tHist_flux_{0}_integral)".format(nuSpec))
  exec("mHist_flux_{0}_integral.SetName(\"integratedFlux_{0}\")".format(nuSpec))
  
  ## Populate error bands
  for systName,universePrefix,nUniverses in XS_SYSTS + DETECTOR_SYSTS + G4_SYSTS + OTHER_SYSTS:
    exec("mHist_flux_{0}_integral.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(nuSpec))

  for systName,universePrefix,nUniverses in FLUX_SYSTS:
    exec("mHist_flux_{0}_integral.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(nuSpec))
    for i in range(nUniverses):
      # Calculate flux integral for this universe
      exec("flux_integral_{1}_{2} = tHist_flux_{0}_{1}_{2}.Integral()".format(nuSpec,systName,i))
      # Populate corresponding hist in error band with bin content
      for j in range(1,nBins_analysis+1):
        exec("mHist_flux_{0}_integral.GetVertErrorBand(systName).GetHist(i).SetBinContent(j,flux_integral_{1}_{2})".format(nuSpec,systName,i))

  ## For some reason, applying the scale factor separately to each universe's integrated flux before filling the error bands went horribly wrong
  exec("mHist_flux_{0}_integral.Scale(scale_factor)".format(nuSpec))
  
  exec("writeHist(mHist_flux_{0}_integral,outFile)".format(nuSpec))

## Construct total flux
mHist_flux = mHist_flux_numu.Clone("flux")
mHist_flux.Add(mHist_flux_numubar)
mHist_flux.Add(mHist_flux_nue)
mHist_flux.Add(mHist_flux_nuebar)

writeHist(mHist_flux,outFile)

mHist_flux_integral = mHist_flux_numu_integral.Clone("integratedFlux")
mHist_flux_integral.Add(mHist_flux_numubar_integral)
mHist_flux_integral.Add(mHist_flux_nue_integral)
mHist_flux_integral.Add(mHist_flux_nuebar_integral)

writeHist(mHist_flux_integral,outFile)

#############################################################################################################
### Assemble MnvH1Ds for Data and BNB EXT Background (which is also data) ###################################
#############################################################################################################

## Pull out CV hists
for sigDef in ["2g1p","2g0p"]:
  exec("tHist_data_selected_{0} = cvFile_{0}.Get(\"nu_uBooNE_{0}_Data\")".format(sigDef))
  exec("tHist_BNB_ext_{0} = cvFile_{0}.Get(\"nu_uBooNE_{0}_BNBext\")".format(sigDef))

for dataDef in ["data_selected","BNB_ext"]:


  ## Add together 2g1p and 2g0p hists
  exec("tHist_{0}_2gnp = tHist_{0}_2g0p.Clone(\"tHist_{0}_2gnp\")".format(dataDef))
  exec("tHist_{0}_2gnp.Scale(POT_scaling)".format(dataDef))
  exec("tHist_{0}_2gnp.Add(tHist_{0}_2g1p)".format(dataDef))

  for sigDef in ["2g0p","2g1p","2gnp"]:

    ## Create MnvH1D from TH1D
    exec("mHist_{0}_{1} = ROOT.PlotUtils.MnvH1D(tHist_{0}_{1})".format(dataDef,sigDef))
    exec("mHist_{0}_{1}.SetName(\"{0}_{1}\")".format(dataDef,sigDef))

    ## Populate error bands
    for systName,universePrefix,nUniverses in XS_SYSTS + FLUX_SYSTS + DETECTOR_SYSTS + G4_SYSTS + OTHER_SYSTS:
      exec("mHist_{0}_{1}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(dataDef,sigDef))

    exec("writeHist(mHist_{0}_{1},outFile)".format(dataDef,sigDef))

#############################################################################################################
### Assemble xsec component MnvHnDs for 2g1p, 2g0p ##########################################################
#############################################################################################################

for sigDef in ["2g1p","2g0p"]:
  
  ## OMIT EXCL: for sigDefexcl in ['inclusive', 'exclusive']:
  for sigDefexcl in ['inclusive']:
    #############################################################################################################
    ### Construct Efficiency Denominator MnvH1D #################################################################
    #############################################################################################################

    ## CV
    # Pull out the TH1D
    exec("tHist_effDenom_{0}_{1}_CV = effDenomFile_{0}_{1}.Get(\"Flux_XS_CV_Dir/nu_uBooNE_AllNCPi0_Signal\")".format(sigDef,sigDefexcl))
    # Copy this into an MnvH1D (no systs yet)
    exec("mHist_effDenom_{0}_{1} = ROOT.PlotUtils.MnvH1D(tHist_effDenom_{0}_{1}_CV)".format(sigDef,sigDefexcl))
    # Rename new hist object
    exec("mHist_effDenom_{0}_{1}.SetName(\"effDenom_{0}_{1}\")".format(sigDef,sigDefexcl))

    ## Loop over cross section and flux systematics
    for systName,universePrefix,nUniverses in XS_SYSTS + FLUX_SYSTS:

      # Create the appropriate error band in the MnvH1D
      exec("mHist_effDenom_{0}_{1}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(sigDef,sigDefexcl))

      # Loop over universes in this category of systematic
      for i in range(nUniverses):
        # Pull out the relevant TH1D and map from TH1D universe numbering (starting at 1) to MnvH1D universe number (starting at 0)
        exec("tHist_effDenom_{0}_{1}_{2}_{3} = effDenomFile_{0}_{1}.Get(\"Flux_XS_{2}_Dir/nu_uBooNE_AllNCPi0_{4}_{5}_Signal\")".format(sigDef,sigDefexcl,systName,i,universePrefix,i+1))
        # Pull out content of relevant bin
        for j in range(1,nBins_analysis+1):
          exec("binVal = tHist_effDenom_{0}_{1}_{2}_{3}.GetBinContent(j)".format(sigDef,sigDefexcl,systName,i)) 
          # Populate corresponding hist in error band with bin content
          exec("mHist_effDenom_{0}_{1}.GetVertErrorBand(systName).GetHist(i).SetBinContent(j,binVal)".format(sigDef,sigDefexcl))


    ## Loop over detector, GEANT4, and other systematics
    for systName,universePrefix,nUniverses in DETECTOR_SYSTS + G4_SYSTS + OTHER_SYSTS:

      # Create the appropriate error band in the MnvH1D
      exec("mHist_effDenom_{0}_{1}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(sigDef,sigDefexcl))

      ## Detector variations don't actually affect the true event rate, so leave it at that

    exec("writeHist(mHist_effDenom_{0}_{1},outFile)".format(sigDef,sigDefexcl))

    #############################################################################################################
    ### Construct Efficiency Numerator & Background MnvH1D ##########################################################
    #############################################################################################################

    ### This is the actual CV
    #########################

    ## Pull out CV hists, which are separate for coherent and non-coherent
    exec("tHist_NCPi0Coh_{0} = cvFile_{0}.Get(\"nu_uBooNE_{0}_NCPi0Coh\")".format(sigDef))
    exec("tHist_NCPi0NotCoh_{0} = cvFile_{0}.Get(\"nu_uBooNE_{0}_NCPi0NotCoh\")".format(sigDef))
    
    ## Add together coherent and non-coherent
    exec("tHist_effNum_{0}_CV = tHist_NCPi0Coh_{0}.Clone(\"tHist_effNum_{0}_CV\")".format(sigDef))
    exec("tHist_effNum_{0}_CV.Add(tHist_NCPi0NotCoh_{0})".format(sigDef))

    ## Pull out CV hists, which are separate for various background components 
    for bkgDef in ["NCDelta","CC1Pi0","BNBOther","Nue","Dirt","OTPCExtra"]:
      exec("tHist_{0}_{1} = cvFile_{1}.Get(\"nu_uBooNE_{1}_{0}\")".format(bkgDef,sigDef))

    ## Add together various background components
    exec("tHist_background_{0}_CV = tHist_NCDelta_{0}.Clone(\"tHist_background_{0}_CV\")".format(sigDef))
    for bkgDef in ["CC1Pi0","BNBOther","Nue","Dirt","OTPCExtra"]:
      exec("tHist_background_{1}_CV.Add(tHist_{0}_{1})".format(bkgDef,sigDef))
    
    for histCat_label in [("effNum","Signal"),("background",'Bkgd')]:
      ## Pull out the value and save as a scalar
      ## This is the actual CV in each bin for the analysis
      binVal_trueCV = []
      for i in range(1,nBins_analysis+1):
        exec("binVal_trueCV.append(tHist_{0}_{1}_CV.GetBinContent(i))".format(histCat_label[0],sigDef))

      # Copy this into an MnvH1D (no systs yet)
      exec("mHist_{0}_{1}_{2} = ROOT.PlotUtils.MnvH1D(tHist_{0}_{1}_CV)".format(histCat_label[0],sigDef,sigDefexcl))
      # Rename new hist object
      exec("mHist_{0}_{1}_{2}.SetName(\"{0}_{1}_{2}\")".format(histCat_label[0],sigDef,sigDefexcl))

      ### This is the CV that will be used for internal comparisons against the reweightable systematics
      #########################

      exec("tHist_{0}_{1}_fakeCV_XS_GENIE = effNumFile_{1}_{2}.Get(\"Flux_XS_CV_Dir/nu_uBooNE_{1}_{3}\")".format(histCat_label[0],sigDef,sigDefexcl,histCat_label[1]))
      binVal_fakeCV_XS_GENIE = []
      for i in range(1,nBins_analysis+1):
        exec("binVal_fakeCV_XS_GENIE.append(tHist_{0}_{1}_fakeCV_XS_GENIE.GetBinContent(i))".format(histCat_label[0],sigDef))       

      ## Loop over cross section and flux systematics
      for systName,universePrefix,nUniverses in XS_SYSTS + FLUX_SYSTS:

        # Create the appropriate error band in the MnvH1D
        exec("mHist_{0}_{1}_{2}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(histCat_label[0],sigDef,sigDefexcl))
        # Loop over universes in this category of systematic
        for i in range(nUniverses):
          # Pull out the relevant TH1D and map from TH1D universe numbering (starting at 1) to MnvH1D universe number (starting at 0)
          exec("tHist_{0}_{1}_{2}_{3} = effNumFile_{1}_{4}.Get(\"Flux_XS_{2}_Dir/nu_uBooNE_{1}_{5}_{6}_{7}\")".format(histCat_label[0],sigDef,systName,i,sigDefexcl,universePrefix,i+1,histCat_label[1]))
          # Pull out content of relevant bin
          for j in range(1,nBins_analysis+1):
            exec("binVal_variation = tHist_{0}_{1}_{2}_{3}.GetBinContent(j)".format(histCat_label[0],sigDef,systName,i))
            # Calculate appropriate shift (applied fractionally for many universe evaluation)
            fracShift = (binVal_variation-binVal_fakeCV_XS_GENIE[j-1])/binVal_fakeCV_XS_GENIE[j-1]
            absShift = fracShift*binVal_trueCV[j-1]

            # Populate corresponding hist in error band with bin content
            exec("mHist_{0}_{1}_{2}.GetVertErrorBand(systName).GetHist(i).SetBinContent(j,binVal_trueCV[j-1]+absShift)".format(histCat_label[0],sigDef,sigDefexcl))

      ## Loop over detector systematics
      for systName,universePrefix,nUniverses in DETECTOR_SYSTS:

        # Create the appropriate error band in the MnvH1D
        exec("mHist_{0}_{1}_{2}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(histCat_label[0],sigDef,sigDefexcl))

        # Pull out the TH1D corresponding to the CV specifically generated to go with the detector system variation
        exec("tHist_{0}_{1}_{2}_fakeCV = effNumFile_{1}_{3}.Get(\"{2}_CV_Dir/nu_uBooNE_{1}_{4}\")".format(histCat_label[0],sigDef,systName,sigDefexcl,histCat_label[1]))
        # Pull out the TH1D corresponding to the detector system variation
        exec("tHist_{0}_{1}_{2}_variation = effNumFile_{1}_{3}.Get(\"{2}_{2}_Dir/nu_uBooNE_{1}_universe_1_{4}\")".format(histCat_label[0],sigDef,systName,sigDefexcl,histCat_label[1]))

        # Pull out content of relevant bin from each distribution
        for i in range(1,nBins_analysis+1):
          # Pull out content of relevant bin from each distribution
          exec("binVal_fakeCV_local = tHist_{0}_{1}_{2}_fakeCV.GetBinContent(i)".format(histCat_label[0],sigDef,systName))
          exec("binVal_variation = tHist_{0}_{1}_{2}_variation.GetBinContent(i)".format(histCat_label[0],sigDef,systName))

          # Calculate appropriate shift (applied absolutely for +/- 1 sigma shifts)
          fracShift = (binVal_variation-binVal_fakeCV_local)/binVal_fakeCV_local
          absShift = fracShift*binVal_trueCV[i-1]

          # Populate corresponding hist in error band with bin content
          exec("mHist_{0}_{1}_{2}.GetVertErrorBand(systName).GetHist(0).SetBinContent(i,binVal_trueCV[i-1]-absShift)".format(histCat_label[0],sigDef,sigDefexcl))
          exec("mHist_{0}_{1}_{2}.GetVertErrorBand(systName).GetHist(1).SetBinContent(i,binVal_trueCV[i-1]+absShift)".format(histCat_label[0],sigDef,sigDefexcl))

      ## Loop over GEANT4 systematics
      for systName,universePrefix,nUniverses in G4_SYSTS:

        # Create the appropriate error band in the MnvH1D
        exec("mHist_{0}_{1}_{2}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(histCat_label[0],sigDef,sigDefexcl))

        # Pull out the TH1D corresponding to the CV specifically generated to go with the GEANT4 variation
        exec("tHist_{0}_{1}_{2}_fakeCV = g4File_{1}_{3}.Get(\"GEANT4_CV_Dir/nu_uBooNE_{1}_{4}\")".format(histCat_label[0],sigDef,systName,sigDefexcl,histCat_label[1]))
        binVal_fakeCV_G4 = [] 
        for i in range(1,nBins_analysis+1):
          exec("binVal_fakeCV_G4.append(tHist_{0}_{1}_{2}_fakeCV.GetBinContent(i))".format(histCat_label[0],sigDef,systName))        

        # Loop over universes in this category of systematic
        for i in range(nUniverses):
          # Pull out the TH1D corresponding to the GEANT4 variation
          exec("tHist_{0}_{1}_{2}_{3} = g4File_{1}_{4}.Get(\"GEANT4_reinteractions_{2}_Geant4_Dir/nu_uBooNE_{1}_universe_{5}_{6}\")".format(histCat_label[0],sigDef,systName,i,sigDefexcl,i+1,histCat_label[1]))

          # Pull out content of relevant bin from each distribution
          for j in range(1,nBins_analysis+1):
            exec("binVal_variation = tHist_{0}_{1}_{2}_{3}.GetBinContent(j)".format(histCat_label[0],sigDef,systName,i))

            # Calculate appropriate shift (applied absolutely for +/- 1 sigma shifts)
            fracShift = (binVal_variation-binVal_fakeCV_G4[j-1])/binVal_fakeCV_G4[j-1]
            absShift = fracShift*binVal_trueCV[j-1]

            # Populate corresponding hist in error band with bin content
            exec("mHist_{0}_{1}_{2}.GetVertErrorBand(systName).GetHist(i).SetBinContent(j,binVal_trueCV[j-1]+absShift)".format(histCat_label[0],sigDef,sigDefexcl))

      ## Loop over other systematics
      for systName,universePrefix,nUniverses in OTHER_SYSTS:

        # Create the appropriate error band in the MnvH1D
        exec("mHist_{0}_{1}_{2}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(histCat_label[0],sigDef,sigDefexcl))

      exec("writeHist(mHist_{0}_{1}_{2},outFile)".format(histCat_label[0],sigDef,sigDefexcl))
    
#############################################################################################################
### Derive xsec component MnvHnDs for 2gnp ##################################################################
#############################################################################################################

for histCat in ["effNum","background"]:

  exec("mHist_{0}_2gnp_inclusive = mHist_{0}_2g0p_inclusive.Clone(\"{0}_2gnp_inclusive\")".format(histCat))
  exec("mHist_{0}_2gnp_inclusive.Scale(POT_scaling)".format(histCat))
  exec("mHist_{0}_2gnp_inclusive.Add(mHist_{0}_2g1p_inclusive)".format(histCat))

  exec("writeHist(mHist_{0}_2gnp_inclusive,outFile)".format(histCat))

## The 2gnp effDenom is just the 2g1p effDenom because the 2g0p sample has been scaled to have equivalent POT to 2g1p
mHist_effDenom_2gnp_inclusive = mHist_effDenom_2g1p_inclusive.Clone("effDenom_2gnp_inclusive")

writeHist(mHist_effDenom_2gnp_inclusive,outFile)

#############################################################################################################
### Loop over 2g1p, 2g0p, 2gnp ##############################################################################
#############################################################################################################

for sigDef in ["2g1p","2g0p","2gnp"]:

  ## OMIT EXCL: for sigDefexcl in ["inclusive","exclusive"]:
  for sigDefexcl in ["inclusive"]:

    #############################################################################################################
    ### Calculate Things ########################################################################################
    #############################################################################################################
    
    if sigDef == "2gnp" and sigDefexcl == "exclusive":
      continue
    
    else:
 
      ## Efficiency
      exec("mHist_eff_{0}_{1} = mHist_effNum_{0}_{1}.Clone(\"eff_{0}_{1}\")".format(sigDef,sigDefexcl))
      exec("mHist_eff_{0}_{1}.Divide(mHist_effNum_{0}_{1},mHist_effDenom_{0}_{1})".format(sigDef,sigDefexcl))

      exec("writeHist(mHist_eff_{0}_{1},outFile)".format(sigDef,sigDefexcl))

      ## Background-subtracted event rate
      exec("mHist_evtRate_{0}_{1} = mHist_data_selected_{0}.Clone(\"evtRate_{0}_{1}\")".format(sigDef,sigDefexcl))
      exec("mHist_evtRate_{0}_{1}.Add(mHist_BNB_ext_{0},-1.)".format(sigDef,sigDefexcl))
      exec("mHist_evtRate_{0}_{1}.Add(mHist_background_{0}_{1},-1.)".format(sigDef,sigDefexcl))

      exec("writeHist(mHist_evtRate_{0}_{1},outFile)".format(sigDef,sigDefexcl))

      ## Write efficiency numerator and background-subtracted event rate as TH1D
      exec("tHist_effNum_{0}_{1} = mHist_effNum_{0}_{1}.GetCVHistoWithError()".format(sigDef,sigDefexcl))
      exec("tHist_effNum_{0}_{1}.SetName(\"tHist_effNum_{0}_{1}\")".format(sigDef,sigDefexcl))
      exec("writeHist(tHist_effNum_{0}_{1},outFile)".format(sigDef,sigDefexcl))

      exec("tHist_evtRate_{0}_{1} = mHist_evtRate_{0}_{1}.GetCVHistoWithError()".format(sigDef,sigDefexcl))
      exec("tHist_evtRate_{0}_{1}.SetName(\"tHist_evtRate_{0}_{1}\")".format(sigDef,sigDefexcl))
      exec("writeHist(tHist_evtRate_{0}_{1},outFile)".format(sigDef,sigDefexcl))

      ## Cross section calculation
      exec("mHist_xSection_{0}_{1} = mHist_evtRate_{0}_{1}.Clone(\"xSection_{0}_{1}\")".format(sigDef,sigDefexcl))
      exec("mHist_xSection_{0}_{1}.Divide(mHist_xSection_{0}_{1},mHist_eff_{0}_{1})".format(sigDef,sigDefexcl))
      exec("mHist_xSection_{0}_{1}.Divide(mHist_xSection_{0}_{1},mHist_flux_integral)".format(sigDef,sigDefexcl))
      exec("mHist_xSection_{0}_{1}.Divide(mHist_xSection_{0}_{1},mHist_POT_{0})".format(sigDef,sigDefexcl)) # Remove units of per POT
      exec("mHist_xSection_{0}_{1}.Divide(mHist_xSection_{0}_{1},mHist_nTargets)".format(sigDef,sigDefexcl))
  
      exec("writeHist(mHist_xSection_{0}_{1},outFile)".format(sigDef,sigDefexcl))

      ## MC signal prediction
      exec("mHist_xSection_mc_{0}_{1} = mHist_effNum_{0}_{1}.Clone(\"xSection_mc_{0}_{1}\")".format(sigDef,sigDefexcl))
      ## We don't want to subtract off the BNB_ext and background prediction because 
      ## the MC signal doesn't include backgrounds
      exec("mHist_xSection_mc_{0}_{1}.Divide(mHist_xSection_mc_{0}_{1},mHist_eff_{0}_{1})".format(sigDef,sigDefexcl))
      exec("mHist_xSection_mc_{0}_{1}.Divide(mHist_xSection_mc_{0}_{1},mHist_flux_integral)".format(sigDef,sigDefexcl))
      exec("mHist_xSection_mc_{0}_{1}.Divide(mHist_xSection_mc_{0}_{1},mHist_POT_{0})".format(sigDef,sigDefexcl)) # Remove units of per POT
      exec("mHist_xSection_mc_{0}_{1}.Divide(mHist_xSection_mc_{0}_{1},mHist_nTargets)".format(sigDef,sigDefexcl))

      ## Pop out all error bands except GENIE from the MC xsection
      for systName,universePrefix,nUniverses in FLUX_SYSTS + DETECTOR_SYSTS + G4_SYSTS + OTHER_SYSTS:
        exec("mHist_xSection_mc_{0}_{1}.PopVertErrorBand(\"{2}\")".format(sigDef,sigDefexcl,systName))

      exec("writeHist(mHist_xSection_mc_{0}_{1},outFile)".format(sigDef,sigDefexcl))

      #############################################################################################################
      ### Extract Covariance Matrix for Unfolding #################################################################
      #############################################################################################################
      
      exec("tMat_cov_evtRate_{0}_{1} = mHist_evtRate_{0}_{1}.GetTotalErrorMatrix()".format(sigDef,sigDefexcl))
      exec("tHist2D_cov_evtRate_{0}_{1} = ROOT.TH2D(\"tHist2D_cov_evtRate_{0}_{1}\", \"\", nBins_analysis, array(\"d\",bin_lowEdges), nBins_analysis, array(\"d\",bin_lowEdges))".format(sigDef,sigDefexcl))   
      
      for i in range(1,nBins_analysis+1): 
        for j in range(1,nBins_analysis+1):
          exec("cov = ROOT.TH2D(tMat_cov_evtRate_{0}_{1}).GetBinContent(i+1,j+1)".format(sigDef,sigDefexcl))
          exec("tHist2D_cov_evtRate_{0}_{1}.SetBinContent(i,j,cov)".format(sigDef,sigDefexcl))

      exec("tHist2D_cov_evtRate_{0}_{1}.SetName(\"tHist2D_cov_evtRate_{0}_{1}\")".format(sigDef,sigDefexcl))   
      exec("writeHist(tHist2D_cov_evtRate_{0}_{1}, outFile)".format(sigDef,sigDefexcl))


#############################################################################################################
### Response Matrix #########################################################################################
#############################################################################################################

## TO-DO: Add path to response matrix
responseIdentity = tHist2D_cov_evtRate_2gnp_inclusive.Clone("responseIdentity")
for i in range(1,nBins_analysis+1):
  for j in range(1,nBins_analysis+1):
    if i == j:
      responseIdentity.SetBinContent(i,j,1.)
    else:
      responseIdentity.SetBinContent(i,j,0.)

for sigDef in ["2g1p","2g0p","2gnp"]:
  exec("tHist2D_response_{0} = responseIdentity.Clone(\"tHist2D_response_{0}\")".format(sigDef))
  exec("tHist2D_response_{0}.SetName(\"tHist2D_response_{0}\")".format(sigDef))
  exec("writeHist(tHist2D_response_{0}, outFile)".format(sigDef))


#############################################################################################################
### Close output file; other business #######################################################################
#############################################################################################################
outFile.Close()

