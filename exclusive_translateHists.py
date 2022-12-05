### Script to take TH1Ds evaluated in various systematic universes (provided by Guanqun circa 2021 Jun 8)
### and package them into MnvH1Ds using the MINERvA Analysis Toolkit

import ROOT
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
cvFilePath_2g1p = "/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/XS_calculation_July2021/XS2g1p_MConly_FinalSelection_CV.SBNspec.root"
cvFile_2g1p = ROOT.TFile(cvFilePath_2g1p)
cvFilePath_2g0p = "/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/XS_calculation_July2021/XS2g0p_MConly_FinalSelection_CV.SBNspec.root"
cvFile_2g0p = ROOT.TFile(cvFilePath_2g0p)

## Efficiency Denominators (from earlier stage of processing; no nonreweightable systematics)
#effDenomFilePath_2g1p = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/earlier_stage/Inclusive/2g1p/variation_spectra/SBNfit_variation_spectra_Flux_XS.root"
#effDenomFilePath_2g1p = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/earlier_stage/Exclusive/2g1p/KE_20MeV/variation_spectra/SBNfit_variation_spectra_Flux_XS.root"
effDenomFilePath_2g1p = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/earlier_stage/Exclusive/2g1p/KE_50MeV/variation_spectra/SBNfit_variation_spectra_Flux_XS.root"
effDenomFile_2g1p = ROOT.TFile(effDenomFilePath_2g1p)
#effDenomFilePath_2g0p = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/earlier_stage/Inclusive/2g0p/variation_spectra/SBNfit_variation_spectra_Flux_XS.root"
#effDenomFilePath_2g0p = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/earlier_stage/Exclusive/2g0p/KE_20MeV/variation_spectra/SBNfit_variation_spectra_Flux_XS.root"
effDenomFilePath_2g0p = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/earlier_stage/Exclusive/2g0p/KE_50MeV/variation_spectra/SBNfit_variation_spectra_Flux_XS.root"
effDenomFile_2g0p = ROOT.TFile(effDenomFilePath_2g0p)

## Final stage; flux, XS, Det systematics included
## (histograms of 2g1p and 2g0p are included in this file with naming "nu_uBooNE_2g1p/2g0p_xxx")
effNumFilePath_2g0p = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/final_stage/Exclusive/NCPi0NoVisProton/KE_50MeV/variation_spectra/Merged_SBNfit_variation_spectra_FluxXSDet.root"
effNumFile_2g0p = ROOT.TFile(effNumFilePath_2g0p)
effNumFilePath_2g1p = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/final_stage/Exclusive/NCPi0OneProton/KE_50MeV/variation_spectra/Merged_SBNfit_variation_spectra_FluxXSDet.root"
effNumFile_2g1p = ROOT.TFile(effNumFilePath_2g1p)

g4FilePath_2g0p = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/final_stage/Exclusive/NCPi0NoVisProton/KE_50MeV/variation_spectra/SBNfit_variation_spectra_GEANT4.root"
g4File_2g0p = ROOT.TFile(g4FilePath_2g0p)
g4FilePath_2g1p = "/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/NCpi_cross_section/final_stage/Exclusive/NCPi0OneProton/KE_50MeV/variation_spectra/SBNfit_variation_spectra_GEANT4.root"
g4File_2g1p = ROOT.TFile(g4FilePath_2g1p)

## Output file 
#outputFilePath = "2021-12-09_out_exclusive_20MeV.root"
outputFilePath = "2021-12-09_out_exclusive_50MeV.root"
outFile = ROOT.TFile(outputFilePath,"recreate")

#############################################################################################################
### Systematic Universes ####################################################################################
#############################################################################################################

## List of cross section systematics
XS_SYSTS = [
  ## [syst_name,universe_prefix,n_universes]
  ["All_UBGenie","universe",1000],
  ["AxFFCCQEshape_UBGenie","minmax",2],
  ["DecayAngMEC_UBGenie","minmax",2],
  ["NormCCCOH_UBGenie","minmax",2],
  ["NormNCCOH_UBGenie","minmax",2],
  ["RPA_CCQE_UBGenie","universe",2],
  ["Theta_Delta2Npi_UBGenie","minmax",2],
  ["VecFFCCQEshape_UBGenie","minmax",2],
]

## List of flux systematics
FLUX_SYSTS = [
  ## [syst_name,universe_prefix,n_universes]
  "expskin_FluxUnisim",
  "horncurrent_FluxUnisim",
  "kminus_PrimaryHadronNormalization",
  "kplus_PrimaryHadronFeynmanScaling",
  "kzero_PrimaryHadronSanfordWang",
  "nucleoninexsec_FluxUnisim",
  "nucleonqexsec_FluxUnisim",
  "nucleontotxsec_FluxUnisim",
  "piminus_PrimaryHadronSWCentralSplineVariation",
  "pioninexsec_FluxUnisim",
  "pionqexsec_FluxUnisim",
  "piontotxsec_FluxUnisim",
  "piplus_PrimaryHadronSWCentralSplineVariation"
]

## Easy way to make the program run faster
#nFluxUniverses = 10
nFluxUniverses = 1000

## List of detector systematics
DETECTOR_SYSTS = [
  "AngleXZ",
  "AngleYZ",
  "LYAtt",
  "LYRay",
  "LY",
  "Recom2",
  "SCE",
  "WireX",
  "WireYZ",
]

G4_SYSTS = [
  "piminus",
  "piplus",
  "proton"
]
 
OTHER_SYSTS = [
  "POT_variation",
  "target_variation"
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
tHist_nTargets = ROOT.TH1D('h_nTargets','h_nTargets',1,array('d',[-1,1]))
tHist_nTargets.SetBinContent(1,n_targets)
## Create MnvH1D from TH1D
mHist_nTargets = ROOT.PlotUtils.MnvH1D(tHist_nTargets)
mHist_nTargets.SetName("nTargets")
## Populate error bands
for systName,universePrefix,nUniverses in XS_SYSTS:
  mHist_nTargets.AddVertErrorBandAndFillWithCV(systName,nUniverses)
for systName in FLUX_SYSTS:
  mHist_nTargets.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)
for systName in DETECTOR_SYSTS:
  mHist_nTargets.AddVertErrorBandAndFillWithCV(systName,2)
for systName in G4_SYSTS:
  mHist_nTargets.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)
for systName in OTHER_SYSTS:
  mHist_nTargets.AddVertErrorBandAndFillWithCV(systName,2)

## Steven Gardiner told us to use a 1% variation  
mHist_nTargets.GetVertErrorBand("target_variation").GetHist(0).SetBinContent(1,n_targets*0.99)
mHist_nTargets.GetVertErrorBand("target_variation").GetHist(1).SetBinContent(1,n_targets*1.01)

writeHist(mHist_nTargets,outFile)

#############################################################################################################
### Calculate POT ###########################################################################################
#############################################################################################################

POT_2g1p = 5.8447*10**20
POT_2g0p = 5.8930*10**20

for sigDef in ["2g0p","2g1p"]:
  ## Put scalar into TH1D
  exec("tHist_POT_{0} = ROOT.TH1D('h_POT_{0}','h_POT_{0}',1,array('d',[-1,1]))".format(sigDef))
  exec("tHist_POT_{0}.SetBinContent(1,POT_{0})".format(sigDef))
  ## Create MnvH1D from TH1D
  exec("mHist_POT_{0} = ROOT.PlotUtils.MnvH1D(tHist_POT_{0})".format(sigDef))
  exec("mHist_POT_{0}.SetName(\"POT_{0}\")".format(sigDef))
  ## Populate error bands
  for systName,universePrefix,nUniverses in XS_SYSTS:
    exec("mHist_POT_{0}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(sigDef))
  for systName in FLUX_SYSTS:
    exec("mHist_POT_{0}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(sigDef))
  for systName in DETECTOR_SYSTS:
    exec("mHist_POT_{0}.AddVertErrorBandAndFillWithCV(systName,2)".format(sigDef))
  for systName in G4_SYSTS:
    exec("mHist_POT_{0}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(sigDef))
  for systName in OTHER_SYSTS:
    exec("mHist_POT_{0}.AddVertErrorBandAndFillWithCV(systName,2)".format(sigDef))

  ## Steven Gardiner told us to use a 2% variation  
  exec("mHist_POT_{0}.GetVertErrorBand(\"POT_variation\").GetHist(0).SetBinContent(1,POT_{0}*0.98)".format(sigDef))
  exec("mHist_POT_{0}.GetVertErrorBand(\"POT_variation\").GetHist(1).SetBinContent(1,POT_{0}*1.02)".format(sigDef))

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

  ## Loop over cross section systematics
  for systName,universePrefix,nUniverses in XS_SYSTS:
  
    # Create the appropriate error band in the MnvH1D
    exec("mHist_flux_{0}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(nuSpec))
  
    ## Cross section variations don't affect the flux, so leave it at that
  
  ## Loop over flux systematics
  for systName in FLUX_SYSTS:

    print "systName: {0}".format(systName)   
    # Create the appropriate error band in the MnvH1D
    exec("mHist_flux_{0}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(nuSpec))
  
    # Loop over universes in this category of systematic
    for i in range(nFluxUniverses):
      # Pull out the relevant TH1D
      exec("tHist_flux_{0}_{1}_{2} = fluxFile.Get(\"{0}_ms_{1}/hE{0}_{1}_ms_{2}\")".format(nuSpec,systName,i))
      # Loop over bins of flux
      for j in range(nBins_flux):
        # Pull out content of the relevant bin
        exec("binVal_{1}_{2} = tHist_flux_{0}_{1}_{2}.GetBinContent(j)".format(nuSpec,systName,i))
        # Populate corresponding hist in error band with bin content
        exec("mHist_flux_{0}.GetVertErrorBand(systName).GetHist(i).SetBinContent(j,binVal_{1}_{2})".format(nuSpec,systName,i))
 
  ## Loop over detector systematics
  for systName in DETECTOR_SYSTS:
  
    # Create the appropriate error band in the MnvH1D
    exec("mHist_flux_{0}.AddVertErrorBandAndFillWithCV(systName,2)".format(nuSpec))
  
    ## Detector variations don't affect the flux, so leave it at that
  
  ## Loop over detector systematics
  for systName in G4_SYSTS:
  
    # Create the appropriate error band in the MnvH1D
    exec("mHist_flux_{0}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(nuSpec))
  
    ## Detector variations don't affect the flux, so leave it at that
  ## Loop over other systematics
  for systName in OTHER_SYSTS:
  
    # Create the appropriate error band in the MnvH1D
    exec("mHist_flux_{0}.AddVertErrorBandAndFillWithCV(systName,2)".format(nuSpec))
  
  exec("writeHist(mHist_flux_{0},outFile)".format(nuSpec))
  
  ### Now integrate flux ###########################
  ##################################################
  
  ## Calculate CV flux integral
  exec("flux_integral_CV = tHist_flux_{0}_CV.Integral()".format(nuSpec))
  scale_factor = 1./(4997.*5*10**8)/(256.35*233.) # Units of POT*cm^2; from fluxFileDir/readme.txt

  ## Put scalar integrated flux into TH1D
  # Note that this hist has the binninng of the analysis, not the binning of the flux
  exec("tHist_flux_{0}_integral = ROOT.TH1D('h_flux_integral','h_flux_integral',1,array('d',[-1,1]))".format(nuSpec))
  exec("tHist_flux_{0}_integral.SetBinContent(1,flux_integral_CV)".format(nuSpec))
  
  ## Create MnvH1D from TH1D
  exec("mHist_flux_{0}_integral = ROOT.PlotUtils.MnvH1D(tHist_flux_{0}_integral)".format(nuSpec))
  exec("mHist_flux_{0}_integral.SetName(\"integratedFlux_{0}\")".format(nuSpec))
  
  ## Populate error bands
  for systName,universePrefix,nUniverses in XS_SYSTS:
    exec("mHist_flux_{0}_integral.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(nuSpec))
  #for systName in FLUX_SYSTS:
  for systName in FLUX_SYSTS:
    exec("mHist_flux_{0}_integral.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(nuSpec))
    for i in range(nFluxUniverses):
      # Calculate flux integral for this universe
      exec("flux_integral_{1}_{2} = tHist_flux_{0}_{1}_{2}.Integral()".format(nuSpec,systName,i))
      # Populate corresponding hist in error band with bin content
      exec("mHist_flux_{0}_integral.GetVertErrorBand(systName).GetHist(i).SetBinContent(1,flux_integral_{1}_{2})".format(nuSpec,systName,i))

  ## For some reason, applying the scale factor separately to each universe's integrated flux before filling the error bands went horribly wrong
  exec("mHist_flux_{0}_integral.Scale(scale_factor)".format(nuSpec))
      
  for systName in DETECTOR_SYSTS:
    exec("mHist_flux_{0}_integral.AddVertErrorBandAndFillWithCV(systName,2)".format(nuSpec))
  
  for systName in G4_SYSTS:
    exec("mHist_flux_{0}_integral.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(nuSpec))
  
  for systName in OTHER_SYSTS:
    exec("mHist_flux_{0}_integral.AddVertErrorBandAndFillWithCV(systName,2)".format(nuSpec))
  
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

## Print out integrated flux
integrated_flux = mHist_flux_integral.GetCVHistoWithError()
integrated_flux_cv = integrated_flux.GetBinContent(1)
integrated_flux_err = integrated_flux.GetBinError(1)
print "--integrated flux-- CV: {0}\terr: {1}".format(integrated_flux_cv,integrated_flux_err)

#############################################################################################################
### Assemble MnvH1Ds for Data and BNB EXT Background (which is also data) ###################################
#############################################################################################################

## Pull out CV hists
for sigDef in ["2g1p","2g0p"]:
  exec("tHist_data_selected_{0} = cvFile_{0}.Get(\"nu_uBooNE_{0}_Data\")".format(sigDef))
  exec("tHist_BNB_ext_{0} = cvFile_{0}.Get(\"nu_uBooNE_{0}_BNBext\")".format(sigDef))

for dataDef in ["data_selected","BNB_ext"]:

  for sigDef in ["2g0p","2g1p"]:

    ## Create MnvH1D from TH1D
    exec("mHist_{0}_{1} = ROOT.PlotUtils.MnvH1D(tHist_{0}_{1})".format(dataDef,sigDef))
    exec("mHist_{0}_{1}.SetName(\"{0}_{1}\")".format(dataDef,sigDef))

    ## Populate error bands
    for systName,universePrefix,nUniverses in XS_SYSTS:
      exec("mHist_{0}_{1}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(dataDef,sigDef))
    for systName in FLUX_SYSTS:
      exec("mHist_{0}_{1}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(dataDef,sigDef))
    for systName in DETECTOR_SYSTS:
      exec("mHist_{0}_{1}.AddVertErrorBandAndFillWithCV(systName,2)".format(dataDef,sigDef))
    for systName in G4_SYSTS:
      exec("mHist_{0}_{1}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(dataDef,sigDef))
    for systName in OTHER_SYSTS:
      exec("mHist_{0}_{1}.AddVertErrorBandAndFillWithCV(systName,2)".format(dataDef,sigDef))

    exec("writeHist(mHist_{0}_{1},outFile)".format(dataDef,sigDef))

#############################################################################################################
### Assemble xsec component MnvHnDs for 2g1p, 2g0p ##########################################################
#############################################################################################################

for sigDef in ["2g1p","2g0p"]:

  #############################################################################################################
  ### Construct Efficiency Denominator MnvH1D #################################################################
  #############################################################################################################

  ## CV
  # Pull out the TH1D
  exec("tHist_effDenom_{0}_CV = effDenomFile_{0}.Get(\"Flux_XS_CV_Dir/nu_uBooNE_AllNCPi0_Signal\")".format(sigDef))
  # Copy this into an MnvH1D (no systs yet)
  exec("mHist_effDenom_{0} = ROOT.PlotUtils.MnvH1D(tHist_effDenom_{0}_CV)".format(sigDef))
  # Rename new hist object
  exec("mHist_effDenom_{0}.SetName(\"effDenom_{0}\")".format(sigDef))

  ## Loop over cross section systematics
  for systName,universePrefix,nUniverses in XS_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_effDenom_{0}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(sigDef))

    # Loop over universes in this category of systematic
    for i in range(nUniverses):
      # Pull out the relevant TH1D and map from TH1D universe numbering (starting at 1) to MnvH1D universe number (starting at 0)
      exec("tHist_effDenom_{0}_{1}_{2} = effDenomFile_{0}.Get(\"Flux_XS_{1}_Dir/nu_uBooNE_AllNCPi0_{3}_{4}_Signal\")".format(sigDef,systName,i,universePrefix,i+1))
      # Pull out content of the one relevant bin
      exec("binVal = tHist_effDenom_{0}_{1}_{2}.GetBinContent(1)".format(sigDef,systName,i))
      # Populate corresponding hist in error band with bin content
      exec("mHist_effDenom_{0}.GetVertErrorBand(systName).GetHist(i).SetBinContent(1,binVal)".format(sigDef))

  ## Loop over flux systematics
  for systName in FLUX_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_effDenom_{0}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(sigDef))

    # Loop over universes in this category of systematic
    for i in range(nFluxUniverses):
      # Pull out the relevant TH1D and map from TH1D universe numbering (starting at 1) to MnvH1D universe number (starting at 0)
      exec("tHist_effDenom_{0}_{1}_{2} = effDenomFile_{0}.Get(\"Flux_XS_{1}_Dir/nu_uBooNE_AllNCPi0_universe_{3}_Signal\")".format(sigDef,systName,i,i+1))
      # Pull out content of the one relevant bin
      exec("binVal = tHist_effDenom_{0}_{1}_{2}.GetBinContent(1)".format(sigDef,systName,i))
      # Populate corresponding hist in error band with bin content
      exec("mHist_effDenom_{0}.GetVertErrorBand(systName).GetHist(i).SetBinContent(1,binVal)".format(sigDef))

  ## Loop over detector systematics
  for systName in DETECTOR_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_effDenom_{0}.AddVertErrorBandAndFillWithCV(systName,2)".format(sigDef))

    ## Detector variations don't actually affect the true event rate, so leave it at that

  for systName in G4_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_effDenom_{0}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(sigDef))

  for systName in OTHER_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_effDenom_{0}.AddVertErrorBandAndFillWithCV(systName,2)".format(sigDef))

  exec("writeHist(mHist_effDenom_{0},outFile)".format(sigDef))

  #############################################################################################################
  ### Construct Efficiency Numerator MnvH1D ###################################################################
  #############################################################################################################

  ### This is the actual CV
  #########################

  ## Pull out CV hists, which are separate for coherent and non-coherent
  exec("tHist_NCPi0Coh_{0} = cvFile_{0}.Get(\"nu_uBooNE_{0}_NCPi0Coh\")".format(sigDef))
  exec("tHist_NCPi0NotCoh_{0} = cvFile_{0}.Get(\"nu_uBooNE_{0}_NCPi0NotCoh\")".format(sigDef))
  
  ## Add together coherent and non-coherent
  exec("tHist_effNum_{0}_CV = tHist_NCPi0Coh_{0}.Clone(\"tHist_effNum_{0}_CV\")".format(sigDef))
  exec("tHist_effNum_{0}_CV.Add(tHist_NCPi0NotCoh_{0})".format(sigDef))

  ## Pull out the value and save as a scalar
  ## This is the actual CV in this bin for the analysis
  exec("binVal_trueCV = tHist_effNum_{0}_CV.GetBinContent(1)".format(sigDef))

  # Copy this into an MnvH1D (no systs yet)
  exec("mHist_effNum_{0} = ROOT.PlotUtils.MnvH1D(tHist_effNum_{0}_CV)".format(sigDef))
  # Rename new hist object
  exec("mHist_effNum_{0}.SetName(\"effNum_{0}\")".format(sigDef))

  ### This is the CV that will be used for internal comparisons against the reweightable systematics
  #########################

  exec("tHist_effNum_{0}_fakeCV_XS_GENIE = effNumFile_{0}.Get(\"Flux_XS_CV_Dir/nu_uBooNE_{0}_Signal\")".format(sigDef))
  exec("binVal_fakeCV_XS_GENIE = tHist_effNum_{0}_fakeCV_XS_GENIE.GetBinContent(1)".format(sigDef)) 

  ## Loop over cross section systematics
  for systName,universePrefix,nUniverses in XS_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_effNum_{0}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(sigDef))

    # Loop over universes in this category of systematic
    for i in range(nUniverses):
      # Pull out the relevant TH1D and map from TH1D universe numbering (starting at 1) to MnvH1D universe number (starting at 0)
      exec("tHist_effNum_{0}_{1}_{2} = effNumFile_{0}.Get(\"Flux_XS_{1}_Dir/nu_uBooNE_{0}_{3}_{4}_Signal\")".format(sigDef,systName,i,universePrefix,i+1))
      # Pull out content of the one relevant bin
      exec("binVal_variation = tHist_effNum_{0}_{1}_{2}.GetBinContent(1)".format(sigDef,systName,i))

      # Calculate appropriate shift (applied fractionally for many universe evaluation)
      fracShift = (binVal_variation-binVal_fakeCV_XS_GENIE)/binVal_fakeCV_XS_GENIE
      absShift = fracShift*binVal_trueCV

      # Populate corresponding hist in error band with bin content
      exec("mHist_effNum_{0}.GetVertErrorBand(systName).GetHist(i).SetBinContent(1,binVal_trueCV+absShift)".format(sigDef))

  ## Loop over flux systematics
  for systName in FLUX_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_effNum_{0}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(sigDef))

    # Loop over universes in this category of systematic
    for i in range(nFluxUniverses):
      # Pull out the relevant TH1D and map from TH1D universe numbering (starting at 1) to MnvH1D universe number (starting at 0)
      exec("tHist_effNum_{0}_{1}_{2} = effNumFile_{0}.Get(\"Flux_XS_{1}_Dir/nu_uBooNE_{0}_universe_{3}_Signal\")".format(sigDef,systName,i,i+1))
      # Pull out content of the one relevant bin
      exec("binVal_variation = tHist_effNum_{0}_{1}_{2}.GetBinContent(1)".format(sigDef,systName,i))

      # Calculate appropriate shift (applied fractionally for many universe evaluation)
      fracShift = (binVal_variation-binVal_fakeCV_XS_GENIE)/binVal_fakeCV_XS_GENIE
      absShift = fracShift*binVal_trueCV

      # Populate corresponding hist in error band with bin content
      exec("mHist_effNum_{0}.GetVertErrorBand(systName).GetHist(i).SetBinContent(1,binVal_trueCV+absShift)".format(sigDef))

  ## Loop over detector systematics
  for systName in DETECTOR_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_effNum_{0}.AddVertErrorBandAndFillWithCV(systName,2)".format(sigDef))

    # Pull out the TH1D corresponding to the CV specifically generated to go with the detector system variation
    exec("tHist_effNum_{0}_{1}_fakeCV = effNumFile_{0}.Get(\"{1}_CV_Dir/nu_uBooNE_{0}_Signal\")".format(sigDef,systName))
    # Pull out the TH1D corresponding to the detector system variation
    exec("tHist_effNum_{0}_{1}_variation = effNumFile_{0}.Get(\"{1}_{1}_Dir/nu_uBooNE_{0}_universe_1_Signal\")".format(sigDef,systName))

    # Pull out content of the one relevant bin from each distribution
    exec("binVal_fakeCV_local = tHist_effNum_{0}_{1}_fakeCV.GetBinContent(1)".format(sigDef,systName,i))
    exec("binVal_variation = tHist_effNum_{0}_{1}_variation.GetBinContent(1)".format(sigDef,systName,i))

    # Calculate appropriate shift (applied absolutely for +/- 1 sigma shifts
    fracShift = (binVal_variation-binVal_fakeCV_local)/binVal_fakeCV_local
    absShift = fracShift*binVal_trueCV

    # Populate corresponding hist in error band with bin content
    exec("mHist_effNum_{0}.GetVertErrorBand(systName).GetHist(0).SetBinContent(1,binVal_trueCV-absShift)".format(sigDef))
    exec("mHist_effNum_{0}.GetVertErrorBand(systName).GetHist(1).SetBinContent(1,binVal_trueCV+absShift)".format(sigDef))

  ## Loop over GEANT4 systematics
  for systName in G4_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_effNum_{0}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(sigDef))

    # Pull out the TH1D corresponding to the CV specifically generated to go with the GEANT4 variation
    exec("tHist_effNum_{0}_{1}_fakeCV = g4File_{0}.Get(\"GEANT4_CV_Dir/nu_uBooNE_{0}_Signal\")".format(sigDef,systName))
    exec("binVal_fakeCV_G4 = tHist_effNum_{0}_{1}_fakeCV.GetBinContent(1)".format(sigDef,systName))

    # Loop over universes in this category of systematic
    for i in range(nFluxUniverses):
      # Pull out the TH1D corresponding to the GEANT4 variation
      exec("tHist_effNum_{0}_{1}_{2} = g4File_{0}.Get(\"GEANT4_reinteractions_{1}_Geant4_Dir/nu_uBooNE_{0}_universe_{3}_Signal\")".format(sigDef,systName,i,i+1))

      # Pull out content of the one relevant bin from each distribution
      exec("binVal_variation = tHist_effNum_{0}_{1}_{2}.GetBinContent(1)".format(sigDef,systName,i))

      # Calculate appropriate shift (applied absolutely for +/- 1 sigma shifts
      fracShift = (binVal_variation-binVal_fakeCV_G4)/binVal_fakeCV_G4
      absShift = fracShift*binVal_trueCV

      # Populate corresponding hist in error band with bin content
      exec("mHist_effNum_{0}.GetVertErrorBand(systName).GetHist(i).SetBinContent(1,binVal_trueCV+absShift)".format(sigDef))

  ## Loop over other systematics
  for systName in OTHER_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_effNum_{0}.AddVertErrorBandAndFillWithCV(systName,2)".format(sigDef))

  exec("writeHist(mHist_effNum_{0},outFile)".format(sigDef))

  #############################################################################################################
  ### Construct Background MnvH1D #############################################################################
  #############################################################################################################

  ### This is the actual CV
  #########################

  ## Pull out CV hists, which are separate for various background components 
  for bkgDef in ["NCDelta","CC1Pi0","BNBOther","Nue","Dirt","OTPCExtra"]:
    exec("tHist_{0}_{1} = cvFile_{1}.Get(\"nu_uBooNE_{1}_{0}\")".format(bkgDef,sigDef))
  
  ## Add together various background components
  exec("tHist_background_{0}_CV = tHist_NCDelta_{0}.Clone(\"tHist_background_{0}_CV\")".format(sigDef))
  for bkgDef in ["CC1Pi0","BNBOther","Nue","Dirt","OTPCExtra"]:
    exec("tHist_background_{1}_CV.Add(tHist_{0}_{1})".format(bkgDef,sigDef))

  ## Pull out the value and save as a scalar
  exec("binVal_trueCV = tHist_background_{0}_CV.GetBinContent(1)".format(sigDef,systName,i))

  # Copy this into an MnvH1D (no systs yet)
  exec("mHist_background_{0} = ROOT.PlotUtils.MnvH1D(tHist_background_{0}_CV)".format(sigDef))
  # Rename new hist object
  exec("mHist_background_{0}.SetName(\"background_{0}\")".format(sigDef))

  ### This is the CV that will be used for internal comparisons against the reweightable systematics
  #########################

  exec("tHist_background_{0}_fakeCV_XS_GENIE = effNumFile_{0}.Get(\"Flux_XS_CV_Dir/nu_uBooNE_{0}_Bkgd\")".format(sigDef))
  exec("binVal_fakeCV_XS_GENIE = tHist_background_{0}_fakeCV_XS_GENIE.GetBinContent(1)".format(sigDef)) 

  ## TO-DO The below code is identical to the above chunk for the signal hists. These could instead 
  ## be implemented together with a loop over "Bkgd", "Signal"

  ## Loop over cross section systematics
  for systName,universePrefix,nUniverses in XS_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_background_{0}.AddVertErrorBandAndFillWithCV(systName,nUniverses)".format(sigDef))

    # Loop over universes in this category of systematic
    for i in range(nUniverses):
      # Pull out the relevant TH1D and map from TH1D universe numbering (starting at 1) to MnvH1D universe number (starting at 0)
      exec("tHist_background_{0}_{1}_{2} = effNumFile_{0}.Get(\"Flux_XS_{1}_Dir/nu_uBooNE_{0}_{3}_{4}_Bkgd\")".format(sigDef,systName,i,universePrefix,i+1))
      # Pull out content of the one relevant bin
      exec("binVal_variation = tHist_background_{0}_{1}_{2}.GetBinContent(1)".format(sigDef,systName,i))

      # Calculate appropriate shift (applied fractionally for many universe evaluation)
      fracShift = (binVal_variation-binVal_fakeCV_XS_GENIE)/binVal_fakeCV_XS_GENIE
      absShift = fracShift*binVal_trueCV

      # Populate corresponding hist in error band with bin content
      exec("mHist_background_{0}.GetVertErrorBand(systName).GetHist(i).SetBinContent(1,binVal_trueCV+absShift)".format(sigDef))

  ## Loop over flux systematics
  for systName in FLUX_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_background_{0}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(sigDef))

    # Loop over universes in this category of systematic
    for i in range(nFluxUniverses):
      # Pull out the relevant TH1D and map from TH1D universe numbering (starting at 1) to MnvH1D universe number (starting at 0)
      exec("tHist_background_{0}_{1}_{2} = effNumFile_{0}.Get(\"Flux_XS_{1}_Dir/nu_uBooNE_{0}_universe_{3}_Bkgd\")".format(sigDef,systName,i,i+1))
      # Pull out content of the one relevant bin
      exec("binVal_variation = tHist_background_{0}_{1}_{2}.GetBinContent(1)".format(sigDef,systName,i))

      # Calculate appropriate shift (applied fractionally for many universe evaluation)
      fracShift = (binVal_variation-binVal_fakeCV_XS_GENIE)/binVal_fakeCV_XS_GENIE
      absShift = fracShift*binVal_trueCV

      # Populate corresponding hist in error band with bin content
      exec("mHist_background_{0}.GetVertErrorBand(systName).GetHist(i).SetBinContent(1,binVal_trueCV+absShift)".format(sigDef))

  ## Loop over detector systematics
  for systName in DETECTOR_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_background_{0}.AddVertErrorBandAndFillWithCV(systName,2)".format(sigDef))

    # Pull out the TH1D corresponding to the CV specifically generated to go with the detector system variation
    exec("tHist_background_{0}_{1}_fakeCV = effNumFile_{0}.Get(\"{1}_CV_Dir/nu_uBooNE_{0}_Bkgd\")".format(sigDef,systName))
    # Pull out the TH1D corresponding to the detector system variation
    exec("tHist_background_{0}_{1}_variation = effNumFile_{0}.Get(\"{1}_{1}_Dir/nu_uBooNE_{0}_universe_1_Bkgd\")".format(sigDef,systName))

    # Pull out content of the one relevant bin from each distribution
    exec("binVal_fakeCV_local = tHist_background_{0}_{1}_fakeCV.GetBinContent(1)".format(sigDef,systName,i))
    exec("binVal_variation = tHist_background_{0}_{1}_variation.GetBinContent(1)".format(sigDef,systName,i))

    # Calculate appropriate shift
    fracShift = (binVal_variation-binVal_fakeCV_local)/binVal_fakeCV_local
    absShift = fracShift*binVal_trueCV

    # Populate corresponding hist in error band with bin content
    exec("mHist_background_{0}.GetVertErrorBand(systName).GetHist(0).SetBinContent(1,binVal_trueCV-absShift)".format(sigDef))
    exec("mHist_background_{0}.GetVertErrorBand(systName).GetHist(1).SetBinContent(1,binVal_trueCV+absShift)".format(sigDef))

  ## Loop over GEANT4 systematics
  for systName in G4_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_background_{0}.AddVertErrorBandAndFillWithCV(systName,nFluxUniverses)".format(sigDef))

    # Pull out the TH1D corresponding to the CV specifically generated to go with the GEANT4 variation
    exec("tHist_background_{0}_{1}_fakeCV = g4File_{0}.Get(\"GEANT4_CV_Dir/nu_uBooNE_{0}_Bkgd\")".format(sigDef,systName))
    exec("binVal_fakeCV_G4 = tHist_background_{0}_{1}_fakeCV.GetBinContent(1)".format(sigDef,systName))

    # Loop over universes in this category of systematic
    for i in range(nFluxUniverses):
      # Pull out the TH1D corresponding to the GEANT4 variation
      exec("tHist_background_{0}_{1}_{2} = g4File_{0}.Get(\"GEANT4_reinteractions_{1}_Geant4_Dir/nu_uBooNE_{0}_universe_{3}_Bkgd\")".format(sigDef,systName,i,i+1))

      # Pull out content of the one relevant bin from each distribution
      exec("binVal_variation = tHist_background_{0}_{1}_{2}.GetBinContent(1)".format(sigDef,systName,i))

      # Calculate appropriate shift (applied absolutely for +/- 1 sigma shifts
      fracShift = (binVal_variation-binVal_fakeCV_G4)/binVal_fakeCV_G4
      absShift = fracShift*binVal_trueCV

      # Populate corresponding hist in error band with bin content
      exec("mHist_background_{0}.GetVertErrorBand(systName).GetHist(i).SetBinContent(1,binVal_trueCV+absShift)".format(sigDef))

  ## Loop over other systematics
  for systName in OTHER_SYSTS:

    # Create the appropriate error band in the MnvH1D
    exec("mHist_background_{0}.AddVertErrorBandAndFillWithCV(systName,2)".format(sigDef))

  exec("writeHist(mHist_background_{0},outFile)".format(sigDef))

#############################################################################################################
### Loop over 2g1p, 2g0p ####################################################################################
#############################################################################################################

for sigDef in ["2g1p","2g0p"]:

  #############################################################################################################
  ### Calculate Things ########################################################################################
  #############################################################################################################

  ## Efficiency
  exec("mHist_eff_{0} = mHist_effNum_{0}.Clone(\"eff_{0}\")".format(sigDef))
  exec("mHist_eff_{0}.Divide(mHist_effNum_{0},mHist_effDenom_{0})".format(sigDef))

  exec("writeHist(mHist_eff_{0},outFile)".format(sigDef))

  ## Background-subtracted event rate
  exec("mHist_evtRate_{0} = mHist_data_selected_{0}.Clone(\"evtRate_{0}\")".format(sigDef))
  exec("mHist_evtRate_{0}.Add(mHist_BNB_ext_{0},-1.)".format(sigDef))
  exec("mHist_evtRate_{0}.Add(mHist_background_{0},-1.)".format(sigDef))

  exec("writeHist(mHist_evtRate_{0},outFile)".format(sigDef))

  ## Cross section calculation
  exec("mHist_xSection_{0} = mHist_evtRate_{0}.Clone(\"xSection_{0}\")".format(sigDef))
  exec("mHist_xSection_{0}.Divide(mHist_xSection_{0},mHist_eff_{0})".format(sigDef))
  exec("mHist_xSection_{0}.Divide(mHist_xSection_{0},mHist_flux_integral)".format(sigDef))
  exec("mHist_xSection_{0}.Divide(mHist_xSection_{0},mHist_POT_{0})".format(sigDef)) # Remove units of per POT
  exec("mHist_xSection_{0}.Divide(mHist_xSection_{0},mHist_nTargets)".format(sigDef))
  
  exec("writeHist(mHist_xSection_{0},outFile)".format(sigDef))

  ## MC signal prediction
  exec("mHist_xSection_mc_{0} = mHist_effNum_{0}.Clone(\"xSection_mc_{0}\")".format(sigDef))
  ## We don't want to subtract off the BNB_ext and background prediction because 
  ## the MC signal doesn't include backgrounds
  exec("mHist_xSection_mc_{0}.Divide(mHist_xSection_mc_{0},mHist_eff_{0})".format(sigDef))
  exec("mHist_xSection_mc_{0}.Divide(mHist_xSection_mc_{0},mHist_flux_integral)".format(sigDef))
  exec("mHist_xSection_mc_{0}.Divide(mHist_xSection_mc_{0},mHist_POT_{0})".format(sigDef)) # Remove units of per POT
  exec("mHist_xSection_mc_{0}.Divide(mHist_xSection_mc_{0},mHist_nTargets)".format(sigDef))

  ## Pop out all error bands except GENIE from the MC xsection
  for systName in FLUX_SYSTS + DETECTOR_SYSTS + G4_SYSTS + OTHER_SYSTS:
    exec("mHist_xSection_mc_{0}.PopVertErrorBand(\"{1}\")".format(sigDef,systName))

  exec("writeHist(mHist_xSection_mc_{0},outFile)".format(sigDef))

#############################################################################################################
### Close output file; other business #######################################################################
#############################################################################################################
outFile.Close()
