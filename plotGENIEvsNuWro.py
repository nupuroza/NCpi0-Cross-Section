import ROOT
import argparse
import os
from customHistAndPlotMethods import makeEnv_TCanvas

## Plots GENIE and NuWro signal, background, and combined 2g1p and 2g0p reco distributions with associated errors.

## Default output directory for Leon. Feel free to change it locally. 
out_dir = "/uboone/data/users/ltong/gLEE/NCPi0/GENIEvsNuWro"

## Use first argument from command line as out_dir if provided.
parser = argparse.ArgumentParser(description='Script to plot GENIE and NuWro signal, background, and combined 2g1p and 2g0p reco distributions with associated errors')
parser.add_argument('out_dir', help='Path to ouput directory', type=str,nargs='?')
p = parser.parse_args()
if p.out_dir > 0:
  out_dir = p.out_dir

## Create output directory if it doesn't exist
if not os.path.isdir(out_dir):
  print "Making output directory {0}".format(out_dir)
  os.system( "mkdir %s" % out_dir )

## Set batch mode
ROOT.gROOT.SetBatch()

## This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

## GENIE and NuWro paths and files
fGENIE_2g1p_exclusive_path = "/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/MajorMerge_GGE_mark/working_dir/ToTH1D/VersionNextGen_SamePOT_Aug2023/variation_spectra/SBNfit_variation_spectra_exclusive_2g1p.root"
fGENIE_2g0p_exclusive_path = "/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/MajorMerge_GGE_mark/working_dir/ToTH1D/VersionNextGen_SamePOT_Aug2023/variation_spectra/SBNfit_variation_spectra_exclusive_2g0p.root"
fNuWro_path = "/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/MajorMerge_GGE_mark/working_dir/ToTH1D/NuWro_FakeData_Generation_NEW/NuWro_Apr2023_v6_CV.SBNspec.root"
fGENIE_2g1p_exclusive = ROOT.TFile(fGENIE_2g1p_exclusive_path, "read")
fGENIE_2g0p_exclusive = ROOT.TFile(fGENIE_2g0p_exclusive_path, "read")
fNuWro = ROOT.TFile(fNuWro_path, "read")

## Retrieve GENIE and NuWro reco signal, background, and fake data histograms. GENIE "fakedata" is simply the sum of its reco signal and background histograms
GENIE_2g1p_exclusive_signal = fGENIE_2g1p_exclusive.Get("exclusive_2g1p_CV_Dir/Sys2g1p_numerator_reco_Signal")
GENIE_2g0p_exclusive_signal = fGENIE_2g0p_exclusive.Get("exclusive_2g0p_CV_Dir/Sys2g0p_numerator_reco_Signal")
GENIE_2g1p_exclusive_background = fGENIE_2g1p_exclusive.Get("exclusive_2g1p_CV_Dir/Sys2g1p_numerator_reco_Bkgd")
GENIE_2g0p_exclusive_background = fGENIE_2g0p_exclusive.Get("exclusive_2g0p_CV_Dir/Sys2g0p_numerator_reco_Bkgd")
GENIE_2g1p_exclusive_fakedata = GENIE_2g1p_exclusive_signal + GENIE_2g1p_exclusive_background
GENIE_2g0p_exclusive_fakedata = GENIE_2g0p_exclusive_signal + GENIE_2g0p_exclusive_background
NuWro_2g1p_exclusive_signal = fNuWro.Get("nu_uBooNE_breakdown_2g1psig")
NuWro_2g0p_exclusive_signal = fNuWro.Get("nu_uBooNE_breakdown_2g0psig")
NuWro_2g1p_exclusive_background = fNuWro.Get("nu_uBooNE_breakdown_2g1pbkg")
NuWro_2g0p_exclusive_background = fNuWro.Get("nu_uBooNE_breakdown_2g0pbkg")
NuWro_2g1p_exclusive_fakedata = fNuWro.Get("nu_uBooNE_fakedata_2g1p")
NuWro_2g0p_exclusive_fakedata = fNuWro.Get("nu_uBooNE_fakedata_2g0p")

## Loop over all plots to be made
for np in ["2g1p", "2g0p"]:
    for datatype in ["signal", "background", "fakedata"]:
        with makeEnv_TCanvas(out_dir + "/{0}_exclusive_{1}.png".format(np, datatype)) as canvas:

            ## Clone GENIE histograms with different variable names to the same variable within the for loop for conciseness of code
            exec("local_GENIE_exclusive = GENIE_{0}_exclusive_{1}.Clone(\"local_Sys{0}_numerator_reco_{1}\")".format(np, datatype))

            ## Set GENIE histogram Draw properties
            local_GENIE_exclusive.SetLineColor(ROOT.kRed)
            local_GENIE_exclusive.SetMarkerColor(ROOT.kRed)
            local_GENIE_exclusive.SetFillColor(ROOT.kRed)
            local_GENIE_exclusive.SetFillStyle(0)
            local_GENIE_exclusive.SetLineWidth(2)
            local_GENIE_exclusive.GetXaxis().SetTitle("Reco #pi^{{0}} momentum (GeV) [{0}]".format(np))
            local_GENIE_exclusive.GetYaxis().SetTitle("# Reco {0} Events".format(datatype))
            local_GENIE_exclusive.SetFillStyle(0)

            ## Get content and error of GENIE bin with max y value
            GENIE_maxybin = local_GENIE_exclusive.GetMaximumBin()
            GENIE_maxy = local_GENIE_exclusive.GetBinContent(GENIE_maxybin)
            GENIE_maxyerror = local_GENIE_exclusive.GetBinError(GENIE_maxybin)

            ## Clone NuWro histograms with different variable names to the same variable within the for loop for conciseness of code
            exec("local_NuWro_exclusive = NuWro_{0}_exclusive_{1}.Clone(\"local_nu_uBooNE_reco_{1}_{0}\")".format(np, datatype))

            ## Scale NuWro histograms to GENIE POT. 2g1p and 2g0p treated separately due to having different POT historically.
            if np == "2g1p":
                local_NuWro_exclusive.Scale(6.868/3.0041393)
            elif np == "2g0p":
                local_NuWro_exclusive.Scale(6.868/3.0041393)
            else:
                print("Only 2g1p and 2g0p configurations are supported.")
                exit(0)

            ## Set NuWro histogram Draw properties. Plotted on same canvas as GENIE histogram, so doesn't need own titles.
            local_NuWro_exclusive.SetLineWidth(2)
            local_NuWro_exclusive.SetMarkerSize(2)
            local_NuWro_exclusive.SetMarkerStyle(20)

            ## Get content and error of NuWro bin with max y value
            NuWro_maxybin = local_NuWro_exclusive.GetMaximumBin()
            NuWro_maxy = local_NuWro_exclusive.GetBinContent(NuWro_maxybin)
            NuWro_maxyerror = local_NuWro_exclusive.GetBinError(NuWro_maxybin)

            ## Choose maximum of GENIE or NuWro maxy + maxyerror and multiply by 1.1 for maximum y axis value of combined histogram. Set minimum to 0.
            actual_maxy = max(GENIE_maxy + GENIE_maxyerror, NuWro_maxy + NuWro_maxyerror)*1.1
            local_GENIE_exclusive.SetMaximum(actual_maxy)
            local_GENIE_exclusive.SetMinimum(0)

            ## Draw GENIE histogram, GENIE error, and NuWro histogram with error. GENIE histogram drawn first drawn without errors; then error is drawn with hashed fill style.
            ## DrawCopy() prevents first drawn local_GENIE_exclusive histogram from changing its fill style once drawn.
            local_GENIE_exclusive.DrawCopy("HIST")
            local_GENIE_exclusive.SetFillStyle(3545)
            local_GENIE_exclusive.Draw("E2 SAME")
            local_NuWro_exclusive.Draw("E1P SAME")

            ## Draw borderless legend.
            legend = ROOT.TLegend(0.555, 0.5, 0.90, 0.7)
            legend.SetBorderSize(0)
            legend.AddEntry(local_GENIE_exclusive, "GENIE")
            legend.AddEntry(local_NuWro_exclusive, "NuWro")
            legend.Draw()
