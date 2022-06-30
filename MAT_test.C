#include "TMatrixD.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "PlotUtils/MnvH1D.h"


#include "WienerSVDUnfolder.hh"
#include "DAgostiniUnfolder.hh"

// Define functions to convert histograms between TH1D/TH2D & TMatrix

// Fill 2D histogram into matrix 
/// If TH2D(i, j) = Matrix(i, j), rowcolumn = kTRUE, else rowcolumn = kFALSE              
void TH2DtoTMatrix(const TH2D* histo, TMatrixD& mat, bool rowcolumn)
{
    for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
        for(Int_t j=0; j<histo->GetNbinsY(); j++)
        {
            if(rowcolumn) mat(i, j) = histo->GetBinContent(i+1, j+1);
            else mat(j, i) = histo->GetBinContent(i+1, j+1);
        }
    }
}

// Fill 1D histogram into matrix 
// If TH1D(i, j) = Matrix(i, j), rowcolumn = kTRUE, else rowcolumn = kFALSE
void TH1DtoTMatrix(const TH1D* histo, TMatrixD& mat)
{
  for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
      mat(i, 0) = histo->GetBinContent(i+1);
    }
}

//Fill matrix to 2D histogram
void TMatrixtoTH2D(const TMatrixD& mat, TH2D* histo)
{
  for(Int_t i=0; i<mat.GetNrows(); i++)
    {
      for(Int_t j=0; j<mat.GetNcols(); j++)
        {
	  histo->SetBinContent(i+1, j+1, mat(i, j));
        }
    }
}

//Fill matrix to 1D histogram
void TMatrixtoTH1D(const TMatrixD& mat, TH1D* histo)
{
  for(Int_t i=0; i<mat->GetNrows(); i++)
    {
      histo->SetBinContent(i+1, mat->operator()(i, 0));
    }
}


using namespace std;

void MAT_test()
{
    TFile* f = new TFile("/uboone/data/users/noza/gLEE/xsection/2022-06-29_out.root", "READ");
    
    // Pull out measured signal MnvH1D from input file
    PlotUtils::MnvH1D *mHist_data_signal = (PlotUtils::MnvH1D*)f->Get("evtRate_2g1p_inclusive");
    TH1D Data_signal = mHist_data_signal->GetCVHistoWithStatError();
    // Extract covariance 
    TMatrixD data_covmat = mHist_data_signal->GetTotalErrorMatrix();
    // Pull out response matrix from input file
    TH2D* Response = (TH2D*)f->Get("response_2g1p_inclusive");
    // Pull out predicted signal MnvH1D from input file
    PlotUtils::MnvH1D *mHist_prior_true_signal = (PlotUtils::MnvH1D*)f->Get("effNum_2g1p_inclusive");  
    TH1D Prior_true_signal = mHist_prior_true_signal->GetCVHistoWithStatError();

    // Store truth #bins/binning
    Int_t n = Prior_true_signal->GetNbinsX();
    Double_t Nuedges[n+1];
    for(int i=0; i<n+1; i++){
        Nuedges[i] = Prior_true_signal->GetBinLowEdge(i+1);
    }
    // Store reco #bins
    Int_t m = Data_signal->GetNbinsX();

    // Construct matrices for input
    TMatrixD data_signal(m,1);
    TMatrixD prior_true_signal(n,1);
    TMatrixD response(m, n);
   
    // Convert inputs into TMatrixD
    TH1DtoTMatrix(Data_signal, data_signal);
    TH1DtoTMatrix(Prior_true_signal, prior_true_signal);
    TH2DtoTMatrix(Response, response, kFALSE);


    // Initialize unfolder
    // constexpr int NUM_DAGOSTINI_ITERATIONS = 6;
    std::unique_ptr< Unfolder > unfolder (
					  new WienerSVDUnfolder( true,
            WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv )
            //new DAgostiniUnfolder( NUM_DAGOSTINI_ITERATIONS ) 
					  );

    // Run unfolder
    auto result = unfolder->unfold( data_signal, data_covmat, response, prior_true_signal );

    //Pull out unfolded signal and covariance from results
    auto unfolded_signal = result.unfolded_signal_.get();
    auto unfolded_covariance = result.cov_matrix_.get();

    // Construct matrices for output
    TMatrixD Unfolded_signal(n,1);
    TMatrixD Unfolded_covariance(n, n);

    // Convert outputs into TH1D/TH2D  
    TMatrixtoTH1D(unfolded_signal, Unfolded_signal);
    TMatrixtoTH2D(unfolded_covariance, Unfolded_covariance);
 
    // Write unfolder output to file
    TFile* file = new TFile("/uboone/data/users/noza/gLEE/xsection/2022-06-30_unfolded.root", "RECREATE"); 

    Unfolded_signal->Write();
    Unfolded_covariance->Write();
    file->Close();

    return;
}

