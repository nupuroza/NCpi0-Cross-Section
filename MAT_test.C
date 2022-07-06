#include "TMatrixD.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "PlotUtils/MnvH1D.h"

#include "WienerSVDUnfolder.hh"
#include "DAgostiniUnfolder.hh"

// Create TMatrixD from TH2D
// If TH2D(i, j) = Matrix(i, j), invert_matrix = kFALSE, else invert_matrix = kTRUE
// Used to allow for an input TH2D that has an opposite vertical/horizontal axis convention              
TMatrixD TH2DtoTMatrixD(const TH2D& input_hist, bool invert_matrix)
{
  Int_t nBinsX = input_hist.GetNbinsX();
  Int_t nBinsY = input_hist.GetNbinsY();
  TMatrixD output_matrix(nBinsX,nBinsY);
  for(Int_t i=0; i<nBinsX; i++)
  {
    for(Int_t j=0; j<nBinsY; j++)
    {
      if(invert_matrix) output_matrix(j, i) = input_hist.GetBinContent(i+1, j+1);
      else          output_matrix(i, j) = input_hist.GetBinContent(i+1, j+1);
    }
  }
  return output_matrix;
}

// Create TMatrixD from TH1D
// Pad TMatrixD to include underflow/overflow bins
TMatrixD TH1DtoTMatrixD(const TH1D& input_hist)
{
  Int_t nBins = input_hist.GetNbinsX()+2;
  TMatrixD output_matrix(nBins,1);
  for(Int_t i=0; i<nBins; i++)
    {
      output_matrix(i, 0) = input_hist.GetBinContent(i);
    }
  return output_matrix;
}

// Create TH2D from TMatrixD
// reference_hist should have the desired binning, but will not
// itself be filled
TH2D TMatrixDtoTH2D(const TMatrixD& input_mat, const TH1D& reference_hist)
{
  Int_t nBins = reference_hist.GetNbinsX();
  Double_t binEdges[nBins+1];
  for(int i=0; i<nBins+1; i++){
      binEdges[i] = reference_hist.GetBinLowEdge(i+1);
  }
  TH2D output_hist("","",nBins,binEdges,nBins,binEdges);
  for(Int_t i=0; i<input_mat.GetNrows(); i++)
    {
      for(Int_t j=0; j<input_mat.GetNcols(); j++)
      {
	      output_hist.SetBinContent(i+1, j+1, input_mat(i, j));
      }
    }
  return output_hist;
}

// Create TH1D from TMatrixD
// reference_hist should have the desired binning, but will not
// itself be filled
TH1D TMatrixDtoTH1D(const TMatrixD& input_mat, const TH1D& reference_hist)
{
  Int_t nBins = reference_hist.GetNbinsX();
  Double_t binEdges[nBins+1];
  for(int i=0; i<nBins+1; i++){
      binEdges[i] = reference_hist.GetBinLowEdge(i+1);
  }
  TH1D output_hist("","",nBins,binEdges);
  for(Int_t i=0; i<input_mat.GetNrows(); i++)
    {
      output_hist.SetBinContent(i+1, input_mat.operator()(i, 0));
    }
  return output_hist;
}

using namespace std;

void MAT_test()
{
    TFile* f = new TFile("/uboone/data/users/noza/gLEE/xsection/2022-06-29_out.root", "READ");
    
    // Pull out measured signal MnvH1D from input file
    PlotUtils::MnvH1D *mHist_data_signal = (PlotUtils::MnvH1D*)f->Get("evtRate_2g1p_inclusive");
    TH1D tHist_data_signal = mHist_data_signal->GetCVHistoWithStatError();
    // Extract covariance 
    TMatrixD tMat_data_covmat = mHist_data_signal->GetTotalErrorMatrix();
    // Pull out response matrix from input file
    TH2D* tHist2D_response = (TH2D*)f->Get("response_2g1p_inclusive");
    // Pull out predicted signal MnvH1D from input file
    PlotUtils::MnvH1D *mHist_prior_true_signal = (PlotUtils::MnvH1D*)f->Get("effNum_2g1p_inclusive");  
    TH1D tHist_prior_true_signal = mHist_prior_true_signal->GetCVHistoWithStatError();

    // Convert inputs into TMatrixD
    TMatrixD tMat_data_signal = TH1DtoTMatrixD(tHist_data_signal);
    TMatrixD tMat_prior_true_signal = TH1DtoTMatrixD(tHist_prior_true_signal);
    TMatrixD tMat_response = TH2DtoTMatrixD(*tHist2D_response, kTRUE);

    // DEBUG
    Int_t data_rows = tMat_data_signal.GetNrows();
    Int_t data_cols = tMat_data_signal.GetNcols();

    std::cout << "data_rows: " << data_rows << "\t" << "data_cols: " << data_cols << std::endl;

    // Initialize unfolder
    // constexpr int NUM_DAGOSTINI_ITERATIONS = 6;
    std::unique_ptr< Unfolder > unfolder (
					  new WienerSVDUnfolder( true,
            WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv )
            //new DAgostiniUnfolder( NUM_DAGOSTINI_ITERATIONS ) 
					  );

    // Run unfolder
    auto result = unfolder->unfold( tMat_data_signal, tMat_data_covmat, tMat_response, tMat_prior_true_signal );

    //Pull out unfolded signal and covariance from results
    auto tMat_unfolded_signal = result.unfolded_signal_.get();
    auto tMat_unfolded_covariance = result.cov_matrix_.get();
    auto tMat_errprop_matrix = result.err_prop_matrix_.get();

    // Convert outputs into TH1D/TH2D  
    TH1D tHist_unfolded_signal = TMatrixDtoTH1D(*tMat_unfolded_signal, tHist_prior_true_signal);
    TH2D tHist2D_unfolded_covariance = TMatrixDtoTH2D(*tMat_unfolded_covariance, tHist_prior_true_signal);

    // Construct MnvH1D
     
    // Write unfolder output to file
    TFile* file = new TFile("/uboone/data/users/noza/gLEE/xsection/2022-06-30_unfolded.root", "RECREATE"); 

    tHist_unfolded_signal.Write();
    tHist2D_unfolded_covariance.Write();
    file->Close();

    return;
}

