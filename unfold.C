#include "TMatrixD.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "PlotUtils/MnvH1D.h"

#include "WienerSVDUnfolder.hh"
#include "DAgostiniUnfolder.hh"

// Create TMatrixD from TH2D
// Pad TMatrixD to include underflow/overflow bins
// If TH2D(i, j) = Matrix(i, j), invert_matrix = kFALSE, else invert_matrix = kTRUE
// Used to allow for an input TH2D that has an opposite vertical/horizontal axis convention              
TMatrixD TH2DtoTMatrixD(const TH2D& input_hist, bool invert_matrix)
{
  Int_t nBinsX = input_hist.GetNbinsX()+2;
  Int_t nBinsY = input_hist.GetNbinsY()+2;
  TMatrixD output_matrix(nBinsX,nBinsY);
  for(Int_t i=0; i<nBinsX; i++)
  {
    for(Int_t j=0; j<nBinsY; j++)
    {
      if(invert_matrix) output_matrix(j, i) = input_hist.GetBinContent(i, j);
      else          output_matrix(i, j) = input_hist.GetBinContent(i, j);
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

void unfold(std::string filePath_in)
{
    TFile* file_in = new TFile(filePath_in.c_str(),"READ");
    file_in->Cp("test_output.root");
    file_in->Close();
    TFile* file_out = new TFile("test_output.root","UPDATE");
 
    // Pull out measured signal MnvH1D from input file
    PlotUtils::MnvH1D *mHist_data_signal_folded = (PlotUtils::MnvH1D*)file_out->Get("evtRate_2g1p_inclusive");
    TH1D tHist_data_signal = mHist_data_signal_folded->GetCVHistoWithStatError();
    // Extract covariance 
    TMatrixD tMat_data_covmat = mHist_data_signal_folded->GetTotalErrorMatrix();
    // Pull out response matrix from input file
    TH2D* tHist2D_response = (TH2D*)file_out->Get("response_2g1p_inclusive");
    // Pull out predicted signal MnvH1D from input file
    PlotUtils::MnvH1D *mHist_prior_true_signal = (PlotUtils::MnvH1D*)file_out->Get("effNum_2g1p_inclusive");  
    TH1D tHist_prior_true_signal = mHist_prior_true_signal->GetCVHistoWithStatError();

    // Convert inputs into TMatrixD
    TMatrixD tMat_data_signal = TH1DtoTMatrixD(tHist_data_signal);
    TMatrixD tMat_prior_true_signal = TH1DtoTMatrixD(tHist_prior_true_signal);
    TMatrixD tMat_response = TH2DtoTMatrixD(*tHist2D_response, kTRUE);


    // DEBUG
    Int_t data_rows = tMat_data_signal.GetNrows();
    Int_t data_cols = tMat_data_signal.GetNcols();

    std::cout << "data_rows: " << data_rows << "\t" << "data_cols: " << data_cols << std::endl;

    Int_t cov_rows = tMat_data_covmat.GetNrows();
    Int_t cov_cols = tMat_data_covmat.GetNcols();

    std::cout << "cov_rows: " << cov_rows << "\t" << "cov_cols: " << cov_cols << std::endl;

    Int_t response_rows = tMat_response.GetNrows();
    Int_t response_cols = tMat_response.GetNcols();

    std::cout << "response_rows: " << response_rows << "\t" << "response_cols: " << response_cols << std::endl;

    Int_t mc_rows = tMat_prior_true_signal.GetNrows();
    Int_t mc_cols = tMat_prior_true_signal.GetNcols();

    std::cout << "mc_rows: " << mc_rows << "\t" << "mc_cols: " << mc_cols << std::endl;


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
    auto tMat_data_signal_unfolded = result.unfolded_signal_.get();
    auto tMat_unfolded_covariance = result.cov_matrix_.get();
    auto tMat_errprop_matrix = result.err_prop_matrix_.get();

    // Convert outputs into TH1D/TH2D  
    TH1D tHist_data_signal_unfolded = TMatrixDtoTH1D(*tMat_data_signal_unfolded, tHist_prior_true_signal);
    TH2D tHist2D_unfolded_covariance = TMatrixDtoTH2D(*tMat_unfolded_covariance, tHist_prior_true_signal);

    // Construct MnvH1D
    PlotUtils::MnvH1D mHist_data_signal_unfolded = PlotUtils::MnvH1D(tHist_data_signal_unfolded);
    mHist_data_signal_unfolded.AddMissingErrorBandsAndFillWithCV(*mHist_data_signal_folded);
    
    file_out->cd(); 
    mHist_data_signal_unfolded.Write();
    tHist2D_unfolded_covariance.Write();
    file_out->Close();

    return;
}

