#include "TMatrixD.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "PlotUtils/MnvH1D.h"

#include "WienerSVDUnfolder.hh"
#include "DAgostiniUnfolder.hh"

// Create TMatrixD from TH2D
// Used to allow for an input TH2D that has an opposite vertical/horizontal axis convention              
TMatrixD TH2DtoTMatrixD(const TH2D& input_hist, bool invert_matrix, bool include_underflow_X, bool include_overflow_X, bool include_underflow_Y, bool include_overflow_Y)
{
  // Number of bins are increased by one unit for each of the underflow and overflow that is included
  Int_t nRows = input_hist.GetNbinsX();
  Int_t nColumns = input_hist.GetNbinsY();
  
  if(include_underflow_X) nRows = nRows+1;
  if(include_overflow_X)  nRows = nRows+1;

  if(include_underflow_Y) nColumns = nColumns+1;
  if(include_overflow_Y)  nColumns = nColumns+1;
   
  TMatrixD output_matrix(nRows,nColumns);
  // The following loops start at zero only if we're including underflow
  Int_t offset_X = include_underflow_X ? 0 : 1;  
  Int_t offset_Y = include_underflow_Y ? 0 : 1;
  for(Int_t i=0; i<nRows; i++)
  {
    for(Int_t j=0; j<nColumns; j++)
    {
      if(invert_matrix) output_matrix(i, j) = input_hist.GetBinContent(j+offset_Y, i+offset_X);
      else              output_matrix(i, j) = input_hist.GetBinContent(i+offset_X, j+offset_Y);
    }
  }
  return output_matrix;
}

// Create TMatrixD from TH1D
TMatrixD TH1DtoTMatrixD(const TH1D& input_hist, bool include_underflow, bool include_overflow)
{
  // Number of bins are increased by one unit for each of the underflow and overflow that is included
  Int_t nRows = input_hist.GetNbinsX();

  if(include_underflow) nRows = nRows+1;
  if(include_overflow)  nRows = nRows+1;

  TMatrixD output_matrix(nRows,1);
  // The following loop starts at zero only if we're including underflow
  Int_t offset = include_underflow ? 0 : 1;  
  for(Int_t i=0; i<nRows; i++)
    {
      output_matrix(i, 0) = input_hist.GetBinContent(i+offset);
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
    TFile* file_in = new TFile((filePath_in+".root").c_str(),"READ");
    file_in->Cp((filePath_in+"_unfolded.root").c_str());
    file_in->Close();
    TFile* file_out = new TFile((filePath_in+"_unfolded.root").c_str(),"UPDATE");
 
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
 
    // Decide whether to include or exclude underflow/overflow
    // Think about a good tolerance to compare the bin content against
    double threshold = 10e-6;

    // Check if underflow is zero for data_signal (reco space)
    double binVal_underflow_reco = tHist_data_signal.GetBinContent(0);
    bool include_underflow_reco = binVal_underflow_reco > threshold ? 1 : 0;

    // Check if overflow is zero for data_signal (reco space)
    Int_t nBins_reco = tHist_data_signal.GetNbinsX();
    double binVal_overflow_reco = tHist_data_signal.GetBinContent(nBins_reco+1);
    bool include_overflow_reco = binVal_overflow_reco > threshold ? 1 : 0;

    // Check if underflow is zero for prior_true_signal (reco space)
    double binVal_underflow_true = tHist_prior_true_signal.GetBinContent(0);
    bool include_underflow_true = binVal_underflow_true > threshold ? 1 : 0;

    // Check if overflow is zero for prior_true_signal (reco space)
    Int_t nBins_true = tHist_prior_true_signal.GetNbinsX();
    double binVal_overflow_true = tHist_prior_true_signal.GetBinContent(nBins_true+1);
    bool include_overflow_true = binVal_overflow_true > threshold ? 1 : 0;

    include_underflow_reco = 0;
    include_overflow_reco = 0;
    include_underflow_true = 0;
    include_overflow_true = 0;   

    // Convert inputs into TMatrixD
    TMatrixD tMat_data_signal = TH1DtoTMatrixD(tHist_data_signal, include_underflow_reco, include_overflow_reco);
    TMatrixD tMat_prior_true_signal = TH1DtoTMatrixD(tHist_prior_true_signal, include_underflow_true, include_overflow_true);
    TMatrixD tMat_response = TH2DtoTMatrixD(*tHist2D_response, false, include_underflow_reco, include_overflow_reco, include_underflow_true, include_overflow_true);

    // If either include_underflow_reco or include_overflow_reco is false, replace tMat_data_covmat with tMat_data_covmat->GetSub(x1,x2,y1,y2);
    // Usage: TMatrixT< Element > GetSub (Int_t row_lwb, Int_t row_upb, Int_t col_lwb, Int_t col_upb, Option_t *option="S") const 
    Int_t lowerBound = include_underflow_reco ? 0 : 1;
    Int_t upperBound = include_overflow_reco ? nBins_reco+1 : nBins_reco;
    TMatrixD tMat_data_covmat_final = tMat_data_covmat.GetSub(lowerBound, upperBound, lowerBound, upperBound);   

    // Print out number of rows and columns in each input matrix
    Int_t data_rows = tMat_data_signal.GetNrows();
    Int_t data_cols = tMat_data_signal.GetNcols();
    std::cout << "data_rows: " << data_rows << "\t" << "data_cols: " << data_cols << std::endl;

    Int_t cov_final_rows = tMat_data_covmat_final.GetNrows();
    Int_t cov_final_cols = tMat_data_covmat_final.GetNcols();
    std::cout << "cov_final_rows: " << cov_final_rows << "\t" << "cov_final_cols: " << cov_final_cols << std::endl;

    Int_t response_rows = tMat_response.GetNrows();
    Int_t response_cols = tMat_response.GetNcols();
    std::cout << "response_rows: " << response_rows << "\t" << "response_cols: " << response_cols << std::endl;

    Int_t mc_rows = tMat_prior_true_signal.GetNrows();
    Int_t mc_cols = tMat_prior_true_signal.GetNcols();
    std::cout << "mc_rows: " << mc_rows << "\t" << "mc_cols: " << mc_cols << std::endl;

    // Initialize unfolder
    //constexpr int NUM_DAGOSTINI_ITERATIONS = 4;
    std::unique_ptr< Unfolder > unfolder (
					  new WienerSVDUnfolder( true,
            WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv )
            //new DAgostiniUnfolder( NUM_DAGOSTINI_ITERATIONS ) 
					  );
    
    // Run unfolder
    auto result = unfolder->unfold( tMat_data_signal, tMat_data_covmat_final, tMat_response, tMat_prior_true_signal );

    //Pull out unfolded signal and covariance from results
    auto tMat_data_signal_unfolded = result.unfolded_signal_.get();
    auto tMat_unfolded_covariance = result.cov_matrix_.get();
    auto tMat_errprop_matrix = result.err_prop_matrix_.get();

    // Convert outputs into TH1D/TH2D  
    TH1D tHist_data_signal_unfolded = TMatrixDtoTH1D(*tMat_data_signal_unfolded, tHist_prior_true_signal);
    TH2D tHist2D_unfolded_covariance = TMatrixDtoTH2D(*tMat_unfolded_covariance, tHist_prior_true_signal);
    TH2D tHist2D_covariance = TMatrixDtoTH2D(tMat_data_covmat_final, tHist_data_signal);

    // Construct MnvH1D
    PlotUtils::MnvH1D mHist_data_signal_unfolded = PlotUtils::MnvH1D(tHist_data_signal_unfolded);
    mHist_data_signal_unfolded.AddMissingErrorBandsAndFillWithCV(*mHist_data_signal_folded);
    
    mHist_data_signal_unfolded.SetName("unfolded_evtRate_2g1p_inclusive");
    mHist_data_signal_unfolded.Write();   
    tHist2D_unfolded_covariance.SetName("unfolded_cov_evtRate_2g1p_inclusive");
    tHist2D_unfolded_covariance.Write();
    tHist2D_covariance.SetName("cov_evtRate_2g1p_inclusive");
    tHist2D_covariance.Write();
    file_out->Close();
    return;
}

