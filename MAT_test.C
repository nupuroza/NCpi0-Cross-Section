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
TMatrixD TH1DtoTMatrixD(const TH1D& input_hist)
{
  Int_t nBins = input_hist.GetNbinsX();
  TMatrixD output_matrix(nBins,1);
  for(Int_t i=0; i<nBins; i++)
    {
      output_matrix(i, 0) = input_hist.GetBinContent(i+1);
    }
  return output_matrix;
}

// Create TH2D from TMatrixD
// reference_hist should have the desired binning, but will not
// itself be filled
TH2D TMatrixDtoTH2D(const TMatrixD& input_mat, const TH2D& reference_hist)
{
  Int_t nBinsX = reference_hist.GetNbinsX();
  Int_t nBinsY = reference_hist.GetNbinsY();
  Double_t binEdgesX[nBinsX+1];
  for(int i=0; i<nBinsX+1; i++){
      binEdgesX[i] = reference_hist.GetXaxis()->GetBinLowEdge(i+1);
  }
  Double_t binEdgesY[nBinsY+1];
  for(int i=0; i<nBinsY+1; i++){
      binEdgesY[i] = reference_hist.GetYaxis()->GetBinLowEdge(i+1);
  }
  TH2D output_hist("","",nBinsX,binEdgesX,nBinsY,binEdgesY);
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
    TH1D *Data_signal = new TH1D(mHist_data_signal->GetCVHistoWithStatError());
    // Extract covariance 
    TMatrixD data_covmat = mHist_data_signal->GetTotalErrorMatrix();
    // Pull out response matrix from input file
    TH2D* Response = (TH2D*)f->Get("response_2g1p_inclusive");
    // Pull out predicted signal MnvH1D from input file
    PlotUtils::MnvH1D *mHist_prior_true_signal = (PlotUtils::MnvH1D*)f->Get("effNum_2g1p_inclusive");  
    TH1D *Prior_true_signal = new TH1D(mHist_prior_true_signal->GetCVHistoWithStatError());

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

