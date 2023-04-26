#include "TMatrixD.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "PlotUtils/MnvH1D.h"

#include "MatrixConverters.h"

#include "../stv-analysis-new/WienerSVDUnfolder.hh"
#include "../stv-analysis-new/DAgostiniUnfolder.hh"

// -----------------------------------------------------
// Specify Unfolder specs
// -----------------------------------------------------
using DCC = DAgostiniUnfolder::ConvergenceCriterion;
using RMT = WienerSVDUnfolder::RegularizationMatrixType;

const RMT MY_REGULARIZATION = RMT::kIdentity;
//const RMT MY_REGULARIZATION = RMT::kFirstDeriv;
//const RMT MY_REGULARIZATION = RMT::kSecondDeriv;
//const int NUM_ITERATIONS = 5;

std::string unfolding_spec = "WSVD-kIdentity";
//std::string unfolding_spec = "WSVD-kFirstDerivative";
//std::string unfolding_spec = "WSVD-kSecondDerivative";
//std::string unfolding_spec = "DAgostini-5-iteration";
std::string sigDef = "2g1p_exclusive";

// -----------------------------------------------------
// Main method
// -----------------------------------------------------
void unfold(std::string filePath_in)
{

  // -----------------------------------------------------
  // Clone input file to ouput file
  // -----------------------------------------------------

  // If input file path includes ".root" suffix, strip that off
  if(filePath_in.substr(filePath_in.length()-5,filePath_in.length())==".root"){
    filePath_in = filePath_in.substr(0,filePath_in.length()-5);
  }

  TFile* file_in = new TFile((filePath_in+".root").c_str(),"READ");
  file_in->Cp((filePath_in+"_unfolded_"+unfolding_spec+".root").c_str());
  file_in->Close();
 
  // -----------------------------------------------------
  // Pull in response matrix from response matrix file
  // -----------------------------------------------------

  // Navigate to response matrix file in same directory as input hist file
  std::size_t last_slash_pos = filePath_in.find_last_of("/");
  std::string fileDir = filePath_in.substr(0,last_slash_pos);
  TFile* file_in_response = new TFile((fileDir+"/response_"+sigDef+"_v3_d22_23.root").c_str(),"READ");
 
  // Copy all contents of response matrix file to output file 
  TFile* file_out = new TFile((filePath_in+"_unfolded_"+unfolding_spec+".root").c_str(),"UPDATE");
  file_in_response->cd();
  TList* keys = gDirectory->GetListOfKeys();

  for(int i=0; i<keys->GetEntries(); i++){
    TKey* key = (TKey*)keys->At(i);
    TObject* obj = key->ReadObj();
    // This is a bit kludgy, and won't work if there is more than one TMatrix in the response matrix file
    // But if the object to copy is a TMatrix, we need to be more explicit about its name, or it won't get
    // copied over.
    if (TString(obj->ClassName()) == "TMatrixT<double>") {
      file_out->cd();
      obj->Write("response_matrix");
    }
    else{
      file_out->cd();
      obj->Write();
    }
  }
  
  // We don't need this file open any more
  file_in_response->Close();

  // Everything needed can now be pulled from file_out
  // -----------------------------------------------------

  // -----------------------------------------------------
  // Pull requisite unfolding ingredients from file
  // -----------------------------------------------------

  // Pull out measured signal MnvH1D from input file
  PlotUtils::MnvH1D *mHist_data_signal_folded = (PlotUtils::MnvH1D*)file_out->Get(("evtRate_"+sigDef).c_str()); 
  TH1D tHist_data_signal = mHist_data_signal_folded->GetCVHistoWithStatError();
  // Extract covariance 
  TMatrixD tMat_data_covmat = mHist_data_signal_folded->GetTotalErrorMatrix();
  // Pull out predicted signal MnvH1D from input file //needs to be in true space; it's a reference true space
  PlotUtils::MnvH1D *mHist_prior_true_signal = (PlotUtils::MnvH1D*)file_out->Get(("effDenom_"+sigDef).c_str());  
  TH1D tHist_prior_true_signal = mHist_prior_true_signal->GetCVHistoWithStatError();

  // Pull out response matrix 
  TMatrixD* tMat_smearcept = (TMatrixD*)file_out->Get("response_matrix");
  //TH1D* tHist_prior_true_signal = (TH1D*)file_out->Get("htrue");

  // -----------------------------------------------------
  // Handle underflow/overflow bins and make sure all
  // inputs are self-consistent 
  // -----------------------------------------------------

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

  // -----------------------------------------------------
  // Convert unfolding ingredients into matrices
  // -----------------------------------------------------

//  // tMat_data_covmat will always include the underflow/overflow from mHist_data_signal_folded
//  // by construction, so we always want to include the underflow/overflow for the other inputs
//  // to the unfolding. Hence all of the `True` arguments below
//  
//  // Convert inputs into TMatrixD
//  TMatrixD tMat_data_signal = TH1DtoTMatrixD(tHist_data_signal, true, true);
//  TMatrixD tMat_prior_true_signal = TH1DtoTMatrixD(tHist_prior_true_signal, true, true);
//  TMatrixD tMat_smearcept = TH2DtoTMatrixD(*tHist2D_response, false, true, true, true, true); // Note the first bool here isn't related to underflow/overflow
//  //TMatrixD tMat_smearcept = TH2DtoTMatrixD(*tHist2D_response, true, true, true, true, true);

  // Convert inputs into TMatrixD
  TMatrixD tMat_data_signal = TH1DtoTMatrixD(tHist_data_signal, include_underflow_reco, include_overflow_reco);
  TMatrixD tMat_prior_true_signal = TH1DtoTMatrixD(tHist_prior_true_signal, include_underflow_true, include_overflow_true);
  //TMatrixD tMat_smearcept = TH2DtoTMatrixD(*tHist2D_response, false, include_underflow_reco, include_overflow_reco, include_underflow_true, include_overflow_true);
  //TMatrixD tMat_smearcept = TH2DtoTMatrixD(*tHist2D_response, true, include_underflow_reco, include_overflow_reco, include_underflow_true, include_overflow_true);

  // If either include_underflow_reco or include_overflow_reco is false, replace tMat_data_covmat with tMat_data_covmat->GetSub(x1,x2,y1,y2);
  // Usage: TMatrixT< Element > GetSub (Int_t row_lwb, Int_t row_upb, Int_t col_lwb, Int_t col_upb, Option_t *option="S") const 
  Int_t lowerBound = include_underflow_reco ? 0 : 1;
  Int_t upperBound = include_overflow_reco ? nBins_reco+1 : nBins_reco;
  TMatrixD tMat_data_covmat_final = tMat_data_covmat.GetSub(lowerBound, upperBound, lowerBound, upperBound);   
  TMatrixD tMat_smearcept_final = tMat_smearcept->GetSub(lowerBound, upperBound, lowerBound, upperBound);   

  // Print out number of rows and columns in each input matrix

  std::cout << "binVal_underflow_reco: " << binVal_underflow_reco << std::endl;
  std::cout << "include_underflow_reco: " << include_underflow_reco << std::endl;
  std::cout << "binVal_overflow_reco: " << binVal_overflow_reco << std::endl;
  std::cout << "include_overflow_reco: " << include_overflow_reco << std::endl;
  std::cout << "binVal_underflow_true: " << binVal_underflow_true << std::endl;
  std::cout << "include_underflow_true: " << include_underflow_true << std::endl;
  std::cout << "binVal_overflow_true: " << binVal_overflow_true << std::endl;
  std::cout << "include_overflow_true: " << include_overflow_true << std::endl;

  Int_t data_rows = tMat_data_signal.GetNrows();
  Int_t data_cols = tMat_data_signal.GetNcols();
  std::cout << "data_rows: " << data_rows << "\t" << "data_cols: " << data_cols << std::endl;

  Int_t cov_final_rows = tMat_data_covmat_final.GetNrows();
  Int_t cov_final_cols = tMat_data_covmat_final.GetNcols();
  std::cout << "cov_final_rows: " << cov_final_rows << "\t" << "cov_final_cols: " << cov_final_cols << std::endl;

  Int_t smearcept_rows = tMat_smearcept_final.GetNrows();
  Int_t smearcept_cols = tMat_smearcept_final.GetNcols();
  std::cout << "smearcept_rows: " << smearcept_rows << "\t" << "smearcept_cols: " << smearcept_cols << std::endl;
  std::cout << "tMat_smearcept_final.Print(): " << std::endl;
  //tMat_smearcept_final.Print();

  Int_t mc_rows = tMat_prior_true_signal.GetNrows();
  Int_t mc_cols = tMat_prior_true_signal.GetNcols();
  std::cout << "mc_rows: " << mc_rows << "\t" << "mc_cols: " << mc_cols << std::endl;

  // -----------------------------------------------------
  // Configure Unfolder and run 
  // -----------------------------------------------------

  // Instantiate an object derived from the Unfolder base class
  WienerSVDUnfolder unfolder( true, MY_REGULARIZATION );
  //WienerSVDUnfolder unfolder( false, MY_REGULARIZATION );

  //DAgostiniUnfolder unfolder( NUM_ITERATIONS );

  //DAgostiniUnfolder unfolder( DCC:FigureOfMerit, 0.025 );

  // Perform the unfolding
  UnfoldedMeasurement result = unfolder.unfold( tMat_data_signal, tMat_data_covmat_final,
    tMat_smearcept_final, tMat_prior_true_signal );

  // -----------------------------------------------------
  // Extract output from Unfolder; write to output file 
  // -----------------------------------------------------

  // Output is struct with these data members
  //std::unique_ptr< TMatrixD > unfolded_signal_;
  //std::unique_ptr< TMatrixD > cov_matrix_;
  //std::unique_ptr< TMatrixD > unfolding_matrix_;
  //std::unique_ptr< TMatrixD > err_prop_matrix_;
  //std::unique_ptr< TMatrixD > add_smear_matrix_;

  //Pull out unfolded signal and covariance from results
  auto tMat_data_signal_unfolded = result.unfolded_signal_.get();
  auto tMat_unfolded_covariance = result.cov_matrix_.get();
  auto tMat_errprop_matrix = result.err_prop_matrix_.get();

  // -----------------------------------------------------
  // Convert unfolder outputs back to hists 
  // -----------------------------------------------------

  // Convert outputs into TH1D/TH2D  
  TH1D tHist_data_signal_unfolded = TMatrixDtoTH1D(*tMat_data_signal_unfolded, tHist_prior_true_signal);
  TH2D tHist2D_unfolded_covariance = TMatrixDtoTH2D(*tMat_unfolded_covariance, tHist_prior_true_signal);
  TH2D tHist2D_covariance = TMatrixDtoTH2D(tMat_data_covmat_final, tHist_data_signal);

  // -----------------------------------------------------
  // 

  // OOF, I think this is a whole piece that still needs to get sorted out... Maybe it's not so bad?
  // Definitely we're not set up yet to do the unfolding in every systematic universe simultaneously,
  // but the covariance that get spit out of the unfolder is correct and the diagonal errors can be given
  // back to the CV hist...
 
  //for(Int_t i=1; i<nBins_true+1; i++)
  //{    
  //double a = tHist_data_signal_unfolded.GetBinError(i);
  //double b = tHist2D_unfolded_covariance.GetBinError(i,i);
  //tHist_data_signal_unfolded.SetBinError(i,b);
  //double c = tHist_data_signal_unfolded.GetBinError(i);
  //std::cout << "tHist_before: " << a << "\t" << "tHist_after: " << c << "\t" << "cov: " << b << std::endl;  
  //}

  // 
  // -----------------------------------------------------

  // Construct MnvH1D
  PlotUtils::MnvH1D mHist_data_signal_unfolded = PlotUtils::MnvH1D(tHist_data_signal_unfolded);
  mHist_data_signal_unfolded.AddMissingErrorBandsAndFillWithCV(*mHist_data_signal_folded);
    
  // -----------------------------------------------------
  // Write unfolded results to output file
  // -----------------------------------------------------

  mHist_data_signal_unfolded.SetName(("unfolded_evtRate_"+sigDef).c_str());
  mHist_data_signal_unfolded.Write();   
  tHist2D_unfolded_covariance.SetName(("unfolded_cov_evtRate_"+sigDef).c_str());
  tHist2D_unfolded_covariance.Write();
  tHist2D_covariance.SetName(("cov_evtRate_"+sigDef).c_str());
  tHist2D_covariance.Write();

  //file_out->Close();
  return;
}

