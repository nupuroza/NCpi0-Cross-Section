#include "TFile.h"
#include "TMatrixD.h"
#include "TH1D.h"
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

// -----------------------------------------------------
// Main method
// -----------------------------------------------------
void example_unfolding_gLEE_inputs() {

  // -----------------------------------------------------
  // Specify input file; clone to ouput file
  // -----------------------------------------------------

  TFile* file_in = new TFile("/uboone/data/users/finer/gLEE/NCPi0/2023-01-26_unfolding-investigation/2023-01-26_out_unfolded.root","READ");
  TFile* file_in2 = new TFile("/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/MajorMerge_GGE_mark/working_dir/ToTH1D/ResponseMaker/NCPi0_2g1p_Response_v2.root","READ");
  file_in->Cp("output_example_unfolding_gLEE-inputs.root");
  file_in->Close();
  TFile* file_out = new TFile("output_example_unfolding_gLEE-inputs.root","UPDATE");

  // -----------------------------------------------------
  // Get inputs for unfolding sorted out
  // -----------------------------------------------------

//*//  // Pull out response matrix and prior_true_signal from input file
//*//  TH2D* tHist_smearcept = (TH2D*)file_out->Get("response_2g1p_inclusive");
//*//  TH1D* tHist_prior_true_signal = (TH1D*)file_out->Get("effDenom_2g1p_inclusive");
//*//
//*//  // Derive data_signal directly from smearcept matrix (for closure test)
//*//  // The projection of the smearcept matrix onto the truth axis is the selection efficiency
//*//  TH1D* tHist_eff = (TH1D*) tHist_smearcept->ProjectionY();
//*//  TH1D* tHist_smearcept_projectionX = (TH1D*) tHist_smearcept->ProjectionX();
//*//  TH1D* tHist_smearcept_projectionY = (TH1D*) tHist_smearcept->ProjectionY();
//*//  // This can be used with the prior_true_signal distribution to create a closure test "data" input
//*//  TH1D* tHist_closure_test_data_signal = (TH1D*) tHist_eff->Clone();
//*//  tHist_closure_test_data_signal->Multiply(tHist_prior_true_signal);
//*//
//*//  // Convert tHnDs into TMatrices
//*//  TMatrixD tMat_smearcept = TH2DtoTMatrixD(*tHist_smearcept, false, true, true, true, true);
//*//  TMatrixD tMat_data_signal = TH1DtoTMatrixD(*tHist_closure_test_data_signal, true, true);
//*//  TMatrixD tMat_prior_true_signal = TH1DtoTMatrixD(*tHist_prior_true_signal, true, true);

//*//  // I think for the purpose of a very basic closure test it shouldn't matter what the covariance matrix is (because we want a basic enough unfolding that the covariance isn't taken into account/perfect closure is guaranteed?
//*//  // But probably for the code to not complain it needs to have the right dimensions, so I'll just pull out the covariance I use for the gLEE analysis... or maybe I could even use the same matrix as what is being used for smearcept??  
//*//  PlotUtils::MnvH1D *mHist_dummy_covariance = (PlotUtils::MnvH1D*)file_out->Get("effNum_2g1p_inclusive");
//*//  TH1D tHist_dummy_covariance = mHist_dummy_covariance->GetCVHistoWithStatError();
//*//  TMatrixD tMat_data_covmat = mHist_dummy_covariance->GetTotalErrorMatrix();

// Test with inputs directly from Mark's file
//------------------------------------------- 
  // Pull out response matrix and prior_true_signal from Mark's file
  TMatrixD* tMat_smearcept = (TMatrixD*)file_in2->Get("response_matrix");
  TH1D* tHist_prior_true_signal = (TH1D*)file_in2->Get("htrue");

  // Convert prior_true_signal into TMatrixD
  TMatrixD tMat_prior_true_signal = TH1DtoTMatrixD(*tHist_prior_true_signal, false, false);

  // Convert response (smearcept) matrix to histogram, so that projection can be made
  TH2D tHist_smearcept = TMatrixDtoTH2D(*tMat_smearcept,*tHist_prior_true_signal);
  TH1D* tHist_smearcept_projectionX = (TH1D*) tHist_smearcept.ProjectionX();
  TH1D* tHist_smearcept_projectionY = (TH1D*) tHist_smearcept.ProjectionY();

  // Derive data_signal directly from smearcept matrix and prior_true_signal (for closure test)
  TMatrixD* tMat_data_signal = new TMatrixD(*tMat_smearcept,TMatrixD::EMatrixCreatorsOp2::kMult,tMat_prior_true_signal);
  // Convert closure test "data" input TMatrixD into TH1D
  TH1D tHist_closure_test_data_signal = TMatrixDtoTH1D(*tMat_data_signal,*tHist_prior_true_signal);

  // The projection of the smearcept matrix onto the truth axis is the selection efficiency
  TH1D* tHist_eff = (TH1D*) tHist_smearcept.ProjectionX();
  //! // This can be used with the prior_true_signal distribution to create a closure test "data" input
  //! TH1D* tHist_closure_test_data_signal = (TH1D*) tHist_eff->Clone();
  //! tHist_closure_test_data_signal->Multiply(tHist_prior_true_signal);

  //! // Convert closure test "data" input tHist into TMatrix
  //! TMatrixD tMat_data_signal = TH1DtoTMatrixD(*tHist_closure_test_data_signal,false,false);

//------------------------------------------- 
// Test with inputs directly from Mark's file

  // -----------------------------------------------------
  // Configure Unfolder and run 
  // -----------------------------------------------------

  // Instantiate an object derived from the Unfolder base class
  //WienerSVDUnfolder unfolder( true, MY_REGULARIZATION );
  //WienerSVDUnfolder unfolder( false, MY_REGULARIZATION );

  //DAgostiniUnfolder unfolder( NUM_ITERATIONS );
  DAgostiniUnfolder unfolder( 1 );

  //DAgostiniUnfolder unfolder( DCC:FigureOfMerit, 0.025 );

  // Perform the unfolding
  //*//UnfoldedMeasurement result = unfolder.unfold( tMat_data_signal, tMat_data_covmat,
  //*//  *tMat_smearcept, tMat_prior_true_signal );
  UnfoldedMeasurement result = unfolder.unfold( *tMat_data_signal, *tMat_smearcept,
    *tMat_smearcept, tMat_prior_true_signal );

  // -----------------------------------------------------
  // Extract output from Unfolder; write to output file 
  // -----------------------------------------------------

  // Output is struct with these data members
  //std::unique_ptr< TMatrixD > unfolded_signal_;
  //std::unique_ptr< TMatrixD > cov_matrix_;
  //std::unique_ptr< TMatrixD > unfolding_matrix_;
  //std::unique_ptr< TMatrixD > err_prop_matrix_;
  //std::unique_ptr< TMatrixD > add_smear_matrix_;

  //Pull out unfolded signal
  auto tMat_data_signal_unfolded = result.unfolded_signal_.get();

  // Convert unfolded signal to TH1D
  TH1D tHist_data_signal_unfolded = TMatrixDtoTH1D(*tMat_data_signal_unfolded,*tHist_prior_true_signal);
 
  // Write folded and unfolded signals and prior true signal for closure test to output file
  tHist_closure_test_data_signal.SetName("closure_test_signal");
  tHist_closure_test_data_signal.Write();
  tHist_prior_true_signal->SetName("closure_test_prior_true_signal");
  tHist_prior_true_signal->Write();
  tHist_smearcept_projectionX->SetName("closure_test_smearcept_projectionX");
  tHist_smearcept_projectionX->Write();
  tHist_smearcept_projectionY->SetName("closure_test_smearcept_projectionY");
  tHist_smearcept_projectionY->Write();
  tHist_eff->SetName("closure_test_eff");
  tHist_eff->Write();
  tHist_data_signal_unfolded.SetName("closure_test_signal_unfolded");
  tHist_data_signal_unfolded.Write();

}
