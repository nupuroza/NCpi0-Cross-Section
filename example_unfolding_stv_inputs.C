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
void example_unfolding_stv_inputs() {

  // -----------------------------------------------------
  // Specify input file; clone to ouput file
  // -----------------------------------------------------

  TFile* file_in = new TFile("../stv-analysis-new/unf_input_dump.root","READ");
  file_in->Cp("output_example_unfolding_stv-inputs.root");
  file_in->Close();
  TFile* file_out = new TFile("output_example_unfolding_stv-inputs.root","UPDATE");

  // -----------------------------------------------------
  // Get inputs for unfolding sorted out
  // -----------------------------------------------------

  // Pull out response matrix, covariance, and prior_true_signal from input file
  TMatrixD* tMat_smearcept = (TMatrixD*)file_out->Get("smearcept");
  TMatrixD* tMat_data_covmat = (TMatrixD*)file_out->Get("data_covmat");
  TMatrixD* tMat_prior_true_signal = (TMatrixD*)file_out->Get("prior_true_signal");

  // Define reference histogram for TMatrixD to THnD conversions
  TH1D tHist_ref = TH1D("tHist_ref","tHist_ref",15,0,14);

  // Convert prior_true_signal to histogram, so that it can be used to calculate closure test "data" input
  TH1D tHist_prior_true_signal = TMatrixDtoTH1D(*tMat_prior_true_signal,tHist_ref);

  // Derive data_signal directly from smearcept matrix and prior_true_signal (for closure test)
  TMatrixD* tMat_data_signal = new TMatrixD(*tMat_smearcept,TMatrixD::EMatrixCreatorsOp2::kMult,*tMat_prior_true_signal);

  // Convert closure test "data" input TMatrixD into TH1D
  TH1D tHist_closure_test_data_signal = TMatrixDtoTH1D(*tMat_data_signal,tHist_ref);

  // -----------------------------------------------------
  // Configure Unfolder and run 
  // -----------------------------------------------------

  // Instantiate an object derived from the Unfolder base class
  //WienerSVDUnfolder unfolder( true, MY_REGULARIZATION );
  WienerSVDUnfolder unfolder( false, MY_REGULARIZATION );

  //DAgostiniUnfolder unfolder( NUM_ITERATIONS );
  //DAgostiniUnfolder unfolder( 1 );

  //DAgostiniUnfolder unfolder( DCC:FigureOfMerit, 0.025 );

  // Perform the unfolding
  UnfoldedMeasurement result = unfolder.unfold( *tMat_data_signal, *tMat_data_covmat,
    *tMat_smearcept, *tMat_prior_true_signal );

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
  TH1D tHist_data_signal_unfolded = TMatrixDtoTH1D(*tMat_data_signal_unfolded,tHist_ref);
 
  // Write folded and unfolded signals and prior true signal for closure test to output file
  tHist_closure_test_data_signal.SetName("closure_test_signal");
  tHist_closure_test_data_signal.Write();
  tHist_prior_true_signal.SetName("closure_test_prior_true_signal");
  tHist_prior_true_signal.Write();
  tHist_data_signal_unfolded.SetName("closure_test_signal_unfolded");
  tHist_data_signal_unfolded.Write();

}
