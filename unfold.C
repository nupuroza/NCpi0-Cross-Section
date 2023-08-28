#include "TMatrixD.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "PlotUtils/MnvH1D.h"

#include "MatrixConverters.h"

#include "../stv-analysis-new/WienerSVDUnfolder.hh"
#include "../stv-analysis-new/DAgostiniUnfolder.hh"


// -----------------------------------------------------
// Bit that should be evaluated separately once for each
// signal definition (sigDef) that we care about
// -----------------------------------------------------
void execute_unfolding(TFile* file_out, std::string sigDef, bool useWienerSVD, bool closureTest, WienerSVDUnfolder::RegularizationMatrixType MY_REGULARIZATION, int NUM_ITERATIONS)
{

  std::cout << "I'm inside execute_unfolding(). I'm running over " << sigDef << std::endl;

  // -----------------------------------------------------
  // Pull requisite unfolding ingredients from file
  // -----------------------------------------------------

  // Pull out measured signal MnvH1D from input file (evtRate is background-subtracted event rate)
  PlotUtils::MnvH1D *mHist_data_signal_folded = (PlotUtils::MnvH1D*)file_out->Get(("evtRate_"+sigDef).c_str()); 
  TH1D tHist_data_signal = mHist_data_signal_folded->GetCVHistoWithStatError();
  // Extract covariance 
  TMatrixD tMat_data_covmat = mHist_data_signal_folded->GetTotalErrorMatrix();
  // Pull out predicted signal MnvH1D from input file // needs to be in true space, happens to also be the efficiency
  // denominator in our xsec extraction
  PlotUtils::MnvH1D *mHist_prior_true_signal = (PlotUtils::MnvH1D*)file_out->Get(("effDenom_"+sigDef).c_str());  
  TH1D tHist_prior_true_signal = mHist_prior_true_signal->GetCVHistoWithStatError();

  // Pull out response matrix; adopt jargon of Unfolder tool, "smearcept" == "smearing + acceptance"
  TMatrixD* tMat_smearcept = (TMatrixD*)file_out->Get(("response_matrix_"+sigDef).c_str());

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

  // Check if underflow is zero for prior_true_signal (true space)
  double binVal_underflow_true = tHist_prior_true_signal.GetBinContent(0);
  bool include_underflow_true = binVal_underflow_true > threshold ? 1 : 0;

  // Check if overflow is zero for prior_true_signal (true space)
  Int_t nBins_true = tHist_prior_true_signal.GetNbinsX();
  double binVal_overflow_true = tHist_prior_true_signal.GetBinContent(nBins_true+1);
  bool include_overflow_true = binVal_overflow_true > threshold ? 1 : 0;

  //-------------------------------------------------------------
  // Deal with reducing the matrix scope first so that if needed
  // the correct reduced smearcept matrix is used in closure test
  //-------------------------------------------------------------

  // If either include_underflow_reco or include_overflow_reco is false, replace tMat_data_covmat with tMat_data_covmat->GetSub(x1,x2,y1,y2);
  // Usage: TMatrixT< Element > GetSub (Int_t row_lwb, Int_t row_upb, Int_t col_lwb, Int_t col_upb, Option_t *option="S") const 
  Int_t lowerBound = include_underflow_reco ? 0 : 1;
  Int_t upperBound = include_overflow_reco ? nBins_reco+1 : nBins_reco;
  TMatrixD tMat_data_covmat_final = tMat_data_covmat.GetSub(lowerBound, upperBound, lowerBound, upperBound);   
  TMatrixD tMat_smearcept_final = tMat_smearcept->GetSub(lowerBound, upperBound, lowerBound, upperBound);   

  // --------------------------------------------------------
  // Then convert other unfolding ingredients into matrices
  // --------------------------------------------------------

  TMatrixD tMat_prior_true_signal = TH1DtoTMatrixD(tHist_prior_true_signal, include_underflow_true, include_overflow_true);
  TMatrixD tMat_data_signal = TH1DtoTMatrixD(tHist_data_signal, include_underflow_reco, include_overflow_reco);

  // --------------------------------------------------------
  // Derive reco-space generator prediction -- needed for
  // closure test but also generally useful to have later
  // --------------------------------------------------------

  // Derive generator prediction in reco space directly from smearcept matrix and prior_true_signal
  TMatrixD tMat_prior_reco_signal = TMatrixD(tMat_smearcept_final,TMatrixD::EMatrixCreatorsOp2::kMult,tMat_prior_true_signal);
  TH1D tHist_prior_reco_signal = TMatrixDtoTH1D(tMat_prior_reco_signal,tHist_data_signal);

  // If this is a closure test, data_signal should be generator prediction in reco space
  if(closureTest){
    //tMat_data_signal = TMatrixD(tMat_smearcept_final,TMatrixD::EMatrixCreatorsOp2::kMult,tMat_prior_true_signal);
    tMat_data_signal = tMat_prior_reco_signal;
  }

  // ----------------------------------------------------------
  // Print out number of rows and columns in each input matrix
  // For debugging purposes; generally useful to keep around
  // ----------------------------------------------------------

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
  std::cout << "folded_signal_rows (data): " << data_rows << "\t" << "folded_signal_cols (data): " << data_cols << std::endl;

  Int_t cov_final_rows = tMat_data_covmat_final.GetNrows();
  Int_t cov_final_cols = tMat_data_covmat_final.GetNcols();
  std::cout << "cov_final_rows (data): " << cov_final_rows << "\t" << "cov_final_cols (data): " << cov_final_cols << std::endl;

  Int_t mc_rows = tMat_prior_true_signal.GetNrows();
  Int_t mc_cols = tMat_prior_true_signal.GetNcols();
  std::cout << "prior_true_signal_rows (mc): " << mc_rows << "\t" << "prior_true_signal_cols (mc): " << mc_cols << std::endl;

  Int_t smearcept_rows = tMat_smearcept_final.GetNrows();
  Int_t smearcept_cols = tMat_smearcept_final.GetNcols();
  std::cout << "smearcept_rows (data:mc): " << smearcept_rows << "\t" << "smearcept_cols (data:mc): " << smearcept_cols << std::endl;

  // -----------------------------------------------------
  // Configure Unfolder and run 
  // -----------------------------------------------------

  // Instantiate an object derived from the Unfolder base class
  Unfolder* unfolder; 
  if(useWienerSVD){
    unfolder = new WienerSVDUnfolder( true, MY_REGULARIZATION );
  } 
  else{
    unfolder = new DAgostiniUnfolder( NUM_ITERATIONS );
  }

  // Perform the unfolding
  UnfoldedMeasurement result = unfolder->unfold( tMat_data_signal, tMat_data_covmat_final,
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
  auto tMat_add_smear_matrix = result.add_smear_matrix_.get();

  // -----------------------------------------------------
  // Calculate smeared true signal distribution
  // -----------------------------------------------------
  
  // Note: TMatrix axes are inverted from TH2D, so for tMat_add_smear_matrix, rows correspond to reco, and columns correspond to true (see unfold.C)
  TMatrixD* tMat_smeared_true_signal = new TMatrixD(*tMat_add_smear_matrix,TMatrixD::EMatrixCreatorsOp2::kMult,tMat_prior_true_signal);

  // -----------------------------------------------------
  // Convert unfolder outputs back to hists 
  // -----------------------------------------------------

  // Convert outputs into TH1D/TH2D  
  TH1D tHist_data_signal_unfolded = TMatrixDtoTH1D(*tMat_data_signal_unfolded, tHist_prior_true_signal);
  TH2D tHist2D_unfolded_covariance = TMatrixDtoTH2D(*tMat_unfolded_covariance, tHist_prior_true_signal);
  TH2D tHist2D_covariance = TMatrixDtoTH2D(tMat_data_covmat_final, tHist_data_signal);
  TH2D tHist2D_add_smear_matrix = TMatrixDtoTH2D(*tMat_add_smear_matrix, tHist_prior_true_signal);
  TH1D tHist_smeared_true_signal = TMatrixDtoTH1D(*tMat_smeared_true_signal, tHist_prior_true_signal);

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

  tHist2D_add_smear_matrix.SetName(("add_smear_matrix_"+sigDef).c_str());
  tHist2D_add_smear_matrix.Write();
  tHist_prior_true_signal.SetName(("prior_true_signal_"+sigDef).c_str());
  tHist_prior_true_signal.Write();
  tHist_prior_reco_signal.SetName(("prior_reco_signal_"+sigDef).c_str());
  tHist_prior_reco_signal.Write();
  tHist_smeared_true_signal.SetName(("smeared_true_signal_"+sigDef).c_str());
  tHist_smeared_true_signal.Write();

  //file_out->Close();
  return;
}

// -----------------------------------------------------
// Main method
// -----------------------------------------------------
void unfold(std::string filePath_in, bool useWienerSVD, std::string unfoldingConfig, bool closureTest)
{

  // -----------------------------------------------------
  // Specify Unfolder specs
  // -----------------------------------------------------

  std::string unfolding_spec = "";
  using RMT = WienerSVDUnfolder::RegularizationMatrixType;
  RMT MY_REGULARIZATION; 
  int NUM_ITERATIONS = 0;
 
  if(useWienerSVD){
    if(unfoldingConfig=="kIdentity"){
      MY_REGULARIZATION = RMT::kIdentity;
    }
    else if(unfoldingConfig=="kFirstDerivative"){    
      MY_REGULARIZATION = RMT::kFirstDeriv;
    }
    else if(unfoldingConfig== "kSecondDerivative"){
      MY_REGULARIZATION = RMT::kSecondDeriv;
    }
    unfolding_spec = ("WSVD-"+unfoldingConfig).c_str();
  }
  // If we're not using Wiener-SVD, then we're using D'Agostini
  else{
    NUM_ITERATIONS = std::stoi(unfoldingConfig);
    unfolding_spec = ("DAgostini-"+unfoldingConfig+"-iteration").c_str();
  } 
 
  // -----------------------------------------------------
  // Clone input file to ouput file
  // -----------------------------------------------------

  // If input file path includes ".root" suffix, strip that off
  if(filePath_in.substr(filePath_in.length()-5,filePath_in.length())==".root"){
    filePath_in = filePath_in.substr(0,filePath_in.length()-5);
  }

  TFile* file_in = new TFile((filePath_in+".root").c_str(),"READ");
  // A little hacky, but the best I figured out for how to handle the
  // if/else and also the fact that ().c_str() returns a const char*
  // -----------------------------------------------------------------
  std::string output_file_path_base = filePath_in+"_unfolded_"+unfolding_spec;
  std::string output_file_path_suffix = closureTest ? "_closureTest.root" : ".root";
  std::string output_file_path = output_file_path_base + output_file_path_suffix;
  std::cout << "DEBUG\toutput_file_path: " << output_file_path << std::endl;

  file_in->Cp(output_file_path.c_str()); // this is the point at which the output file is created
  file_in->Close();
 
  // -----------------------------------------------------
  // Pull in response matrix from response matrix file
  // and copy contents to output file
  // -----------------------------------------------------

  // Navigate to response matrix file in same directory as input hist file
  std::size_t last_slash_pos = filePath_in.find_last_of("/");
  std::string fileDir = filePath_in.substr(0,last_slash_pos);
  TFile* file_in_response = new TFile((fileDir+"/response_matrices_exclusive.root").c_str(),"READ");
 
  // Copy all contents of response matrix file to output file 
  TFile* file_out = new TFile(output_file_path.c_str(),"UPDATE");
  file_in_response->cd();
  TList* keys = gDirectory->GetListOfKeys();

  // Go through the contents of the file one-by-one
  for(int i=0; i<keys->GetEntries(); i++){
    TKey* key = (TKey*)keys->At(i);
    TObject* obj = key->ReadObj();
    // TMatrixD objects are a little trickier, so we find them by their names
    // and then explicitly propagate the object names to the output file
    if (TString(obj->ClassName()) == "TMatrixT<double>") {
      TString objName = key->GetName();
      if (objName == "response_matrix_2g1p_exclusive") {
        file_out->cd();
        obj->Write("response_matrix_2g1p_exclusive");
      }
      else if(objName == "response_matrix_2g0p_exclusive") {
        file_out->cd();
        obj->Write("response_matrix_2g0p_exclusive");
      }
    }
    // The rest of the objects are much more straightforward to deal with
    else{
      file_out->cd();
      obj->Write();
    }
  }
  
  // We don't need this file open any more
  file_in_response->Close();

  // Everything needed can now be pulled from file_out
  // -----------------------------------------------------

  execute_unfolding(file_out,"2g1p_exclusive",useWienerSVD,closureTest,MY_REGULARIZATION,NUM_ITERATIONS);
  execute_unfolding(file_out,"2g0p_exclusive",useWienerSVD,closureTest,MY_REGULARIZATION,NUM_ITERATIONS);

  return;
}

