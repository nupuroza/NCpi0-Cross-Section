#include "TMatrixD.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "TRandom.h"
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

    std::string sigDefnp = sigDef.substr(0, sigDef.length()-10);
    // Pull out measured signal MnvH1D from input file (evtRate is background-subtracted event rate)
    PlotUtils::MnvH1D *mHist_data_signal_folded = (PlotUtils::MnvH1D*)file_out->Get(("evtRate_"+sigDef).c_str()); 
    TH1D tHist_data_signal = mHist_data_signal_folded->GetCVHistoWithStatError();
    // Extract covariance, including MC statistical uncertainty.
    PlotUtils::MnvH1D *mHist_fakedata_mc = (PlotUtils::MnvH1D*)file_out->Get(("fakedata_mc_"+sigDef).c_str());
    TMatrixD tMat_data_covmat = mHist_fakedata_mc->GetTotalErrorMatrix(true);
    // Pull out predicted signal MnvH1D from input file // needs to be in true space, happens to also be the efficiency
    // denominator in our xsec extraction
    PlotUtils::MnvH1D *mHist_prior_true_signal = (PlotUtils::MnvH1D*)file_out->Get(("effDenom_"+sigDef).c_str());  
    TH1D tHist_prior_true_signal = mHist_prior_true_signal->GetCVHistoWithStatError(); 
    // Pull out NuWro truth for smearing.
    TH1D *tHist_nuwro_true_signal = (TH1D*) file_out -> Get(("nu_uBooNE_denom_" + sigDefnp).c_str());
    // Pull out response matrix; adopt jargon of Unfolder tool, "smearcept" == "smearing + acceptance"
    TMatrixD* tMat_smearcept = (TMatrixD*)file_out->Get(("response_matrix_"+sigDef).c_str());
    // Pull out selected data for data statistical uncertainties.
    PlotUtils::MnvH1D *mHist_data_selected = (PlotUtils::MnvH1D*) file_out -> Get(("data_selected_" + sigDefnp).c_str());

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
    Int_t lowerBoundTrue = include_underflow_true ? 0 : 1; 
    Int_t upperBoundTrue = include_overflow_true ? nBins_true+1 : nBins_true;

    TMatrixD tMat_data_covmat_final = tMat_data_covmat.GetSub(lowerBound, upperBound, lowerBound, upperBound);   
    TMatrixD tMat_smearcept_final = tMat_smearcept->GetSub(lowerBound, upperBound, lowerBoundTrue, upperBoundTrue);   

    // --------------------------------------------------------
    // Then convert other unfolding ingredients into matrices
    // --------------------------------------------------------

    TMatrixD tMat_prior_true_signal = TH1DtoTMatrixD(tHist_prior_true_signal, include_underflow_true, include_overflow_true);
    TMatrixD tMat_data_signal = TH1DtoTMatrixD(tHist_data_signal, include_underflow_reco, include_overflow_reco);
    TMatrixD tMat_nuwro_true_signal = TH1DtoTMatrixD(*tHist_nuwro_true_signal, include_underflow_true, include_overflow_true);

    // --------------------------------------------------------
    // Derive reco-space generator prediction -- needed for
    // closure test but also generally useful to have later
    // --------------------------------------------------------

    // Derive generator prediction in reco space directly from smearcept matrix and prior_true_signal
    TMatrixD tMat_folded_true_signal = TMatrixD(tMat_smearcept_final,TMatrixD::EMatrixCreatorsOp2::kMult,tMat_prior_true_signal); 
    TH1D tHist_folded_true_signal = TMatrixDtoTH1D(tMat_folded_true_signal,tHist_data_signal, include_underflow_true);
    PlotUtils::MnvH1D *mHist_prior_reco_signal = (PlotUtils::MnvH1D*)file_out->Get(("effNum_reco_"+sigDef).c_str());  
    TH1D tHist_prior_reco_signal = mHist_prior_reco_signal->GetCVHistoWithStatError();
    // This one always closes. Use in emergency to show closure.
    // TH1D tHist_prior_reco_signal = tHist_folded_true_signal;
    TMatrixD tMat_prior_reco_signal_initial = TMatrixD(nBins_reco + 2, 1, tHist_prior_reco_signal.GetArray());
    TMatrixD tMat_prior_reco_signal = tMat_prior_reco_signal_initial.GetSub(lowerBound, upperBound, 0, 0);

    // --------------------------------------------------------
    // Extract covariance matrix for each vertical error band
    // to be unfolded and used for plotting error breakdown.
    // --------------------------------------------------------
    std::vector<std::string> error_bands = mHist_fakedata_mc -> GetVertErrorBandNames();
    std::vector<TMatrixD> folded_covariance_matrices;
    TMatrixD tMat_folded_data_stat_covariance_matrix = TMatrixD(nBins_reco + 2, nBins_reco + 2);
    TMatrixD tMat_folded_mc_stat_covariance_matrix = TMatrixD(nBins_reco + 2, nBins_reco + 2);
    for(std::string error_band : error_bands)
        folded_covariance_matrices.push_back(mHist_fakedata_mc -> GetVertErrorBand(error_band) -> CalcCovMx());

    // If this is a closure test, data_signal should be generator prediction in reco space
    // Extract data statistical uncertainty using Neyman-Pearson method (see calculateChi2.py).
    if(closureTest){
        tMat_data_signal = tMat_prior_reco_signal;
        for(int i = 0; i < upperBound + include_underflow_reco; i++){
            double data_stat_error = mHist_fakedata_mc -> GetBinContent(i + lowerBound);
            double mc_stat_error = mHist_fakedata_mc -> GetBinError(i + lowerBound);
            tMat_data_covmat_final[i][i] += data_stat_error;
            tMat_folded_data_stat_covariance_matrix[i+lowerBound][i+lowerBound] = data_stat_error;
            tMat_folded_mc_stat_covariance_matrix[i+lowerBound][i+lowerBound] = mc_stat_error;
    }}
    else{
        for(int i = 0; i < upperBound + include_underflow_reco; i++){
            double data_stat_error = 3/(1/mHist_data_selected -> GetBinContent(i + lowerBound) + 2/mHist_fakedata_mc -> GetBinContent(i + lowerBound));
            double mc_stat_error = mHist_fakedata_mc -> GetBinError(i + lowerBound);
            tMat_data_covmat_final[i][i] += data_stat_error;
            tMat_folded_data_stat_covariance_matrix[i+lowerBound][i+lowerBound] = data_stat_error;
            tMat_folded_mc_stat_covariance_matrix[i+lowerBound][i+lowerBound] = mc_stat_error;
    }}
    // Add data and MC statistical uncertainties to vector of folded covariance matrices.
    error_bands.push_back("data_statistical");
    folded_covariance_matrices.push_back(tMat_folded_data_stat_covariance_matrix);
    error_bands.push_back("mc_statistical");
    folded_covariance_matrices.push_back(tMat_folded_mc_stat_covariance_matrix);

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

    std::cout << "data signal" << std::endl;
    tMat_data_signal.Print();
    std::cout << "covmat final" << std::endl;
    tMat_data_covmat_final.Print();
    std::cout << "smearcept final" << std::endl;
    tMat_smearcept_final.Print();
    std::cout << "prior true signal" << std::endl;
    tMat_prior_true_signal.Print();

    // Perform the unfolding
    UnfoldedMeasurement result = unfolder->unfold( tMat_data_signal, tMat_data_covmat_final,
            tMat_smearcept_final, tMat_prior_true_signal );

    std::cout << "Finished unfolding" << std::endl;

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
    auto tMat_unfolding_matrix = result.unfolding_matrix_.get();
    auto tMat_errprop_matrix = result.err_prop_matrix_.get();
    auto tMat_add_smear_matrix = result.add_smear_matrix_.get();

    // Refold the unfolded data
    TMatrixD refolded_true_signal(tMat_smearcept_final, TMatrixD::EMatrixCreatorsOp2::kMult, *tMat_data_signal_unfolded);
    
    // line moved from deleted piece due to relevance to overall code
    TH1D tHist_data_signal_unfolded = TMatrixDtoTH1D(*tMat_data_signal_unfolded, tHist_prior_true_signal, include_underflow_true);

    // -----------------------------------------------------
    // Calculate smeared true signal distribution
    // -----------------------------------------------------

    // Note: TMatrix axes are inverted from TH2D, so for tMat_add_smear_matrix, rows correspond to reco, and columns correspond to true (see unfold.C)
    TMatrixD* tMat_smeared_true_signal = new TMatrixD(*tMat_add_smear_matrix,TMatrixD::EMatrixCreatorsOp2::kMult,tMat_prior_true_signal);
    TMatrixD* tMat_smeared_nuwro_signal = new TMatrixD(*tMat_add_smear_matrix, TMatrixD::EMatrixCreatorsOp2::kMult, tMat_nuwro_true_signal);

    // -----------------------------------------------------
    // Convert unfolder outputs back to hists 
    // -----------------------------------------------------

    // Convert outputs into TH1D/TH2D   
    TH2D tHist2D_unfolded_covariance = TMatrixDtoTH2D(*tMat_unfolded_covariance, tHist_prior_true_signal, include_underflow_true);
    TH2D tHist2D_covariance = TMatrixDtoTH2D(tMat_data_covmat_final, tHist_data_signal, include_underflow_true);
    TH2D tHist2D_errprop_matrix = TMatrixDtoTH2D(*tMat_errprop_matrix, tHist_prior_true_signal, include_underflow_true);
    TH2D tHist2D_add_smear_matrix = TMatrixDtoTH2D(*tMat_add_smear_matrix, tHist_prior_true_signal, include_underflow_true);
    TH1D tHist_smeared_true_signal = TMatrixDtoTH1D(*tMat_smeared_true_signal, tHist_prior_true_signal, include_underflow_true);
    TH1D tHist_smeared_nuwro_signal = TMatrixDtoTH1D(*tMat_smeared_nuwro_signal, tHist_prior_true_signal, include_underflow_true);
    //TH1D tHist_bias = TMatrixDtoTH1D(*bias, tHist_prior_true_signal, include_underflow_true);

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
    TMatrixD *tMat_unfolding_matrix_final = new TMatrixD(nBins_reco + 2, nBins_true + 2);
    for(int i = 0; i < tMat_unfolding_matrix -> GetNrows(); i++)
        for(int j = 0; j< tMat_unfolding_matrix -> GetNcols(); j++)
            tMat_unfolding_matrix_final -> operator()(i + lowerBound, j + lowerBoundTrue) = tMat_unfolding_matrix -> operator()(i, j);
    int matrix_number = 0;
    for(const auto& error_band : error_bands){
        TMatrixD folded_covariance_matrix = folded_covariance_matrices.at(matrix_number);
        TMatrixD unfolded_covariance_matrix_intermediate = TMatrixD(*tMat_unfolding_matrix_final, TMatrixD::EMatrixCreatorsOp2::kMult, folded_covariance_matrix);
        TMatrixD unfolded_covariance_matrix = TMatrixD(unfolded_covariance_matrix_intermediate, TMatrixD::EMatrixCreatorsOp2::kMultTranspose, *tMat_unfolding_matrix_final);
        mHist_data_signal_unfolded.FillSysErrorMatrix(error_band.c_str(), unfolded_covariance_matrix);
        matrix_number++;
    }

    // -----------------------------------------------------
    // Write unfolded results to output file
    // -----------------------------------------------------

    mHist_data_signal_unfolded.SetName(("unfolded_evtRate_"+sigDef).c_str());
    mHist_data_signal_unfolded.Write();   
    tHist2D_unfolded_covariance.SetName(("unfolded_cov_evtRate_"+sigDef).c_str());
    tHist2D_unfolded_covariance.Write();
    tHist2D_covariance.SetName(("cov_evtRate_"+sigDef).c_str());
    tHist2D_covariance.Write();
    tHist2D_errprop_matrix.SetName(("errprop_matrix_"+sigDef).c_str());
    tHist2D_errprop_matrix.Write();
    tHist2D_add_smear_matrix.SetName(("add_smear_matrix_"+sigDef).c_str());
    tHist2D_add_smear_matrix.Write();
    tHist_prior_true_signal.SetName(("prior_true_signal_"+sigDef).c_str());
    tHist_prior_true_signal.Write();
    tHist_prior_reco_signal.SetName(("prior_reco_signal_"+sigDef).c_str());
    tHist_prior_reco_signal.Write();
    tHist_folded_true_signal.SetName(("folded_true_signal_" + sigDef).c_str());
    tHist_folded_true_signal.Write();
    tHist_smeared_true_signal.SetName(("smeared_true_signal_"+sigDef).c_str());
    tHist_smeared_true_signal.Write();
    tHist_smeared_nuwro_signal.SetName(("smeared_nuwro_signal_" + sigDef).c_str());
    tHist_smeared_nuwro_signal.Write();
    //tHist_bias.SetName(("bias_"+sigDef).c_str());
    //tHist_bias.Write();

    //file_out->Close();
    return;
}

// -----------------------------------------------------
// Main method
// -----------------------------------------------------
void unfold(std::string filePath_in, bool useWienerSVD, std::string unfoldingConfig, bool closureTest)
{

    gROOT -> SetBatch(kTRUE);
    gStyle -> SetOptStat(0);
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
    std::size_t second_last_slash_pos = filePath_in.find_last_of("/", last_slash_pos - 1); 
    std::string fileDir = filePath_in.substr(0, second_last_slash_pos); 

    TFile* file_in_response = new TFile((fileDir+"/response_matrices/response_matrices_exclusive.root").c_str(),"READ");


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
            else if(objName == "response_matrix_2gXp_inclusive") {
                file_out->cd();
                obj->Write("response_matrix_2gXp_inclusive");
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

    execute_unfolding(file_out,"2g1p_exclusive",useWienerSVD,closureTest,MY_REGULARIZATION,NUM_ITERATIONS);
    execute_unfolding(file_out,"2g0p_exclusive",useWienerSVD,closureTest,MY_REGULARIZATION,NUM_ITERATIONS);
    execute_unfolding(file_out,"2gXp_inclusive",useWienerSVD,closureTest,MY_REGULARIZATION,NUM_ITERATIONS);


}
