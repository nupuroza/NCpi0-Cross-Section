#include "TMatrixD.h"
#include "../stv-analysis-new/WienerSVDUnfolder.hh"
#include "../stv-analysis-new/DAgostiniUnfolder.hh"

using namespace std;

void unfold_debug()
{
    double array_debug_data [10] = {11.12,90.42,127.7,85.99,78.39,60.29,46.88,21.98,13.39,11.33};
    TMatrixD tMat_debug_data(10,1,array_debug_data);

    double array_debug_cov [100] = {10.89,65.82,96.57,67.32,70.91,47.64,36.73,16.51,7.617,8.618, 
                                    65.82,604.1,880.6,648.1,652.3,465.3,381.6,166.4,73.07,86.44, 
                                    96.57,880.6,1383,1025,1015,707.1,573.9,242.4,109.7,127.8, 
                                    67.32,648.1,1025,812.1,774.9,546,439.1,181.4,81.91,97.83, 
                                    70.91,652.3,1015,774.9,782.2,535.4,430.3,182.3,80.3,97.14, 
                                    47.64,465.3,707.1,546,535.4,406.7,325.3,140.5,62.21,75.18, 
                                    36.73,381.6,573.9,439.1,430.3,325.3,285.5,119.8,54,62.83, 
                                    16.51,166.4,242.4,181.4,182.3,140.5,119.8,57.36,23.74,28.48, 
                                    7.617,73.07,109.7,81.91,80.3,62.21,54,23.74,14.17,12.78, 
                                    8.618,86.44,127.8,97.83,97.14,75.18,62.83,28.48,12.78,16.97};
    TMatrixD tMat_debug_cov(10,10,array_debug_cov);

    double array_debug_response [100] = {0.009166,0.001807,0.0001383,0.0001034,8.559e-05,0,0,0,0.00022,0,
                                         0.005538,0.01284,0.004886,0.001413,0.0004279,0.0009621,0.0002698,0.0002961,0.00022,0.0001544,
                                         0.0003819,0.003586,0.01217,0.005619,0.002824,0.001641,0.001169,0.0002961,0.0008799,0,
                                         0.0003819,0.0005003,0.003572,0.01155,0.005435,0.003396,0.001709,0.001333,0,0.0003088,
                                         0.0001909,0.0001668,0.0003687,0.004171,0.01138,0.005716,0.003058,0.0008884,0.0011,0.0003088,
                                         0,5.559e-05,4.609e-05,0.0006895,0.004322,0.01019,0.005845,0.002369,0.0008799,0.000772,
                                         0,2.78e-05,4.609e-05,0.0002758,0.001027,0.004301,0.007464,0.00533,0.0033,0.0006176,
                                         0,5.559e-05,6.914e-05,0.0001379,0.0001284,0.001132,0.005665,0.007403,0.004839,0.001853,
                                         0,0,2.305e-05,3.447e-05,0.0001712,0.0003962,0.001799,0.00607,0.005939,0.001853,
                                         0,5.559e-05,4.609e-05,0,0,0.000283,8.993e-05,0.001925,0.005059,0.0088};
    TMatrixD tMat_debug_response(10,10,array_debug_response);

    double array_debug_prior_true_signal [10] = {725.7,4986.0,6013.0,4020.0,3238.0,2448.0,1541.0,935.9,630.0,897.6};
    TMatrixD tMat_debug_prior_true_signal(10,1,array_debug_prior_true_signal);

    //// Even simpler debugging test, which turns off the Wiener filter and only performs an inversion 
    // Initialize unfolder
    //constexpr int NUM_DAGOSTINI_ITERATIONS = 3;
    std::unique_ptr< Unfolder > unfolder (
					  new WienerSVDUnfolder( false,
            WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv )
            //new DAgostiniUnfolder( NUM_DAGOSTINI_ITERATIONS ) 
					  );
    
//    // Initialize unfolder
//    //constexpr int NUM_DAGOSTINI_ITERATIONS = 3;
//    std::unique_ptr< Unfolder > unfolder (
//					  new WienerSVDUnfolder( true,
//            WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv )
//            //new DAgostiniUnfolder( NUM_DAGOSTINI_ITERATIONS ) 
//					  );
    
    // Run unfolder
    auto result = unfolder->unfold( tMat_debug_data, tMat_debug_cov, tMat_debug_response, tMat_debug_prior_true_signal );

//     //Pull out unfolded signal and covariance from results
//     auto tMat_data_signal_unfolded = result.unfolded_signal_.get();
//     auto tMat_unfolded_covariance = result.cov_matrix_.get();
//     auto tMat_errprop_matrix = result.err_prop_matrix_.get();

    return;
}

