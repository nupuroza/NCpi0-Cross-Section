//#include "TFile.h"
#include "TH2D.h"
#include "PlotUtils/MnvH1D.h"

using namespace std;

void MAT_test()
{
    TFile* f = new TFile("/uboone/data/users/finer/gLEE/NCPi0/2022-06-28_revampTest/2022-06-28_out.root", "READ");
    
    // Pull out MnvH1D from input file
    PlotUtils::MnvH1D *mHist_true_sig = (PlotUtils::MnvH1D*)f->Get("effNum_2g1p_inclusive");  
    TH1D true_sig = mHist_true_sig->GetCVHistoWithStatError();  
    TH2D* resp = (TH2D*)f->Get("tHist2D_response_2g1p_inclusive");

    return;
}

