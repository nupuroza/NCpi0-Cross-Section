#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
//from gardiner
#include <iomanip>
#include <sstream>

#include "TRandom3.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"

//base class from gardiner
/*
#include "ConstrainedCalculator.hh"
#include "FiducialVolume.hh"
#include "MCC9Unfolder.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"
*/
#include "WienerSVDUnfolder.hh"
#include "DAgostiniUnfolder.hh"

using namespace std;

//***TODO*** keep them in utils later 
//from 2017 WSVD github
void V2H(const TVectorD vec, TH1D* histo)
{
    // Fill vector to histogram,                                                                                                                                                   
    for(Int_t i=0; i<vec.GetNrows(); i++)
    {
        histo->SetBinContent(i+1, vec(i));
    }
}

void H2V(const TH1D* histo, TVectorD& vec)
{
    // Fill 1D histogram into matrix                                                                                                                                               
    for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
        vec(i) = histo->GetBinContent(i+1);
    }
}

void H2M(const TH2F* histo, TMatrixD& mat, bool rowcolumn)
{
    // Fill 2D histogram into matrix                                                                                                                                               
    // If TH2D(i, j) = Matrix(i, j), rowcolumn = kTRUE, else rowcolumn = kFALSE                                                                                                    
    for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
        for(Int_t j=0; j<histo->GetNbinsY(); j++)
        {
            if(rowcolumn) mat(i, j) = histo->GetBinContent(i+1, j+1);
            else mat(j, i) = histo->GetBinContent(i+1, j+1);
        }
    }
}

void H2M1(const TH1F* histo, TMatrixD& mat, bool rowcolumn)
{
  // Fill 2D histogram into matrix                                                                                                                                                                        
  // If TH2D(i, j) = Matrix(i, j), rowcolumn = kTRUE, else rowcolumn = kFALSE                                                                                                                             
  for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
      if(rowcolumn) mat(i, 0) = histo->GetBinContent(i+1);
      //else mat(j, i) = histo->GetBinContent(i+1, j+1);
    }
}



void M2H(const TMatrixD& mat, TH2F* histo)
{
    // Fill matrix to histogram                                                                                                                                                    
  for(Int_t i=0; i<mat.GetNrows(); i++)
    {
      for(Int_t j=0; j<mat.GetNcols(); j++)
        {
	  histo->SetBinContent(i+1, j+1, mat(i, j));
        }
    }
}

void M2H1_spec(TMatrixD& mat, TH1F* histo)
{
  // Fill matrix to histogram                                                                                                                                                                               
  for(Int_t i=0; i<mat.GetNrows(); i++)
    {
      histo->SetBinContent(i+1, mat(i, 0));
    }
}



void M2H_spec(TMatrixD* mat, TH2F* histo)
{
  // Fill matrix to histogram                                                                                                                                                                             
  for(Int_t i=0; i<mat->GetNrows(); i++)
    {
      for(Int_t j=0; j<mat->GetNcols(); j++)
        {
          histo->SetBinContent(i+1, j+1, mat->operator()(i, j));
        }
    }
}

void M2H1(TMatrixD* mat, TH1F* histo)
{
  // Fill matrix to histogram                                                                                                                                                                             
  for(Int_t i=0; i<mat->GetNrows(); i++)
    {
      histo->SetBinContent(i+1, mat->operator()(i, 0));
    }
}

void unfold_final()
{
    TFile* f = new TFile("/uboone/data/users/noza/gLEE/xsection/2022-03-17_out.root", "READ");
    
    // Ingredients - true_signal, data_signal, response, covariance  
    TH1F *true_sig = (TH1F*)f->Get("tHist_effNum_2gnp_inclusive");  
    TH1F *mes_sig = (TH1F*)f->Get("tHist_evtRate_2gnp_inclusive");
    TH2F* resp = (TH2F*)f->Get("tHist2D_response_2gnp");
    TH2F* cov = (TH2F*)f->Get("tHist2D_cov_evtRate_2gnp_inclusive");

    Int_t n = true_sig->GetNbinsX();
    Double_t Nuedges[n+1];
    for(int i=0; i<n+1; i++){
        Nuedges[i] = true_sig->GetBinLowEdge(i+1);
    }

    Int_t m = mes_sig->GetNbinsX();

    // construct matrices (for 1D/2D histogram) for input - m=n=8 (for ccpi0)
    TMatrixD signal(n,1);
    TMatrixD measure(m,1);
    TMatrixD response(m, n);
    TMatrixD covariance(m, m);
   
    //convert ingredients into TMatrixD
    H2M1(true_sig, signal, kTRUE);
    H2M1(mes_sig, measure, kTRUE);
    H2M(resp, response, kTRUE);
    H2M(cov, covariance, kTRUE);

    // construct to record additinal smearing matrix and wiener filter (diagomal matrix) elements. 
    TH2F* smear = new TH2F("smear","Additional Smearing Matirx",n,Nuedges,n,Nuedges);
    TH1F* unf_signal = new TH1F("unf_signal","unfolded spectrum",n,Nuedges);
    TH2F* unf_cov = new TH2F("unf_cov","Unfolded covariance", n, Nuedges, n, Nuedges);

    //**TODO** argc, argv later 
    TFile* file = new TFile("/uboone/data/users/noza/gLEE/xsection/2022-03-17_unfolded.root", "RECREATE"); 

    //Wiener-SVD

    std::unique_ptr< Unfolder > unfolder (
					  new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kIdentity )
					  );
    auto result = unfolder->unfold( measure, covariance, response, signal );
    
    WienerSVDUnfolder* wsvd_ptr = dynamic_cast< WienerSVDUnfolder* >(
								     unfolder.get() );
    TMatrixD A_C_(m,m);
    TMatrixD unfold_sig_(m,1);
    TMatrixD unfold_cov_(m,m);
    
    if ( wsvd_ptr ) {

      //additional smearing matrix
      auto A_C = wsvd_ptr->additional_smearing_matrix();
      A_C_ = A_C;
      M2H(A_C, smear);
    }
    auto unfold_sig = result.unfolded_signal_.get();

    auto unfold_cov = result.cov_matrix_.get();
    
    M2H_spec(unfold_cov, unf_cov);
    
    M2H1(unfold_sig, unf_signal);

    H2M1(unf_signal, unfold_sig_, kTRUE);
    H2M(unf_cov, unfold_cov_, kTRUE);

    for (int b=1;b<9;b++) unf_signal->SetBinError(b,sqrt((*unfold_cov)(b-1,b-1)));

    TH1F* diff = new TH1F("diff","Fractional difference of unfolded signal and true signal model",n, Nuedges);
    for(Int_t i=1; i<=n; i++)
      {
        Double_t s2s = 0;
        Double_t u = unf_signal->GetBinContent(i);
        Double_t t = signal(i-1, 0);
        if(t!=0) s2s = u/t - 1;
        else s2s = 1.; // deal with  t=0
        diff->SetBinContent(i, s2s); // in percentage 
      }
   
    // intrinsic bias (Ac-I)*s_bar formula
    TH1F* bias = new TH1F("bias","intrinsic bias w.r.t. model",n, Nuedges);
    TH1F* bias2 = new TH1F("bias2","intrinsic bias2 w.r.t. unfolded result",n, Nuedges);
    TMatrixD unit(n,n);
    unit.UnitMatrix();
    
    TMatrixD intrinsicbias = (A_C_ - unit)*signal;
    TMatrixD intrinsicbias2 = (A_C_ - unit)*unfold_sig_;
    
    TMatrixD Smeared_true_sig = A_C_ * signal;
    
    TH1F* smeared_true_sig = new TH1F("smeared_true_sig", "smeared true signal", n, Nuedges);
    M2H1_spec(Smeared_true_sig , smeared_true_sig);
    
    for(int i=0; i<n; i++)
      {
        if(signal(i,0)!=0) intrinsicbias(i,0) = intrinsicbias(i,0)/signal(i,0);
        else intrinsicbias(i,0)=0.;
        if(unfold_sig_(i,0)!=0) intrinsicbias2(i,0) = intrinsicbias2(i,0)/unfold_sig_(i,0);
        else intrinsicbias2(i,0)=0.;
      }
    
    M2H1_spec(intrinsicbias, bias);
    M2H1_spec(intrinsicbias2, bias2);

    // diagonal uncertainty
    TH1F* fracError = new TH1F("fracError", "Fractional uncertainty", n, Nuedges);
    TH1F* absError = new TH1F("absError", "absolute uncertainty", n, Nuedges);
    for(int i=1; i<=n; i++)
      {
        fracError->SetBinContent(i, TMath::Sqrt(unfold_cov_(i-1, i-1))/unfold_sig_(i-1, 0));
        absError->SetBinContent(i, TMath::Sqrt(unfold_cov_(i-1, i-1)));
      }
   
    /// MSE
    TH1F* MSE = new TH1F("MSE", "Mean Square Error: variance+bias^2", n, Nuedges);
    TH1F* MSE2 = new TH1F("MSE2", "Mean Square Error: variance", n, Nuedges);
    for(int i=0; i<n; i++)
      {
        MSE->SetBinContent(i+1, TMath::Power(intrinsicbias2(i,0)*unfold_sig_(i,0),2)+unfold_cov_(i,i));
        MSE2->SetBinContent(i+1, unfold_cov_(i,i));
      }

  
    smear->Write();
    unf_signal->Write();
    unf_cov->Write();
    diff->Write();
    bias->Write();
    bias2->Write();
    fracError->Write();
    absError->Write();
    MSE->Write();
    MSE2->Write();
    smeared_true_sig->Write();
    file->Close();

    return 1;
}

