void makeXsecRatioErrorSummary(){

    //******************* MC numbers **********************
    //                              comb    1p      0p
    std::vector<double> xs_genie = {1.68e-38,0.722e-38,0.775e-38}; 
    std::vector<double> xs_genie_scale  = {1.68,0.722,0.775}; 

    std::vector<double> E_genie =  {0.389e-38,0.183e-38,0.160e-28};
    std::vector<double> E_genie_scale =  {0.389,0.183,0.160};

    std::vector<double> uncertainty_total     = {0.19461,0.24458,0.24194};
    std::vector<double> uncertainty_stat      = {0.16024,0.10640,0.11965};
    std::vector<double> uncertainty_detector  = {0.08934,0.08157,0.06775};
    std::vector<double> uncertainty_flux      = {0.00794,0.13408,0.13218};
    std::vector<double> uncertainty_geant4    = {0.01530,0.02327,0.01700};
    std::vector<double> uncertainty_genie     = {0.06258,0.15109,0.14617};
    std::vector<double> uncertainty_pot       = {0.00000,0.02001,0.02001};
    std::vector<double> uncertainty_targets   = {0.00000,0.01000,0.01000};

    // ******************* Just Plotting Below****************8
    TCanvas *c = new TCanvas();//"c","c",100,500);
    c->cd();

    gStyle->SetOptStat(0);

    //************

    TH1D *tHist_uncertainty_total     = new TH1D("total","total",3,0.0,1.0);
    TH1D *tHist_uncertainty_stat      = new TH1D("stat","stat",3,0.0,1.0);
    TH1D *tHist_uncertainty_detector  = new TH1D("detector","detector",3,0.0,1.0);
    TH1D *tHist_uncertainty_flux      = new TH1D("flux","flux",3,0.0,1.0);
    TH1D *tHist_uncertainty_geant4    = new TH1D("geant4","geant4",3,0.0,1.0);
    TH1D *tHist_uncertainty_genie     = new TH1D("genie","genie",3,0.0,1.0);
    TH1D *tHist_uncertainty_pot       = new TH1D("pot","pot",3,0.0,1.0);
    TH1D *tHist_uncertainty_targets   = new TH1D("targets","targets",3,0.0,1.0);

    for(int i=0; i< 3; i++){
        tHist_uncertainty_total->SetBinContent(i+1,uncertainty_total[i]);
        tHist_uncertainty_stat->SetBinContent(i+1,uncertainty_stat[i]);
        tHist_uncertainty_detector->SetBinContent(i+1,uncertainty_detector[i]);
        tHist_uncertainty_flux->SetBinContent(i+1,uncertainty_flux[i]);
        tHist_uncertainty_geant4->SetBinContent(i+1,uncertainty_geant4[i]);
        tHist_uncertainty_genie->SetBinContent(i+1,uncertainty_genie[i]);
        tHist_uncertainty_pot->SetBinContent(i+1,uncertainty_pot[i]);
        tHist_uncertainty_targets->SetBinContent(i+1,uncertainty_targets[i]);
    }

    tHist_uncertainty_total->SetLineColor(kBlack);
    tHist_uncertainty_total->SetLineWidth(3);
    tHist_uncertainty_total->Draw("hist");

    tHist_uncertainty_total->SetMaximum(0.45);
    tHist_uncertainty_total->SetMinimum(0);
    tHist_uncertainty_total->SetTitle("");
    tHist_uncertainty_total->GetYaxis()->SetTitle("Fractional Uncertainty");
    tHist_uncertainty_total->GetYaxis()->SetTitleOffset(0.8);
    tHist_uncertainty_total->GetYaxis()->SetTitleSize(0.05);
    tHist_uncertainty_total->GetXaxis()->SetLabelOffset(999);
    tHist_uncertainty_total->GetXaxis()->SetLabelSize(0);
    tHist_uncertainty_total->GetXaxis()->SetTickLength(0.);

    tHist_uncertainty_stat->SetLineColor(kBlack);
    tHist_uncertainty_stat->SetLineStyle(2);
    tHist_uncertainty_stat->SetLineWidth(3);
    tHist_uncertainty_stat->DrawCopy("hist same");

    tHist_uncertainty_detector->SetLineColor(kOrange+1);
    tHist_uncertainty_detector->SetLineWidth(3);
    tHist_uncertainty_detector->DrawCopy("hist same");

    tHist_uncertainty_flux->SetLineColor(kBlue-7);
    tHist_uncertainty_flux->SetLineWidth(3);
    tHist_uncertainty_flux->DrawCopy("hist same");

    tHist_uncertainty_geant4->SetLineColor(kGreen-3);
    tHist_uncertainty_geant4->SetLineWidth(3);
    tHist_uncertainty_geant4->DrawCopy("hist same");

    tHist_uncertainty_genie->SetLineColor(kRed-7);
    tHist_uncertainty_genie->SetLineWidth(3);
    tHist_uncertainty_genie->DrawCopy("hist same");

    tHist_uncertainty_pot->SetLineColor(kCyan);
    tHist_uncertainty_pot->SetLineWidth(3);
    tHist_uncertainty_pot->DrawCopy("hist same");

    tHist_uncertainty_targets->SetLineColor(kMagenta);
    tHist_uncertainty_targets->SetLineWidth(3);
    tHist_uncertainty_targets->DrawCopy("hist same");

    // Vertical lines separating panels
    TLine *l1 = new TLine(1.0/3.0,0,1.0/3.0,0.45);
    TLine *l2 = new TLine(2.0/3.0,0,2.0/3.0,0.45);
    l1->SetLineColor(kBlack);
    l1->SetLineWidth(3);
    l2->SetLineColor(kBlack);
    l2->SetLineWidth(3);
    l1->Draw();
    l2->Draw();

    TH1D *tHist_dummy = new TH1D("dummy","dummy",3,0.0,1.0);
    tHist_dummy->SetLineColor(kWhite);

//    // Panel labels
//    TLegend *l = new TLegend(0.10,0.00,0.89,0.10);
//    //l->AddEntry(tHist_dummy,"#splitline{#sigma_{NC 1 #pi^{0} + 1 p}/}{#sigma_{NC 1 #pi^{0} + 0 p}}","lp");
//    l->AddEntry(tHist_dummy,"#frac{#sigma_{NC 1 #pi^{0} + 1 p}}{#sigma_{NC 1 #pi^{0} + 0 p}}","lp");
//    l->AddEntry(tHist_dummy,"#splitline{Exclusive}{NC 1 #pi^{0} + 1 proton}","lp");
//    l->AddEntry(tHist_dummy,"#splitline{Exclusive}{NC 1 #pi^{0} + 0 proton}","lp");
//    l->SetNColumns(3);
//    l->SetLineColor(kWhite);
//    l->SetLineWidth(0);
//    l->SetFillStyle(0);
//    l->Draw();

    // Legend to display ratio data point
    TLegend *leg1 = new TLegend(0.13,0.03,0.4,0.07);
    //TLegend *l = new TLegend(0.10,0.00,0.89,0.10);
    leg1->AddEntry(tHist_dummy,"#frac{#sigma_{NC 1 #pi^{0} + 1 p}}{#sigma_{NC 1 #pi^{0} + 0 p}}","lp");
    leg1->SetNColumns(1);
    leg1->SetLineColor(kWhite);
    leg1->SetLineWidth(0);
    leg1->SetFillStyle(0);
    leg1->Draw();

    // Legend to display exclusive measurement data points
    TLegend *leg2 = new TLegend(0.34,0.00,0.875,0.10);
    //TLegend *l = new TLegend(0.10,0.00,0.89,0.10);
    leg2->AddEntry(tHist_dummy,"#splitline{Exclusive}{NC 1 #pi^{0} + 1 proton}","lp");
    leg2->AddEntry(tHist_dummy,"#splitline{Exclusive}{NC 1 #pi^{0} + 0 proton}","lp");
    leg2->SetNColumns(2);
    leg2->SetLineColor(kWhite);
    leg2->SetLineWidth(0);
    leg2->SetFillStyle(0);
    leg2->Draw();

    // Curve labels
    TLegend *lb = new TLegend(0.375,0.60,0.89,0.875);
    lb->AddEntry(tHist_uncertainty_total,"Total","l");
    lb->AddEntry(tHist_uncertainty_stat,"Statistical","l");
    lb->AddEntry(tHist_uncertainty_detector,"Detector","l");
    lb->AddEntry(tHist_uncertainty_flux,"Flux","l");
    lb->AddEntry(tHist_uncertainty_geant4,"Geant4","l");
    lb->AddEntry(tHist_uncertainty_genie,"GENIE","l");
    lb->AddEntry(tHist_uncertainty_pot,"POT","l");
    lb->AddEntry(tHist_uncertainty_targets,"N_{targets}","l");
    lb->SetNColumns(2);
    lb->SetLineColor(kBlack);
    lb->Draw();

    // "MicroBooNE" label
    TLatex latex;
    latex.SetTextSize(0.05);
    latex.SetTextAlign(13);  //align at top
    latex.DrawLatex(.05,0.43,"MicroBooNE");

    // Save plot as image
    c->SaveAs("/uboone/data/users/finer/gLEE/NCPi0/2023-02_ratio-work/errorSummary_xsec_ratio_panels.pdf","pdf");

}
