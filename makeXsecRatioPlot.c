void makeXsecRatioPlot(){

    //----------------------------------------------------
    //------------- XSEC NUMBERS -------------------------
    //----------------------------------------------------

    //***************** 2g1p ***********************
    double xs_2g1p = 0.444e-38;
    double xs_2g1p_scale = xs_2g1p/1e-38;

    double E_2g1p = sqrt(pow(0.098,2)+pow(0.047,2))*1e-38;
    double E_2g1p_scale = E_2g1p/1e-38;

    double Estat_2g1p = 0.047e-38;
    double Estat_2g1p_scale = Estat_2g1p/1e-38;

    std::cout<<" XS for 2g1p "<<xs_2g1p_scale<<" +/- "<<E_2g1p_scale<<" x 10^-38 cm^2/Atom "<<std::endl;
    
    //***************** 2g0p ***********************
    double xs_2g0p = 0.624e-38; 
    double xs_2g0p_scale  = xs_2g0p/1e-38 ;

    double E_2g0p = sqrt(pow(0.131,2)+pow(0.075,2))*1e-38;
    double E_2g0p_scale = E_2g0p/1e-38;

    double Estat_2g0p = 0.075e-38;
    double Estat_2g0p_scale = Estat_2g0p/1e-38;

    std::cout<<" XS for 2g0p "<<xs_2g0p_scale<<" +/- "<<E_2g0p_scale<<" x 10^-38 cm^2/Atom "<<std::endl;

    //***************** 2g1p/2g0p Ratio ***********************
    double xs_ratio = 0.710216920178; 
    double xs_ratio_scale  = xs_ratio;

    double E_ratio = 0.13821;
    double E_ratio_scale = E_ratio;

    double Estat_ratio =  0.11381;
    double Estat_ratio_scale = Estat_ratio;

    std::cout<<" XS for ratio "<<xs_ratio_scale<<" +/- "<<E_ratio_scale<<" with "<<std::endl;

    //******************* MC numbers **********************
    //                              ratio    1p      0p
    std::vector<double> xs_genie = {0.93248,0.722e-38,0.775e-38}; 
    std::vector<double> xs_genie_scale  = {0.93248,0.722,0.775}; 

    std::vector<double> E_genie =  {0.102,0.183e-38,0.160e-38};
    std::vector<double> E_genie_scale =  {0.102,0.183,0.160};

    std::vector<double> xs_genie2_scale = {0.72318/1.0033,0.72318,1.0033};
    std::vector<double> xs_neut_scale = {0.61092/0.76273,0.61092,0.76273};
    std::vector<double> xs_nuwro_scale = {0.58938/0.96381,0.58938,0.96381};

    //----------------------------------------------------
    //------------- PLOTTING BELOW ----------------------
    //----------------------------------------------------

    TCanvas* c0 = new TCanvas("myCanvasName","The Canvas Title",800,600);

    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(7);
    gStyle->SetLineStyleString(10,"80 15 25 15");

    //----------------------------------------------------
    //------------- LEFT PANEL
    //----------------------------------------------------

    TPad *pad1 = new TPad("pad1", "pad1",0.0,0.0,0.33,1.0);
    pad1->SetLeftMargin(0.22);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    pad1->cd();

    TH1D *Dat_MC_ratio = new TH1D("MC_ratio","MC_ratio",1,0.0,1.0);
    TH1D *Dat_genie2_ratio = new TH1D("genie2_ratio","genie2_ratio",1,0.0,1.0);
    TH1D *Dat_nuwro_ratio = new TH1D("nuwro_ratio","nuwro_ratio",1,0.0,1.0);
    TH1D *Dat_neut_ratio = new TH1D("neut_ratio","neut_ratio",1,0.0,1.0);

    for(int i=0; i< 1; i++){
        Dat_MC_ratio->SetBinContent(i+1,xs_genie_scale[i]);
        Dat_MC_ratio->SetBinError(i+1,E_genie_scale[i]); 
        Dat_genie2_ratio->SetBinContent(i+1,xs_genie2_scale[i]);
        Dat_neut_ratio->SetBinContent(i+1,xs_neut_scale[i]);
        Dat_nuwro_ratio->SetBinContent(i+1,xs_nuwro_scale[i]);
    }

    Dat_MC_ratio->SetMaximum(1.75);
    Dat_MC_ratio->SetMinimum(0.25);
    Dat_MC_ratio->SetTitle("");
    Dat_MC_ratio->GetYaxis()->SetTitle("Ratio");
    Dat_MC_ratio->GetYaxis()->SetTitleOffset(0.8);
    Dat_MC_ratio->GetYaxis()->SetTitleSize(0.12);
    Dat_MC_ratio->GetYaxis()->SetLabelOffset(0.02);
    Dat_MC_ratio->GetYaxis()->SetLabelSize(0.07);
    Dat_MC_ratio->GetYaxis()->SetTickLength(0.06);
    Dat_MC_ratio->GetYaxis()->SetNdivisions(605,kTRUE);
    Dat_MC_ratio->GetXaxis()->SetLabelOffset(999);
    Dat_MC_ratio->GetXaxis()->SetLabelSize(0);
    Dat_MC_ratio->GetXaxis()->SetTickLength(0.);

    Dat_MC_ratio->SetLineColor(kRed-7);
    Dat_MC_ratio->SetLineWidth(4);
    Dat_MC_ratio->DrawCopy("hist");
    TH1D * clone_Dat_MC_ratio = (TH1D*)Dat_MC_ratio->Clone("datmcratioclone");

    Dat_MC_ratio->SetFillColor(kRed-7);
    Dat_MC_ratio->SetFillStyle(3554);
    Dat_MC_ratio->Draw("E2 same");
    
    clone_Dat_MC_ratio->DrawCopy("hist same");

    Dat_genie2_ratio->SetLineColor(kOrange+1);
    Dat_genie2_ratio->SetLineStyle(2);
    Dat_genie2_ratio->SetLineWidth(3);
    Dat_genie2_ratio->DrawCopy("hist same");

    Dat_neut_ratio->SetLineColor(kBlue-7);
    Dat_neut_ratio->SetLineStyle(9);
    Dat_neut_ratio->SetLineWidth(3);
    Dat_neut_ratio->DrawCopy("hist same");

    //!! gStyle->SetLineStyleString(10,"80 15 25 15");

    Dat_nuwro_ratio->SetLineColor(kGreen-3);
    Dat_nuwro_ratio->SetLineStyle(10);
    Dat_nuwro_ratio->SetLineWidth(3);
    Dat_nuwro_ratio->DrawCopy("hist same");

    //************

    TH1D *Dat_ratio = new TH1D("dat","dat",100,0.0,1.0);
    Dat_ratio->SetBinContent(51,xs_ratio_scale);
    Dat_ratio->SetBinError(51,E_ratio_scale); //15
    Dat_ratio->SetLineColor(kBlack);
    Dat_ratio->SetMarkerColor(kBlack);
    Dat_ratio->SetLineWidth(2);
    Dat_ratio->SetMarkerStyle(20);
    Dat_ratio->DrawCopy("E1p same");

    Dat_ratio->SetBinError(51,Estat_ratio_scale);
    Dat_ratio->Draw("E1p same");

    //************

    // Legend to display ratio data point
    TLegend *leg1 = new TLegend(0.33,0.12,0.85,0.2);
    leg1->AddEntry(Dat_ratio,"#frac{#sigma_{NC 1 #pi^{0} + 1 proton}}{#sigma_{NC 1 #pi^{0} + 0 proton}}","lp");
    leg1->SetNColumns(1);
    leg1->SetLineColor(kWhite);
    leg1->SetLineWidth(0);
    leg1->SetFillStyle(0);
    leg1->Draw();

    // "MicroBooNE" label
    TLatex latex;
    latex.SetTextSize(0.1);
    latex.SetTextAlign(13);  //align at top
    latex.DrawLatex(.05,1.7,"MicroBooNE");

    c0->cd();

    //----------------------------------------------------
    //------------- RIGHT PANELS
    //----------------------------------------------------

    TPad *pad2 = new TPad("pad2", "pad2",0.33,0.0,1.0,1.0);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();

    TH1D *Dat_MC = new TH1D("MC","MC",2,0.0,1.0);
    TH1D *Dat_genie2 = new TH1D("genie2","genie2",2,0.0,1.0);
    TH1D *Dat_nuwro = new TH1D("nuwro","nuwro",2,0.0,1.0);
    TH1D *Dat_neut = new TH1D("neut","neut",2,0.0,1.0);

    for(int i=1; i< xs_genie_scale.size(); i++){
        Dat_MC->SetBinContent(i,xs_genie_scale[i]);
        Dat_MC->SetBinError(i,E_genie_scale[i]); 
        Dat_genie2->SetBinContent(i,xs_genie2_scale[i]);
        Dat_neut->SetBinContent(i,xs_neut_scale[i]);
        Dat_nuwro->SetBinContent(i,xs_nuwro_scale[i]);
    }

    Dat_MC->SetMaximum(2.5);
    Dat_MC->SetMinimum(0);
    Dat_MC->SetTitle("");
    Dat_MC->GetYaxis()->SetTitle("#sigma_{NC 1 #pi^{0}} [10^{-38} cm^{2}/Atom]");
    Dat_MC->GetYaxis()->SetTitleOffset(1.2);
    Dat_MC->GetYaxis()->SetTitleSize(0.05);
    Dat_MC->GetXaxis()->SetLabelOffset(999);
    Dat_MC->GetXaxis()->SetLabelSize(0);
    Dat_MC->GetXaxis()->SetTickLength(0.);

    Dat_MC->SetLineColor(kRed-7);
    Dat_MC->SetLineWidth(4);
    Dat_MC->DrawCopy("hist");
    TH1D * clone_Dat_MC = (TH1D*)Dat_MC->Clone("datmcclone");

    Dat_MC->SetFillColor(kRed-7);
    Dat_MC->SetFillStyle(3554);
    Dat_MC->Draw("E2 same");
    
    clone_Dat_MC->DrawCopy("hist same");

    Dat_genie2->SetLineColor(kOrange+1);
    Dat_genie2->SetLineStyle(2);
    Dat_genie2->SetLineWidth(3);
    Dat_genie2->DrawCopy("hist same");

    Dat_neut->SetLineColor(kBlue-7);
    Dat_neut->SetLineStyle(9);
    Dat_neut->SetLineWidth(3);
    Dat_neut->DrawCopy("hist same");

    Dat_nuwro->SetLineColor(kGreen-3);
    Dat_nuwro->SetLineStyle(10);
    Dat_nuwro->SetLineWidth(3);
    Dat_nuwro->DrawCopy("hist same");

    //************

    TH1D *Dat_2g1p = new TH1D("dat2g1p","dat2g1p",200,0.0,1.0);
    Dat_2g1p->SetBinContent(51,xs_2g1p_scale);
    Dat_2g1p->SetBinError(51,E_2g1p_scale); //15
    Dat_2g1p->SetLineColor(kBlack); //kGreen-3
    Dat_2g1p->SetMarkerColor(kBlack);
    Dat_2g1p->SetLineWidth(2);
    Dat_2g1p->SetMarkerStyle(20);
    Dat_2g1p->DrawCopy("E1p same");

    Dat_2g1p->SetBinError(51,Estat_2g1p_scale);
    Dat_2g1p->Draw("E1p same");

    //************

    TH1D *Dat_2g0p = new TH1D("dat2g0p","dat2g0p",200,0.0,1.0);
    Dat_2g0p->SetBinContent(151,xs_2g0p_scale);
    Dat_2g0p->SetBinError(151,E_2g0p_scale); //15
    Dat_2g0p->SetLineColor(kBlack);//kBlue -7
    Dat_2g0p->SetMarkerColor(kBlack);
    Dat_2g0p->SetLineWidth(2);
    Dat_2g0p->SetMarkerStyle(20);
    Dat_2g0p->DrawCopy("E1p same");

    Dat_2g0p->SetBinError(151,Estat_2g0p_scale);
    Dat_2g0p->Draw("E1p same");

    TLine *l1 = new TLine(1.0/2.0,0,1.0/2.0,2.5);
    l1->SetLineColor(kBlack);
    l1->SetLineWidth(3);
    l1->Draw();

    // Legend to display exclusive measurement data points
    TLegend *leg2 = new TLegend(0.15,0.06,0.885,0.23);
    leg2->AddEntry(Dat_2g1p,"#splitline{Exclusive}{NC 1 #pi^{0} + 1 proton}","lp");
    leg2->AddEntry(Dat_2g0p,"#splitline{Exclusive}{NC 1 #pi^{0} + 0 proton}","lp");
    leg2->SetNColumns(2);
    leg2->SetLineColor(kWhite);
    leg2->SetLineWidth(0);
    leg2->SetFillStyle(0);
    leg2->Draw();

    TLegend *lb = new TLegend(0.2,0.55,0.85,0.85);
    lb->AddEntry(Dat_ratio, "Runs 1-3", "ep");
    lb->AddEntry(Dat_MC,"GENIE v3.0.6 (G18_10a_02_11a)","lf");
    lb->AddEntry(Dat_genie2,"GENIE v2.12.10","l");
    lb->AddEntry(Dat_neut,"NEUT v5.4.0.1","l");
    lb->AddEntry(Dat_nuwro,"NuWro v19.02.1","l");
    lb->SetNColumns(1);
    lb->SetLineColor(kBlack);
    lb->Draw();

    c0->Draw();

    //----------------------------------------------------
    //------------- WRITE OUT PLOT 
    //----------------------------------------------------

    // Save plot as image
    c0->SaveAs("/uboone/data/users/finer/gLEE/NCPi0/2023-02_ratio-work/xsec_ratio_panels.pdf","pdf");

}
