void makeXsecSummaryPlot(){

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


    //***************** Inclusive Combined ***********************
    double xs_comb = 1.243e-38; 
    double xs_comb_scale  = xs_comb/1e-38 ;

    double E_comb =  sqrt(pow(0.185,2)+pow(0.076,2))*1e-38;//1.9332595587e-39;
    double E_comb_scale = E_comb/1e-38;

    double Estat_comb =  0.076e-38;
    double Estat_comb_scale = Estat_comb/1e-38;

    std::cout<<" XS for comb "<<xs_comb_scale<<" +/- "<<E_comb_scale<<" x 10^-38 cm^2/Atom with "<<std::endl;


    //******************* MC numbers **********************
    //                              comb    1p      0p
    std::vector<double> xs_genie = {1.6781555267e-38, 0.722186330878e-39, 0.774476432059-38}; 
    std::vector<double> xs_genie_scale = {1.6781555267, 0.722186330878, 0.774476432059}; 

    std::vector<double> E_genie =  {0.39057457925e-38, 0.184474803085e-38,0.162209481224e-38};
    std::vector<double> E_genie_scale =  {0.39057457925, 0.184474803085,0.162209481224};
    //                              comb    1p      0p
    std::vector<double> xs_genie2_scale = {1.899,0.72318,1.0033};
    std::vector<double> xs_neut_scale = {1.5003,0.61092,0.76273};
    std::vector<double> xs_nuwro_scale = {1.6456,0.58938,0.96381};
    std::vector<double> xs_gibuu_scale = {1.6944,0.73226,0.77836};

    std::cout<<"NUIS: "<<(xs_gibuu_scale[0]-xs_genie_scale[0])/(xs_gibuu_scale[0])*100.0<<"\n";
    std::cout<<"NUIS1p: "<<(xs_gibuu_scale[1]-xs_genie_scale[1])/(xs_gibuu_scale[1])*100.0<<"\n";
    std::cout<<"NUIS0p: "<<(xs_gibuu_scale[2]-xs_genie_scale[2])/(xs_gibuu_scale[2])*100.0<<"\n";
    // std::cout<<" XS for genie "<<xs_genie_scale<<" +/- "<<E_genie_scale<<" x 10^-38 cm^2/Atom "<<std::endl;

    // Chi^2 nonsence

    std::vector<std::string> nams = {"Inc","2g1p","2g0p"};
    std::vector<double> data = {xs_comb_scale,xs_2g1p_scale,xs_2g0p_scale};
    std::vector<double> errdata = {E_comb_scale,E_2g1p_scale,E_2g0p_scale};
    for(int i=0; i< 3; i++){

        double genie3_chi = pow(xs_genie_scale[i]-data[i],2)/pow(errdata[i],2);
        double genie2_chi = pow(xs_genie2_scale[i]-data[i],2)/pow(errdata[i],2);
        double neut_chi = pow(xs_neut_scale[i]-data[i],2)/pow(errdata[i],2);
        double nuwro_chi = pow(xs_nuwro_scale[i]-data[i],2)/pow(errdata[i],2);

        std::cout<<" ------- "<<nams[i]<<" ----- "<<std::endl;
        std::cout<<"Genie3 "<<genie3_chi<<" sqrt(chi): "<<sqrt(genie3_chi)<<std::endl;  
        std::cout<<"Bits "<<xs_genie_scale[i]<<" "<<data[i]<<" "<<errdata[i]<<std::endl;
        std::cout<<"Genie2 "<<genie2_chi<<" sqrt(chi): "<<sqrt(genie2_chi)<<std::endl;
        std::cout<<"Neut "<<neut_chi<<" sqrt(chi): "<<sqrt(neut_chi)<<std::endl;
        std::cout<<"Nuwro "<<nuwro_chi<<" sqrt(chi): "<<sqrt(nuwro_chi)<<std::endl;

    }




    // ******************* Just Plotting Below****************8
    TCanvas *c = new TCanvas();//"c","c",100,500);
    c->cd();

    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(7);
    //gStyle->SetErrorX(0.00001);

    TH1D *Dat_comb = new TH1D("dat","dat",300,0.0,1.0);
    Dat_comb->SetBinContent(51,xs_comb_scale);
    Dat_comb->SetBinError(51,E_comb_scale); //15
    Dat_comb->SetLineColor(kBlack);
    Dat_comb->SetMarkerColor(kBlack);
    Dat_comb->SetLineWidth(2);
    Dat_comb->SetMarkerStyle(20);

    Dat_comb->Draw("E1p");
    Dat_comb->SetMaximum(2.5);
    Dat_comb->SetMinimum(0);
    //Dat_comb->SetAxisRange(0.4,0.6,"X");
    Dat_comb->SetTitle("");
    Dat_comb->GetYaxis()->SetTitle("#sigma_{NC 1 #pi^{0}} [10^{-38} cm^{2}/Atom]");
    Dat_comb->GetYaxis()->SetTitleOffset(0.8);
    Dat_comb->GetYaxis()->SetTitleSize(0.05);
    Dat_comb->GetXaxis()->SetLabelOffset(999);
    Dat_comb->GetXaxis()->SetLabelSize(0);
    Dat_comb->GetXaxis()->SetTickLength(0.);

    Dat_comb->DrawCopy("E1p same");

    Dat_comb->SetBinError(51,Estat_comb_scale);
    Dat_comb->Draw("E1p same");

    //************

    TH1D *Dat_MC = new TH1D("MC","MC",3,0.0,1.0);

    for(int i=0; i< xs_genie_scale.size(); i++){
        Dat_MC->SetBinContent(i+1,xs_genie_scale[i]);
        Dat_MC->SetBinError(i+1,E_genie_scale[i]); 
    }

    Dat_MC->SetLineColor(kRed-7);
    Dat_MC->SetLineWidth(4);
    Dat_MC->DrawCopy("hist same");
    TH1D * clone_Dat_MC = (TH1D*)Dat_MC->Clone("datmcclone");

    Dat_MC->SetFillColor(kRed-7);
    Dat_MC->SetFillStyle(3554);
    Dat_MC->Draw("E2 same");

    clone_Dat_MC->DrawCopy("hist same");

    TH1D *Dat_genie2 = new TH1D("genie2","genie2",3,0.0,1.0);
    TH1D *Dat_nuwro = new TH1D("nuwro","nuwro",3,0.0,1.0);
    TH1D *Dat_gibuu = new TH1D("gibuu","gibuu",3,0.0,1.0);
    TH1D *Dat_neut = new TH1D("neut","neut",3,0.0,1.0);

    for(int i=0; i< xs_genie_scale.size(); i++){
        Dat_genie2->SetBinContent(i+1,xs_genie2_scale[i]);
        Dat_neut->SetBinContent(i+1,xs_neut_scale[i]);
        Dat_nuwro->SetBinContent(i+1,xs_nuwro_scale[i]);
        Dat_gibuu->SetBinContent(i+1,xs_gibuu_scale[i]);
    }

    Dat_genie2->SetLineColor(kOrange+1);
    Dat_genie2->SetLineStyle(2);
    Dat_genie2->SetLineWidth(4);
    Dat_genie2->DrawCopy("hist same");

    Dat_neut->SetLineColor(kBlue-7);
    Dat_neut->SetLineStyle(9);
    Dat_neut->SetLineWidth(3);
    Dat_neut->DrawCopy("hist same");

    gStyle->SetLineStyleString(10,"80 15 25 15");

    Dat_nuwro->SetLineColor(kGreen-3);
    Dat_nuwro->SetLineStyle(10);
    Dat_nuwro->SetLineWidth(3);
    Dat_nuwro->DrawCopy("hist same");

    Dat_gibuu->SetLineColor(kMagenta);
    Dat_gibuu->SetLineStyle(9);
    Dat_gibuu->SetLineWidth(4);
    //Dat_gibuu->DrawCopy("hist same");

    //************

    TH1D *Dat_2g0p = new TH1D("dat2g0p","dat2g0p",300,0.0,1);
    Dat_2g0p->SetBinContent(251,xs_2g0p_scale);
    Dat_2g0p->SetBinError(251,E_2g0p_scale); //15
    Dat_2g0p->SetLineColor(kBlack);//kBlue -7
    Dat_2g0p->SetMarkerColor(kBlack);
    Dat_2g0p->SetLineWidth(2);
    Dat_2g0p->SetMarkerStyle(20);
    Dat_2g0p->DrawCopy("E1p same");

    Dat_2g0p->SetBinError(251,Estat_2g0p_scale);
    Dat_2g0p->DrawClone("E1p same");

    //************

    TH1D *Dat_2g1p = new TH1D("dat2g1p","dat2g1p",300,0.0,1.0);
    Dat_2g1p->SetBinContent(151,xs_2g1p_scale);
    Dat_2g1p->SetBinError(151,E_2g1p_scale); //15
    Dat_2g1p->SetLineColor(kBlack); //kGreen-3
    Dat_2g1p->SetMarkerColor(kBlack);
    Dat_2g1p->SetLineWidth(2);
    Dat_2g1p->SetMarkerStyle(20);
    Dat_2g1p->DrawCopy("E1p same");

    Dat_2g1p->SetBinError(151,Estat_2g1p_scale);
    Dat_2g1p->Draw("E1p same");

    //************

    TLine *l1 = new TLine(1.0/3.0,0,1.0/3.0,2.5);
    TLine *l2 = new TLine(2.0/3.0,0,2.0/3.0,2.5);
    l1->SetLineColor(kBlack);
    l1->SetLineWidth(3);
    l2->SetLineColor(kBlack);
    l2->SetLineWidth(3);
    l1->Draw();
    l2->Draw();

    TLegend *l = new TLegend(0.13,0.06,0.89,0.23);
    l->AddEntry(Dat_comb,"#splitline{Semi-Inclusive}{NC 1 #pi^{0}           }","lp");
    l->AddEntry(Dat_2g1p,"#splitline{Exclusive}{NC 1 #pi^{0} + 1 proton}","lp");
    l->AddEntry(Dat_2g0p,"#splitline{Exclusive}{NC 1 #pi^{0} + 0 proton}","lp");
    l->SetNColumns(3);
    l->SetLineColor(kWhite);
    l->SetLineWidth(0);
    l->SetFillStyle(0);
    l->Draw();

    TLegend *lb = new TLegend(0.45,0.55,0.85,0.85);
    lb->AddEntry(Dat_MC,"GENIE v3.0.6 (G18_10a_02_11a)","lf");
    lb->AddEntry(Dat_genie2,"GENIE v2.12.10","l");
    lb->AddEntry(Dat_neut,"NEUT v5.4.0.1","l");
    lb->AddEntry(Dat_nuwro,"NuWro v19.02.1","l");
    //lb->AddEntry(Dat_gibuu,"NUISENCE Genie v3","l");
    lb->SetNColumns(1);
    lb->SetLineColor(kBlack);
    lb->Draw();

    TLatex latex;
    latex.SetTextSize(0.05);
    latex.SetTextAlign(13);  //align at top
    latex.DrawLatex(.05,2.4,"MicroBooNE");

    c->SaveAs("/uboone/data/users/finer/gLEE/NCPi0/2023-02_ratio-work/xSec_summary_panels.pdf","pdf");


}

