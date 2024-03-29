
// Method to get gLEE tree loaded for proper looping-over
// Originally from /uboone/app/users/markrl/useful_scripts/plothelper.h
TTree* loadgLEE(std::string filename, std::string int_dir){

    TFile*f = new TFile((filename).c_str(),"read");
    TTree*v = (TTree*)f->Get((int_dir+"/vertex_tree").c_str());
    TTree*s = (TTree*)f->Get((int_dir+"/simple_tree").c_str());
    TTree*e = (TTree*)f->Get((int_dir+"/eventweight_tree").c_str());
    v->AddFriend(s);
    v->AddFriend(e);

    return v;
}

// Method that creates response matrix to be used in 1D NCpi0 xsec extraction
void ResponseMaker(std::string outDir){

    // -----------------------------------------------------
    // Interface with gLEE tuple 
    // -----------------------------------------------------

    // Despite being labeled as 2g1p, because this file is from an earlier stage of the selection, it is inclusive of 2g0p and can be used to derive both response 2g0p and 2g1p matrices
    TTree * v = (TTree*)loadgLEE("/pnfs/uboone/persistent/users/markross/Jan2022_gLEE_files/NCPi0CrossSection/2g1p_v4/sbnfit_2g1p_NextGen_v4_stage_-1_ext_Denom_NCPi0_CutFromBNB_Run123_v50.5.root","singlephoton");

    std::cout<<v->GetEntries()<<std::endl;

    // Prescription for calculating reconstructed pion momentum
    // In hive framework, this variable is aliased as "reco_pion_momentum"
    // Pmom_pi = sqrt(Ppix^2+Ppiy^2+Ppiz^2);
    // Ppix = Pgamma1_x+Pgamma2_x;
    // Pgamma_x = energy*dir_x;
    std::string reco_var =
    "sqrt(\
      ( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_implied_dirx[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_implied_dirx[(i_shr[1])] )\
     *( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_implied_dirx[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_implied_dirx[(i_shr[1])] )\
      +\
      ( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_implied_diry[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_implied_diry[(i_shr[1])] )\
     *( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_implied_diry[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_implied_diry[(i_shr[1])] )\
      +\
      ( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_implied_dirz[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_implied_dirz[(i_shr[1])] )\
     *( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_implied_dirz[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_implied_dirz[(i_shr[1])] )\
    )/1000";
 
    // Prescription for calculating truth pion momentum (this is already stored in a branch)
    std::string true_var = "mctruth_exiting_pi0_mom";

    TTreeFormula * reco_form = new TTreeFormula("recof",reco_var.c_str(), v);
    TTreeFormula * true_form = new TTreeFormula("truef",true_var.c_str(), v);

    // Prescription for calculating event weight
    // -----------------------------------------

    // Bit that's common to 2g1p and 2g0p
    std::string additional_weight_common =
    "(\
      (simple_pot_weight*6.7873e20/4.9669582e+21)\
     *(m_flash_optfltr_pe_beam >20 && m_flash_optfltr_pe_veto < 20)\
     *(MCFlux_NuPosX > (0.0-1.55) && MCFlux_NuPosX < (256.35-1.55) && MCFlux_NuPosY > (-116.5+0.97) && MCFlux_NuPosY < (116.5+0.97) && MCFlux_NuPosZ > 0.0+0.1 && MCFlux_NuPosZ < 1036.8+0.1)\
     *( (run_number >= 4952 && run_number <= 7770)*0.943100 + ( (run_number >= 8317 && run_number <=  13696) || (run_number >= 13697 && run_number <= 14116) || (run_number >= 14117 && run_number <= 18960) )*1.020139 )\
    )";

    // Bit that's unique to 2g1p
    std::string additional_weight_2g1p = additional_weight_common + "*(Sum$(mctruth_exiting_proton_energy-0.93827 > 0.05)==1)";

    // Bit that's unique to 2g0p
    std::string additional_weight_2g0p = additional_weight_common + "*(Sum$(mctruth_exiting_proton_energy-0.93827 > 0.05)==0)";

    // Prescription for determining which events 
    // pass selection cuts
    // -----------------------------------------
    TTreeFormula * pass_form_2g1p = new TTreeFormula("pass","(simple_2g1p_NextGen_v4COSMIC_mva>0.894 && simple_2g1p_NextGen_v4BNB_mva >0.737)", v);
    //TTreeFormula * pass_form_2g0p = new TTreeFormula("pass","(simple_2g0p_NextGen_v4COSMIC_mva>0.944 && simple_2g0p_NextGen_v4BNB_mva >0.731)", v);
    TTreeFormula * pass_form_2g0p = new TTreeFormula("pass","(simple_2g1p_NextGen_v4COSMIC_mva>0.894 && simple_2g1p_NextGen_v4BNB_mva >0.737)", v); // Placeholder

    // Prescription for normalizing selection 
    // and implementing some cuts
    // -----------------------------------------
    TTreeFormula * norm_form_2g1p = new TTreeFormula("norm",(additional_weight_2g1p).c_str(), v);
    TTreeFormula * norm_form_2g0p = new TTreeFormula("norm",(additional_weight_2g0p).c_str(), v);

    // -----------------------------------------------------
    // Construct response matrix 
    // -----------------------------------------------------

    // Specify reco and true analysis binning schemes
    std::vector<double> reco_bins = {0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675,  0.9};
    std::vector<double> true_bins = {0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675,  0.9};

    // Create htrue, hreco, and response hists
    TH1D * htrue_2g1p = new TH1D("htrue_2g1p","htrue_2g1p",true_bins.size()-1, &true_bins[0]);
    TH1D * htrue_2g0p = new TH1D("htrue_2g0p","htrue_2g0p",true_bins.size()-1, &true_bins[0]);
    TH1D * hreco_2g1p = new TH1D("hreco_2g1p","hreco_2g1p",reco_bins.size()-1, &reco_bins[0]);
    TH1D * hreco_2g0p = new TH1D("hreco_2g0p","hreco_2g0p",reco_bins.size()-1, &reco_bins[0]);
    TH2D * resp_2g1p = new TH2D("Response_2g1p","Response_2g1p",reco_bins.size()-1, &reco_bins[0], true_bins.size()-1, &true_bins[0]);
    TH2D * resp_2g0p = new TH2D("Response_2g0p","Response_2g0p",reco_bins.size()-1, &reco_bins[0], true_bins.size()-1, &true_bins[0]);

    // Create response TMatrix (large enough to include underflow/underflow from corresponding TH2D)
    TMatrixD mat_2g1p(reco_bins.size()+1,true_bins.size()+1);
    mat_2g1p.Zero();
    TMatrixD mat_2g0p(reco_bins.size()+1,true_bins.size()+1);
    mat_2g0p.Zero();

    // Loop through event tree; fill htrue, hreco, and response hists
    for(int i=0; i<v->GetEntries();i++){

        v->GetEntry(i);

        // Technically not necessary, but Mark says that this avoids some inconsistent obscure ROOT bug
        reco_form->GetNdata();
        true_form->GetNdata();
        pass_form_2g1p->GetNdata();
        pass_form_2g0p->GetNdata();
        norm_form_2g1p->GetNdata();
        norm_form_2g0p->GetNdata();

        double r = reco_form->EvalInstance();
        double t = true_form->EvalInstance();
        double p_2g1p = pass_form_2g1p->EvalInstance();
        double p_2g0p = pass_form_2g0p->EvalInstance();
        double w_2g1p = norm_form_2g1p->EvalInstance();
        double w_2g0p = norm_form_2g0p->EvalInstance();

        if(i%20000==0)std::cout<<i<<" "<<v->GetEntries()<<std::endl;
    //    std::cout<<r<<"\t\t"<<t<<"\t\t"<<p<<std::endl;

        if(p_2g1p){
            resp_2g1p->Fill(r,t,w_2g1p);
            htrue_2g1p->Fill(t,w_2g1p);
            hreco_2g1p->Fill(r,w_2g1p);
        }
        else if(p_2g0p){
            resp_2g0p->Fill(r,t,w_2g0p);
            htrue_2g0p->Fill(t,w_2g0p);
            hreco_2g0p->Fill(r,w_2g0p);
        }
        else{
            //resp->Fill(-999,t,w);
            htrue_2g1p->Fill(t,w_2g1p);
            htrue_2g0p->Fill(t,w_2g0p);
        }

    }
    std::cout << "Finished main loop; constructed htrue, hreco, and response hists for both 2g1p and 2g0p signal definitions." << std::endl;

    for(int i=0;i< hreco_2g1p->GetNbinsX()+2; i++){ // Both loops need to go to nbins+1 to capture the overflow
        for(int a=0; a<htrue_2g1p->GetNbinsX()+2; a++){ // It doesn't matter whether we loop over 2g1p or 2g0p to get the bin counts
            
            //2g1p  
            double eff_2g1p = resp_2g1p->GetBinContent(i,a)/htrue_2g1p->GetBinContent(a);
            resp_2g1p->SetBinContent(i,a, eff_2g1p); // Resets content of Response TH2D so there is consistency
            mat_2g1p(i,a) = resp_2g1p->GetBinContent(i,a); 

            //2g0p  
            double eff_2g0p = resp_2g0p->GetBinContent(i,a)/htrue_2g0p->GetBinContent(a);
            resp_2g0p->SetBinContent(i,a, eff_2g0p); // Resets content of Response TH2D so there is consistency
            mat_2g0p(i,a) = resp_2g0p->GetBinContent(i,a); 

        }
    }
    std::cout << "Constructed response TMatrices." << std::endl;
    
    // -----------------------------------------------------
    // Construct output file 
    // -----------------------------------------------------

    // Create output file
    TFile *fout = new TFile((outDir+"/response_matrices_exclusive.root").c_str(),"RECREATE");
    fout->cd();

    // Write out response matrix TH2Ds
    std::cout << "Writing response matrix TH2Ds to output file" << std::endl;
    resp_2g1p->GetXaxis()->SetTitle("Reco Var");
    resp_2g1p->GetYaxis()->SetTitle("True Var");
    resp_2g1p->Write();
    resp_2g0p->GetXaxis()->SetTitle("Reco Var");
    resp_2g0p->GetYaxis()->SetTitle("True Var");
    resp_2g0p->Write();

    // Write out response matrix TMatrices
    std::cout << "Writing response matrix TMatrices to output file" << std::endl;
    mat_2g1p.Write("response_matrix");
    mat_2g0p.Write("response_matrix");
    
    // Write out htrue and hreco
    std::cout << "Writing htrue and hreco TH1Ds to output file" << std::endl;
    htrue_2g1p->GetXaxis()->SetTitle("True Var");
    hreco_2g1p->GetXaxis()->SetTitle("Reco Var");
    htrue_2g1p->Write();
    hreco_2g1p->Write();
    htrue_2g0p->GetXaxis()->SetTitle("True Var");
    hreco_2g0p->GetXaxis()->SetTitle("Reco Var");
    htrue_2g0p->Write();
    hreco_2g0p->Write();

    // Close output file
    fout->Close();
    return ;

}
