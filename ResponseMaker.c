
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
    std::string additional_weight=
    "(\
      (simple_pot_weight*6.7873e20/4.9669582e+21)\
     *(m_flash_optfltr_pe_beam >20 && m_flash_optfltr_pe_veto < 20)\
     *(MCFlux_NuPosX > (0.0-1.55) && MCFlux_NuPosX < (256.35-1.55) && MCFlux_NuPosY > (-116.5+0.97) && MCFlux_NuPosY < (116.5+0.97) && MCFlux_NuPosZ > 0.0+0.1 && MCFlux_NuPosZ < 1036.8+0.1)\
     *( (run_number >= 4952 && run_number <= 7770)*0.943100 + ( (run_number >= 8317 && run_number <=  13696) || (run_number >= 13697 && run_number <= 14116) || (run_number >= 14117 && run_number <= 18960) )*1.020139 )\
     *(Sum$(mctruth_exiting_proton_energy-0.93827 > 0.05)==1)\
    )";

    // Prescription for determining which events pass selection cuts
    // This would change for 2g0p
    TTreeFormula * pass_form = new TTreeFormula("pass","(simple_2g1p_NextGen_v4COSMIC_mva>0.894 && simple_2g1p_NextGen_v4BNB_mva >0.737)", v);

    // Prescription for normalizing selection and implementing some cuts
    TTreeFormula * norm_form = new TTreeFormula("norm",(additional_weight).c_str(), v);

    // -----------------------------------------------------
    // Construct response matrix 
    // -----------------------------------------------------

    // Specify reco and true analysis binning schemes
    std::vector<double> reco_bins = {0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675,  0.9};
    std::vector<double> true_bins = {0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675,  0.9};

    // Create htrue, hreco, and response hists
    TH1D * htrue = new TH1D("htrue","htrue",true_bins.size()-1, &true_bins[0]);
    TH1D * hreco = new TH1D("hreco","hreco",reco_bins.size()-1, &reco_bins[0]);
    TH2D * resp = new TH2D("Response","Response",reco_bins.size()-1, &reco_bins[0], true_bins.size()-1, &true_bins[0]);

    // Create response TMatrix (large enough to include underflow/underflow from corresponding TH2D)
    TMatrixD mat(reco_bins.size()+1,true_bins.size()+1);
    mat.Zero();

    // Loop through event tree; fill htrue, hreco, and response hists
    for(int i=0; i<v->GetEntries();i++){

        v->GetEntry(i);

        // Technically not necessary, but Mark says that this avoids some inconsistent obscure ROOT bug
        reco_form->GetNdata();
        true_form->GetNdata();
        pass_form->GetNdata();
        norm_form->GetNdata();

        double r = reco_form->EvalInstance();
        double t = true_form->EvalInstance();
        double p = pass_form->EvalInstance();
        double w = norm_form->EvalInstance();

        if(i%20000==0)std::cout<<i<<" "<<v->GetEntries()<<std::endl;
    //    std::cout<<r<<"\t\t"<<t<<"\t\t"<<p<<std::endl;

        if(p){
            resp->Fill(r,t,w);
            htrue->Fill(t,w);
            hreco->Fill(r,w);
        }
        else{
            //resp->Fill(-999,t,w);
            htrue->Fill(t,w);
        }

    }
    std::cout << "Finished main loop; constructed htrue, hreco, and response hists." << std::endl;

    for(int i=0;i< hreco->GetNbinsX()+2; i++){ // Both loops need to go to nbins+1 to capture the overflow
        for(int a=0; a<htrue->GetNbinsX()+2; a++){
         
            double eff = resp->GetBinContent(i,a)/htrue->GetBinContent(a);
            resp->SetBinContent(i,a, eff); // Resets content of Response TH2D so there is consistency
            mat(i,a) = resp->GetBinContent(i,a); 

        }
    }
    std::cout << "Constructed response TMatrix." << std::endl;
    
    // -----------------------------------------------------
    // Construct output file 
    // -----------------------------------------------------

    // Create output file
    TFile *fout = new TFile((outDir+"/response_2g1p_exclusive_v3_d22_23.root").c_str(),"RECREATE");
    fout->cd();

    // Write out response matrix TH2D
    resp->GetXaxis()->SetTitle("Reco Var");
    resp->GetYaxis()->SetTitle("True Var");
    std::cout << "Writing response matrix TH2D to output file" << std::endl;
    resp->Write();

    // Write out response matrix TMatrix
    std::cout << "Writing response matrix TMatrix to output file" << std::endl;
    mat.Write("response_matrix");
    
    // Write out htrue and hreco
    std::cout << "Writing htrue and hreco TH1Ds to output file" << std::endl;
    htrue->GetXaxis()->SetTitle("True Var");
    hreco->GetXaxis()->SetTitle("Reco Var");
    htrue->Write();
    hreco->Write();

    // Close output file
    fout->Close();
    return ;

}
