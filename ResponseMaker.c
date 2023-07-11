
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

    TTree * v_2g1p = (TTree*)loadgLEE("/pnfs/uboone/persistent/users/markross/Jan2022_gLEE_files/NCPi0CrossSection/2g1p_v4/sbnfit_2g1p_NextGen_v4_stage_-1_ext_Denom_NCPi0_CutFromBNB_Run123_v50.5.root","singlephoton");
    TTree * v_2g0p = (TTree*)loadgLEE("/pnfs/uboone/persistent/users/markross/Jan2022_gLEE_files/NCPi0CrossSection/2g0p_v4/sbnfit_2g0p_NextGen_v4_stage_-1_ext_Denom_NCPi0_CutFromBNB_Run123_v50.5.root","singlephoton");

    std::cout<<"Number of entries in 2g1p tuple: "<<v_2g1p->GetEntries()<<std::endl;
    std::cout<<"Number of entries in 2g0p tuple: "<<v_2g0p->GetEntries()<<std::endl;

    if (v_2g1p->GetEntries()!=v_2g0p->GetEntries()) {
        std::cout << "The number of entries in the 2g1p and 2g0p input tuples are unequal and parts of this script won't work correctly, so I'm exiting..." << std::endl;
        exit(1); // Terminate with a non-zero exit status
    }

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

    TTreeFormula * reco_form_2g1p = new TTreeFormula("recof",reco_var.c_str(), v_2g1p);
    TTreeFormula * reco_form_2g0p = new TTreeFormula("recof",reco_var.c_str(), v_2g0p);
    TTreeFormula * true_form_2g1p = new TTreeFormula("truef",true_var.c_str(), v_2g1p);
    TTreeFormula * true_form_2g0p = new TTreeFormula("truef",true_var.c_str(), v_2g0p);

    // Prescription for calculating event weight
    // -----------------------------------------

    // Bit that's common to 2g1p and 2g0p
    // 1. Scaling simulation POT to data POT. Simulation POT in denominator here. Data POT is in additional_weight.
    // 2. Only select events with 20 photons during beam window and fewer than 20 photons in veto time window.
    // 3. Remove events outside detector region.
    // 4. Scale Monte Carlo runs to number of events captured in actual data runs.
    std::string additional_weight_common =
    "(\
      (simple_pot_weight/4.9669582e+21)\
     *(m_flash_optfltr_pe_beam >20 && m_flash_optfltr_pe_veto < 20)\
     *(MCFlux_NuPosX > (0.0-1.55) && MCFlux_NuPosX < (256.35-1.55) && MCFlux_NuPosY > (-116.5+0.97) && MCFlux_NuPosY < (116.5+0.97) && MCFlux_NuPosZ > 0.0+0.1 && MCFlux_NuPosZ < 1036.8+0.1)\
     *( (run_number >= 4952 && run_number <= 7770)*0.943100 + ( (run_number >= 8317 && run_number <=  13696) || (run_number >= 13697 && run_number <= 14116) || (run_number >= 14117 && run_number <= 18960) )*1.020139 )\
    )";

    // Bit that's unique to 2g1p
    std::string additional_weight_2g1p = additional_weight_common + "*(Sum$(mctruth_exiting_proton_energy-0.93827 > 0.05)==1)*6.7873e20";

    // Bit that's unique to 2g0p
    std::string additional_weight_2g0p = additional_weight_common + "*(Sum$(mctruth_exiting_proton_energy-0.93827 > 0.05)==0)*5.8930e20";

    // Prescription for determining which events 
    // pass selection cuts
    // -----------------------------------------
    TTreeFormula * pass_form_2g1p = new TTreeFormula("pass","(simple_2g1p_NextGen_v4COSMIC_mva>0.894 && simple_2g1p_NextGen_v4BNB_mva >0.737)", v_2g1p);
    TTreeFormula * pass_form_2g0p = new TTreeFormula("pass","(simple_2g0p_NextGen_v4COSMIC_mva>0.944 && simple_2g0p_NextGen_v4BNB_mva >0.731)", v_2g0p);

    // Prescription for normalizing selection 
    // and implementing some cuts
    // -----------------------------------------
    TTreeFormula * norm_form_2g1p = new TTreeFormula("norm",(additional_weight_2g1p).c_str(), v_2g1p);
    TTreeFormula * norm_form_2g0p = new TTreeFormula("norm",(additional_weight_2g0p).c_str(), v_2g0p);

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

    // Number of true and misids for each reco id type. Unknown true ids are when norm_form_2g1p and norm_form_2g0p both evaluate (as w_2g1p and w_2g0p) to 0.
    int N_true_2g1p = 0;
    int N_false_2g1p = 0;
    int N_unknown_2g1p = 0;
    int N_true_2g0p = 0;
    int N_false_2g0p = 0;
    int N_unknown_2g0p = 0;

    // Loop through event tree; fill htrue, hreco, and response hists
    for(int i=0; i<v_2g1p->GetEntries();i++){

        v_2g1p->GetEntry(i);
        v_2g0p->GetEntry(i);

        // Technically not necessary, but Mark says that this avoids some inconsistent obscure ROOT bug
        reco_form_2g1p->GetNdata();
        true_form_2g1p->GetNdata();
        pass_form_2g1p->GetNdata();
        norm_form_2g1p->GetNdata();

        reco_form_2g0p->GetNdata();
        true_form_2g0p->GetNdata();
        pass_form_2g0p->GetNdata();
        norm_form_2g0p->GetNdata();

        double r_2g1p = reco_form_2g1p->EvalInstance();
        double t_2g1p = true_form_2g1p->EvalInstance();
        double p_2g1p = pass_form_2g1p->EvalInstance();
        double w_2g1p = norm_form_2g1p->EvalInstance();

        double r_2g0p = reco_form_2g0p->EvalInstance();
        double t_2g0p = true_form_2g0p->EvalInstance();
        double p_2g0p = pass_form_2g0p->EvalInstance();
        double w_2g0p = norm_form_2g0p->EvalInstance();

        if(i%20000==0)std::cout<<i<<" "<<v_2g1p->GetEntries()<<std::endl;
    //    std::cout<<r<<"\t\t"<<t<<"\t\t"<<p<<std::endl;

        if(p_2g1p){
            resp_2g1p->Fill(r_2g1p,t_2g1p,w_2g1p);
            hreco_2g1p->Fill(r_2g1p,w_2g1p);
	    if(w_2g1p != 0){
		    if(w_2g0p != 0){
			    printf("Both w_2g1p and w_2g0p are nonzero! Exiting.\n");
			    exit(0);
		    }
		    N_true_2g1p++;
	    }
	    else if(w_2g0p != 0)
		    N_false_2g1p++;
	    else{
		    N_unknown_2g1p++;
	    }
        }
        else if(p_2g0p){
            resp_2g0p->Fill(r_2g0p,t_2g0p,w_2g0p);
            hreco_2g0p->Fill(r_2g0p,w_2g0p);
	    if(w_2g0p != 0){
		    if(w_2g1p != 0){
			    printf("Both w_2g0p and w_2g1p are nonzero! Exiting.\n");
			    exit(0);
		    }
		    N_true_2g0p++;
	    }
	    else if(w_2g1p != 0)
		    N_false_2g0p++;
	    else
		    N_unknown_2g0p++;
        }

	htrue_2g1p->Fill(t_2g1p, w_2g1p);
	htrue_2g0p->Fill(t_2g0p, w_2g0p);

    }
    std::cout << "Finished main loop; constructed htrue, hreco, and response hists for both 2g1p and 2g0p signal definitions." << std::endl;
    printf("Correct 2g1p Identifications: %d\n", N_true_2g1p);
    printf("False 2g1p Indentifications: %d\n", N_false_2g1p);
    printf("Unknown 2g1p Identifications: %d\n", N_unknown_2g1p);
    printf("Correct 2g0p Identifications: %d\n", N_true_2g0p);
    printf("False 2g0p Identifications: %d\n", N_false_2g0p);
    printf("Unknown 2g0p Identifiactions: %d\n", N_unknown_2g0p);

    for(int i=0;i< hreco_2g1p->GetNbinsX()+2; i++){ // Both loops need to go to nbins+2 to capture the underflow and overflow
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
    mat_2g1p.Write("response_matrix_2g1p_exclusive");
    mat_2g0p.Write("response_matrix_2g0p_exclusive");
    
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
