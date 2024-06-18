
// Method to get gLEE tree loaded for proper looping-over
// Originally from /exp/uboone/app/users/markrl/useful_scripts/plothelper.h
TTree* loadgLEE(std::string filename, std::string int_dir){

    TFile*f = new TFile((filename).c_str(),"read");
    TTree*v = (TTree*)f->Get((int_dir+"/vertex_tree").c_str());
    TTree*s = (TTree*)f->Get((int_dir+"/simple_tree").c_str());
    TTree*e = (TTree*)f->Get((int_dir+"/eventweight_tree").c_str());
    TTree*t = (TTree*)f->Get((int_dir+"/true_eventweight_tree").c_str());
    v->AddFriend(s);
    v->AddFriend(e);
    v->AddFriend(t);

    return v;
}

// Method that creates response matrix and updated cross section systematic universes to be used in 1D NCpi0 xsec extraction
void ResponseMaker(std::string outDir){

    // -----------------------------------------------------
    // Interface with gLEE tuple 
    // -----------------------------------------------------

    TTree *v_2g1p = (TTree*)loadgLEE("/mnt/morrigan/NCPi0_XS_data/sbnfit_2g1p_NextGen_v4_stage_-1_ext_Denom_NCPi0_CutFromBNB_Run123_v50.5.root","singlephoton");
    TTree *v_2g0p = (TTree*)loadgLEE("/mnt/morrigan/NCPi0_XS_data/sbnfit_2g0p_NextGen_v4_stage_-1_ext_Denom_NCPi0_CutFromBNB_Run123_v50.5.root","singlephoton");
 
    // Copy variation spectra files from SBNfit to own directory in order to update cross section universes to only account for differences in response.
    TFile *spectrain_2g1p = new TFile("/mnt/morrigan/NCPi0_XS_data/SBNfit_variation_spectra_exclusive_2g1p.root", "read");
    TFile *spectrain_2g0p = new TFile("/mnt/morrigan/NCPi0_XS_data/SBNfit_variation_spectra_exclusive_2g0p.root", "read");
    std::string spectraout_2g1p_path = "/app/users/crbergner/data/variation_spectra/SBNfit_variation_spectra_exclusive_2g1p.root";
    std::string spectraout_2g0p_path = "/app/users/crbergner/data/variation_spectra/SBNfit_variation_spectra_exclusive_2g0p.root";
    spectrain_2g1p -> Cp(spectraout_2g1p_path.c_str());
    spectrain_2g0p -> Cp(spectraout_2g0p_path.c_str());
    spectrain_2g1p -> Close();
    spectrain_2g0p -> Close();
    TFile *spectraout_2g1p = new TFile(spectraout_2g1p_path.c_str(), "UPDATE");
    TFile *spectraout_2g0p = new TFile(spectraout_2g0p_path.c_str(), "UPDATE");

    std::cout<<"Number of entries in 2g1p tuple: "<<v_2g1p->GetEntries()<<std::endl;
    std::cout<<"Number of entries in 2g0p tuple: "<<v_2g0p->GetEntries()<<std::endl;

    if (v_2g1p->GetEntries()!=v_2g0p->GetEntries()) {
        std::cout << "The number of entries in the 2g1p and 2g0p input tuples are unequal and parts of this script won't work correctly, so I'm exiting..." << std::endl;
        exit(1); // Terminate with a non-zero exit status
    }

    // Generate dictionaries
    gInterpreter->GenerateDictionary("std::map<std::string,std::vector<double>>", "vector;map;string");
    gInterpreter->GenerateDictionary("std::map<std::string, int>", "map;string");

    // Prescription for calculating reconstructed pion momentum
    // In hive framework, this variable is aliased as "reco_pion_momentum"
    // Pmom_pi = sqrt(Ppix^2+Ppiy^2+Ppiz^2);
    // Ppix = Pgamma1_x+Pgamma2_x;
    // Pgamma_x = energy*dir_x;
    // implied_dir is based on reco vertex and is used when a hadronic track is present (2g1p)
    std::string reco_var_2g1p =
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
    
    std::string reco_var_2g0p =
        "sqrt(\
        ( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_dirx[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_dirx[(i_shr[1])] )\
        *( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_dirx[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_dirx[(i_shr[1])] )\
        +\
        ( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_diry[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_diry[(i_shr[1])] )\
        *( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_diry[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_diry[(i_shr[1])] )\
        +\
        ( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_dirz[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_dirz[(i_shr[1])] )\
        *( (1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486) * reco_shower_dirz[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486) * reco_shower_dirz[(i_shr[1])] )\
        )/1000";

    // Prescription for calculating truth pion momentum (this is already stored in a branch)
    std::string true_var = "mctruth_exiting_pi0_mom";

    // Implement kinematic variable definitions as TTreeFormula for each TTree.
    TTreeFormula * reco_form_2g1p = new TTreeFormula("recof",reco_var_2g1p.c_str(), v_2g1p);
    TTreeFormula * reco_form_2g0p = new TTreeFormula("recof",reco_var_2g0p.c_str(), v_2g0p);
    TTreeFormula * true_form_2g1p = new TTreeFormula("truef",true_var.c_str(), v_2g1p);
    TTreeFormula * true_form_2g0p = new TTreeFormula("truef",true_var.c_str(), v_2g0p);

    // Prescription for calculating event weight
    // -----------------------------------------

    // Bit that's common to 2g1p and 2g0p
    // 1. Scaling simulation POT to data POT. Simulation POT in denominator here. Data POT is in additional_weight.
    // 2. Remove events outside detector region.
    // 3. Scale Monte Carlo runs to number of events captured in actual data runs.
    // 4. Remove events with pi+/-
    // 5. Remove events with exiting photons
    std::string additional_weight_common =
        "(\
        (simple_pot_weight*6.868e20/4.9669582e+21)\
        *(MCFlux_NuPosX > (0.0-1.55) && MCFlux_NuPosX < (256.35-1.55) && MCFlux_NuPosY > (-116.5+0.97) && MCFlux_NuPosY < (116.5+0.97) && MCFlux_NuPosZ > 0.0+0.1 && MCFlux_NuPosZ < 1036.8+0.1)\
        *( (run_number >= 4952 && run_number <= 7770)*0.943100 + ( (run_number >= 8317 && run_number <=  13696) || (run_number >= 13697 && run_number <= 14116) || (run_number >= 14117 && run_number <= 18960) )*1.020139 )\
        *(mctruth_num_exiting_pipm == 0)\
        *(mctruth_num_exiting_photons == 0)\
        )";

    // Bit that's unique to 2g1p
    std::string additional_weight_2g1p = additional_weight_common + "*(Sum$(mctruth_exiting_proton_energy-0.93827 > 0.05)==1)";

    // Bit that's unique to 2g0p
    std::string additional_weight_2g0p = additional_weight_common + "*(Sum$(mctruth_exiting_proton_energy-0.93827 > 0.05)==0)";

    // Prescription for determining which events 
    // pass selection cuts
    // -----------------------------------------
    // 1. Only select events with 20 photons during beam window and fewer than 20 photons in veto time window.
    // 2. Combined topological, preselection, and BDT selection cuts.
    TTreeFormula * pass_form_2g1p = new TTreeFormula("pass","(m_flash_optfltr_pe_beam >20 && m_flash_optfltr_pe_veto < 20) && (simple_2g1p_NextGen_v4COSMIC_mva>0.894 && simple_2g1p_NextGen_v4BNB_mva >0.737)", v_2g1p);
    TTreeFormula * pass_form_2g0p = new TTreeFormula("pass","(m_flash_optfltr_pe_beam >20 && m_flash_optfltr_pe_veto < 20) && (simple_2g0p_NextGen_v4COSMIC_mva>0.944 && simple_2g0p_NextGen_v4BNB_mva >0.731)", v_2g0p);

    // Prescription for determining which events
    // pass signal cuts
    // -----------------------------------------
    TTreeFormula * norm_form_2g1p = new TTreeFormula("norm",(additional_weight_2g1p).c_str(), v_2g1p);
    TTreeFormula * norm_form_2g0p = new TTreeFormula("norm",(additional_weight_2g0p).c_str(), v_2g0p);

    // Load maps to universe weights
    // -----------------------------------------
    std::map<std::string, std::vector<double>> *univ_weight_2g1p = 0;
    std::map<std::string, std::vector<double>> *univ_weight_2g0p = 0;
    v_2g1p -> SetBranchAddress("mcweight", &univ_weight_2g1p);
    v_2g0p -> SetBranchAddress("mcweight", &univ_weight_2g0p);

    // -----------------------------------------------------
    // Construct response matrix 
    // -----------------------------------------------------

//this area is under Cricket construction be warned !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // Specify reco and true analysis binning schemes
    // std::vector<double> reco_bins = {0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675,  0.9}; 
    // std::vector<double> true_bins = {0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675,  0.9}; 
    
    // int nbins_reco = reco_bins.size() - 1;
    // int nbins_true = true_bins.size() - 1;

    //fill vectors with getLowBin values instead of manual numbers 3
    //set number of bins by variation spectra 4
    // make histogram 5
    //debug 6
    //call Leon 7


/*1 -> use a pointer to get the number of bins */
    //how to load in the correct file for the numBins variable?
    TH1D *trueHistTest = (TH1D*) spectraout_2g1p -> Get("exclusive_2g1p_CV_Dir/Sys2g1p_numerator_truth_Signal");
    TH1D *recoHistTest = (TH1D*) spectraout_2g1p -> Get("exclusive_2g1p_CV_Dir/Sys2g1p_numerator_reco_Signal");

    int nbins_true = trueHistTest->GetNbinsX() - 1; //(get bins with the objects)
    int nbins_reco = recoHistTest->GetNbinsX() - 1;

/*2 -> use the number of bins to make the size of an empty array */
    std::vector<double> reco_bins(nbins_true + 1);
    std::vector<double> true_bins(nbins_reco + 1);

/*3 -> fill array with lowBin values through for loop */
    for(int i=1; i<(nbins_true+2); i++){ //true
        true_bins.at(i-1) = recoHistTest->GetBinLowEdge(i);
        std::cout << true_bins.at(i-1) << std::endl;
    }
    for(int i=1; i<(nbins_reco+2); i++){ //reco
        reco_bins.at(i-1) = trueHistTest->GetBinLowEdge(i); 
    }

    

//end Cricket construction !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   std::cout << "I am here." << std::endl;
    // Create htrue, hreco, and response hists
   PlotUtils::MnvH1D * htrue_2g1p = new PlotUtils::MnvH1D("htrue_2g1p","htrue_2g1p",nbins_true, &true_bins[0]);
   PlotUtils::MnvH1D * htrue_2g0p = new PlotUtils::MnvH1D("htrue_2g0p","htrue_2g0p",nbins_true, &true_bins[0]);
   PlotUtils::MnvH1D * hreco_2g1p = new PlotUtils::MnvH1D("hreco_2g1p","hreco_2g1p",nbins_reco, &reco_bins[0]);
   PlotUtils::MnvH1D * hreco_2g0p = new PlotUtils::MnvH1D("hreco_2g0p","hreco_2g0p",nbins_reco, &reco_bins[0]);
   PlotUtils::MnvH2D * resp_2g1p = new PlotUtils::MnvH2D("Response_2g1p","Response_2g1p",nbins_reco, &reco_bins[0], nbins_true, &true_bins[0]);
   PlotUtils::MnvH2D * resp_2g0p = new PlotUtils::MnvH2D("Response_2g0p","Response_2g0p",nbins_reco, &reco_bins[0], nbins_true, &true_bins[0]);

    // Create vertical error bands for cross section systematic universes that need to be updated.
    // Vertical error bands store all the universes of a given category.
    int nMultiverses = 1000;
    int nMinMaxUniverses = 2;
    std::map<std::string, int> XS_SYSTS = {{"All_UBGenie", nMultiverses}, {"AxFFCCQEshape_UBGenie", nMinMaxUniverses}, {"DecayAngMEC_UBGenie", nMinMaxUniverses}, {"NormCCCOH_UBGenie", nMinMaxUniverses}, {"NormNCCOH_UBGenie", nMinMaxUniverses}, {"RPA_CCQE_UBGenie", nMinMaxUniverses}, {"Theta_Delta2Npi_UBGenie", nMinMaxUniverses}, {"VecFFCCQEshape_UBGenie", nMinMaxUniverses}};
    for(auto syst : XS_SYSTS){
        htrue_2g1p -> AddVertErrorBand(syst.first, syst.second);
        htrue_2g0p -> AddVertErrorBand(syst.first, syst.second);
        hreco_2g1p -> AddVertErrorBand(syst.first, syst.second);
        hreco_2g0p -> AddVertErrorBand(syst.first, syst.second);
        resp_2g1p -> AddVertErrorBand(syst.first, syst.second);
        resp_2g0p -> AddVertErrorBand(syst.first, syst.second);
    }

    // Create response TMatrix (large enough to include underflow/underflow from corresponding TH2D)
    TMatrixD mat_2g1p(nbins_reco + 2, nbins_true + 2);
    mat_2g1p.Zero();
    TMatrixD mat_2g0p(nbins_reco + 2, nbins_true + 2);
    mat_2g0p.Zero();

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

        // Evaluate TTreeFormulas for each event.
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

        // Fill 2D response histogram and reco histogram if event passes selection cut with weight determined by additional_weight (truth level signal weight).
        // Background events that pass selection cut will be filled with weight of 0.
        if(p_2g1p){
            resp_2g1p->Fill(t_2g1p, r_2g1p, w_2g1p);
            hreco_2g1p->Fill(r_2g1p, w_2g1p);
        }
        else if(p_2g0p){
            resp_2g0p->Fill(t_2g0p,r_2g0p, w_2g0p);
            hreco_2g0p->Fill(r_2g0p, w_2g0p);
        }
            
        // Always fill true histograms with true momentum value and signal weight (0 if not signal).
        htrue_2g1p->Fill(t_2g1p, w_2g1p);
        htrue_2g0p->Fill(t_2g0p, w_2g0p);

        // Fill vertical error bands with same events using corresponding universe weights. This gives systematic uncertainties through reweighting.
        // The MnvVertErrorBand code seems like it should fill the CV universe histogram as well, but that is not the case through testing,
        // so the CV histograms are filled separately in the preceding code.
        for(auto syst : XS_SYSTS){
            if(p_2g1p){
                resp_2g1p->FillVertErrorBand(syst.first, t_2g1p, r_2g1p, (*univ_weight_2g1p)[syst.first], w_2g1p);
                hreco_2g1p->FillVertErrorBand(syst.first, r_2g1p, (*univ_weight_2g1p)[syst.first], w_2g1p);
            }
            else if(p_2g0p){
                resp_2g0p->FillVertErrorBand(syst.first, t_2g0p,r_2g0p, (*univ_weight_2g0p)[syst.first], w_2g0p);
                hreco_2g0p->FillVertErrorBand(syst.first, r_2g0p, (*univ_weight_2g0p)[syst.first], w_2g0p);
            }    
            htrue_2g1p->FillVertErrorBand(syst.first, t_2g1p, (*univ_weight_2g1p)[syst.first], w_2g1p);
            htrue_2g0p->FillVertErrorBand(syst.first, t_2g0p, (*univ_weight_2g0p)[syst.first], w_2g0p);
    }}
    std::cout << "Finished main loop; constructed htrue, hreco, and response hists for both 2g1p and 2g0p signal definitions." << std::endl;
    
    // Divide each response bin by corresponding true bin. Only applies to CV hists.
    for(int i = 0; i < nbins_reco + 2; i++){ // Both loops need to go to nbins+2 to capture the underflow and overflow
        for(int a = 0; a < nbins_true + 2; a++){ // It doesn't matter whether we loop over 2g1p or 2g0p to get the bin counts
            double eff_2g1p = resp_2g1p->GetBinContent(a, i)/htrue_2g1p->GetBinContent(a);
            if(isnan(eff_2g1p))
                eff_2g1p = 0;
            resp_2g1p -> SetBinContent(a, i, eff_2g1p);
            double eff_2g0p = resp_2g0p->GetBinContent(a, i)/htrue_2g0p->GetBinContent(a);
            if(isnan(eff_2g0p))
                eff_2g0p = 0;
            resp_2g0p -> SetBinContent(a, i, eff_2g0p);
        }}

    // ---------------------------------
    // Update Cross Section Universes
    // ---------------------------------

    // Pull out CV truth distributions from SBNfit variation spectra files. The updated truth universes should be constant, removing uncertainties on the true cross section, which is what is being calculated.
    // To calculate updated reco universes, multiply CV truth by universe response matrices. The updated universes will only account for effects of cross section systematics on efficiency and smearing of the CV truth distribution.
    TH1D *htrue_resp_2g1p = (TH1D*) spectraout_2g1p -> Get("exclusive_2g1p_CV_Dir/Sys2g1p_denominator_truth_Signal");
    TH1D *htrue_resp_2g0p = (TH1D*) spectraout_2g0p -> Get("exclusive_2g0p_CV_Dir/Sys2g0p_denominator_truth_Signal");

    // Loop through vertical error band
    for(std::string vert_error_band_name : resp_2g1p -> GetVertErrorBandNames()){
        // Load directory containing the corresponding systematic universes in the new variation spectra files.
        TDirectory *systcat_2g1p_dir = spectraout_2g1p -> GetDirectory(("exclusive_2g1p_" + vert_error_band_name + "_Dir").c_str());
        TDirectory *systcat_2g0p_dir = spectraout_2g0p -> GetDirectory(("exclusive_2g0p_" + vert_error_band_name + "_Dir").c_str());
        // Loop through universes within vertical error band. All histograms should have the same number of universes in each error band, so using resp_2g1p is arbitrary.
        for(int histnum = 0; histnum < resp_2g1p -> GetVertErrorBand(vert_error_band_name) -> GetNHists(); histnum++){
            // Pull out the corresponding universe's histograms.
            TH2D *resp_2g1p_hist = resp_2g1p -> GetVertErrorBand(vert_error_band_name) -> GetHist(histnum);
            TH1D *true_2g1p_hist = htrue_2g1p -> GetVertErrorBand(vert_error_band_name) -> GetHist(histnum);
            TH2D *resp_2g0p_hist = resp_2g0p -> GetVertErrorBand(vert_error_band_name) -> GetHist(histnum);
            TH1D *true_2g0p_hist = htrue_2g0p -> GetVertErrorBand(vert_error_band_name) -> GetHist(histnum);
            // Divide each response bin by corresponding truth bin.
            for(int i = 0; i < nbins_reco + 2; i++){ // Both loops need to go to nbins+2 to capture the underflow and overflow
                for(int a = 0; a < nbins_true + 2; a++){ // It doesn't matter whether we loop over 2g1p or 2g0p to get the bin counts
                    double eff_2g1p = resp_2g1p_hist->GetBinContent(a, i)/true_2g1p_hist->GetBinContent(a);
                    if(isnan(eff_2g1p))
                        eff_2g1p = 0;
                    resp_2g1p_hist -> SetBinContent(a, i, eff_2g1p);
                    double eff_2g0p = resp_2g0p_hist->GetBinContent(a, i)/true_2g0p_hist->GetBinContent(a);
                    if(isnan(eff_2g0p))
                        eff_2g0p = 0;
                    resp_2g0p_hist -> SetBinContent(a, i, eff_2g0p);
                }}
            // Convert histograms to matrices and multiply to obtain new reco universes.
            // htrue_(2g1p/2g0p)_resp -> GetArray() has one more element than true_(2g1p/2g0p)_mat (an extra overflow bin), but it will simply be ignored, which is the desired behavior.
            TMatrixD *resp_2g1p_mat = new TMatrixD(nbins_reco + 2, nbins_true + 2, resp_2g1p_hist -> GetArray());
            TMatrixD *resp_2g0p_mat = new TMatrixD(nbins_reco + 2, nbins_true + 2, resp_2g0p_hist -> GetArray());
            TMatrixD *true_2g1p_mat = new TMatrixD(nbins_true + 2, 1, htrue_resp_2g1p -> GetArray());
            TMatrixD *true_2g0p_mat = new TMatrixD(nbins_true + 2, 1, htrue_resp_2g0p -> GetArray());
            TMatrixD *reco_resp_2g1p_mat = new TMatrixD(*resp_2g1p_mat, TMatrixD::kMult, *true_2g1p_mat);
            TMatrixD *reco_resp_2g0p_mat = new TMatrixD(*resp_2g0p_mat, TMatrixD::kMult, *true_2g0p_mat);
            // Prescribe names of histograms to be replaced.
            std::string hreco_resp_2g1p_name;
            std::string hreco_resp_2g0p_name;
            std::string htrue_resp_2g1p_name;
            std::string htrue_resp_2g0p_name;
            std::string heffNum_truth_resp_2g1p_name;
            std::string heffNum_truth_resp_2g0p_name;
            if(vert_error_band_name == "All_UBGenie" || vert_error_band_name == "RPA_CCQE_UBGenie"){
                hreco_resp_2g1p_name = "Sys2g1p_numerator_reco_universe_" + std::to_string(histnum + 1) + "_Signal";
                hreco_resp_2g0p_name = "Sys2g0p_numerator_reco_universe_" + std::to_string(histnum + 1) + "_Signal";
                htrue_resp_2g1p_name = "Sys2g1p_denominator_truth_universe_" + std::to_string(histnum + 1) + "_Signal";
                htrue_resp_2g0p_name = "Sys2g0p_denominator_truth_universe_" + std::to_string(histnum + 1) + "_Signal";
                heffNum_truth_resp_2g1p_name = "Sys2g1p_numerator_truth_universe_" + std::to_string(histnum + 1) + "_Signal";
                heffNum_truth_resp_2g0p_name = "Sys2g0p_numerator_truth_universe_" + std::to_string(histnum + 1) + "_Signal";
            }
            else{
                hreco_resp_2g1p_name = "Sys2g1p_numerator_reco_minmax_" + std::to_string(histnum + 1) + "_Signal";
                hreco_resp_2g0p_name = "Sys2g0p_numerator_reco_minmax_" + std::to_string(histnum + 1) + "_Signal";
                htrue_resp_2g1p_name = "Sys2g1p_denominator_truth_minmax_" + std::to_string(histnum + 1) + "_Signal";
                htrue_resp_2g0p_name = "Sys2g0p_denominator_truth_minmax_" + std::to_string(histnum + 1) + "_Signal";
                heffNum_truth_resp_2g1p_name = "Sys2g1p_numerator_truth_minmax_" + std::to_string(histnum + 1) + "_Signal";
                heffNum_truth_resp_2g0p_name = "Sys2g0p_numerator_truth_minmax_" + std::to_string(histnum + 1) + "_Signal";
            }
            // Convert new reco distributions back to histograms for writing.
            TH1D *hreco_resp_2g1p = new TH1D(hreco_resp_2g1p_name.c_str(), hreco_resp_2g1p_name.c_str(), nbins_reco, &reco_bins[0]);
            TH1D *hreco_resp_2g0p = new TH1D(hreco_resp_2g0p_name.c_str(), hreco_resp_2g0p_name.c_str(), nbins_reco, &reco_bins[0]);
            for(int bin = 0; bin < nbins_reco + 2; bin++){
                hreco_resp_2g1p -> SetBinContent(bin, (*reco_resp_2g1p_mat)(bin, 0));
                hreco_resp_2g0p -> SetBinContent(bin, (*reco_resp_2g0p_mat)(bin, 0));
            }
            // Get the new efficiency numerator truth distributions for use in calculating the efficiency.
            TH1D *heffNum_truth_resp_2g1p = resp_2g1p_hist -> ProjectionX(heffNum_truth_resp_2g1p_name.c_str());
            TH1D *heffNum_truth_resp_2g0p = resp_2g0p_hist -> ProjectionX(heffNum_truth_resp_2g0p_name.c_str());
            for(int bin = 0; bin < nbins_true + 2; bin++){
                heffNum_truth_resp_2g1p -> SetBinContent(bin, heffNum_truth_resp_2g1p -> GetBinContent(bin)*htrue_resp_2g1p -> GetBinContent(bin));
                heffNum_truth_resp_2g0p -> SetBinContent(bin, heffNum_truth_resp_2g0p -> GetBinContent(bin)*htrue_resp_2g0p -> GetBinContent(bin));
            }
            // Cd to the correct directory and overwrite existing distributions.
            systcat_2g1p_dir -> cd();
            hreco_resp_2g1p -> Write(hreco_resp_2g1p_name.c_str(), TObject::kWriteDelete);
            htrue_resp_2g1p -> Write(htrue_resp_2g1p_name.c_str(), TObject::kWriteDelete);
            heffNum_truth_resp_2g1p -> Write(heffNum_truth_resp_2g1p_name.c_str(), TObject::kWriteDelete);
            systcat_2g0p_dir -> cd();
            hreco_resp_2g0p -> Write(hreco_resp_2g0p_name.c_str(), TObject::kWriteDelete);
            htrue_resp_2g0p -> Write(htrue_resp_2g0p_name.c_str(), TObject::kWriteDelete);
            heffNum_truth_resp_2g0p -> Write(heffNum_truth_resp_2g0p_name.c_str(), TObject::kWriteDelete);
        }
        systcat_2g1p_dir -> Write();
        systcat_2g0p_dir -> Write();
    }
    

    std::cout << "Finished dividing." << std::endl;

    // Fill response matrix with the same values.
    for(int i=0;i< hreco_2g1p->GetNbinsX()+2; i++){ // Both loops need to go to nbins+2 to capture the underflow and overflow
        for(int a=0; a<htrue_2g1p->GetNbinsX()+2; a++){ // It doesn't matter whether we loop over 2g1p or 2g0p to get the bin counts
            mat_2g1p(i,a) = resp_2g1p->GetBinContent(a,i); // NOTE: TMatrix uses (row, column) for indices, so while (i, a) corresponds to (xbin, ybin) in resp_2g1p, it is reversed in mat_2g1p.
            mat_2g0p(i,a) = resp_2g0p->GetBinContent(a,i); // Also, matrices have their "origin" at the top left, while 2D histograms have theirs at the bottom left.
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
    spectraout_2g1p -> Close();
    spectraout_2g0p -> Close();
    return ;

}
