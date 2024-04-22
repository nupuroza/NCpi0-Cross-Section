#include "plothelper.h"

// Returns TH1 of var with specified cuts and binning using the TTree::Draw method
// tin is the input TTree
// var is the definition of the x variable in the 1D histogram being created
// cuts is the definition of event weight (0 of if definitions not met; simple_pot_weight otherwise)
// nam is internal name of TH1 object. Just make sure they don't overlap one another.
TH1 * getTH1(TTree * tin, std::string var, std::string cuts, std::string nam){

    tin->Draw((var+">>+"+nam).c_str() , ("("+cuts+ ")").c_str(),"goff");
    TH1* th1 = (TH1*)gDirectory->Get(nam.c_str());
    th1->SetLineWidth(1);
    th1->SetStats(0);
    th1->SetDirectory(0);

    return th1;
}

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

// Breaks down selected GENIE and (to be implemented) NuWro events into signal and background components.
// Default out_dir is for Leon's own convenience. Feel free to change it locally.
void SigBkgdBreakdown(std::string out_dir = "/exp/uboone/data/users/ltong/gLEE/NCPi0/GENIEvsNuWro/"){
    gROOT -> SetBatch(kTRUE);

    // Load input files
    std::string filename_spectra = "/exp/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/MajorMerge_GGE_mark/working_dir/ToTH1D/VersionNextGen_SamePOT_Aug2023/variation_spectra/SBNfit_variation_spectra_exclusive_2g1p.root";
    std::string filename_2g1p = "/pnfs/uboone/persistent/users/markross/Jan2022_gLEE_files/NCPi0CrossSection/UpdatedPOT/2g1p_v4/sbnfit_2g1p_NextGen_v4_stage_3_AllMC.root";
    std::string filename_2g0p = "/pnfs/uboone/persistent/users/markross/Jan2022_gLEE_files/NCPi0CrossSection/UpdatedPOT/2g0p_v4/sbnfit_2g0p_NextGen_v4_stage_3_AllMC.root";
    std::string cosmics_2g1p = "/pnfs/uboone/persistent/users/markross/Jan2022_gLEE_files/NCPi0CrossSection/UpdatedPOT/2g1p_v4/sbnfit_2g1p_NextGen_v4_stage_3_BNBext.root";
    std::string cosmics_2g0p = "/pnfs/uboone/persistent/users/markross/Jan2022_gLEE_files/NCPi0CrossSection/UpdatedPOT/2g0p_v4/sbnfit_2g0p_NextGen_v4_stage_3_BNBext.root";

    // Load TH1 from variation spectra file to be used for binning
    TFile * fsp = new TFile(filename_spectra.c_str(), "read");
    TH1 * hsp = (TH1*) fsp -> Get("exclusive_2g1p_CV_Dir/Sys2g1p_numerator_reco_Signal");
    TH1 * href = (TH1*) hsp -> Clone("href");
    href -> SetTitle("");
    int nBins = href -> GetNbinsX();
    for(int i = 1; i < nBins + 1; i++){
        href -> SetBinContent(i,0.);
        href -> SetBinError(i,0.);
    }
    // Contains signal and background Monte Carlo chosen by our BDT
    TTree * tmc_2g1p = (TTree*)loadgLEE(filename_2g1p.c_str(), "singlephoton");
    TTree * tmc_2g0p = (TTree*)loadgLEE(filename_2g0p.c_str(), "singlephoton");
    // Contains cosmic data
    TTree * tcos_2g1p = (TTree*)loadgLEE(cosmics_2g1p.c_str(), "singlephoton");
    TTree * tcos_2g0p = (TTree*)loadgLEE(cosmics_2g0p.c_str(), "singlephoton");

    // Reco momentum definition
    std::string variable_2g1p = "sqrt(((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_implied_dirx[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_implied_dirx[(i_shr[1])])*((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_implied_dirx[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_implied_dirx[(i_shr[1])]) + ((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_implied_diry[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_implied_diry[(i_shr[1])])*((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_implied_diry[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_implied_diry[(i_shr[1])]) + ((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_implied_dirz[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_implied_dirz[(i_shr[1])])*((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_implied_dirz[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_implied_dirz[(i_shr[1])]))/1000";
    std::string variable_2g0p = "sqrt(((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_dirx[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_dirx[(i_shr[1])])*((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_dirx[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_dirx[(i_shr[1])]) + ((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_diry[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_diry[(i_shr[1])])*((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_diry[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_diry[(i_shr[1])]) + ((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_dirz[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_dirz[(i_shr[1])])*((1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486)*reco_shower_dirz[i_shr[0]] + (1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486)*reco_shower_dirz[(i_shr[1])]))/1000";
    std::string x_label = "Reco #pi^{0} momentum [GeV]";

    // Names of signal categories (used in legend) and their definitions
    std::vector<std::string> nams_sig = {"Quasi-Elastic", "Deep Inelastic", "Resonant #Delta(1232)", "Other Resonant", "Coherent", "Other"};
    std::vector<std::string> cutdef_sig = {
        "GTruth_Gscatter == 1",
        "GTruth_Gscatter == 3",
        "GTruth_Gscatter == 4 && GTruth_ResNum == 0",
        "GTruth_Gscatter == 4 && GTruth_ResNum != 0",
        "GTruth_Gscatter == 5"
    };
    // Push back "Other" signal that is any signal that doesn't meet the definition of any of the categories listed.
    std::string other_sig= "1";
    for(auto &s: cutdef_sig){ other_sig+="&& !("+s+")"; }
    cutdef_sig.push_back(other_sig);

    // Definition of 2g1p and 2g0p signal. Fiducial volume listed separately, as it is also used in defining background categories.
    std::string fiducial_volume = "(MCFlux_NuPosX > (0.0-1.55) && MCFlux_NuPosX < (256.35-1.55) && MCFlux_NuPosY > (-116.5+0.97) && MCFlux_NuPosY < (116.5+0.97) && MCFlux_NuPosZ > 0.0+0.1 && MCFlux_NuPosZ < 1036.8+0.1)";
    std::string sig_def =  "(mctruth_cc_or_nc==1 && mctruth_num_exiting_pi0==1 && mctruth_num_exiting_pipm==0 && mctruth_num_exiting_photons == 0 && " + fiducial_volume;
    std::string sig_def_2g1p = sig_def + " && (Sum$(mctruth_exiting_proton_energy-0.93827 > 0.05)==1))";
    std::string sig_def_2g0p = sig_def + " && (Sum$(mctruth_exiting_proton_energy-0.93827 > 0.05)==0))";	

    // Color used for each signal category
    std::vector<int> cols_sig = {kMagenta, kOrange - 7, kCyan + 1, kGreen, kOrange, kGray};

    // Names of background categories (used in legend) and their definitions. Must be careful to avoid overlapping background definitions.
    // Our fiducial volume includes the entire active volume of the TPC, so events outside the fiducial volume are labeled "Out-of-TPC". All events outside the fiducial volume are classified as such regardless of their interaction type.
    // A few NC delta radiative events also result in pi+/- production. These are negligible but excluded from the pi+/- backgrounds.
    // 2g0p signal misids are events matching the 2g0p signal definition in the selected 2g1p data and vice versa.
    std::vector<std::string> nams_bkgd = {"Out-of-TPC", "CC #nu_{#mu} 1#pi^{0}", "NC #pi^{+}/#pi^{-} 1#pi^{0}", "Other NC #pi^{+}/#pi^{-}", "NC >1#pi^{0}", "NC #eta#rightarrow#gamma#gamma"};
    std::vector<std::string> cutdef_bkgd = {
        "!" + fiducial_volume,
        "mctruth_cc_or_nc == 0 && mctruth_nu_pdg == 14 && mctruth_num_exiting_pi0 == 1 &&" + fiducial_volume,
        "mctruth_cc_or_nc == 1 && mctruth_num_exiting_pipm > 0 && mctruth_num_exiting_pi0 == 1 &&" + fiducial_volume,
        "mctruth_cc_or_nc == 1 && mctruth_num_exiting_pipm > 0 && mctruth_num_exiting_pi0 != 1 &&" + fiducial_volume,
        "mctruth_cc_or_nc == 1 && mctruth_num_exiting_pi0 > 1 && mctruth_num_exiting_pipm == 0 &&" + fiducial_volume,
        "mctruth_cc_or_nc == 1 && Sum$(mctruth_daughters_pdg == 221) > 0 && Sum$(mctruth_daughters_pdg == 22 && mctruth_daughters_status_code == 1) == 2 &&" + fiducial_volume
    };
    std::vector<std::string> nams_bkgd_2g1p = {"NC 1#pi^{0} 2g0p Signal Misid.", "BNB Other"};
    nams_bkgd_2g1p.insert(nams_bkgd_2g1p.begin(), nams_bkgd.begin(), nams_bkgd.end());
    std::vector<std::string> cutdef_bkgd_2g1p = cutdef_bkgd;
    cutdef_bkgd_2g1p.push_back(sig_def_2g0p);
    std::string other_bkgd_2g1p = "1";
    for(auto &b: cutdef_bkgd_2g1p){ other_bkgd_2g1p+="&& !("+b+")"; };
    cutdef_bkgd_2g1p.push_back(other_bkgd_2g1p);
    std::vector<std::string> nams_bkgd_2g0p = {"NC 1#pi^{0} 2g1p Signal Misid.", "BNB Other"};
    nams_bkgd_2g0p.insert(nams_bkgd_2g0p.begin(), nams_bkgd.begin(), nams_bkgd.end());
    std::vector<std::string> cutdef_bkgd_2g0p = cutdef_bkgd;
    cutdef_bkgd_2g0p.push_back(sig_def_2g1p);
    std::string other_bkgd_2g0p = "1";
    for(auto &b: cutdef_bkgd_2g0p){ other_bkgd_2g0p+="&& !("+b+")"; };
    cutdef_bkgd_2g0p.push_back(other_bkgd_2g0p);

    // Definition of 2g1p and 2g0p background. Anything that does not meet their signal definition is defined as background.
    std::string bkgd_def_2g1p = "!" + sig_def_2g1p;
    std::string bkgd_def_2g0p = "!" + sig_def_2g0p;

    // Color used for each background category
    std::vector<int> cols_bkgd = {kOrange - 7, kBlue, kCyan, kMagenta, kYellow - 9, kCyan + 3, kRed, kGray};

    // Names of pi0 parent particle categories (used in legend) and their definitions. Must be careful to avoid overlapping definitions.
    // Events in a category must originate from the corresponding neutrino interaction and have corresponding particles from those interactions selected as the 2 reconstructed pi0 showers.
    // All events outside the fiducial volume are classified as dirt regardless of their interaction type or pi0 parent particles.
    // Cosmic contamination refers to events in a slice with a real neutrino interaction in which one or more of the selected pi0 showers were actually cosmics.
    std::vector<std::string> nams_parent = {"Out-of-TPC", "CC #nu_{#mu} 1#pi^{0}", "NC #pi^{+}/#pi^{-} 1#pi^{0}", "NC >1#pi^{0}", "#pi^{0} from Rescattering", "Cosmic Contamination", "#eta#rightarrow#gamma#gamma", "Multiple Backgrounds"};
    std::vector<std::string> cutdef_parent = {
        "!" + fiducial_volume,
        "mctruth_cc_or_nc == 0 && mctruth_nu_pdg == 14 && mctruth_num_exiting_pi0 == 1 && Sum$(sim_shower_parent_pdg == 111) == 2 &&" + fiducial_volume,
        "mctruth_cc_or_nc == 1 && mctruth_num_exiting_pipm > 0 && mctruth_num_exiting_pi0 == 1 && Sum$(sim_shower_parent_pdg == 111) == 2 &&" + fiducial_volume,
        "mctruth_cc_or_nc == 1 && mctruth_num_exiting_pi0 > 1 && Sum$(sim_shower_parent_pdg == 111) == 2 &&" + fiducial_volume,
        "mctruth_num_exiting_pi0 == 0 && Sum$(sim_shower_parent_pdg == 111) == 2 &&" + fiducial_volume,
        "Sum$(sim_shower_matched == 0) > 1 &&" + fiducial_volume,
        "Sum$(mctruth_daughters_pdg == 221) > 0 && Sum$(mctruth_daughters_pdg == 22 && mctruth_daughters_status_code == 1) == 2 && Sum$(sim_shower_pdg == 22) == 2 && " + fiducial_volume,
        "sim_shower_parent_pdg[0] != sim_shower_parent_pdg[1] &&" + fiducial_volume
    };
    std::vector<std::string> nams_parent_2g1p = {"NC 1#pi^{0} 2g0p Signal Misid.", "BNB Other"};
    nams_parent_2g1p.insert(nams_parent_2g1p.begin(), nams_parent.begin(), nams_parent.end());
    std::vector<std::string> cutdef_parent_2g1p = cutdef_parent;
    cutdef_parent_2g1p.push_back(sig_def_2g0p + "&& Sum$(sim_shower_parent_pdg == 111) == 2");
    std::string other_parent_2g1p = "1";
    for(auto &b: cutdef_parent_2g1p){ other_parent_2g1p+="&& !("+b+")"; };
    cutdef_parent_2g1p.push_back(other_parent_2g1p);
    std::vector<std::string> nams_parent_2g0p = {"NC 1#pi^{0} 2g1p Signal Misid.", "BNB Other"};
    nams_parent_2g0p.insert(nams_parent_2g0p.begin(), nams_parent.begin(), nams_parent.end());
    std::vector<std::string> cutdef_parent_2g0p = cutdef_parent;
    cutdef_parent_2g0p.push_back(sig_def_2g1p + "&& Sum$(sim_shower_parent_pdg == 111) == 2");
    std::string other_parent_2g0p = "1";
    for(auto &b: cutdef_parent_2g0p){ other_parent_2g0p+="&& !("+b+")"; };
    cutdef_parent_2g0p.push_back(other_parent_2g0p);

    // Color used for each parent particle category
    std::vector<int> cols_parent = {kOrange - 7, kBlue, kCyan, kYellow - 9, kViolet, kGreen + 1, kCyan + 3, kOrange, kRed, kGray};

    // Move cosmic contamination background to the end. The name is actually a placeholder, as cosmic contamination will be combined with cosmic data into a single category labeled "cosmics".
    cutdef_parent_2g1p.push_back(cutdef_parent[5]);
    cutdef_parent_2g1p.erase(cutdef_parent_2g1p.begin() + 5);
    cutdef_parent_2g0p.push_back(cutdef_parent[5]);
    cutdef_parent_2g0p.erase(cutdef_parent_2g0p.begin() + 5);
    nams_parent_2g1p.push_back(nams_parent[5]);
    nams_parent_2g1p.erase(nams_parent_2g1p.begin() + 5);
    nams_parent_2g0p.push_back(nams_parent[5]);
    nams_parent_2g0p.erase(nams_parent_2g0p.begin() + 5);
    cols_parent.push_back(cols_parent[5]);
    cols_parent.erase(cols_parent.begin() + 5);

    // Create THStack and legend objects for signal plots
    THStack *stk_2g1p_sig = new THStack();
    THStack *stk_2g0p_sig = new THStack();
    TLegend *l_sig_2g1p = new TLegend(0.11,0.69,0.89,0.89);
    l_sig_2g1p->SetLineWidth(0);
    l_sig_2g1p->SetFillStyle(0);
    l_sig_2g1p->SetLineColor(kWhite);
    l_sig_2g1p->SetNColumns(2);
    TLegend *l_sig_2g0p = new TLegend(0.11,0.69,0.89,0.89);
    l_sig_2g0p->SetLineWidth(0);
    l_sig_2g0p->SetFillStyle(0);
    l_sig_2g0p->SetLineColor(kWhite);
    l_sig_2g0p->SetNColumns(2);

    // Create pt object for public plots.
    TPaveText *pt = new TPaveText(0.5, 0.64, 0.89, 0.69, "NDC");
    pt -> SetBorderSize(0);
    pt -> SetFillColorAlpha(0, 0);
    pt -> AddText("MicroBooNE Simulation In Progress");

    // Loop through signal categories for 2g1p and 2g0p and fill each histogram into corresponding THStack object.
    // Integral of category across spectrum calculated and included in legend.
    for(int i=0; i< nams_sig.size(); i++){
        std::cout<<i<<" cutdef "<<cutdef_sig[i]<<std::endl;
        TH1D* thtmp_2g1p = (TH1D*) href -> Clone(("2g1p_sig_" + std::to_string(i)).c_str());
        TH1D* thtmp_2g0p = (TH1D*) href -> Clone(("2g0p_sig_" + std::to_string(i)).c_str());
        //Multiply signal definition by category definition to get category definition within signal; then scale by POT weight.
        thtmp_2g1p = (TH1D*)getTH1(tmc_2g1p,variable_2g1p, "simple_pot_weight*("+cutdef_sig[i]+")*"+sig_def_2g1p,"2g1p_sig_" + std::to_string(i));
        thtmp_2g0p = (TH1D*)getTH1(tmc_2g0p,variable_2g0p, "simple_pot_weight*("+cutdef_sig[i]+")*"+sig_def_2g0p,"2g0p_sig_" + std::to_string(i));
        double total_2g1p = thtmp_2g1p -> Integral();
        double total_2g0p = thtmp_2g0p -> Integral();
        thtmp_2g1p->SetFillColor(cols_sig[i]);
        thtmp_2g1p->SetLineColor(kBlack);
        thtmp_2g1p->SetLineWidth(1);
        thtmp_2g0p->SetFillColor(cols_sig[i]);
        thtmp_2g0p->SetLineColor(kBlack);
        thtmp_2g0p->SetLineWidth(1);
        stk_2g1p_sig->Add(thtmp_2g1p);
        stk_2g0p_sig->Add(thtmp_2g0p);
        // to_string_prec from plothelper.h
        l_sig_2g1p->AddEntry(thtmp_2g1p,(nams_sig[i] + " (" + to_string_prec(total_2g1p, 1) + ")").c_str(),"f");
        l_sig_2g0p->AddEntry(thtmp_2g0p,(nams_sig[i] + " (" + to_string_prec(total_2g0p, 1) + ")").c_str(),"f");
    }

    // Draw and save the signal figures.
    TCanvas *c = new TCanvas();
    stk_2g1p_sig->Draw("hist");
    stk_2g1p_sig->SetMaximum(stk_2g1p_sig->GetMaximum()*1.4);
    stk_2g1p_sig->GetXaxis()->SetTitle(x_label.c_str());
    stk_2g1p_sig->GetYaxis()->SetTitle("Events");
    stk_2g1p_sig->SetTitle("GENIE 2g1p Exclusive Signal Breakdown");
    l_sig_2g1p->Draw();
    pt -> Draw();
    c -> Print((out_dir + "/2g1p_exclusive_signal_breakdown.png").c_str());
    c -> Clear();
    stk_2g0p_sig->Draw("hist");
    stk_2g0p_sig->SetMaximum(stk_2g0p_sig->GetMaximum()*1.4);
    stk_2g0p_sig->GetXaxis()->SetTitle(x_label.c_str());
    stk_2g0p_sig->GetYaxis()->SetTitle("Events");
    stk_2g0p_sig->SetTitle("GENIE 2g0p Exclusive Signal Breakdwon");
    l_sig_2g0p->Draw();
    pt -> Draw();
    c -> Print((out_dir + "/2g0p_exclusive_signal_breakdown.png").c_str());

    // Repeat for background categories
    THStack *stk_2g1p_bkgd = new THStack();
    THStack *stk_2g0p_bkgd = new THStack();

    TLegend *l_bkgd_2g1p = new TLegend(0.11,0.69,0.89,0.89);
    l_bkgd_2g1p->SetLineWidth(0);
    l_bkgd_2g1p->SetFillStyle(0);
    l_bkgd_2g1p->SetLineColor(kWhite);
    l_bkgd_2g1p->SetNColumns(2);
    
    TLegend *l_bkgd_2g0p = new TLegend(0.11,0.69,0.89,0.89);
    l_bkgd_2g0p->SetLineWidth(0);
    l_bkgd_2g0p->SetFillStyle(0);
    l_bkgd_2g0p->SetLineColor(kWhite);
    l_bkgd_2g0p->SetNColumns(2);
    
    for(int i=0; i< nams_bkgd_2g1p.size(); i++){
        std::cout<<i<<" cutdef "<<cutdef_bkgd_2g1p[i]<<std::endl;
        TH1D* thtmp_2g1p = (TH1D*) href -> Clone(("2g1p_bkgd_" + std::to_string(i)).c_str());
        TH1D* thtmp_2g0p = (TH1D*) href -> Clone(("2g0p_bkgd_" + std::to_string(i)).c_str());
        thtmp_2g1p = (TH1D*)getTH1(tmc_2g1p,variable_2g1p, "simple_pot_weight*("+cutdef_bkgd_2g1p[i]+")*"+bkgd_def_2g1p,"2g1p_bkgd_" + std::to_string(i));
        thtmp_2g0p = (TH1D*)getTH1(tmc_2g0p,variable_2g0p, "simple_pot_weight*("+cutdef_bkgd_2g0p[i]+")*"+bkgd_def_2g0p,"2g0p_bkgd_" + std::to_string(i));
        double total_2g1p = thtmp_2g1p -> Integral();
        double total_2g0p = thtmp_2g0p -> Integral();
        thtmp_2g1p->SetFillColor(cols_bkgd[i]);
        thtmp_2g1p->SetLineColor(kBlack);
        thtmp_2g1p->SetLineWidth(1);
        thtmp_2g0p->SetFillColor(cols_bkgd[i]);
        thtmp_2g0p->SetLineColor(kBlack);
        thtmp_2g0p->SetLineWidth(1);
        stk_2g1p_bkgd->Add(thtmp_2g1p);
        stk_2g0p_bkgd->Add(thtmp_2g0p);
        l_bkgd_2g1p->AddEntry(thtmp_2g1p,(nams_bkgd_2g1p[i] + " (" + to_string_prec(total_2g1p, 1) + ")").c_str(),"f");
        l_bkgd_2g0p->AddEntry(thtmp_2g0p,(nams_bkgd_2g0p[i] + " (" + to_string_prec(total_2g0p, 1) + ")").c_str(),"f");
    }
    
    // All cosmics are backgrounds and get added separately.
    std::cout << nams_bkgd_2g1p.size() << " Cosmics" << std::endl;
    TH1D* thtmp_2g1p = (TH1D*) href -> Clone("2g1p_cosmics");
    TH1D* thtmp_2g0p = (TH1D*) href -> Clone("2g0p_cosmics");
    thtmp_2g1p = (TH1D*)getTH1(tcos_2g1p, variable_2g1p, "simple_pot_weight", "2g1p_cosmics");
    thtmp_2g0p = (TH1D*)getTH1(tcos_2g0p, variable_2g0p, "simple_pot_weight", "2g0p_cosmics");
    double total_2g1p = thtmp_2g1p -> Integral();
    double total_2g0p = thtmp_2g0p -> Integral();
    thtmp_2g1p->SetFillColor(kGreen);
    thtmp_2g1p->SetFillStyle(3244);
    thtmp_2g1p->SetLineColor(kBlack);
    thtmp_2g1p->SetLineWidth(1);
    thtmp_2g0p->SetFillColor(kGreen);
    thtmp_2g0p->SetFillStyle(3244);
    thtmp_2g0p->SetLineColor(kBlack);
    thtmp_2g0p->SetLineWidth(1);
    stk_2g1p_bkgd->Add(thtmp_2g1p);
    stk_2g0p_bkgd->Add(thtmp_2g0p);
    l_bkgd_2g1p->AddEntry(thtmp_2g1p,("Cosmics (" + to_string_prec(total_2g1p, 1) + ")").c_str(),"f");
    l_bkgd_2g0p->AddEntry(thtmp_2g0p,("Cosmics (" + to_string_prec(total_2g0p, 1) + ")").c_str(),"f");
    
    stk_2g1p_bkgd->Draw("hist");
    stk_2g1p_bkgd->SetMaximum(stk_2g1p_bkgd->GetMaximum()*1.4);
    stk_2g1p_bkgd->GetXaxis()->SetTitle(x_label.c_str());
    stk_2g1p_bkgd->GetYaxis()->SetTitle("Events");
    stk_2g1p_bkgd->SetTitle("GENIE 2g1p Exclusive Background Neutrino Interaction");
    l_bkgd_2g1p->Draw();
    pt -> Draw();
    c -> Print((out_dir + "/2g1p_exclusive_background_breakdown.png").c_str());
    c -> Clear();
    stk_2g0p_bkgd->Draw("hist");
    stk_2g0p_bkgd->SetMaximum(stk_2g0p_bkgd->GetMaximum()*1.4);
    stk_2g0p_bkgd->GetXaxis()->SetTitle(x_label.c_str());
    stk_2g0p_bkgd->GetYaxis()->SetTitle("Events");
    stk_2g0p_bkgd->SetTitle("GENIE 2g0p Exclusive Background Neutrino Interaction");
    l_bkgd_2g0p->Draw();
    pt -> Draw();
    c -> Print((out_dir + "/2g0p_exclusive_background_breakdown.png").c_str());

    // Repeat for parent particle categories
    THStack *stk_2g1p_parent = new THStack();
    THStack *stk_2g0p_parent = new THStack();

    TLegend *l_parent_2g1p = new TLegend(0.11,0.69,0.89,0.89);
    l_parent_2g1p->SetLineWidth(0);
    l_parent_2g1p->SetFillStyle(0);
    l_parent_2g1p->SetLineColor(kWhite);
    l_parent_2g1p->SetNColumns(2);
    
    TLegend *l_parent_2g0p = new TLegend(0.11,0.69,0.89,0.89);
    l_parent_2g0p->SetLineWidth(0);
    l_parent_2g0p->SetFillStyle(0);
    l_parent_2g0p->SetLineColor(kWhite);
    l_parent_2g0p->SetNColumns(2);

    TH1D* thcos_2g1p = new TH1D();
    TH1D* thcos_2g0p = new TH1D();
    
    for(int i=0; i< nams_parent_2g1p.size(); i++){
        std::cout<<i<<" cutdef "<<cutdef_parent_2g1p[i]<<std::endl;
        thtmp_2g1p = (TH1D*) href -> Clone(("2g1p_parent_" + std::to_string(i)).c_str());
        thtmp_2g0p = (TH1D*) href -> Clone(("2g0p_parent_" + std::to_string(i)).c_str());
        thtmp_2g1p = (TH1D*)getTH1(tmc_2g1p,variable_2g1p, "simple_pot_weight*("+cutdef_parent_2g1p[i]+")*"+bkgd_def_2g1p,"2g1p_parent_" + std::to_string(i));
        thtmp_2g0p = (TH1D*)getTH1(tmc_2g0p,variable_2g0p, "simple_pot_weight*("+cutdef_parent_2g0p[i]+")*"+bkgd_def_2g0p,"2g0p_parent_" + std::to_string(i));
        double total_2g1p = thtmp_2g1p -> Integral();
        double total_2g0p = thtmp_2g0p -> Integral();
        thtmp_2g1p->SetFillColor(cols_parent[i]);
        thtmp_2g1p->SetLineColor(kBlack);
        thtmp_2g1p->SetLineWidth(1);
        thtmp_2g0p->SetFillColor(cols_parent[i]);
        thtmp_2g0p->SetLineColor(kBlack);
        thtmp_2g0p->SetLineWidth(1);
        // Add most histograms directly to the stack. Cosmic contamination is saved and combined with cosmic data before being added.
        if(i != nams_parent_2g1p.size() - 1){
            stk_2g1p_parent->Add(thtmp_2g1p);
            stk_2g0p_parent->Add(thtmp_2g0p);
            l_parent_2g1p->AddEntry(thtmp_2g1p,(nams_parent_2g1p[i] + " (" + to_string_prec(total_2g1p, 1) + ")").c_str(),"f");
            l_parent_2g0p->AddEntry(thtmp_2g0p,(nams_parent_2g0p[i] + " (" + to_string_prec(total_2g0p, 1) + ")").c_str(),"f");
        }
        else{
            thcos_2g1p = thtmp_2g1p;
            thcos_2g0p = thtmp_2g0p;
        }
    }
    
    // All cosmics are backgrounds and get added separately.
    std::cout << nams_parent_2g1p.size() << " Cosmics" << std::endl;
    thtmp_2g1p = (TH1D*) href -> Clone("2g1p_cosmics");
    thtmp_2g0p = (TH1D*) href -> Clone("2g0p_cosmics");
    thtmp_2g1p = (TH1D*)getTH1(tcos_2g1p, variable_2g1p, "simple_pot_weight", "2g1p_cosmics");
    thtmp_2g0p = (TH1D*)getTH1(tcos_2g0p, variable_2g0p, "simple_pot_weight", "2g0p_cosmics");
    thcos_2g1p -> Add(thtmp_2g1p);
    thcos_2g0p -> Add(thtmp_2g0p);
    total_2g1p = thcos_2g1p -> Integral();
    total_2g0p = thcos_2g0p -> Integral();
    thcos_2g1p->SetFillColor(kGreen);
    thcos_2g1p->SetFillStyle(3244);
    thcos_2g1p->SetLineColor(kBlack);
    thcos_2g1p->SetLineWidth(1);
    thcos_2g0p->SetFillColor(kGreen);
    thcos_2g0p->SetFillStyle(3244);
    thcos_2g0p->SetLineColor(kBlack);
    thcos_2g0p->SetLineWidth(1);
    stk_2g1p_parent->Add(thcos_2g1p);
    stk_2g0p_parent->Add(thcos_2g0p);
    l_parent_2g1p->AddEntry(thcos_2g1p,("Cosmics (" + to_string_prec(total_2g1p, 1) + ")").c_str(),"f");
    l_parent_2g0p->AddEntry(thcos_2g0p,("Cosmics (" + to_string_prec(total_2g0p, 1) + ")").c_str(),"f");
    
    stk_2g1p_parent->Draw("hist");
    stk_2g1p_parent->SetMaximum(stk_2g1p_parent->GetMaximum()*1.4);
    stk_2g1p_parent->GetXaxis()->SetTitle(x_label.c_str());
    stk_2g1p_parent->GetYaxis()->SetTitle("Events");
    stk_2g1p_parent->SetTitle("GENIE 2g1p Exclusive Background Failure Modes");
    l_parent_2g1p->Draw();
    pt -> Draw();
    c -> Print((out_dir + "/2g1p_exclusive_background_reco_showers_parent.png").c_str());
    c -> Clear();
    stk_2g0p_parent->Draw("hist");
    stk_2g0p_parent->SetMaximum(stk_2g0p_parent->GetMaximum()*1.4);
    stk_2g0p_parent->GetXaxis()->SetTitle(x_label.c_str());
    stk_2g0p_parent->GetYaxis()->SetTitle("Events");
    stk_2g0p_parent->SetTitle("GENIE 2g0p Exclusive Background Failure Modes");
    l_parent_2g0p->Draw();
    pt -> Draw();
    c -> Print((out_dir + "/2g0p_exclusive_background_reco_showers_parent.png").c_str());

    // Turn off batch mode
    gROOT -> SetBatch(kFALSE);
}
