#include<TLorentzVector.h>
#include<vector>
#include<iostream>
#include<string>
#include<fstream>
#include<filesystem>
#include<math.h>
#include<map>

#include"HeaderFiles/LFuncs.h"
#include"HeaderFiles/BackEnd.h"
#include"HeaderFiles/DataBase.h"
#include"HeaderFiles/DataStructures.h"
#include"HeaderFiles/SelCutAlgo.h"
#include"HeaderFiles/SelectionCuts.h"

using std::cout;
using std::string;
using std::vector;
using std::map;

gStyle->SetOptStat(111111);
ROOT::EnableImplicitMT();

// structs


// funcs


void DataBase::NoCut(DataStructs::FrameAndData&fd)
{
  auto at_least_two = [](ROOT::VecOps::RVec<Float_t> value){
    return value.size() >= 2;
  };
  auto print = [](ROOT::VecOps::RVec<Float_t> mu_pt, ROOT::VecOps::RVec<Float_t> mu_phi, ROOT::VecOps::RVec<int> mu_status, ROOT::VecOps::RVec<Float_t> el_pt, ROOT::VecOps::RVec<Float_t> el_phi, ROOT::VecOps::RVec<int> el_status){
    cout <<"Muon PT: " << mu_pt << " Phi: " << mu_phi << " status: " << mu_status << "Electron PT: " << el_pt << " Phi: " << el_phi <<" status: " << el_status << "\n";  
  };
  auto print_el = [](ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<int> status){
    cout <<"Electron: " << pt  << " status: " << status <<  "\n";  
  };
    auto jetbnum = [](int jetBnum){
    return jetBnum ==2;
  };
  // fd.Filter(SC::pt_g_10_cut, {"GenElectron_PT"}, "GenElectronG10PT");
  // fd.Filter(SC::pt_g_10_cut, {"GenMuon_PT"}, "GenMuonG10PT");
  // fd.node.Foreach(print, {"GenMuon_PT", "GenMuon_Phi", "GenMuon_Status", "GenElectron_PT", "GenElectron_Phi", "GenElectron_Status"});
  // fd.Filter(jetbnum, {"Jet_BNum"}, "2BJets");
  // fd.node.Foreach(print_el, {"GenElectron_PT", "GenElectron_Status"});
  // fd.Filter(at_least_two, {"Jet_TauTagPT"}, "at_least_two");
}

void DataBase::Analyse(DataStructs::FrameAndData& fd)
{
  fd.Define("NeutrinoPT2", LFuncs::get_col_neutrinopt2, {"MissingET.MET", "Jet_TauTagPhi", "MissingET.Phi"});
  fd.Define("NeutrinoPT1", LFuncs::get_col_neutrinopt1, {"MissingET.MET", "Jet_TauTagPhi", "MissingET.Phi", "NeutrinoPT2"});
  fd.Define("Jet_TauTagNeutrinoPT", LFuncs::add_col_pt, {"Jet_TauTagPT", "NeutrinoPT1", "NeutrinoPT2"});
  fd.Define("TauTagJetDeltaPhi", LFuncs::get_delta_phi, {"Jet_TauTagPhi"});
  fd.Define("AngleBetweenMET", LFuncs::met_jet_ang, {"MissingET.Phi", "Jet_TauTagPhi"});
  // fd.Define("AngleBetweenMET", wraparound_fix, {"Jet_DTauTagPhi", "Jet_DTauTagPT", "MissingET.Phi", "MissingET.MET"});
  fd.Define("TauJetInvMass", LFuncs::inv_mass_ml, {"Jet_TauTagPT", "Jet_TauTagEta", "Jet_TauTagPhi"});
  fd.Define("TauJetInvMassWithNeutrino", LFuncs::inv_mass_ml, {"Jet_TauTagNeutrinoPT", "Jet_TauTagPhi", "Jet_TauTagEta"});
  
  //FILTERS
  fd.Filter(SC::opp_sign, {"Jet_TauTagCharge"}, "OSCut");
  fd.Filter(SC::met_angle_diff_fine, {"AngleBetweenMET"}, "METbetweenJets2E-6");
  // fd.Filter(SC::delta_phi_1p2, {"TauTagJetDeltaPhi"}, "DeltaPhi1p2");

}
void DataBase::TruthAnalysis(DataStructs::FrameAndData& fd){
  // truth analysis was all about matching the generator level particles to the jets and calculating the relevant quantities such
  // as the Delta R the Delta Phi and inv mass with and without neutrinos

}
void DataBase::GeneratorLevelTauAnalysis(DataStructs::FrameAndData& fd){
  // generator level tau analysis was purely on focused on the generator information of the taus and no reco level variables
  fd.Define("GenTauInvMass", LFuncs::inv_mass_ml, {"Tau_PT", "Tau_Phi", "Tau_Eta"});
  cout << fd.frame.GetColumnType("Particle.Status") << " <<<<<THE TYPE OF STATUS\n";
  //fd.node.Foreach([](ROOT::VecOps::RVec<Float_t> tau_pt, ROOT::VecOps::RVec<Float_t> tau_phi, ROOT::VecOps::RVec<Float_t> tau_eta, ROOT::VecOps::RVec<Float_t> tau_status){cout<<"TauPT: "<<tau_pt<< " TauEta: " << tau_eta << " TauPhi: " << tau_phi << " TauStatus: "<< tau_status << " \n";}, {"Tau_PT", "Tau_Phi", "Tau_Eta", "Tau_Status"});

}
void DataBase::SelectionCutAnalysis(DataStructs::FrameAndData& fd){
  // selection cut was all about restricting the inv mass region and then figuring out what the optimal selection cut was
  // Need to restrict IMass range to 80-100
  // then generate columns with the number of events for pt of specific cut.
  // for loop for range and for loop for all events.  thats selcutalgo i think
  // then need to save the example thing and save the new csv file
  // get_integral(fd, column_name, selection_cut, sel_cut_values, resolution)
  //get invmass


  // fd.Define("NeutrinoPT2", LFuncs::get_co/ // fd.Filter(filter_open_met, {"MissingET.Eta", "Jet_TauTagEta"});
  // PrintOff Truth Neutrino PT and reco Neutrino PT
  auto print_neut= [](ROOT::VecOps::RVec<Float_t> neutrino_pt, ROOT::VecOps::RVec<int> status){
    cout <<"tau_pt: " << neutrino_pt  << " status: " << status <<  "\n";  
  };
  auto jetbnum = [](int jetBnum){
    return jetBnum ==2;
  };
  auto get_transverse = [](ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> phi){
    return sqrt(pow(1.777,2) + pow(pt[0] * cos(phi[0]) + pt[1] * cos(phi[1]),2) + pow(pt[0] * sin(phi[0]) + pt[1] * sin(phi[1]), 2) );
  };
  auto get_transverse_mass = [](RV<Float_t> lep_pt, RV<Float_t> met, RV<Float_t> lep_phi, RV<Float_t> met_phi){
    return sqrt(2 * lep_pt[0] * met[0] * (1- cos(LFuncs::get_delta_phi_special(lep_phi[0], met_phi[0]) ) ));
  };
  fd.Define("NeutrinoPT2", LFuncs::get_col_neutrinopt2, {"MissingET.MET", "Jet_DTauTagPhi", "MissingET.Phi"});
  fd.Define("NeutrinoPT1", LFuncs::get_col_neutrinopt1, {"MissingET.MET", "Jet_DTauTagPhi", "MissingET.Phi", "NeutrinoPT2"});
  fd.Define("Jet_DTauTagNeutrinoPT", LFuncs::add_col_pt, {"Jet_DTauTagPT", "NeutrinoPT1", "NeutrinoPT2"});
  fd.Define("TruthMatchDeltaPhi", LFuncs::get_delta_phi, {"Jet_DTauTagPhi"});
  fd.Define("AngleBetweenMET", LFuncs::met_jet_ang, {"MissingET.Phi", "Jet_DTauTagPhi"});
  // fd.Define("AngleBetweenMET", wraparound_fix, {"Jet_DTauTagPhi", "Jet_DTauTagPT", "MissingET.Phi", "MissingET.MET"});
  fd.Define("TauJetInvMass", LFuncs::inv_mass_ml, {"Jet_DTauTagPT", "Jet_DTauTagEta", "Jet_DTauTagPhi"});
  fd.Define("TauJetInvMassWithNeutrino", LFuncs::inv_mass_ml, {"Jet_DTauTagNeutrinoPT", "Jet_DTauTagPhi", "Jet_DTauTagEta"});
  fd.Define("GenTauInvMass", LFuncs::inv_mass_ml, {"Tau_PT", "Tau_Phi", "Tau_Eta"});

  // Truth Match to Truth MET
  // fd.Define("TruthNeutrinoPT2", LFuncs::get_col_neutrinopt2, {"TruthMissingET_MET", "Jet_DTauTagPhi", "TruthMissingET_Phi"});
  // fd.Define("TruthNeutrinoPT1", LFuncs::get_col_neutrinopt1, {"TruthMissingET_MET", "Jet_DTauTagPhi", "TruthMissingET_Phi", "TruthNeutrinoPT2"});
  // fd.Define("Jet_DTauTagTruthNeutrinoPT", LFuncs::add_col_pt, {"Jet_DTauTagPT", "TruthNeutrinoPT1", "TruthNeutrinoPT2"});
  // fd.Define("AngleBetweenTruthMET", LFuncs::met_jet_ang, {"TruthMissingET_Phi", "Jet_DTauTagPhi"});
  // fd.Define("TauJetInvMassWithTruthNeutrino", LFuncs::inv_mass_ml, {"Jet_DTauTagTruthNeutrinoPT", "Jet_DTauTagPhi", "Jet_DTauTagEta"});
  // fd.Define("DifferenceBetweenTruthRecoMET", LFuncs::get_difference_between_vectors, {"TruthMissingET_MET", "MissingET.MET"});
  // fd.Define("NormalisedDifferenceBetweenTruthRecoMET", LFuncs::get_normalised_by_v1_difference_between_vectors, {"TruthMissingET_MET", "MissingET.MET"});
  // fd.Define("TauTagDeltaR", LFuncs::get_DeltaR, {"Jet_TauTagPhi", "Jet_TauTagEta"});
  // TRUTH JET DEFINITIONS.
  // fd.Define("TruthMatchDeltaPhiTRUTHJET", LFuncs::get_delta_phi, {"Jet_TruthTauMatchPT", "Jet_TruthTauMatchPhi"});
  // fd.Define("AngleBetweenMETTRUTHJET", LFuncs::met_jet_ang, {"MissingET.Phi", "Jet_TruthTauMatchPhi"});
  // fd.Define("Jet_TruthTauMatchNeutrinoPT", LFuncs::add_col_pt, {"Jet_TruthTauMatchPT", "NeutrinoPT1", "NeutrinoPT2"});
  // fd.Define("TauJetInvMassTRUTHJET", LFuncs::inv_mass_ml, {"Jet_TruthTauMatchPT", "Jet_TruthTauMatchEta", "Jet_TruthTauMatchPhi"});
  // fd.Define("TauJetInvMassWithNeutrinoTRUTHJET", LFuncs::inv_mass_ml, {"Jet_TruthTauMatchNeutrinoPT", "Jet_TruthTauMatchPhi", "Jet_TruthTauMatchEta"});


  
  // fd.Filter(SC::met_angle_diff, {"AngleBetweenTruthMET"}, "TruthMETbetweenJets0p025");
  // fd.Filter(SC::met_angle_diff_fine, {"AngleBetweenTruthMET"}, "TruthMETbetweenJets2E-6");

  // fd.Filter(SC::met_angle_diff, {"AngleBetweenMET"}, "METbetweenJets0p025");
  // fd.Filter(SC::size_0, {"Electron.PT", "Muon.PT"}, "NoLep");
  // fd.Filter(jetbnum, {"Jet_BNum"}, "2BJets");
  fd.Filter(SC::opp_sign, {"Jet_DTauTagCharge"}, "OSCut");
  fd.Filter(SC::jet_number_cut, {"Jet_Num"}, "g7jet"); 
  fd.Filter(SC::met_angle_diff_fine, {"AngleBetweenMET"}, "METbetweenJets2E-6");
  
  // fd.Filter(SC::norm_truth_reco_diff_cut, {"NormalisedDifferenceBetweenTruthRecoMET"}, "METNormDiff0p3");
  fd.Filter(SC::deltaRcut0p3, {"TruthMatchedDeltaR"}, "DeltaRCut0p3");
  // fd.Filter(SC::truth_and_fake_cut, {"TruthMatchedDeltaR"}, "1Truth1Fake");
  // fd.Filter(SC::all_fake, {"TruthMatchedDeltaR"}, "allFake");
  // fd.Filter(SC::pt_tau_cut, {"Jet_DTauTagPT"}, "jetptg20");
  // fd.Filter(SC::gen_inv_mass_g20, {"GenTauInvMass"}, "GenTauInvMassG20");
  // fd.Filter(SC::deltaRcut0p2, {"DeltaRJetTauSel"}, "DeltaRCut0p2");
  // fd.node.Foreach(print_neut, {"Tau_PT", "Tau_Status"});

  // fd.Filter(SC::delta_phi_1, {"TruthMatchDeltaPhi"}, "DeltaPhi1");
  // fd.Filter(SC::delta_phi_1p2, {"TruthMatchDeltaPhi"}, "DeltaPhi1p2");
  // fd.Filter(SC::delta_phi_1p6, {"TruthMatchDeltaPhi"} , "DeltaPhi1p6");
  // fd.Filter(SC::delta_phi_2, {"TruthMatchDeltaPhi"}, "DeltaPhi2");
  // fd.frame.Foreach([](ROOT::VecOps::RVec<Float_t> gen_tau_pt, ROOT::VecOps::RVec<Float_t> reco_tau_pt){cout<<gen_tau_pt<<" "<<reco_tau_pt<<"\n";}, {"Tau_PT","Jet_DTauTagNeutrinoPT"});
  // vector<vector<Integral<Float_t>>> inout = SCAlgo::get_integral<Float_t>(fd, "TruthMatchDeltaPhi", "TauJetInvMassWithNeutrino", SC::y_lt_x_inout, {0, 3.2}, 10);
  // BackEnd::serialise_integral(inout[0], "DeltaPhiIn" + fd.save_string);
  // BackEnd::serialise_integral(inout[1], "DeltaPhiOut" + fd.save_string);
}

void DataBase::SemiLeptonicAnalysis(DataStructs::FrameAndData& fd){
  // fd.Filter(SC::semi_leptonic_cut, {"Jet_DTauTagCharge", "GenLep_Charge"}, "SemiLeptonicCut");
  auto size_1 = [](RV<int> charge){
    cout<< charge.size();
    return charge.size() == 1;
  };
  auto get_semi_inv_mass = [](RV<Float_t> jet_pt, RV<Float_t> jet_phi, RV<Float_t> jet_eta, RV<Float_t> lep_pt, RV<Float_t> lep_phi, RV<Float_t> lep_eta){
    return sqrt(2* jet_pt[0] * lep_pt[0] * (cosh(jet_eta[0] - lep_eta[0]) - cos(jet_phi[0] - lep_phi[0])));
  };
  auto pt_iso = [](RV<Float_t> iso_var){
    return iso_var[0] < 3;
  };
  auto iso_var_cut = [](RV<Float_t> iso_var){
    return iso_var[0] < 0.08;
  };
  auto get_transverse_mass = [](RV<Float_t> lep_pt, RV<Float_t> met, RV<Float_t> lep_phi, RV<Float_t> met_phi){
    return sqrt(2 * lep_pt[0] * met[0] * (1- cos(LFuncs::get_delta_phi_special(lep_phi[0], met_phi[0]) ) ));
  };
  auto m_t_cut = [](Float_t mt){
    return mt < 60;
  };
  auto m_t_cut_5 = [](Float_t mt){
    return mt < 50;
  };
  auto make_decay_products = [](RV<Float_t> p1, RV<Float_t> p2){
    RV<Float_t> res{}; res.push_back(p1[0]); res.push_back(p2[0]);
    return res;
  };
  auto lep_pt_cut = [](RV<Float_t> pt){return pt[0] > 25;};
  auto bjet_cut = [](RV<Float_t> pt){
    return pt.size() >= 1;
  };
  auto get_omega = [](RV<Float_t> met_phi, RV<Float_t> dp_phi, Float_t metjetang){
    Float_t min{100};
    int index{0};
    int i{0};
    for(auto phi : dp_phi){
      Float_t dphi{LFuncs::get_delta_phi_special(met_phi[0], phi)};
      if(dphi < min){min = dphi; index = i;}
      i++;
    }
    Float_t small_omega{min/LFuncs::get_delta_phi_special(dp_phi[0], dp_phi[1])};
    Float_t omega{};
    Float_t in_jet{2E-6};
    // jettautag is index == 0, lep index == 1
    if(metjetang < in_jet && index == 0){omega= small_omega;}
    else if(metjetang < in_jet && index == 1){omega = 1-small_omega;}
    else if(metjetang > in_jet && index == 0){omega = -small_omega;}
    else if(metjetang > in_jet && index == 1){omega = small_omega + 1;}

    return omega;
  };
  auto omega_cut = [](Float_t omega){
    return -0.2 < omega && omega < 1.2;
  };
  // fd.Filter(size_1, {"Jet_DTauTagCharge"}, "1TMTauJet");
  fd.Filter(SC::semileptonic_cut, {"Jet_TauTagCharge", "Lep_Charge"}, "SemiLeptonicCut");
  // fd.Filter(isolation_cut, {"Lep_IsolationVar"}, "LepIsolationCut0p1");
  fd.Define("SemiLepInvMass", get_semi_inv_mass, {"Jet_TauTagPT", "Jet_TauTagPhi", "Jet_TauTagEta", "Lep_PT", "Lep_Phi", "Lep_Eta"});
  fd.Define("TransverseMass", get_transverse_mass, {"Lep_PT", "MissingET.MET", "Lep_Phi", "MissingET.Phi"});
  fd.Define("DecProd_PT", make_decay_products, {"Jet_TauTagPT", "Lep_PT"});
  fd.Define("DecProd_Phi", make_decay_products, {"Jet_TauTagPhi", "Lep_Phi"});
  fd.Define("DecProd_Eta", make_decay_products, {"Jet_TauTagEta", "Lep_Eta"});
  fd.Define("NeutrinoPT2", LFuncs::get_col_neutrinopt2, {"MissingET.MET", "DecProd_Phi", "MissingET.Phi"});
  fd.Define("NeutrinoPT1", LFuncs::get_col_neutrinopt1, {"MissingET.MET", "DecProd_Phi", "MissingET.Phi", "NeutrinoPT2"});
  fd.Define("DecProd_ColPT", LFuncs::add_col_pt, {"DecProd_PT", "NeutrinoPT1", "NeutrinoPT2"});
  fd.Define("SemiLepInvMassWithNeutrino", LFuncs::inv_mass_ml, {"DecProd_ColPT", "DecProd_PT", "DecProd_Eta"});
  fd.Define("AngleBetweenMET", LFuncs::met_jet_ang, {"MissingET.Phi", "DecProd_Phi"});
  fd.Define("SemiLepDeltaPhi", LFuncs::get_delta_phi, {"DecProd_Phi"});
  fd.Define("Omega", get_omega, {"MissingET.Phi", "DecProd_Phi", "AngleBetweenMET"});
  
  

  // filters
  fd.Filter(m_t_cut, {"TransverseMass"}, "TransverseMassL60");
  // fd.Filter(m_t_cut_5, {"TransverseMass"}, "TransverseMassL50");
  // fd.Filter(pt_iso, {"Lep_SumPT"}, "SumPTL3");
  fd.Filter(iso_var_cut, {"Lep_IsolationVar"}, "IsoVarL0p08");
  fd.Filter(omega_cut, {"Omega"}, "OmegabtwnM0p2a1p2");
  // fd.Filter(SC::met_angle_diff_fine, {"AngleBetweenMET"}, "METbetweenJets2E-6");
  // fd.Filter(SC::met_angle_diff_less_fine, {"AngleBetweenMET"}, "MetBetweenJets3E-6");
  // fd.Filter(SC::met_g_30, {"MissingET.MET"}, "METG30");
  fd.Filter(lep_pt_cut, {"Lep_PT"}, "LepPTG25");
  // fd.Filter(SC::delta_phi_1, {"SemiLepDeltaPhi"}, "DeltaPhiL1");
  fd.Filter(SC::delta_phi_1p6, {"SemiLepDeltaPhi"}, "DeltaPhiL1p6");
  // fd.Filter(SC::delta_phi_2, {"SemiLepDeltaPhi"}, "DeltaPhiL2");
  // fd.Filter(SC::jet_number_cut, {"Jet_Num"}, "JetNumG5");
  fd.Filter(bjet_cut, {"Jet_BTagPT"}, "GE1BJet");
  // fd.Filter()
}

void DataBase::EAnalysis(DataStructs::FrameAndData& fd){
  // comparitive analysis to the general analysis to see how the electrons behave for a comparison with rhys and with the taus.

}

void DataBase::MomentumTest(DataStructs::FrameAndData& fd){
  fd.Define("GenTauInvMass", LFuncs::inv_mass_ml, {"Tau_PT", "Tau_Phi", "Tau_Eta"});
}

void analysis_procedure(string analysis_mode, string mode, string weighting){
  
  vector<string> files = BackEnd::load_files(analysis_mode);
  cout << "Analysis Mode: " << analysis_mode << "\nNumber of files Loaded: " <<files.size()<<"\n";
  string outputFileName = DataBase::database[analysis_mode].fname;
  // HistInfo: {name of hist, x axis title, Column name, nbins, lbound, ubound, units, offset by .5}

  // Analysis code
  cout<<"Creating RDF\n";
  FrameAndData fd{FrameAndData(mode, files)};
  cout<<"PreProcessColumns'\n";
  BackEnd::pre_process_columns(fd);
  BackEnd::initial_cut(fd, DataBase::database[analysis_mode].sel_cuts);
  // Analyse
  cout<<"Analysing\n";
  DataBase::database[analysis_mode].a_func(fd);
  cout<<"Plotting\n";
  BackEnd::plot(outputFileName + fd.save_string, fd, DataBase::database[analysis_mode].hists, weighting);
  BackEnd::save_hist_info(outputFileName, DataBase::database[analysis_mode].hists);
}



void GeneralAnalyser()
{
  // Loading Info
  string mode{"Delphes"};
  vector<string> analysis_modes{"SemiLeptnoicAnalysis"};
  string weighting{"Raw"};
  for(auto analysis_mode : analysis_modes){
    cout<<"\n"<<analysis_mode<<" Analysis\n\n";
    analysis_procedure(analysis_mode, mode, weighting);
  }
}
