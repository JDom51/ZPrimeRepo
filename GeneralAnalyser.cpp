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
  fd.Define("TruthNeutrinoPT2", LFuncs::get_col_neutrinopt2, {"TruthMissingET_MET", "Jet_DTauTagPhi", "TruthMissingET_Phi"});
  fd.Define("TruthNeutrinoPT1", LFuncs::get_col_neutrinopt1, {"TruthMissingET_MET", "Jet_DTauTagPhi", "TruthMissingET_Phi", "TruthNeutrinoPT2"});
  fd.Define("Jet_DTauTagTruthNeutrinoPT", LFuncs::add_col_pt, {"Jet_DTauTagPT", "TruthNeutrinoPT1", "TruthNeutrinoPT2"});
  fd.Define("AngleBetweenTruthMET", LFuncs::met_jet_ang, {"TruthMissingET_Phi", "Jet_DTauTagPhi"});
  fd.Define("TauJetInvMassWithTruthNeutrino", LFuncs::inv_mass_ml, {"Jet_DTauTagTruthNeutrinoPT", "Jet_DTauTagPhi", "Jet_DTauTagEta"});
  fd.Define("DifferenceBetweenTruthRecoMET", LFuncs::get_difference_between_vectors, {"TruthMissingET_MET", "MissingET.MET"});
  fd.Define("NormalisedDifferenceBetweenTruthRecoMET", LFuncs::get_normalised_by_v1_difference_between_vectors, {"TruthMissingET_MET", "MissingET.MET"});
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
  // fd.Filter(SC::jet_number_cut, {"Jet_Num"}, "g4jet"); 
  fd.Filter(SC::met_angle_diff_fine, {"AngleBetweenMET"}, "METbetweenJets2E-6");
  
  // fd.Filter(SC::norm_truth_reco_diff_cut, {"NormalisedDifferenceBetweenTruthRecoMET"}, "METNormDiff0p3");
  fd.Filter(SC::deltaRcut0p3, {"DeltaRJetTauSel"}, "DeltaRCut0p3");
  // fd.Filter(SC::truth_and_fake_cut, {"DeltaRJetTauSel"}, "1Truth1Fake");
  // fd.Filter(SC::all_fake, {"DeltaRJetTauSel"}, "allFake");
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

void DataBase::EAnalysis(DataStructs::FrameAndData& fd){
  // comparitive analysis to the general analysis to see how the electrons behave for a comparison with rhys and with the taus.

}

void DataBase::MomentumTest(DataStructs::FrameAndData& fd){
  fd.Define("GenTauInvMass", LFuncs::inv_mass_ml, {"Tau_PT", "Tau_Phi", "Tau_Eta"});
}

void analysis_procedure(string analysis_mode, string mode){
  vector<string> files = BackEnd::load_files(analysis_mode);
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
  BackEnd::plot(outputFileName + fd.save_string, fd, DataBase::database[analysis_mode].hists);
  BackEnd::save_hist_info(outputFileName, DataBase::database[analysis_mode].hists);
}

void GeneralAnalyser()
{
  // Loading Info
  string mode{"Delphes"};
  vector<string> analysis_modes{"GenLevel"};
  for(auto analysis_mode : analysis_modes){
    cout<<"\n"<<analysis_mode<<" Analysis\n\n";
    analysis_procedure(analysis_mode, mode);
  }
}
