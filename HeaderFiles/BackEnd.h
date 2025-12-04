#ifndef BACKEND
#define BACKEND
#include<sstream>
#include<filesystem>
#include<fstream>
#include<iostream>
#include<memory>
#include <future> // for parallel execution


#include"DataStructures.h"
#include"LFuncs.h"
#include"DataBase.h"

using DataStructs::FrameAndData;
using DataStructs::SelectionCut;
using DataStructs::HistInfo;
using std::cout;


namespace BackEnd
{
  vector<string> load_files(string& analysis_mode)
  {
    vector<string> files{};
    string file{};
    // std::ifstream open_files{"OpenFiles.txt"};
    // string open_files{}
    // if(std::filesystem::exists("/gluster/data/atlas/jdombrowski/DATA/" + analysis_mode + ".root")){
      // exists = true;
      // return {"/gluster/data/atlas/jdombrowski/DATA/" + analysis_mode + ".root"};
    // }


    // while(getline(open_files, file))
    for(auto file :  DataBase::database[analysis_mode].ifiles)
    {
      size_t delimiter{file.find(",")};
      cout<<"/gluster/data/atlas/" + file.substr(0, delimiter) + "/DATA/" + file.substr(delimiter + 1) + ".root\n";
      files.push_back("/gluster/data/atlas/" + file.substr(0, delimiter) + "/DATA/" + file.substr(delimiter + 1) + ".root");
    }
    return files;
  }
  void initial_cut(FrameAndData& fd, vector<DataStructs::SelectionCut>& selcuts){
    for(auto sel_cut : selcuts){
      fd.Filter(sel_cut.selection_cut, sel_cut.variables, sel_cut.variables[0]);
    }
  }
  // void save_histograms(string outputFileName, vector<ROOT::RDF::RResultPtr<TH1D>>& histograms, vector<HistInfo>& hist_info)
  // {
    // std::stringstream directory{}; directory << "/gluster/data/atlas/jdombrowski/Histograms/" << outputFileName;
    // Save histograms
    // std::filesystem::create_directory(directory.str());
    // map<int, std::unique_ptr<TCanvas>> cs{};
    // for(int i{0}; i < histograms.size(); i++)
    // {
      // cs[i] = std::make_unique<TCanvas>();
      // cout<<"Drawing...\n";
      // std::stringstream x_title{}; x_title << hist_info[i].x_axis << " (" << hist_info[i].units << ")";
      // histograms[i]->GetXaxis()->SetTitle(x_title.str().c_str()); histograms[i]->GetYaxis()->SetTitle("Events");
      // histograms[i]->GetXaxis()->CenterTitle(); histograms[i]->GetYaxis()->CenterTitle();
      // histograms[i]->Draw();
    // }
    // cout<<"Saving...\n";
    // for(int i{0}; i < histograms.size(); i++){
      // std::stringstream file_name{}; file_name << directory.str() << "/" << hist_info[i].name << ".png";
      // cs[i]->Print(file_name.str().c_str());
    // }
  // }
  // gptr  code 

void save_histograms(
    const std::string& outputFileName,
    std::vector<ROOT::RDF::RResultPtr<TH1D>>& histograms,
    std::vector<HistInfo>& hist_info)
{
    // Enable multi-threading (all available cores)

    // Output ROOT file
    std::filesystem::path dir = "/gluster/data/atlas/jdombrowski/Histograms";
    std::filesystem::create_directories(dir);
    std::string file_path = (dir / (outputFileName + ".root")).string();
    TFile out_file(file_path.c_str(), "RECREATE");

    // Materialize all histograms (compute them in parallel)
    for (auto& hptr : histograms) {
        hptr.GetValue(); // triggers computation
    }

    // Write all histograms to file
    for (size_t i = 0; i < histograms.size(); ++i) {
        TH1D* h = histograms[i].GetPtr(); // get pointer
        h->SetName(hist_info[i].name.c_str()); // optional rename
        h->GetXaxis()->SetTitle((hist_info[i].x_axis + " (" + hist_info[i].units + ")").c_str());
        h->GetYaxis()->SetTitle("Events");
        // h->Scale(6.24 * pow(10, -5)); // Scale to rhys' ATLAS lum 36.6fb-1 epsilon 0.028, cross section = 29fb
        // h->Scale(2.2 * pow(10, -3)); // Scale to my cross section and predicted ATLAS lum 36.6fb-1 and epsilon 0.95, cross section 29fb.
        h->Write(); // write to ROOT file
    }

    out_file.Close();

    std::cout << "Saved " << histograms.size() << " histograms to " << file_path << std::endl;
}


  void save_hist_info(string outputFileName, const vector<HistInfo>& histograms){
    ROOT::EnableImplicitMT();
    std::stringstream complete_file_name{}; complete_file_name << "/gluster/data/atlas/jdombrowski/HistogramInfo/" << outputFileName << ".txt";
    std::ofstream histograms_file{complete_file_name.str()};
    std::stringstream save_stringstream{};
    save_stringstream << "{";
    for(auto h : histograms)
    {
      save_stringstream << "{" << h.name << ", " << h.x_axis << ", " << h.columns[0] << ", " << h.nbins << ", " << h.lbound << ", " << h.ubound << ", " << h.units << ", " << h.take_point_five << "}, ";
    }
    string save_string{save_stringstream.str()};
    save_string[save_string.size()-2] = '}';
    histograms_file << save_string;
    histograms_file.close();
  }

  void plot(string outputFileName, FrameAndData& fd, vector<HistInfo>& his, string weighting){
    vector<ROOT::RDF::RResultPtr<TH1D>> histograms(his.size());
    for(int i{0}; i < his.size(); i++){
      cout<<"Histo1D\n";
      if(his[i].mode == 0){
        histograms[i] = fd.Histo1D("Hist", his[i].name, his[i].nbins, his[i].lbound - 0.5 * his[i].take_point_five, his[i].ubound - 0.5 * his[i].take_point_five, his[i].columns[0]);
      }
    }
    for(int i{0}; i < his.size(); i++){
      histograms[i]->Scale(DataBase::weighting.at(weighting));
    }
    save_histograms(weighting + "_" + outputFileName, histograms, his);

  }

  void pre_process_columns(FrameAndData& fd){
    // Constants To Pass
    fd.Define("PARTICLEONE", LFuncs::create_one, {"Particle.Mass"});
    fd.Define("PARTICLEZERO", LFuncs::create_zero, {"Particle.Mass"});
    fd.Define("PARTICLENEGATIVEONE", LFuncs::create_minus_one, {"Particle.Mass"});
    // Tagging
    fd.Define("Jet_BTagIndicies", LFuncs::get_indicies_is_one, {"Jet.BTag"});
    fd.Define("Jet_TauTagIndicies", LFuncs::get_indicies_is_one, {"Jet.TauTag"});
    // fd.Define("GenJet_TauTagIndicies", LFuncs::get_gen_indicies, {"Jet.TauTag", "GenJet.TauTag"});
    // fd.Define("PositiveParticleIndicies", LFuncs::get_indicies_int, {"Particle.Charge", "PARTICLEONE"});
    // fd.Define("NeutralParticleIndicies", LFuncs::get_indicies_int, {"Particle.Charge", "PARTICLEZERO"});
    // fd.Define("NegativeParticleIndicies", LFuncs::get_indicies_int, {"Particle.Charge", "PARTICLENEGATIVEONE"});
    // Jets
    fd.Define("Jet_Num", LFuncs::get_size<Float_t>, {"Jet.PT"});
    // BTagged Jets
    fd.Define("Jet_BTagMass", LFuncs::use_indicies, {"Jet.Mass", "Jet_BTagIndicies"});
    fd.Define("Jet_BTagPT", LFuncs::use_indicies, {"Jet.PT", "Jet_BTagIndicies"});
    fd.Define("Jet_BTagEta", LFuncs::use_indicies, {"Jet.Eta", "Jet_BTagIndicies"});
    fd.Define("Jet_BTagPhi", LFuncs::use_indicies, {"Jet.Phi", "Jet_BTagIndicies"});
    fd.Define("Jet_BNum", LFuncs::get_size<unsigned int>, {"Jet_BTagIndicies"});
    // TauTagged Jets
    fd.Define("Jet_TauTagMass", LFuncs::use_indicies, {"Jet.Mass", "Jet_TauTagIndicies"});
    fd.Define("Jet_TauTagPT", LFuncs::use_indicies, {"Jet.PT", "Jet_TauTagIndicies"});
    fd.Define("Jet_TauTagEta", LFuncs::use_indicies, {"Jet.Eta", "Jet_TauTagIndicies"});
    fd.Define("Jet_TauTagPhi", LFuncs::use_indicies, {"Jet.Phi", "Jet_TauTagIndicies"});
    fd.Define("Jet_TauTagCharge", LFuncs::use_indicies_int, {"Jet.Charge", "Jet_TauTagIndicies"});
    fd.Define("Jet_TauNum", LFuncs::get_size<unsigned int>, {"Jet_TauTagIndicies"});
    fd.Define("Jet_TauTagPx", LFuncs::get_px, {"Jet_TauTagPT", "Jet_TauTagPhi"});
    fd.Define("Jet_TauTagPy", LFuncs::get_py, {"Jet_TauTagPT", "Jet_TauTagPhi"});
    fd.Define("Jet_TauTagPz", LFuncs::get_pz, {"Jet_TauTagPT", "Jet_TauTagEta"});
    // Not B or Tau Jets
    fd.Define("Jet_NotBTauTagMass", LFuncs::use_neither_indicies, {"Jet.Mass", "Jet_BTagIndicies", "Jet_TauTagIndicies"});
    fd.Define("Jet_NotBTauTagPT", LFuncs::use_neither_indicies, {"Jet.PT", "Jet_BTagIndicies", "Jet_TauTagIndicies"});
    fd.Define("Jet_NotBTauTagEta", LFuncs::use_neither_indicies, {"Jet.Eta", "Jet_BTagIndicies", "Jet_TauTagIndicies"});
    fd.Define("Jet_NotBTauTagPhi", LFuncs::use_neither_indicies, {"Jet.Phi", "Jet_BTagIndicies", "Jet_TauTagIndicies"});
    fd.Define("Jet_NotBTauNum", LFuncs::get_neither_size, {"Jet.TauTag", "Jet_BTagIndicies", "Jet_TauTagIndicies"});
    // GenJet TauTag // THIS DOESNT WORK AS GENJETTAUTAG INDICIES DONT WORK 
    // fd.Define("GenJet_TauTagMass", LFuncs::use_indicies, {"Jet.Mass", "GenJet_TauTagIndicies"});
    // fd.Define("GenJet_TauTagPT", LFuncs::use_indicies, {"Jet.PT", "GenJet_TauTagIndicies"});
    // fd.Define("GenJet_TauTagEta", LFuncs::use_indicies, {"Jet.Eta", "GenJet_TauTagIndicies"});
    // fd.Define("GenJet_TauTagPhi", LFuncs::use_indicies, {"Jet.Phi", "GenJet_TauTagIndicies"});
    // fd.Define("GenJet_TauNum", LFuncs::get_size<unsigned int>, {"GenJet_TauTagIndicies"});
    // fd.Define("GenJet_TauTagPx", LFuncs::get_px, {"GenJet_TauTagPT", "GenJet_TauTagPhi"});
    // fd.Define("GenJet_TauTagPy", LFuncs::get_py, {"GenJet_TauTagPT", "GenJet_TauTagPhi"});
    // fd.Define("GenJet_TauTagPz", LFuncs::get_pz, {"GenJet_TauTagPT", "GenJet_TauTagEta"});
    
    // Gen Particle Indicies
    fd.Define("Particle_TotInvMass", LFuncs::inv_mass_xyz, {"Particle.E", "Particle.Px", "Particle.Py", "Particle.Pz"});

    fd.Define("Particle_ElectronIndicies", LFuncs::get_electron_indicies, {"Particle.PID", "Particle.Status"});
    
    fd.Define("Particle_MuonIndicies", LFuncs::get_muon_indicies, {"Particle.PID", "Particle.Status"});
    

    fd.Define("Particle_TauIndicies_uncut", LFuncs::get_tau_indicies, {"Particle.PID", "Particle.Status", "Particle.PT"});
    fd.Define("Intermediate_TauIndicies", LFuncs::get_intermediate_tau_indicies, {"Particle_TauIndicies_uncut", "Particle.PT"});
    fd.Define("Particle_TauIndicies", LFuncs::get_cut_tau_indicies, {"Particle_TauIndicies_uncut", "Intermediate_TauIndicies"});
    fd.Define("Particle_TauNeutrinoIndicies", LFuncs::get_tau_neutrino_indicies, {"Particle.PID", "Particle.Status", "Intermediate_TauIndicies"});
    // - Tau
    fd.Define("Tau_PT", LFuncs::use_indicies, {"Particle.PT", "Particle_TauIndicies"});
    fd.Define("Tau_Eta", LFuncs::use_indicies, {"Particle.Eta", "Particle_TauIndicies"});
    fd.Define("Tau_Phi", LFuncs::use_indicies, {"Particle.Phi", "Particle_TauIndicies"});
    fd.Define("Tau_Energy", LFuncs::use_indicies, {"Particle.E", "Particle_TauIndicies"});
    fd.Define("Tau_Mass", LFuncs::use_indicies, {"Particle.Mass", "Particle_TauIndicies"});
    fd.Define("Tau_PID", LFuncs::use_indicies_int, {"Particle.PID", "Particle_TauIndicies"});
    fd.Define("Tau_Status", LFuncs::use_indicies_int, {"Particle.Status", "Particle_TauIndicies"});
    fd.Define("Tau_IsPU", LFuncs::use_indicies_int, {"Particle.IsPU", "Particle_TauIndicies"});
    fd.Define("Tau_Alpha", LFuncs::get_alpha, {"Tau_Eta"});
    fd.Define("GenTauDeltaR", LFuncs::calc_delta_r, {"Tau_Phi", "Tau_Eta"});
    fd.Define("Tau_Num", LFuncs::get_size<Float_t>, {"Tau_PT"});
    // Tau Neutrino
    fd.Define("TauNeutrino_PT", LFuncs::use_indicies, {"Particle.PT", "Particle_TauNeutrinoIndicies"});
    fd.Define("TauNeutrino_Eta", LFuncs::use_indicies, {"Particle.Eta", "Particle_TauNeutrinoIndicies"});
    fd.Define("TauNeutrino_Phi", LFuncs::use_indicies, {"Particle.Phi", "Particle_TauNeutrinoIndicies"});
    fd.Define("TauNeutrino_Energy", LFuncs::use_indicies, {"Particle.E", "Particle_TauNeutrinoIndicies"});
    fd.Define("TauNeutrino_Mass", LFuncs::use_indicies, {"Particle.Mass", "Particle_TauNeutrinoIndicies"});
    fd.Define("TauNeutrino_Status", LFuncs::use_indicies_int, {"Particle.Status", "Particle_TauNeutrinoIndicies"});
    fd.Define("TauNeutrino_IsPU", LFuncs::use_indicies_int, {"Particle.IsPU", "Particle_TauNeutrinoIndicies"});
    fd.Define("TauNeutrino_Alpha", LFuncs::get_alpha, {"TauNeutrino_Eta"});
    //fd.Define("GenTauDeltaR", LFuncs::get_DeltaR, {"Tau_Phi", "Tau_Eta"});
    fd.Define("TauNeutrino_Num", LFuncs::get_size<Float_t>, {"TauNeutrino_PT"});
    // TruthMET
    fd.Define("TruthMissingET_MET", LFuncs::get_truth_met_met, {"TauNeutrino_PT", "TauNeutrino_Phi"});
    fd.Define("TruthMissingET_Phi", LFuncs::get_truth_met_phi, {"TauNeutrino_PT", "TauNeutrino_Phi"});
    // Delta R Tau
    fd.Define("DeltaRJetTau", LFuncs::get_delta_r, {"Jet_TauTagPhi", "Jet_TauTagEta", "Tau_Phi", "Tau_Eta"});
    fd.Define("TruthMatchTau", LFuncs::get_truth_match, {"DeltaRJetTau"});
    fd.Define("TruthMatchIndicies", LFuncs::get_truth_match_indicies, {"TruthMatchTau"});
    fd.Define("TruthMatchedDeltaR", LFuncs::get_truth_match_indicies_dr, {"TruthMatchIndicies", "DeltaRJetTau"});
    fd.Define("DeltaR_Num", LFuncs::get_size<unsigned int>, {"TruthMatchIndicies"});

    // fd.Define("NTruthMatchIndicies", LFuncs::get_non_truth_match_indicies, {"TruthMatchTau"});
    // fd.Define("NTruthMatchedDeltaR", LFuncs::get_truth_match_indicies_dr, {"NTruthMatchIndicies", "DeltaRJetTau"});
    // Delta R Tau to Generic Jet
    // fd.Define("DeltaRJetTau1TRUTHJET", LFuncs::get_delta_r_1, {"Tau_Eta", "Jet.Eta", "Tau_Phi", "Jet.Phi"});
    // fd.Define("DeltaRJetTau2TRUTHJET", LFuncs::get_delta_r_2, {"Tau_Eta", "Jet.Eta", "Tau_Phi", "Jet.Phi"});
    // fd.Define("DeltaRIndiciesTRUTHJET", LFuncs::get_delta_r_indicies, {"DeltaRJetTau1TRUTHJET", "DeltaRJetTau2TRUTHJET"});
    // fd.Define("DeltaRJetTauSelTRUTHJET", LFuncs::use_delr_indicies, {"DeltaRJetTau1TRUTHJET", "DeltaRJetTau2TRUTHJET", "TruthMatchIndiciesTRUTHJET"});
    // fd.Define("DeltaR_NumTRUTHJET", LFuncs::get_size<unsigned int>, {"TruthMatchIndiciesTRUTHJET"});


    // Selected TauTag Jets DeltaR
    fd.Define("Jet_DTauTagMass", LFuncs::use_indicies, {"Jet_TauTagMass", "TruthMatchIndicies"});
    fd.Define("Jet_DTauTagPT", LFuncs::use_indicies, {"Jet_TauTagPT", "TruthMatchIndicies"});
    fd.Define("Jet_DTauTagEta", LFuncs::use_indicies, {"Jet_TauTagEta", "TruthMatchIndicies"});
    fd.Define("Jet_DTauTagPhi", LFuncs::use_indicies, {"Jet_TauTagPhi", "TruthMatchIndicies"});
    fd.Define("Jet_DTauTagCharge", LFuncs::use_indicies_int, {"Jet_TauTagCharge", "TruthMatchIndicies"});
    fd.Define("Jet_DTauNum", LFuncs::get_size<unsigned int>, {"TruthMatchIndicies"});
    fd.Define("Jet_DTauTagPx", LFuncs::get_px, {"Jet_DTauTagPT", "Jet_DTauTagPhi"});
    fd.Define("Jet_DTauTagPy", LFuncs::get_py, {"Jet_DTauTagPT", "Jet_DTauTagPhi"});
    fd.Define("Jet_DTauTagPz", LFuncs::get_pz, {"Jet_DTauTagPT", "Jet_DTauTagEta"});
    fd.Define("Jet_DTauTagAlpha", LFuncs::get_alpha, {"Jet_DTauTagEta"}); // Alpha is tan-1 (pt/pz)
    fd.Define("TauJetDeltaR", LFuncs::calc_delta_r, {"Jet_DTauTagPhi", "Jet_DTauTagEta"});
    // fd.Define("Jet_DTauTagMT", LFuncs::get_transverse_mass, {"Jet_DTauTag"})
    // Select Non TruthTauTag Jets DeltaR
    // fd.Define("Jet_NDTauTagMass", LFuncs::use_indicies, {"Jet_TauTagMass", "NTruthMatchIndicies"});
    // fd.Define("Jet_NDTauTagPT", LFuncs::use_indicies, {"Jet_TauTagPT", "NTruthMatchIndicies"});
    // fd.Define("Jet_NDTauTagEta", LFuncs::use_indicies, {"Jet_TauTagEta", "NTruthMatchIndicies"});
    // fd.Define("Jet_NDTauTagPhi", LFuncs::use_indicies, {"Jet_TauTagPhi", "NTruthMatchIndicies"});
    // fd.Define("Jet_NDTauNum", LFuncs::get_size<unsigned int>, {"NTruthMatchIndicies"});
    // fd.Define("Jet_NDTauTagPx", LFuncs::get_px, {"Jet_NDTauTagPT", "Jet_NDTauTagPhi"});
    // fd.Define("Jet_NDTauTagPy", LFuncs::get_py, {"Jet_NDTauTagPT", "Jet_NDTauTagPhi"});
    // fd.Define("Jet_NDTauTagPz", LFuncs::get_pz, {"Jet_NDTauTagPT", "Jet_NDTauTagEta"});
    // fd.Define("Jet_NDTauTagAlpha", LFuncs::get_alpha, {"Jet_NDTauTagEta"}); // Alpha is tan-1 (pt/pz)
    // fd.Define("NTauJetDeltaR", LFuncs::calc_delta_r, {"Jet_NDTauTagPhi", "Jet_NDTauTagEta"});
    // Truth Matched JET to Tau.
    // fd.Define("Jet_TruthTauMatchMass", LFuncs::use_indicies, {"Jet.Mass", "TruthMatchIndiciesTRUTHJET"});
    // fd.Define("Jet_TruthTauMatchPT", LFuncs::use_indicies, {"Jet.PT", "TruthMatchIndiciesTRUTHJET"});
    // fd.Define("Jet_TruthTauMatchEta", LFuncs::use_indicies, {"Jet.Eta", "TruthMatchIndiciesTRUTHJET"});
    // fd.Define("Jet_TruthTauMatchPhi", LFuncs::use_indicies, {"Jet.Phi", "TruthMatchIndiciesTRUTHJET"});
    // fd.Define("Jet_TruthTauMatch", LFuncs::get_size<unsigned int>, {"TruthMatchIndiciesTRUTHJET"});
    // fd.Define("Jet_TruthTauMatchPx", LFuncs::get_px, {"Jet_TruthTauMatchPT", "Jet_TruthTauMatchPhi"});
    // fd.Define("Jet_TruthTauMatchPy", LFuncs::get_py, {"Jet_TruthTauMatchPT", "Jet_TruthTauMatchPhi"});
    // fd.Define("Jet_TruthTauMatchPz", LFuncs::get_pz, {"Jet_TruthTauMatchPT", "Jet_TruthTauMatchEta"});
    
    // GenElectron
    fd.Define("GenElectron_PT",  LFuncs::use_indicies, {"Particle.PT", "Particle_ElectronIndicies"});
    fd.Define("GenElectron_Eta", LFuncs::use_indicies, {"Particle.Eta", "Particle_ElectronIndicies"});
    fd.Define("GenElectron_Phi", LFuncs::use_indicies, {"Particle.Phi", "Particle_ElectronIndicies"});
    fd.Define("GenElectron_Energy", LFuncs::use_indicies, {"Particle.E", "Particle_ElectronIndicies"});
    fd.Define("GenElectron_Charge", LFuncs::use_indicies_int, {"Particle.Charge", "Particle_ElectronIndicies"});
    fd.Define("GenElectron_Mass", LFuncs::use_indicies, {"Particle.Mass", "Particle_ElectronIndicies"});
    fd.Define("GenElectron_PID", LFuncs::use_indicies_int, {"Particle.PID", "Particle_ElectronIndicies"});
    fd.Define("GenElectron_Status", LFuncs::use_indicies_int, {"Particle.Status", "Particle_ElectronIndicies"});
    fd.Define("GenElectron_IsPU", LFuncs::use_indicies_int, {"Particle.IsPU", "Particle_ElectronIndicies"});
    // fd.Define("GenElectron_SumPt", LFuncs::use_indicies, {"Particle.SumPt", "Particle_Electron_indicies"});
    // fd.Define("GenElectron_SumPtCharged", LFuncs::use_indicies, {"Particle.SumPtCharged", "Particle_Electron_indicies"});
    // fd.Define("GenElectron_SumPtNeutral", LFuncs::use_indicies, {"Particle.SumPtNeutral", "Particle_Electron_indicies"});
    fd.Define("GenElectron_Num", LFuncs::get_size<Float_t>, {"GenElectron_PT"});
    // Gen Muon
    fd.Define("GenMuon_PT",  LFuncs::use_indicies, {"Particle.PT", "Particle_MuonIndicies"});
    fd.Define("GenMuon_Eta", LFuncs::use_indicies, {"Particle.Eta", "Particle_MuonIndicies"});
    fd.Define("GenMuon_Phi", LFuncs::use_indicies, {"Particle.Phi", "Particle_MuonIndicies"});
    fd.Define("GenMuon_Energy", LFuncs::use_indicies, {"Particle.E", "Particle_MuonIndicies"});
    fd.Define("GenMuon_Charge", LFuncs::use_indicies_int, {"Particle.Charge", "Particle_MuonIndicies"});
    fd.Define("GenMuon_Mass", LFuncs::use_indicies, {"Particle.Mass", "Particle_MuonIndicies"});
    fd.Define("GenMuon_PID", LFuncs::use_indicies_int, {"Particle.PID", "Particle_MuonIndicies"});
    fd.Define("GenMuon_Status", LFuncs::use_indicies_int, {"Particle.Status", "Particle_MuonIndicies"});
    fd.Define("GenMuon_IsPU", LFuncs::use_indicies_int, {"Particle.IsPU", "Particle_MuonIndicies"});
    // fd.Define("GenMuon_SumPt", LFuncs::use_indicies, {"Particle.SumPt", "Particle_Muon_indicies"});
    // fd.Define("GenMuon_SumPtCharged", LFuncs::use_indicies, {"Particle.SumPtCharged", "Particle_Muon_indicies"});
    // fd.Define("GenMuon_SumPtNeutral", LFuncs::use_indicies, {"Particle.SumPtNeutral", "Particle_Muon_indicies"});
    fd.Define("GenMuon_Num", LFuncs::get_size<Float_t>, {"GenMuon_PT"});

    // Indicies for pt cut variables.
    fd.Define("Particle_ElectronIndicies_PTCUT", LFuncs::get_g10_pt_indicies, {"GenElectron_PT", "GenElectron_Eta"});
    fd.Define("Particle_MuonIndicies_PTCUT", LFuncs::get_g10_pt_indicies, {"GenMuon_PT", "GenMuon_Eta"});
    // GenElectron cut
    fd.Define("GenElectronCut_PT",  LFuncs::use_indicies, {"GenElectron_PT", "Particle_ElectronIndicies_PTCUT"});
    fd.Define("GenElectronCut_Eta", LFuncs::use_indicies, {"GenElectron_Eta", "Particle_ElectronIndicies_PTCUT"});
    fd.Define("GenElectronCut_Phi", LFuncs::use_indicies, {"GenElectron_Phi", "Particle_ElectronIndicies_PTCUT"});
    fd.Define("GenElectronCut_Energy", LFuncs::use_indicies, {"GenElectron_Energy", "Particle_ElectronIndicies_PTCUT"});
    fd.Define("GenElectronCut_Charge", LFuncs::use_indicies_int, {"GenElectron_Charge", "Particle_ElectronIndicies_PTCUT"});
    fd.Define("GenElectronCut_Mass", LFuncs::use_indicies, {"GenElectron_Mass", "Particle_ElectronIndicies_PTCUT"});
    fd.Define("GenElectronCut_PID", LFuncs::use_indicies_int, {"Particle.PID", "Particle_MuonIndicies"});
    fd.Define("GenElectronCut_Status", LFuncs::use_indicies_int, {"GenElectron_Status", "Particle_ElectronIndicies_PTCUT"});
    fd.Define("GenElectronCut_IsPU", LFuncs::use_indicies_int, {"GenElectron_IsPU", "Particle_ElectronIndicies_PTCUT"});
    // fd.Define("GenElectronCut_SumPt", LFuncs::use_indicies, {"GenElectron_SumPt", "Particle_Electron_indicies_PTCUT"});
    // fd.Define("GenElectronCut_SumPtCharged", LFuncs::use_indicies, {"GenElectron_SumPtCharged", "Particle_Electron_indicies_PTCUT"});
    // fd.Define("GenElectronCut_SumPtNeutral", LFuncs::use_indicies, {"GenElectron_SumPtNeutral", "Particle_Electron_indicies_PTCUT"});
    fd.Define("GenElectronCut_Num", LFuncs::get_size<Float_t>, {"GenElectronCut_PT"});    
    // GenMuon cut
    fd.Define("GenMuonCut_PT",  LFuncs::use_indicies, {"GenMuon_PT", "Particle_MuonIndicies_PTCUT"});
    fd.Define("GenMuonCut_Eta", LFuncs::use_indicies, {"GenMuon_Eta", "Particle_MuonIndicies_PTCUT"});
    fd.Define("GenMuonCut_Phi", LFuncs::use_indicies, {"GenMuon_Phi", "Particle_MuonIndicies_PTCUT"});
    fd.Define("GenMuonCut_Energy", LFuncs::use_indicies, {"GenMuon_Energy", "Particle_MuonIndicies_PTCUT"});
    fd.Define("GenMuonCut_Charge", LFuncs::use_indicies_int, {"GenMuon_Charge", "Particle_MuonIndicies_PTCUT"});
    fd.Define("GenMuonCut_PID", LFuncs::use_indicies_int, {"Particle.PID", "Particle_MuonIndicies"});
    fd.Define("GenMuonCut_Mass", LFuncs::use_indicies, {"GenMuon_Mass", "Particle_MuonIndicies_PTCUT"});
    fd.Define("GenMuonCut_Status", LFuncs::use_indicies_int, {"GenMuon_Status", "Particle_MuonIndicies_PTCUT"});
    fd.Define("GenMuonCut_IsPU", LFuncs::use_indicies_int, {"GenMuon_IsPU", "Particle_MuonIndicies_PTCUT"});
    // fd.Define("GenMuonCut_SumPt", LFuncs::use_indicies, {"GenMuon_SumPt", "Particle_Muon_indicies_PTCUT"});
    // fd.Define("GenMuonCut_SumPtCharged", LFuncs::use_indicies, {"GenMuon_SumPtCharged", "Particle_Muon_indicies_PTCUT"});
    // fd.Define("GenMuonCut_SumPtNeutral", LFuncs::use_indicies, {"GenMuon_SumPtNeutral", "Particle_Muon_indicies_PTCUT"});
    fd.Define("GenMuonCut_Num", LFuncs::get_size<Float_t>, {"GenMuonCut_PT"});
    
    // GenLEPTON
    
    fd.Define("GenLep_PT",  LFuncs::combine, {"GenElectronCut_PT", "GenMuonCut_PT"});
    fd.Define("GenLep_Eta", LFuncs::combine, {"GenElectronCut_Eta", "GenMuonCut_Eta"});
    fd.Define("GenLep_Phi", LFuncs::combine, {"GenElectronCut_Phi", "GenMuonCut_Phi"});
    fd.Define("GenLep_Energy", LFuncs::combine, {"GenElectronCut_Energy", "GenMuonCut_Energy"});
    fd.Define("GenLep_Charge", LFuncs::combine_int, {"GenElectronCut_Charge", "GenMuonCut_Charge"});
    fd.Define("GenLep_Mass", LFuncs::combine, {"GenElectronCut_Mass", "GenMuonCut_Mass"});
    fd.Define("GenLep_PID", LFuncs::combine_int, {"GenElectronCut_PID", "GenMuonCut_PID"});
    fd.Define("GenLep_Status", LFuncs::combine_int, {"GenElectronCut_Status", "GenMuonCut_Status"});
    fd.Define("GenLep_IsPU", LFuncs::combine_int, {"GenElectronCut_IsPU", "GenMuonCut_IsPU"});
    fd.Define("GenLep_Num", LFuncs::get_size<Float_t>, {"GenLep_PT"});

    // Lep
    fd.Define("Lep_PT",  LFuncs::combine, {"Electron.PT", "Muon.PT"});
    fd.Define("Lep_Eta", LFuncs::combine, {"Electron.Eta", "Muon.Eta"});
    fd.Define("Lep_Phi", LFuncs::combine, {"Electron.Phi", "Muon.Phi"});
    // fd.Define("Lep_Energy", LFuncs::combine, {"Electron.Energy", "Muon.Energy"});
    fd.Define("Lep_Charge", LFuncs::combine_int, {"Electron.Charge", "Muon.Charge"});
    // fd.Define("Lep_Mass", LFuncs::combine, {"Electron.Mass", "Muon.Mass"});
    fd.Define("Lep_IsolationVar", LFuncs::combine, {"Electron.IsolationVar", "Muon.IsolationVar"});
    fd.Define("Lep_SumPT", LFuncs::combine, {"Electron.SumPt", "Muon.SumPt"});
    // fd.Define("Lep_PID", LFuncs::combine_int, {"Electron.PID", "Muon.PID"});
    // fd.Define("Lep_Status", LFuncs::combine_int, {"Electron.Status", "Muon.Status"});
    // fd.Define("Lep_IsPU", LFuncs::combine_int, {"Electron.IsPU", "Muon.IsPU"});
    fd.Define("Lep_Num", LFuncs::get_size<Float_t>, {"Lep_PT"});
    
    // - Positive
    // fd.Define("PositiveParticleMass", LFuncs::use_indicies, {"Particle.Mass", "PositiveParticleIndicies"});
    // fd.Define("PositiveParticlePx", LFuncs::use_indicies, {"Particle.Px", "PositiveParticleIndicies"});
    // fd.Define("PositiveParticlePy", LFuncs::use_indicies, {"Particle.Py", "PositiveParticleIndicies"});
    // fd.Define("PositiveParticlePz", LFuncs::use_indicies, {"Particle.Pz", "PositiveParticleIndicies"});
    // fd.Define("PositiveParticleP", LFuncs::use_indicies, {"Particle.P", "PositiveParticleIndicies"});
    // fd.Define("PositiveParticlePT", LFuncs::use_indicies, {"Particle.PT", "PositiveParticleIndicies"});
    // fd.Define("PositiveParticleEta", LFuncs::use_indicies, {"Particle.Eta", "PositiveParticleIndicies"});
    // fd.Define("PositiveParticlePhi", LFuncs::use_indicies, {"Particle.Phi", "PositiveParticleIndicies"});
    // fd.Define("PositiveParticleRapidity", LFuncs::use_indicies, {"Particle.Rapidity", "PositiveParticleIndicies"});
    // - Negative
    // fd.Define("NegativeParticleMass", LFuncs::use_indicies, {"Particle.Mass", "NegativeParticleIndicies"});
    // fd.Define("NegativeParticlePx", LFuncs::use_indicies, {"Particle.Px", "NegativeParticleIndicies"});
    // fd.Define("NegativeParticlePy", LFuncs::use_indicies, {"Particle.Py", "NegativeParticleIndicies"});
    // fd.Define("NegativeParticlePz", LFuncs::use_indicies, {"Particle.Pz", "NegativeParticleIndicies"});
    // fd.Define("NegativeParticleP", LFuncs::use_indicies, {"Particle.P", "NegativeParticleIndicies"});
    // fd.Define("NegativeParticlePT", LFuncs::use_indicies, {"Particle.PT", "NegativeParticleIndicies"});
    // fd.Define("NegativeParticleEta", LFuncs::use_indicies, {"Particle.Eta", "NegativeParticleIndicies"});
    // fd.Define("NegativeParticlePhi", LFuncs::use_indicies, {"Particle.Phi", "NegativeParticleIndicies"});
    // fd.Define("NegativeParticleRapidity", LFuncs::use_indicies, {"Particle.Rapidity", "NegativeParticleIndicies"});
    // - Neutral
    // fd.Define("NeutralParticleMass", LFuncs::use_indicies, {"Particle.Mass", "NeutralParticleIndicies"});
    // fd.Define("NeutralParticlePx", LFuncs::use_indicies, {"Particle.Px", "NeutralParticleIndicies"});
    // fd.Define("NeutralParticlePy", LFuncs::use_indicies, {"Particle.Py", "NeutralParticleIndicies"});
    // fd.Define("NeutralParticlePz", LFuncs::use_indicies, {"Particle.Pz", "NeutralParticleIndicies"});
    // fd.Define("NeutralParticleP", LFuncs::use_indicies, {"Particle.P", "NeutralParticleIndicies"});
    // fd.Define("NeutralParticlePT", LFuncs::use_indicies, {"Particle.PT", "NeutralParticleIndicies"});
    // fd.Define("NeutralParticleEta", LFuncs::use_indicies, {"Particle.Eta", "NeutralParticleIndicies"});
    // fd.Define("NeutralParticlePhi", LFuncs::use_indicies, {"Particle.Phi", "NeutralParticleIndicies"});
    // fd.Define("NeutralParticleRapidity", LFuncs::use_indicies, {"Particle.Rapidity", "NeutralParticleIndicies"});
    // MET
    fd.Define("MET_Px", LFuncs::get_px, {"MissingET.MET", "MissingET.Phi"});
    fd.Define("MET_Py", LFuncs::get_py, {"MissingET.MET", "MissingET.Phi"});
    fd.Define("MET_Pz", LFuncs::get_pz, {"MissingET.MET", "MissingET.Eta"});
    // Electron
    fd.Define("Electron_Num", LFuncs::get_size<Float_t>, {"Electron.PT"});
    // Muon
    fd.Define("Muon_Num", LFuncs::get_size<Float_t>, {"Muon.PT"});
  }
  template<typename T> void serialise_integral(vector<DataStructs::Integral<T>>& integral, string fname){
    std::stringstream ss{}; ss<<"Value,Integral\n";
    for(auto event : integral){
      ss<<event.value<<","<<event.integral<<"\n";
    }
    cout<<ss.str();
    std::ofstream outfile{"/gluster/data/atlas/jdombrowski/SelectionCutIntegrals/" + fname + ".txt"};
    outfile<<ss.str();
    outfile.close();
  }
};

#endif