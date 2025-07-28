#include<TLorentzVector.h>
#include<vector>
#include<iostream>
#include<string>
#include<fstream>

#include"LFuncs.h"
using std::cout;
using std::string;
using std::vector;

gStyle->SetOptStat(111111);

// structs

struct FrameAndData
{
  ROOT::RDataFrame frame;
  ROOT::RDF::RNode node;
  FrameAndData(string mode, vector<string> files) : frame{mode, files}, node{frame} {}
  ROOT::RDF::RResultPtr<TH1D> Histo1D(string hist_id, string hist_name, int nbins, double lbound, double ubound, string plotting_column)
  {
    return node.Histo1D({hist_id.c_str(), hist_name.c_str(), nbins, lbound, ubound}, plotting_column);
  }
  template<typename T>
  void Define(string column_name, T lambda_func, vector<string> var_names){
    node = node.Define(column_name, lambda_func, var_names);
  }
  template<typename Type> void Filter(Type func, vector<string> column_names)
  {
    node = node.Filter(func, column_names);
  }
};

struct HistInfo
{
  string name;
  string x_axis;
  string column;
  int nbins;
  double lbound;
  double ubound;
  string units;
  bool take_point_five;
  HistInfo(string nam, string x_ax, string m_column, int bins, double lowbound, double upbound, string unit) : name{nam}, x_axis{x_ax}, column{m_column}, nbins{bins}, lbound{lowbound}, ubound{upbound}, units{unit}, take_point_five{false} {}
  HistInfo(string nam, string x_ax, string m_column, int bins, double lowbound, double upbound, string unit, bool m_take_point_five) : name{nam}, x_axis{x_ax},  column{m_column}, nbins{bins}, lbound{lowbound}, ubound{upbound}, units{unit}, take_point_five{m_take_point_five} {}
};


// funcs

vector<string> load_files()
{
  vector<string> files{};
  string file{};
  std::ifstream open_files{"OpenFiles.txt"};
  while(getline(open_files, file))
  {
    files.push_back("/gluster/data/atlas/jdombrowski/DATA/" + file);
  }
  return files;
}

void plot(FrameAndData& fd, vector<HistInfo> his)
{
  vector<ROOT::RDF::RResultPtr<TH1D>> histograms{};
  for(auto hi : his)
  {
    if(hi.take_point_five){hi.lbound -= 0.5; hi.ubound -= 0.5;}
    histograms.emplace_back(fd.Histo1D("Hist", hi.name, hi.nbins, hi.lbound, hi.ubound, hi.column));
  }
  for(int i{0}; i < his.size(); i++)
  {
    TCanvas c = TCanvas();
    std::stringstream x_title{}; x_title << his[i].x_axis << " (" << his[i].units << ")";
    histograms[i]->GetXaxis()->SetTitle(x_title.str().c_str()); histograms[i]->GetYaxis()->SetTitle("Events");
    histograms[i]->GetXaxis()->CenterTitle(); histograms[i]->GetYaxis()->CenterTitle();
    histograms[i]->Draw();
    std::stringstream save_name{}; save_name << "/gluster/data/atlas/jdombrowski/Histograms/" << his[i].name << ".png";
    c.Print(save_name.str().c_str());
  }
}

void pre_process_columns(FrameAndData& fd)
{
  // Constants To Pass
  fd.Define("PARTICLEONE", LFuncs::create_one, {"Particle.Mass"});
  fd.Define("PARTICLEZERO", LFuncs::create_zero, {"Particle.Mass"});
  fd.Define("PARTICLENEGATIVEONE", LFuncs::create_minus_one, {"Particle.Mass"});
  // Tagging
  fd.Define("Jet_BTagIndicies", LFuncs::get_indicies_is_one, {"Jet.BTag"});
  fd.Define("Jet_TauTagIndicies", LFuncs::get_indicies_is_one, {"Jet.TauTag"});
  fd.Define("PositiveParticleIndicies", LFuncs::get_indicies_int, {"Particle.Charge", "PARTICLEONE"});
  fd.Define("NeutralParticleIndicies", LFuncs::get_indicies_int, {"Particle.Charge", "PARTICLEZERO"});
  fd.Define("NegativeParticleIndicies", LFuncs::get_indicies_int, {"Particle.Charge", "PARTICLENEGATIVEONE"});
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
  fd.Define("Jet_TauNum", LFuncs::get_size<unsigned int>, {"Jet_TauTagIndicies"});
  fd.Define("Jet_TauTagPx", LFuncs::get_px, {"Jet.PT", "Jet.Phi"});
  fd.Define("Jet_TauTagPy", LFuncs::get_py, {"Jet.PT", "Jet.Phi"});
  fd.Define("Jet_TauTagPz", LFuncs::get_pz, {"Jet.PT", "Jet.Eta"});
  // Not B or Tau Jets
  fd.Define("Jet_NotBTauTagMass", LFuncs::use_neither_indicies, {"Jet.Mass", "Jet_BTagIndicies", "Jet_TauTagIndicies"});
  fd.Define("Jet_NotBTauTagPT", LFuncs::use_neither_indicies, {"Jet.PT", "Jet_BTagIndicies", "Jet_TauTagIndicies"});
  fd.Define("Jet_NotBTauTagEta", LFuncs::use_neither_indicies, {"Jet.Eta", "Jet_BTagIndicies", "Jet_TauTagIndicies"});
  fd.Define("Jet_NotBTauTagPhi", LFuncs::use_neither_indicies, {"Jet.Phi", "Jet_BTagIndicies", "Jet_TauTagIndicies"});
  fd.Define("Jet_NotBTauNum", LFuncs::get_neither_size, {"Jet.TauTag", "Jet_BTagIndicies", "Jet_TauTagIndicies"});
  // Particle
  fd.Define("Particle_TotInvMass", LFuncs::inv_mass_xyz, {"Particle.E", "Particle.Px", "Particle.Py", "Particle.Pz"});
  // - Positive
  fd.Define("PositiveParticleMass", LFuncs::use_indicies, {"Particle.Mass", "PositiveParticleIndicies"});
  fd.Define("PositiveParticlePx", LFuncs::use_indicies, {"Particle.Px", "PositiveParticleIndicies"});
  fd.Define("PositiveParticlePy", LFuncs::use_indicies, {"Particle.Py", "PositiveParticleIndicies"});
  fd.Define("PositiveParticlePz", LFuncs::use_indicies, {"Particle.Pz", "PositiveParticleIndicies"});
  fd.Define("PositiveParticleP", LFuncs::use_indicies, {"Particle.P", "PositiveParticleIndicies"});
  fd.Define("PositiveParticlePT", LFuncs::use_indicies, {"Particle.PT", "PositiveParticleIndicies"});
  fd.Define("PositiveParticleEta", LFuncs::use_indicies, {"Particle.Eta", "PositiveParticleIndicies"});
  fd.Define("PositiveParticlePhi", LFuncs::use_indicies, {"Particle.Phi", "PositiveParticleIndicies"});
  fd.Define("PositiveParticleRapidity", LFuncs::use_indicies, {"Particle.Rapidity", "PositiveParticleIndicies"});
  // - Negative
  fd.Define("NegativeParticleMass", LFuncs::use_indicies, {"Particle.Mass", "NegativeParticleIndicies"});
  fd.Define("NegativeParticlePx", LFuncs::use_indicies, {"Particle.Px", "NegativeParticleIndicies"});
  fd.Define("NegativeParticlePy", LFuncs::use_indicies, {"Particle.Py", "NegativeParticleIndicies"});
  fd.Define("NegativeParticlePz", LFuncs::use_indicies, {"Particle.Pz", "NegativeParticleIndicies"});
  fd.Define("NegativeParticleP", LFuncs::use_indicies, {"Particle.P", "NegativeParticleIndicies"});
  fd.Define("NegativeParticlePT", LFuncs::use_indicies, {"Particle.PT", "NegativeParticleIndicies"});
  fd.Define("NegativeParticleEta", LFuncs::use_indicies, {"Particle.Eta", "NegativeParticleIndicies"});
  fd.Define("NegativeParticlePhi", LFuncs::use_indicies, {"Particle.Phi", "NegativeParticleIndicies"});
  fd.Define("NegativeParticleRapidity", LFuncs::use_indicies, {"Particle.Rapidity", "NegativeParticleIndicies"});
  // - Neutral
  fd.Define("NeutralParticleMass", LFuncs::use_indicies, {"Particle.Mass", "NeutralParticleIndicies"});
  fd.Define("NeutralParticlePx", LFuncs::use_indicies, {"Particle.Px", "NeutralParticleIndicies"});
  fd.Define("NeutralParticlePy", LFuncs::use_indicies, {"Particle.Py", "NeutralParticleIndicies"});
  fd.Define("NeutralParticlePz", LFuncs::use_indicies, {"Particle.Pz", "NeutralParticleIndicies"});
  fd.Define("NeutralParticleP", LFuncs::use_indicies, {"Particle.P", "NeutralParticleIndicies"});
  fd.Define("NeutralParticlePT", LFuncs::use_indicies, {"Particle.PT", "NeutralParticleIndicies"});
  fd.Define("NeutralParticleEta", LFuncs::use_indicies, {"Particle.Eta", "NeutralParticleIndicies"});
  fd.Define("NeutralParticlePhi", LFuncs::use_indicies, {"Particle.Phi", "NeutralParticleIndicies"});
  fd.Define("NeutralParticleRapidity", LFuncs::use_indicies, {"Particle.Rapidity", "NeutralParticleIndicies"});
  // MET
  fd.Define("MET_Px", LFuncs::get_px, {"MissingET.MET", "MissingET.Phi"});
  fd.Define("MET_Py", LFuncs::get_py, {"MissingET.MET", "MissingET.Phi"});
  fd.Define("MET_Pz", LFuncs::get_pz, {"MissingET.MET", "MissingET.Eta"});
  // Electron
  fd.Define("Electron_Num", LFuncs::get_size<Float_t>, {"Electron.PT"});
  // Muon
  fd.Define("Muon_Num", LFuncs::get_size<Float_t>, {"Muon.PT"});
}

void Analyse(FrameAndData& fd)
{
  cout << fd.node.GetColumnType("MissingET.MET") << "\n";
  auto size_of = [](ROOT::VecOps::RVec<unsigned int> tautag_indicies){return tautag_indicies.size() == 2;};
  fd.Filter(size_of, {"Jet_TauTagIndicies"});
  auto get_pt2 = [](ROOT::VecOps::RVec<Float_t> met, ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> phi_e)
  {
    //ROOT::VecOps::RVec<Float_t> pt2{};
    //pt2.push_back(met[0] * ( sin(phi_e[0]) - cos(phi_e) * tan(phi[0]) ) / ( sin(phi[1]) - cos(phi[1]) * tan(phi_e[0]) ) );
    return met[0] * ( sin(phi_e[0]) - cos(phi_e) * tan(phi[0]) ) / ( sin(phi[1]) - cos(phi[1]) * tan(phi_e[0]) );
  };
  auto get_pt1 = [](ROOT::VecOps::RVec<Float_t> met, ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> phi_e, ROOT::VecOps::RVec<Float_t> pt2)
  {
    //ROOT::VecOps::RVec<Float_t> pt1{};
    //pt1.push_back((met[0] * cos( phi_e[0] ) - pt2[0] * cos( phi[1] ) ) / ( cos( phi[0] ) ) );
    return (met[0] * cos( phi_e[0] ) - pt2[0] * cos( phi[1] ) ) / ( cos( phi[0] ) );
  };
  auto get_met = [](ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> phi_e, Float_t pt1, ROOT::VecOps::RVec<Float_t> pt2)
  {
    return ( pt1 * cos(phi[0]) + pt2[0] * cos(phi[1]) ) / cos(phi_e);
  };
  auto get_inv_mass = [](ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> eta)
  {
    return sqrt(2 * pt[0] * pt[1] * ( cosh(eta[0] - eta[1]) - cos(phi[0] - phi[1]) ) );
  };
  auto add = [](ROOT::VecOps::RVec<Float_t> pt, Float_t pt1, ROOT::VecOps::RVec<Float_t> pt2)
  {
    pt[0] += pt1;
    pt[1] += pt2[0];
    return pt;
  };
  auto filter_open_met = [](ROOT::VecOps::RVec<Float_t> met_eta, ROOT::VecOps::RVec<Float_t> tau_eta)
  {
    return abs(met_eta[0]) > abs(tau_eta[0]) && abs(met_eta[0]) > abs(tau_eta[1]);

  };
  fd.Define("Neutrino_PT2", get_pt2, {"MissingET.MET", "Jet_TauTagPhi", "MissingET.Phi"});
  fd.Define("Neutrino_PT1", get_pt1, {"MissingET.MET", "Jet_TauTagPhi", "MissingET.Phi", "Neutrino_PT2"});
  //fd.Define("MET_TEST", get_met, {"Jet_TauTagPhi", "MissingET.Phi", "Neutrino_PT1", "Neutrino_PT2"});
  fd.Define("Jet_TauTagNeutrinoPT", add, {"Jet_TauTagPT", "Neutrino_PT1", "Neutrino_PT2"});
  fd.Define("TauJetInvMass", get_inv_mass, {"Jet_TauTagPT", "Jet_TauTagPhi", "Jet_TauTagEta"});
  fd.Define("TauJetInvMassWithNeutrino", get_inv_mass, {"Jet_TauTagNeutrinoPT", "Jet_TauTagPhi", "Jet_TauTagEta"});
  fd.Filter(filter_open_met, {"MissingET.Eta", "Jet_TauTagEta"});
}

// main

void GeneralAnalyser()
{
  // Loading Info
  string mode{"Delphes"};
  vector<string> files = load_files();
  vector<string> plotting_columns{};
  vector<HistInfo> histograms{
  {"PT of Tau Jets", "PT", "Jet_TauTagPT", 400, 0, 200, "GeV"},
  {"Eta of Tau Jets", "Eta", "Jet_TauTagEta", 400, -4,4, ""},
  {"Phi of Tau Jets", "Phi", "Jet_TauTagPhi", 400, -4, 4, "Rad"},
  {"Missing Transverse Energy", "PT", "MissingET.MET", 400, 0, 200, "GeV"},
  {"Missing Transvere Energy Phi", "Phi", "MissingET.Phi", 400, -4, 4, "Rad"},
  {"Missing Transverse Energy Eta", "Eta", "MissingET.Eta", 400, -4, 4, ""},
  {"METPx", "Px", "MET_Px", 200, 0, 200, "GeV"},
  {"METPy", "Py", "MET_Py", 200, 0, 200, "GeV"},
  {"METPz", "Pz", "MET_Pz", 200, 0, 200, "GeV"},
  {"TauJetPx", "Px", "Jet_TauTagPx", 400, 0, 200, "GeV"},
  {"TauJetPy", "Py", "Jet_TauTagPy", 400, 0, 200, "GeV"},
  {"TauJetPz", "Pz", "Jet_TauTagPz", 400, 0, 200, "GeV"},
  // {"MET TEST", "PT", "MET_TEST", 400, 0, 200, "GeV"},
   {"Inv Mass without Neutrinos", "IMass", "TauJetInvMass", 100, 0, 200, "GeV"},
   {"Inv Mass with Neutrinos", "IMass", "TauJetInvMassWithNeutrino", 100, 0, 200, "GeV"}
  // {"Transverse Momentum of Electrons", "PT", "Electron.PT", 100, 0, 200, "GeV"},
  // {"Transverse Momentum of Muons", "PT", "Muon.PT", 100, 0, 200, "GeV"},
  // {"Number of Electrons", "Num", "Electron_Num", 10, 0, 10, "", true},
  // {"Pseudorapidity of Electrons", "Eta", "Electron.Eta", 40, -4, 4, ""},
  // {"Phi of Electrons", "Phi", "Electron.Phi", 40, -4, 4, ""},
  // {"Number of Muons", "Num", "Muon_Num", 10, 0, 10, "", true},
  // {"Pseudorapidity of Muons", "Eta", "Muon.Eta", 40, -4, 4, ""},
  // {"Phi of Muon", "Phi", "Muon.Phi", 40, -4, 4, ""}
  };

  // Analysis code
  FrameAndData fd{FrameAndData(mode, files)};
  pre_process_columns(fd);
  Analyse(fd);
  plot(fd, histograms);
}
