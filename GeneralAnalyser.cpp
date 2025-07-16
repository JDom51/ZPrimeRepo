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
    std::stringstream save_name{}; save_name << "Histograms/" << his[i].name << ".png";
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

}

void Analyse(FrameAndData& fd)
{

}

// main

void GeneralAnalyser()
{
  // Loading Info
  string mode{"Delphes"};
  vector<string> files = load_files();
  vector<string> plotting_columns{};
  vector<HistInfo> histograms{
    {"Transverse Momentum of Electrons", "PT", "Electron.PT", 200, 0, 200, "GeV", true},
    {"Transverse Momentum of Muons", "PT", "Muon.PT", 200, 0, 200, "GeV"},
  };

  // Analysis code
  FrameAndData fd{FrameAndData(mode, files)};
  pre_process_columns(fd);
  Analyse(fd);
  plot(fd, histograms);
}
