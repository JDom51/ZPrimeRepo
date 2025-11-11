#ifndef DATA_BASE
#define DATA_BASE
#include"DataStructures.h"
#include"SelectionCuts.h"

using DataStructs::FrameAndData;
using DataStructs::HistInfo;
using SC::size_2;

namespace DataBase{
  void NoCut(FrameAndData& fd);
  void Analyse(FrameAndData& fd);
  void TruthAnalysis(FrameAndData& fd);
  void GeneratorLevelTauAnalysis(FrameAndData& fd);
  void SelectionCutAnalysis(FrameAndData& fd);
  void EAnalysis(FrameAndData& fd);
  typedef void (*AnalysisFunc)(FrameAndData& fd);


  map<string, vector<string>> loadable_files{
    {"NoCut", {"jdombrowski,pptata_exh_run_01", "jdombrowski,pptata_exh_run_02",
  "jdombrowski,pptata_exh_run_03", "jdombrowski,pptata_exh_run_04",
  "jdombrowski,pptata_exh_run_05","jdombrowski,pptata_exh_run_06",
  "jdombrowski,pptata_exh_run_07","jdombrowski,pptata_exh_run_08",
  "jdombrowski,pptata_exh_run_09","jdombrowski,pptata_exh_run_10"}},
  {"RecoLevel", {"jdombrowski,pptata_exh_run_01", "jdombrowski,pptata_exh_run_02",
  "jdombrowski,pptata_exh_run_03", "jdombrowski,pptata_exh_run_04",
  "jdombrowski,pptata_exh_run_05","jdombrowski,pptata_exh_run_06",
  "jdombrowski,pptata_exh_run_07","jdombrowski,pptata_exh_run_08",
  "jdombrowski,pptata_exh_run_09","jdombrowski,pptata_exh_run_10"} },
    {"Gen+RecoLevel", {"jdombrowski,pptata_exh_run_01", "jdombrowski,pptata_exh_run_02",
  "jdombrowski,pptata_exh_run_03", "jdombrowski,pptata_exh_run_04",
  "jdombrowski,pptata_exh_run_05","jdombrowski,pptata_exh_run_06",
  "jdombrowski,pptata_exh_run_07","jdombrowski,pptata_exh_run_08",
  "jdombrowski,pptata_exh_run_09","jdombrowski,pptata_exh_run_10"}},
   {"SelectionCut", {"jdombrowski,pptata_exh_run_01", "jdombrowski,pptata_exh_run_02",
  "jdombrowski,pptata_exh_run_03", "jdombrowski,pptata_exh_run_04",
  "jdombrowski,pptata_exh_run_05","jdombrowski,pptata_exh_run_06",
  "jdombrowski,pptata_exh_run_07","jdombrowski,pptata_exh_run_08",
  "jdombrowski,pptata_exh_run_09","jdombrowski,pptata_exh_run_10"}},
    {"GenLevel", {"jdombrowski,pptata_exh_run_01", "jdombrowski,pptata_exh_run_02",
  "jdombrowski,pptata_exh_run_03", "jdombrowski,pptata_exh_run_04",
  "jdombrowski,pptata_exh_run_05","jdombrowski,pptata_exh_run_06",
  "jdombrowski,pptata_exh_run_07","jdombrowski,pptata_exh_run_08",
  "jdombrowski,pptata_exh_run_09","jdombrowski,pptata_exh_run_10"}},
    {"EE", {"jdombrowski,ppzttee_run_01" ,"jdombrowski,ppzttee_run_02","jdombrowski,ppzttee_run_03",
      "jdombrowski,ppzttee_run_04","jdombrowski,ppzttee_run_05","jdombrowski,ppzttee_run_06",
      "jdombrowski,ppzttee_run_07","jdombrowski,ppzttee_run_08","jdombrowski,ppzttee_run_09",
      "jdombrowski,ppzttee_run_10"}}
  };
  map<string, AnalysisFunc> analysis_funcs{
    {"NoCut", *NoCut},
    {"RecoLevel", *Analyse},
    {"Gen+RecoLevel", *TruthAnalysis},
    {"GenLevel", *GeneratorLevelTauAnalysis},
    {"SelectionCut", *SelectionCutAnalysis},
    {"EE", *EAnalysis}
  };

  std::map<string, vector<DataStructs::SelectionCut>> broad_cuts{
    {"NoCut", {}},
    {"RecoLevel", {{size_2, {"Jet_TauTagIndicies"}}}},
    {"Gen+RecoLevel", {{size_2, {"DeltaRIndicies"}}}},
    {"SelectionCut", {{size_2, {"DeltaRIndicies"}}}},
    {"GenLevel", {}},
    {"EE", {}}
  };


  map<string, string> fnames{
    {"NoCut", "NoCut"},
    {"RecoLevel", "Ztata Reco"},
    {"Gen+RecoLevel", "Ztata Gen+Reco"},
    {"SelectionCut", "Selcut"},
    {"GenLevel", "Ztata Gen"},
    {"EE", "Zee, n_e=2"}
  };


  map<string, vector<HistInfo>> histograms{
    {"NoCut", {{"Number of Jets", "Number", "Jet_Num", 10, 0, 10, "", true},
              {"MissingET", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
              {"Number Of BJets", "Number", "Jet_BNum", 10, 0, 10, "", true},
              {"TauJetPT", "PT", "Jet_TauTagPT", 200, 0, 200, "GeV"},
              {"TauPT", "PT", "Tau_PT", 200, 0, 200, "GeV"},
              {"Number of Tau Jets", "Number", "Jet_TauNum", 10, 0, 10, "", true}}},
    {"RecoLevel", {{"MET", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
                   {"Inv Mass without Neutrinos", "IMass", "TauJetInvMass", 200, 0, 200, "GeV"},
                   {"Neutrino PT1", "PT", "Neutrino_PT1", 200, 0, 200, "GeV"},
                   {"Neutrino PT2", "PT", "Neutrino_PT2", 200, 0, 200, "GeV"},
                   {"Inv Mass with Neutrinos", "IMass", "TauJetInvMassWithNeutrino", 200, 0, 200, "GeV"},
                   {"Number of Jets", "Num", "Jet_Num", 10, 0, 10, "", true},
                   {"Number of B Jets", "Num", "Jet_BNum", 10, 0, 10, "", true},
                  {"Angle diff between met and jets", "$\\Delta \\phi$", "AngleBetweenMET", 200, 0, 4, "rad"}}},

    {"Gen+RecoLevel", {
      // {"Number of Jets", "Jet Num", "Jet_Num", 7, 0, 7, "", true},
      {"MET of Tau", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
      {"Delta R of Tau and Selected Tau Jets", "Delta R", "DeltaRJetTauSel", 400, 0, 0.1, ""},
      {"Delta Phi of Tau Leptons", "\\Delta\\Phi", "DeltaPhiGenTau", 100, -0.5, 4, "rad"},
      {"Delta Phi of Truth Matched Jets", "\\Delta\\Phi", "DeltaPhiJetMatchTau", 100, -0.5, 4, "rad"},
      {"Invariant Mass of the truth Taus system0-50", "IMass", "TruthTauInvMass", 100, 0, 50, "GeV"},
      {"Invariant Mass of the truth Taus system0-200", "IMass", "TruthTauInvMass", 100, 0, 200, "GeV"},
      {"PT of Tau Leptons", "Tau PT", "Tau_PT", 800, -200, 200, "GeV"},
      {"Neutrino pt 1", "Neutrino Pt", "Neutrino_PT1", 800, -200, 200, "GeV"},
      {"Neutrino pt 2", "Neutrino Pt", "Neutrino_PT2", 800, -200, 200, "GeV"},
      {"Inv Mass without Neutrinos", "IMass", "TauJetInvMass", 50, 0, 200, "GeV"},
      {"Inv Mass with Neutrinos", "IMass", "TauJetInvMassWithNeutrino", 50, 0, 200, "GeV"},
      {"Number of B tag jets", "Number", "Jet_BNum", 10, 0, 10, "", true},
      {"Angle Difference between MET", "Angle", "AngleBetweenMET", 100, 0, 0.00001, "rad"},
      {"Delta Phi of Tau Leptons", "|Angle Difference|", "DeltaPhiGenTau", 1000, -0.5, 3.5, "rad"},
      {"Total Tau System PT", "PT", "TotalTauPT", 1000,-1000, 1000, "GeV"}
    }},
    {"SelectionCut",{
{"MET of Tau", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
      {"Delta R of Tau and Selected Tau Jets", "Delta R", "DeltaRJetTauSel", 400, 0, 1, ""},
      // {"Delta Phi of Tau Leptons", "\\Delta\\Phi", "DeltaPhiGenTau", 100, -0.5, 4, "rad"},
      {"Delta Phi of Truth Matched Jets", "\\Delta\\Phi", "TruthMatchDeltaPhi", 100, -0.5, 4, "rad"},
      // {"Invariant Mass of the truth Taus system0-50", "IMass", "TruthTauInvMass", 100, 0, 50, "GeV"},
      // {"Invariant Mass of the truth Taus system0-200", "IMass", "TruthTauInvMass", 100, 0, 200, "GeV"},
      {"PT of Tau Leptons", "Tau PT", "Tau_PT", 800, -200, 200, "GeV"},
      {"Neutrino pt 1", "Neutrino Pt", "NeutrinoPT1", 800, -200, 200, "GeV"},
      {"Neutrino pt 2", "Neutrino Pt", "NeutrinoPT2", 800, -200, 200, "GeV"},
      {"Inv Mass without Neutrinos", "IMass", "TauJetInvMass", 50, 0, 200, "GeV   "},
      {"Inv Mass with Neutrinos", "IMass", "TauJetInvMassWithNeutrino", 50, 0, 200, "GeV"},
      {"Number of B tag jets", "Number", "Jet_BNum", 10, 0, 10, "", true},
      {"Angle Difference between MET", "Angle", "AngleBetweenMET", 100, 0, 0.0001, "rad"},
      {"Angle Difference between MET ZOOMOUT", "Angle", "AngleBetweenMET", 1000, 0, 2, "rad"},

      // {"Delta Phi of Tau Leptons", "|Angle Difference|", "DeltaPhiGenTau", 1000, -0.5, 3.5, "rad"},
      // {"Total Tau System PT", "PT", "TotalTauPT", 1000,-1000, 1000, "GeV"}

    }},
    {
      "EE",{ 
      {"Inv Mass of ee system", "IMass", "eeInvMass", 200, 0, 200, "GeV"},
      {"MET", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
      {"eePhi", "Phi", "Electron.Phi", 100, -4, 4, "rad"},
      {"eeDeltaPhi", "Phi", "eeDeltaPhi", 100, -0.5, 4, "rad"},
      {"Total Jet PT", "PT", "Jet_TotPT", 200, 0, 200, "GeV"},
      {"Total ee PT", "PT", "ee_TotPT", 200, 0, 200, "GeV"},
      {"Number of Jets", "Number", "Jet_Num", 10, 0, 10, "", true},
      {"Number of b Jets", "Number", "Jet_BNum", 10, 0, 10, "", true},
      {"Tot PT Difference (jet - ee)", "PT", "PtDiff", 200, -1700, 700, "GeV"}
      }}
  };
};

#endif