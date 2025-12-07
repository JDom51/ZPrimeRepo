#ifndef DATA_BASE
#define DATA_BASE
#include"DataStructures.h"
#include"SelectionCuts.h"

using DataStructs::FrameAndData;
using DataStructs::HistInfo;
using SC::size_2;
using SC::size_0;
using SC::size_0_ui;

namespace DataBase{
  void NoCut(FrameAndData& fd);
  void Analyse(FrameAndData& fd);
  void TruthAnalysis(FrameAndData& fd);
  void GeneratorLevelTauAnalysis(FrameAndData& fd);
  void SelectionCutAnalysis(FrameAndData& fd);
  void EAnalysis(FrameAndData& fd);
  void MomentumTest(FrameAndData& fd);
  void SemiLeptonicAnalysis(DataStructs::FrameAndData& fd);
  void ppzllAnalysis(DataStructs::FrameAndData& fd);
  void ppztataAnalysis(DataStructs::FrameAndData& fd);
  // map<string, DataStructs::GenSelectionCut> sel_cuts{
  //   {"DeltaPhi2", SC::delta_phi_2}
  // };

  map<string, string> pre_path{
    {"cluster", ""},
    {"not", "~/Noether"}
  };
  map<string, float_t> weighting{
    {"Raw", 1},
    {"OpenData", 0.0028},
    {"run2", 0.011}
  };
  map<string, vector<string>> data_sets{
    {"pptata_exh", {"jdombrowski,pptata_exh_run_01", "jdombrowski,pptata_exh_run_02",
                      "jdombrowski,pptata_exh_run_03", "jdombrowski,pptata_exh_run_04",
                      "jdombrowski,pptata_exh_run_05","jdombrowski,pptata_exh_run_06",
                      "jdombrowski,pptata_exh_run_07","jdombrowski,pptata_exh_run_08",
                      "jdombrowski,pptata_exh_run_09","jdombrowski,pptata_exh_run_10",
                      "jdombrowski,pptata_exh_run_12","jdombrowski,pptata_exh_run_13",
                      "jdombrowski,pptata_exh_run_14","jdombrowski,pptata_exh_run_15",
                      "jdombrowski,pptata_exh_run_16","jdombrowski,pptata_exh_run_17",
                      "jdombrowski,pptata_exh_run_18","jdombrowski,pptata_exh_run_19",
                      "jdombrowski,pptata_exh_run_21","jdombrowski,pptata_exh_run_22",
                      "jdombrowski,pptata_exh_run_23","jdombrowski,pptata_exh_run_24",
                      "jdombrowski,pptata_exh_run_25"}},
    {"pptttata_exh", {"jdombrowski,pptttata_exh_run_01", "jdombrowski,pptttata_exh_run_02",
                      "jdombrowski,pptttata_exh_run_03", "jdombrowski,pptttata_exh_run_04",
                      "jdombrowski,pptttata_exh_run_05","jdombrowski,pptttata_exh_run_06",
                      "jdombrowski,pptttata_exh_run_07","jdombrowski,pptttata_exh_run_08",
                      "jdombrowski,pptttata_exh_run_09","jdombrowski,pptttata_exh_run_10"}},
    {"zonly_had_decay", {"jdombrowski,hadronic_top_decay_run_run_01", "jdombrowski,hadronic_top_decay_run_run_02",
                           "jdombrowski,hadronic_top_decay_run_run_03", "jdombrowski,hadronic_top_decay_run_run_04",
                           "jdombrowski,hadronic_top_decay_run_run_05","jdombrowski,hadronic_top_decay_run_run_06",
                           "jdombrowski,hadronic_top_decay_run_run_07","jdombrowski,hadronic_top_decay_run_run_08",
                           "jdombrowski,hadronic_top_decay_run_run_09","jdombrowski,hadronic_top_decay_run_run_10"}},
    {"ppzttee" , {"jdombrowski,ppzttee_run_01" ,"jdombrowski,ppzttee_run_02","jdombrowski,ppzttee_run_03",
            "jdombrowski,ppzttee_run_04","jdombrowski,ppzttee_run_05","jdombrowski,ppzttee_run_06",
            "jdombrowski,ppzttee_run_07","jdombrowski,ppzttee_run_08","jdombrowski,ppzttee_run_09",
            "jdombrowski,ppzttee_run_10"}},
    {"ppzee", {"jdombrowski,ppzee_run_01"}},
    {"ppzmumu",{"jdombrowski,ppzmumu_run_01"}},
    {"ppztata",{"jdombrowski,ppztata_run_01"}},
    {"ppzLL", {"jdombrowski,ppzee_run_01", "jdombrowski,ppzmumu_run_01","jdombrowski,ppztata_run_01"}},
    {"ppttee_exh", {"jdombrowski,ppttee_exh"}},
    {"ppttmumu_exh", {"jdombrowski,ppttmumu_exh"}},
  };
  // pre defined histograms
  map<string, vector<DataStructs::HistInfo>> pdhists{
    {"ppzLL", {
        {"Number of Jets", "Number", "Jet_Num", 20, 0, 20, "", true},
        {"Number Of BJets", "Number", "Jet_BNum", 10, 0, 10, "", true},
        {"MissingET", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
        {"TauJetPT", "PT", "Jet_TauTagPT", 50, 0, 200, "GeV"},
        {"JetPT", "PT", "Jet.PT", 50, 0, 200, "GeV"},
        {"TauPT", "PT", "Tau_PT", 50, 0, 200, "GeV"},
        {"Number of Gen Taus in event", "Number", "Tau_Num", 10, 0, 10, "", true},
        {"Electron PT", "PT", "Electron.PT", 50, 0, 200, "GeV"},
        {"Muon PT", "PT", "Muon.PT", 50, 0, 200, "GeV"},
        {"Electron Number", "Num", "Electron_size", 10, 0, 10, "", true},
        {"Muon Number", "Num", "Muon_size", 10, 0, 10, "", true},
        {"Gen Electron Number", "Num", "GenElectron_Num", 10, 0, 10, "", true},
        {"Gen Muon Number", "Num",  "GenMuon_Num", 10, 0, 10, "", true},
        {"GenElectron_Status", "Status", "GenElectron_Status", 100, 0, 100, "", true},
        {"GenMuon_Status", "Status", "GenMuon_Status", 100, 0, 100, "", true},
        {"GenElectron_IsPU", "PU", "GenElectron_IsPU", 3, 0, 3, "", true},
        {"GenMuon_IsPU", "PU", "GenMuon_IsPU", 3, 0, 3, "", true},
        {"GenMuonCut_Num", "Number", "GenMuonCut_Num", 10, 0, 10, "", true},
        {"GenElectronCut_Num", "Number", "GenElectronCut_Num", 10, 0, 10 , "", true},
        {"GenMuon_Eta", "\\Eta", "GenMuon_Eta", 100, -5, 5, ""},
        {"GenMuonCut_Eta", "\\Eta", "GenMuonCut_Eta", 100, -5, 5, ""},
        {"GenElectron_Eta", "\\Eta", "GenElectron_Eta", 100, -5, 5, ""},
        {"GenElectronCut_Eta", "\\Eta", "GenElectronCut_Eta", 100, -5, 5, ""},
        {"RecoMuonEta", "\\Eta", "Muon.Eta", 100, -5, 5, ""},
        {"RecoElectronEta", "\\Eta", "Electron.Eta", 100, -5, 5, ""},
        {"ElectronRecoSumPt", "PT", "Electron.SumPt", 50, 0, 200, "GeV"},
        {"ElectronRecoSumPtCharged", "PT", "Electron.SumPtCharged", 50, 0, 200, "GeV"},
        {"ElectronRecoSumPtNeutral", "PT", "Electron.SumPtCharged", 50, 0, 200, "GeV"},
        {"MuonRecoSumPt", "PT", "Muon.SumPt", 50, 0, 200, "GeV"},
        {"MuonRecoSumPtCharged", "PT", "Muon.SumPtCharged", 50, 0, 200, "GeV"},
        {"MuonRecoSumPtNeutral", "PT", "Muon.SumPtNeutral", 50, 0, 200, "GeV"},
        {"GenInvMass", "IMass", "GenInvMass", 40, 0, 200, "GeV"},
        {"A", "A", "AngleBetweenMET", 100, 0, 0.0002, "rad"},
        {"\\Delta\\phi between decay products", "\\Delta\\phi", "DeltaPhi", 50, 0, 4, "rad", true},
      }
    }
  };


    // dataset, a_func, broad_cuts, fname, Histograms
  map<string, DataStructs::DataBase> database{
    {"NoCut",
            {data_sets["pptttata_exh"], *NoCut, {}, "NoCut",
                {
                    {"Number of Jets", "Number", "Jet_Num", 20, 0, 20, "", true},
                    {"GenLep_Num", "Number", "GenLep_Num", 20, 0, 20, "", true},
                    {"Number of DTauTag", "Number", "Jet_DTauNum", 10, 0, 10, "", true},
                    {"MissingET", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
                    {"Number Of BJets", "Number", "Jet_BNum", 10, 0, 10, "", true},
                    {"TauJetPT", "PT", "Jet_TauTagPT", 200, 0, 200, "GeV"},
                    {"TauPT", "PT", "Tau_PT", 200, 0, 200, "GeV"},
                    {"Number of Tau Jets", "Number", "Jet_TauNum", 10, 0, 10, "", true},
                    {"Number of DeltaR indicies", "Number", "DeltaR_Num", 3, 0, 3, "", true},
                    {"DeltaR between Gen Taus", "\\Delta R", "GenTauDeltaR", 500, 0, 5, ""},
                    {"Delta R between obs taus", "\\Delta R", "TauJetDeltaR", 500, 0, 5, ""},
                    {"Number of Gen Taus in event", "Number", "Tau_Num", 10, 0, 10, "", true},
                    // {"Generator taus number", "Number", "Tau_Num", 10, 0, 10, "", true},
                    {"Electron PT", "PT", "Electron.PT", 200, 0, 200, "GeV"},
                    {"Muon PT", "PT", "Muon.PT", 200, 0, 200, "GeV"},
                    {"Electron Number", "Num", "Electron_size", 10, 0, 10, "", true},
                    {"Muon Number", "Num", "Muon_size", 10, 0, 10, "", true},
                    {"Gen Electron Number", "Num", "GenElectron_Num", 10, 0, 10, "", true},
                    {"Gen Muon Number", "Num",  "GenMuon_Num", 10, 0, 10, "", true},
                    {"GenElectron_Status", "Status", "GenElectron_Status", 100, 0, 100, "", true},
                    {"GenMuon_Status", "Status", "GenMuon_Status", 100, 0, 100, "", true},
                    {"GenElectron_IsPU", "PU", "GenElectron_IsPU", 3, 0, 3, "", true},
                    {"GenMuon_IsPU", "PU", "GenMuon_IsPU", 3, 0, 3, "", true},
                    {"GenMuonCut_Num", "Number", "GenMuonCut_Num", 10, 0, 10, "", true},
                    {"GenElectronCut_Num", "Number", "GenElectronCut_Num", 10, 0, 10 , "", true},
                    {"GenMuon_Eta", "\\Eta", "GenMuon_Eta", 100, -5, 5, ""},
                    {"GenMuonCut_Eta", "\\Eta", "GenMuonCut_Eta", 100, -5, 5, ""},
                    {"GenElectron_Eta", "\\Eta", "GenElectron_Eta", 100, -5, 5, ""},
                    {"GenElectronCut_Eta", "\\Eta", "GenElectronCut_Eta", 100, -5, 5, ""},
                    {"RecoMuonEta", "\\Eta", "Muon.Eta", 100, -5, 5, ""},
                    {"RecoElectronEta", "\\Eta", "Electron.Eta", 100, -5, 5, ""},
                    {"ElectronRecoSumPt", "PT", "Electron.SumPt", 200, 0, 200, "GeV"},
                    {"ElectronRecoSumPtCharged", "PT", "Electron.SumPtCharged", 200, 0, 200, "GeV"},
                    {"ElectronRecoSumPtNeutral", "PT", "Electron.SumPtCharged", 200, 0, 200, "GeV"},
                    {"MuonRecoSumPt", "PT", "Muon.SumPt", 200, 0, 200, "GeV"},
                    {"MuonRecoSumPtCharged", "PT", "Muon.SumPtCharged", 200, 0, 200, "GeV"},
                    {"MuonRecoSumPtNeutral", "PT", "Muon.SumPtNeutral", 200, 0, 200, "GeV"},
                    // {"ElectronGenSumPt", "PT", "GenElectron_SumPt", 200, 0, 200, "GeV"},
                    // {"ElectronGenSumPtCharged", "PT", "GenElectron_SumPtCharged", 200, 0, 200, "GeV"},
                    // {"ElectronGenSumPtNeutral", "PT", "GenElectron_SumPtCharged", 200, 0, 200, "GeV"},
                    // {"MuonGenSumPt", "PT", "GenMuon_SumPt", 200, 0, 200, "GeV"},
                    // {"MuonGenSumPtCharged", "PT", "GenMuon_SumPtCharged", 200, 0, 200, "GeV"},
                    // {"MuonGenSumPtNeutral", "PT", "GenMuon_SumPtNeutral", 200, 0, 200, "GeV"},


                }
            } 
    },
    {"RecoLevel",
            {data_sets["pptttata_exh"], *Analyse, {{size_2, {"Jet_TauTagIndicies"}}}, "Ztata Reco",
                {
                    {"MET", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
                    {"Delta Phi Between two tau tagged jets", "\\Delta \\phi", "TauTagJetDeltaPhi", 100, 0, 4, "rad"},
                    {"Inv Mass without Neutrinos", "IMass", "TauJetInvMass", 200, 0, 200, "GeV"},
                    {"Neutrino PT1", "PT", "NeutrinoPT1", 200, 0, 200, "GeV"},
                    {"Neutrino PT2", "PT", "NeutrinoPT2", 200, 0, 200, "GeV"},
                    {"Inv Mass with Neutrinos", "IMass", "TauJetInvMassWithNeutrino", 200, 0, 200, "GeV"},
                    {"Inv Mass with Neuitrinos 100 bins", "IMass", "TauJetInvMassWithNeutrino", 100, 0, 200, "GeV"},
                    {"Inv Mass with Neuitrinos 50 bins", "IMass", "TauJetInvMassWithNeutrino", 50, 0, 200, "GeV"},
                    {"Inv Mass with Neuitrinos 20 bins", "IMass", "TauJetInvMassWithNeutrino", 20, 0, 200, "GeV"},
                    {"Number of Jets", "Num", "Jet_Num", 20, 0, 20, "", true},
                    {"Number of B Jets", "Num", "Jet_BNum", 10, 0, 10, "", true},
                    {"Angle diff between met and jets", "$\\Delta \\phi$", "AngleBetweenMET", 200, 0, 4, "rad"},
                        // {"Delta R between obs taus", "\\Delta R", "TauJetDeltaR", 500, 0, 5, ""}
                }
            }
    },
    {"Gen+RecoLevel",
            {data_sets["pptttata_exh"], *TruthAnalysis, {{size_2, {"TruthMatchIndicies"}}}, "Ztata      Gen+Reco",
                {
                    {"Number of Jets", "Number", "Jet_Num", 20, 0, 20, "", true},
                    {"MissingET", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
                    {"Number Of BJets", "Number", "Jet_BNum", 10, 0, 10, "", true},
                    {"TauJetPT", "PT", "Jet_TauTagPT", 200, 0, 200, "GeV"},
                    {"TauPT", "PT", "Tau_PT", 200, 0, 200, "GeV"},
                    {"Number of Tau Jets", "Number", "Jet_TauNum", 10, 0, 10, "", true},
                    {"Number of DeltaR indicies", "Number", "DeltaR_Num", 3, 0, 3, "", true},
                    {"DeltaR between Gen Taus", "\\Delta R", "GenTauDeltaR", 500, 0, 5, ""},
                    {"Delta R between obs taus", "\\Delta R", "TauJetDeltaR", 500, 0, 5, ""},
                    {"Number of Gen Taus in event", "Number", "Tau_Num", 10, 0, 10, "", true},
                    {"Generator taus number", "Number", "Tau_Num", 10, 0, 10, "", true},
                    {"Electron PT", "PT", "Electron.PT", 200, 0, 200, "GeV"},
                    {"Muon PT", "PT", "Muon.PT", 200, 0, 200, "GeV"},
                    {"Electron Number", "Num", "Electron_size", 10, 0, 10, "", true},
                    {"Muon Number", "Num", "Muon_size", 10, 0, 10, "", true},
                    {"Gen Electron Number", "Num", "GenElectron_Num", 10, 0, 10, "", true},
                    {"Gen Muon Number", "Num",  "GenMuon_Num", 10, 0, 10, "", true},
                    {"GenElectron_Status", "Status", "GenElectron_Status", 100, 0, 100, "", true},
                    {"GenMuon_Status", "Status", "GenMuon_Status", 100, 0, 100, "", true},
                    {"GenElectron_IsPU", "PU", "GenElectron_IsPU", 3, 0, 3, "", true},
                    {"GenMuon_IsPU", "PU", "GenMuon_IsPU", 3, 0, 3, "", true},
                    {"GenMuonCut_Num", "Number", "GenMuonCut_Num", 10, 0, 10, "", true},
                    {"GenElectronCut_Num", "Number", "GenElectronCut_Num", 10, 0, 10 , "", true},
                    {"GenMuon_Eta", "\\Eta", "GenMuon_Eta", 100, -5, 5, ""},
                    {"GenMuonCut_Eta", "\\Eta", "GenMuonCut_Eta", 100, -5, 5, ""},
                    {"GenElectron_Eta", "\\Eta", "GenElectron_Eta", 100, -5, 5, ""},
                    {"GenElectronCut_Eta", "\\Eta", "GenElectronCut_Eta", 100, -5, 5, ""},
                    {"RecoMuonEta", "\\Eta", "Muon.Eta", 100, -5, 5, ""},
                    {"RecoElectronEta", "\\Eta", "Electron.Eta", 100, -5, 5, ""},
                    {"ElectronRecoSumPt", "PT", "Electron.SumPt", 200, 0, 200, "GeV"},
                    {"ElectronRecoSumPtCharged", "PT", "Electron.SumPtCharged", 200, 0, 200, "GeV"},
                    {"ElectronRecoSumPtNeutral", "PT", "Electron.SumPtCharged", 200, 0, 200, "GeV"},
                    {"MuonRecoSumPt", "PT", "Muon.SumPt", 200, 0, 200, "GeV"},
                    {"MuonRecoSumPtCharged", "PT", "Muon.SumPtCharged", 200, 0, 200, "GeV"},
                    {"MuonRecoSumPtNeutral", "PT", "Muon.SumPtNeutral", 200, 0, 200, "GeV"},
                    // {"ElectronGenSumPt", "PT", "GenElectron_SumPt", 200, 0, 200, "GeV"},
                    // {"ElectronGenSumPtCharged", "PT", "GenElectron_SumPtCharged", 200, 0, 200, "GeV"},
                    // {"ElectronGenSumPtNeutral", "PT", "GenElectron_SumPtCharged", 200, 0, 200, "GeV"},
                    // {"MuonGenSumPt", "PT", "GenMuon_SumPt", 200, 0, 200, "GeV"},
                    // {"MuonGenSumPtCharged", "PT", "GenMuon_SumPtCharged", 200, 0, 200, "GeV"},
                    // {"MuonGenSumPtNeutral", "PT", "GenMuon_SumPtNeutral", 200, 0, 200, "GeV"},
                }
            }
    },
    {"GenLevel",
            {data_sets["pptttata_exh"], *GeneratorLevelTauAnalysis, {{size_2, {"Particle_TauIndicies"}}}, "Ztata Gen",
                {
                    {"Inv Mass of Generator System", "IMass", "GenTauInvMass", 200, 0, 200, "GeV"}
                }
            }
    },
    {"SemiLeptnoicAnalysis",
            {data_sets["pptttata_exh"], *SemiLeptonicAnalysis, {}, "SemiLeptonic",
              {
                // {"GenMuonCut Number", "Number", "MuonCut_Num", 10, 0, 10, "", true},
                // {"GenElectronCut Number", "Number", "Electron_Num", 10, 0, 10, "", true},
                {"Lep_Num", "Number", "Lep_Num", 5, 0, 5, "", true},
                {"Jet_TauTag Number", "Number", "Jet_TauNum", 10, 0, 10, "", true},
                {"Lep PT", "PT", "Lep_PT", 200, 0, 200, "GeV"},
                {"Jet_TauTagPT", "PT", "Jet_TauTagPT", 200, 0, 200, "GeV"},
                {"Jet_TauTagEta", "\\eta", "Jet_TauTagEta", 200, -5, 5, ""},
                {"Lep_Eta", "\\eta", "Lep_Eta", 200, -5, 5, ""},
                {"delta_phi_between_lep and tau", "\\Delta \\phi", "SemiLepDeltaPhi", 100, 0, 4, "rad", true},
                {"SemiLep invmass without neutrino pt", "IMass", "SemiLepInvMass",  15, 0, 200, "GeV"},
                {"Lep_IsolationVar", "IsolationVarr", "Lep_IsolationVar", 10, 0, 1, ""},
                {"TransverseMass", "M_T", "TransverseMass", 50, 0, 200, "GeV"},
                {"Semi Lep Inv Mass with Neutrino PT", "IMass", "SemiLepInvMassWithNeutrino", 20, 0, 200, "GeV"},
                {"Angle Difference between MET", "Angle", "AngleBetweenMET", 100, 0, 0.0001, "rad"},
                {"Angle Difference between MET ZOOMOUT", "Angle", "AngleBetweenMET", 1000, 0, 2, "rad"},
                {"Neutrino pt 1", "Neutrino Pt", "NeutrinoPT1", 800, -200, 200, "GeV"},
                {"Neutrino pt 2", "Neutrino Pt", "NeutrinoPT2", 800, -200, 200, "GeV"},
                {"Sum PT plot for lepton", "PT", "Lep_SumPT", 200, 0, 100, "GeV"},
                {"MET", "\\not{E_T}", "MissingET.MET", 200, 0, 200, "GeV"},
                {"Jet_BNum", "Number", "Jet_BNum", 10, 0, 10, "", true}, 
                {"Jet_Num", "Number", "Jet_Num", 20 , 0, 20, "", true},
                {"MET", "MET", "MissingET.MET", 100, 0, 200, "GeV"},
                {"Omega", "\\Omega", "Omega", 100, -2, 2, ""},

              }
            }
    },
    {"SelectionCut",
            {data_sets["pptttata_exh"], *SelectionCutAnalysis, {{size_2, {"TruthMatchIndicies"}}, {size_2, {"Particle_TauIndicies"}}}, "Selcut",
                {
                    {"MET of Tau", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
                    // {"Transverse Mass", "M_T", "Jet_TransverseMass", 200, 0, 200, "GeV"},
                    // {"Transverse Mass", "M_T", ""}
                    // {"Delta R of Tau and Selected Tau Jets", "Delta R", "DeltaRJetTauSel", 400, 0, 1, ""},
                    // {"Delta Phi of Tau Leptons", "\\Delta\\Phi", "DeltaPhiGenTau", 100, -0.5, 4, "rad"},
                    {"Delta R of Truth Matching", "\\Delta R", "TruthMatchedDeltaR", 200, 0, 2, ""},
                    // {"Delta R of Non Truth Matching", "\\Delta R", "NTruthMatchedDeltaR", 200, 0, 2, ""},
                    {"Number of Truth Matched Taus", "Number", "Jet_DTauNum", 10, 0, 10, "", true},
                    {"Delta Phi of Truth Matched Jets", "\\Delta\\Phi", "TruthMatchDeltaPhi", 100, -0.5, 4, "rad"},
                    // {"Invariant Mass of the truth Taus system0-50", "IMass", "TruthTauInvMass", 100, 0, 50, "GeV"},
                    // {"Invariant Mass of the truth Taus system0-200", "IMass", "TruthTauInvMass", 100, 0, 200, "GeV"},
                    {"PT of Tau Leptons", "Tau PT", "Tau_PT", 200, 0, 400, "GeV"},
                    {"Neutrino pt 1", "Neutrino Pt", "NeutrinoPT1", 200, -200, 200, "GeV"},
                    {"Neutrino pt 2", "Neutrino Pt", "NeutrinoPT2", 200, -200, 200, "GeV"},
                    {"Inv Mass without Neutrinos", "IMass", "TauJetInvMass", 20, 0, 200, "GeV"},
                    {"Inv Mass with Neutrinos", "IMass", "TauJetInvMassWithNeutrino", 20, 0, 200, "GeV"},
                    {"Number of B tag jets", "Number", "Jet_BNum", 10, 0, 10, "", true},
                    {"Angle Difference between MET", "Angle", "AngleBetweenMET", 100, 0, 0.0001, "rad"},
                    {"Angle Difference between MET ZOOMOUT", "Angle", "AngleBetweenMET", 100, 0, 2, "rad"},
                    {"Inv Mass of Generator System", "IMass", "GenTauInvMass", 200, 0, 200, "GeV"},
                    // {"Inv Mass of TruthMatchGenericJet without Neutrinos", "IMass", "TauJetInvMassTRUTHJET", 50, 0, 200, "GeV   "},
                    // {"Inv Mass of TruthMatchGenericJet with Neutrinos", "IMass", "TauJetInvMassWithNeutrinoTRUTHJET", 50, 0, 200, "GeV"},
                    // {"Angle Difference between MET TRUTHGENERICJET", "Angle", "AngleBetweenMETTRUTHJET", 100, 0, 0.0001, "rad"},
                    // {"Delta Phi of Truth Matched Jets TRUTHGENERICJET", "\\Delta\\Phi", "TruthMatchDeltaPhiTRUTHJET", 100, -0.5, 4, "rad"},
                    {"Tau Obs PT", "PT", "Jet_DTauTagNeutrinoPT", 300, 0, 400, "GeV"},
                    {"Alpha of Gen Tau", "\\alpha", "Tau_Alpha", 100, -4, 4, "rad"},
                    {"Alpha of Truth Matched Tau", "\\alpha", "Jet_DTauTagAlpha", 100, -4, 4, "rad"},
                    {"DeltaR between Gen Taus", "\\Delta R", "GenTauDeltaR", 200, 0, 5, ""},
                    {"Delta R between obs taus", "\\Delta R", "TauJetDeltaR", 200, 0, 5, ""},
                    // {"PT of truthNeutrino1", "PT", "TruthNeutrinoPT2", 300, -100, 200, "GeV"},
                    // {"Pt of truthNeutrino2", "Pt", "TruthNeutrinoPT1", 300, -100, 200, "GeV"},
                    // {"PT of tau jet with truth neutrino included", "PT", "Jet_DTauTagTruthNeutrinoPT", 200, 0, 200, "GeV"},
                    // {"angle difference between the truth met and the jets", "Angle Difference", "AngleBetweenTruthMET", 1000, 0, 1, "rad"},
                    // {"Invariant Mass of Taus with Truth Neutrino PT included", "IMass", "TauJetInvMassWithTruthNeutrino", 100, 0, 200, "GeV/c^2"},
                    // {"Delta R between tau tag jets", "\\Delta R", ""}
                    // {"Difference between the Truth and Reco MET", "\\not{E_T\\mathrm{truth}} - \\not{E_T\\mathrm{reco}}", "DifferenceBetweenTruthRecoMET", 100, -100, 100, "GeV"},
                    // {"NormalisedDifferenceBetweenTruthRecoMET", "\\frac{\\not{E_T\\mathrm{truth}} - \\not{E_T\\mathrm{reco}}}{\\not{E_{T\\mathrm{truth}}}}", "NormalisedDifferenceBetweenTruthRecoMET", 100, -3, 3, ""},
                    // {"Jet_NDTauTagPT", "PT", "Jet_NDTauTagPT", 200, 0, 200, "GeV"},
                    // {"Jet_NDTauTagphi", "\\phi", "Jet_NDTauTagPhi", 200, 0, 200, "rad"},
                    // {"Jet_NDTauTagEta", "\\eta", "Jet_NDTauTagEta", 200, 0, 200, ""},
                    // {"Jet_NDTauTagMass", "IMASS", "Jet_NDTauTagMass", 200, 0, 200, "GeV"},
                    // {"Jet_DTauTagMass", "IMASS", "Jet_DTauTagMass", 200, 0, 200, "GeV"},
                    {"ElectronRecoSumPt", "PT", "Electron.SumPt", 200, 0, 200, "GeV"},
                    {"ElectronRecoSumPtCharged", "PT", "Electron.SumPtCharged", 200, 0, 200, "GeV"},
                    {"ElectronRecoSumPtNeutral", "PT", "Electron.SumPtCharged", 200, 0, 200, "GeV"},
                    {"MuonRecoSumPt", "PT", "Muon.SumPt", 200, 0, 200, "GeV"},
                    {"MuonRecoSumPtCharged", "PT", "Muon.SumPtCharged", 200, 0, 200, "GeV"},
                    {"MuonRecoSumPtNeutral", "PT", "Muon.SumPtNeutral", 200, 0, 200, "GeV"},
                    {"Number of Jets", "Number", "Jet_Num", 20, 0, 20, "", true},

                    // {"Delta Phi of Tau Leptons", "|Angle Difference|", "DeltaPhiGenTau", 1000, -0.5, 3.5, "rad"},
                    // {"Total Tau System PT", "PT", "TotalTauPT", 1000,-1000, 1000, "GeV"}
                }
            }
    },

    {"SelectionCutTRUTHJET",
            {data_sets["pptttata_exh"], *SelectionCutAnalysis, {}, "SelcutTRUTHJET",
                {
                    {"MET of Tau", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
                    // {"Delta R of Tau and Selected Tau Jets", "Delta R", "DeltaRJetTauSel", 400, 0, 1, ""},
                    // {"Delta Phi of Tau Leptons", "\\Delta\\Phi", "DeltaPhiGenTau", 100, -0.5, 4, "rad"},
                    {"Delta Phi of Truth Matched Jets", "\\Delta\\Phi", "TruthMatchDeltaPhi", 100, -0.5, 4, "rad"},
                    // {"Invariant Mass of the truth Taus system0-50", "IMass", "TruthTauInvMass", 100, 0, 50, "GeV"},
                    // {"Invariant Mass of the truth Taus system0-200", "IMass", "TruthTauInvMass", 100, 0, 200, "GeV"},
                    {"PT of Tau Leptons", "Tau PT", "Tau_PT", 800, 0, 400, "GeV"},
                    {"Neutrino pt 1", "Neutrino Pt", "NeutrinoPT1", 800, -200, 200, "GeV"},
                    {"Neutrino pt 2", "Neutrino Pt", "NeutrinoPT2", 800, -200, 200, "GeV"},
                    {"Inv Mass without Neutrinos", "IMass", "TauJetInvMass", 50, 0, 200, "GeV   "},
                    {"Inv Mass with Neutrinos", "IMass", "TauJetInvMassWithNeutrino", 50, 0, 200, "GeV"},
                    {"Number of B tag jets", "Number", "Jet_BNum", 10, 0, 10, "", true},
                    {"Angle Difference between MET", "Angle", "AngleBetweenMET", 100, 0, 0.0001, "rad"},
                    {"Angle Difference between MET ZOOMOUT", "Angle", "AngleBetweenMET", 1000, 0, 2, "rad"},
                    {"Inv Mass of Generator System", "IMass", "GenTauInvMass", 200, 0, 200, "GeV"},
                    {"Inv Mass of TruthMatchGenericJet without Neutrinos", "IMass", "TauJetInvMassTRUTHJET", 50, 0, 200, "GeV   "},
                    {"Inv Mass of TruthMatchGenericJet with Neutrinos", "IMass", "TauJetInvMassWithNeutrinoTRUTHJET", 50, 0, 200, "GeV"},
                    {"Angle Difference between MET TRUTHGENERICJET", "Angle", "AngleBetweenMETTRUTHJET", 100, 0, 0.0001, "rad"},
                    {"Delta Phi of Truth Matched Jets TRUTHGENERICJET", "\\Delta\\Phi", "TruthMatchDeltaPhiTRUTHJET", 100, -0.5, 4, "rad"},
                    {"Tau Obs PT", "PT", "Jet_DTauTagNeutrinoPT", 800, 0, 400, "GeV"}


                    // {"Delta Phi of Tau Leptons", "|Angle Difference|", "DeltaPhiGenTau", 1000, -0.5, 3.5, "rad"},
                    // {"Total Tau System PT", "PT", "TotalTauPT", 1000,-1000, 1000, "GeV"}
                }
            }
    },
    {"SelectionCutZOnly",
            {data_sets["zonly_had_decay"], *SelectionCutAnalysis,  {{size_2, {"DeltaRIndicies"}}}, "SelCutZOnly",
                {
                    {"MET of Tau", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
                    // {"Delta R of Tau and Selected Tau Jets", "Delta R", "DeltaRJetTauSel", 400, 0, 1, ""},
                    // {"Delta Phi of Tau Leptons", "\\Delta\\Phi", "DeltaPhiGenTau", 100, -0.5, 4, "rad"},
                    {"Delta Phi of Truth Matched Jets", "\\Delta\\Phi", "TruthMatchDeltaPhi", 100, -0.5, 4, "rad"},
                    // {"Invariant Mass of the truth Taus system0-50", "IMass", "TruthTauInvMass", 100, 0, 50, "GeV"},
                    // {"Invariant Mass of the truth Taus system0-200", "IMass", "TruthTauInvMass", 100, 0, 200, "GeV"},
                    {"PT of Tau Leptons", "Tau PT", "Tau_PT", 800, 0, 400, "GeV"},
                    {"Neutrino pt 1", "Neutrino Pt", "NeutrinoPT1", 800, -200, 200, "GeV"},
                    {"Neutrino pt 2", "Neutrino Pt", "NeutrinoPT2", 800, -200, 200, "GeV"},
                    {"Inv Mass without Neutrinos", "IMass", "TauJetInvMass", 50, 0, 200, "GeV   "},
                    {"Inv Mass with Neutrinos", "IMass", "TauJetInvMassWithNeutrino", 50, 0, 200, "GeV"},
                    {"Number of B tag jets", "Number", "Jet_BNum", 10, 0, 10, "", true},
                    {"Angle Difference between MET", "Angle", "AngleBetweenMET", 100, 0, 0.0001, "rad"},
                    {"Angle Difference between MET ZOOMOUT", "Angle", "AngleBetweenMET", 1000, 0, 2, "rad"},
                    {"Inv Mass of Generator System", "IMass", "GenTauInvMass", 200, 0, 200, "GeV"},
                    // {"Inv Mass of TruthMatchGenericJet without Neutrinos", "IMass", "TauJetInvMassTRUTHJaaaET", 50, 0, 200, "GeV   "},
                    // {"Inv Mass of TruthMatchGenericJet with Neutrinos", "IMass", "TauJetInvMassWithNeutrinoTRUTHJET", 50, 0, 200, "GeV"},
                    // {"Angle Difference between MET TRUTHGENERICJET", "Angle", "AngleBetweenMETTRUTHJET", 100, 0, 0.0001, "rad"},
                    // {"Delta Phi of Truth Matched Jets TRUTHGENERICJET", "\\Delta\\Phi", "TruthMatchDeltaPhiTRUTHJET", 100, -0.5, 4, "rad"},
                    {"Tau Obs PT", "PT", "Jet_DTauTagNeutrinoPT", 800, 0, 400, "GeV"},
                    //  {"Alpha of Gen Tau", "\\alpha", "Tau_Alpha", 100, -4, 4, "rad"},
                    //  {"Alpha of Truth Matched Tau", "\\alpha", "Jet_DTauTagAlpha", 100, -4, 4, "rad"},
                    {"DeltaR between Gen Taus", "\\Delta R", "GenTauDeltaR", 500, 0, 5, ""},
                    {"Delta R between obs taus", "\\Delta R", "TauJetDeltaR", 500, 0, 5, ""}


                    // {"Delta Phi of Tau Leptons", "|Angle Difference|", "DeltaPhiGenTau", 1000, -0.5, 3.5, "rad"},
                    // {"Total Tau System PT", "PT", "TotalTauPT", 1000,-1000, 1000, "GeV"}

                }
            }
    },
    {"EE",
            {data_sets["ppzttee"], *EAnalysis, {}, "Zee, n_e=2",
                {
                    {"Inv Mass of ee system", "IMass", "eeInvMass", 200, 0, 200, "GeV"},
                    {"MET", "MET", "MissingET.MET", 200, 0, 200, "GeV"},
                    {"eePhi", "Phi", "Electron.Phi", 100, -4, 4, "rad"},
                    {"eeDeltaPhi", "Phi", "eeDeltaPhi", 100, -0.5, 4, "rad"},
                    {"Total Jet PT", "PT", "Jet_TotPT", 200, 0, 200, "GeV"},
                    {"Total ee PT", "PT", "ee_TotPT", 200, 0, 200, "GeV"},
                    {"Number of Jets", "Number", "Jet_Num", 20, 0, 20, "", true},
                    {"Number of b Jets", "Number", "Jet_BNum", 10, 0, 10, "", true},
                    {"Tot PT Difference (jet - ee)", "PT", "PtDiff", 200, -1700, 700, "GeV"}
                }
            }
    },
    {"eeatata",
            {{"jdombrowski,eeatata_run_01"}, *MomentumTest, {}, "eeatata",
                {
                    {"Tau PT", "PT", "Tau_PT", 400, 0, 800, "GeV"},
                    {"Tau Eta", "\\Eta", "Tau_Eta", 100, -5, 5, ""},
                    {"Tau Phi", "\\phi", "Tau_Phi", 100, -4, 4, "rad"},
                    {"Inv Mass of Tau", "IMass", "GenTauInvMass", 200, 0, 200, "GeV"}
                }
            }
    },
    {"eeztata",
            {{"jdombrowski,eeztata_run_01"}, *MomentumTest, {}, "eeztata",
                {
                    {"Tau PT", "PT", "Tau_PT", 400, 0, 800, "GeV"},
                    {"Tau Eta", "\\Eta", "Tau_Eta", 100, -5, 5, ""},
                    {"Tau Phi", "\\phi", "Tau_Phi", 100, -4, 4, "rad"},
                    {"Inv Mass of Tau", "IMass", "GenTauInvMass", 200, 0, 200, "GeV"}
                }
            }
    },
    {"ppzee", {data_sets.at("ppzee"), *ppzllAnalysis, {{size_2, {"Particle_ElectronIndicies_PTCUT"}}, {size_0_ui, {"Particle_MuonIndicies_PTCUT"}}}, "ppzee", pdhists.at("ppzLL")}},
    {"ppzmumu", {data_sets.at("ppzmumu"), *ppzllAnalysis, {{size_2, {"Particle_MuonIndicies_PTCUT"}}, {size_0_ui, {"Particle_ElectronIndicies_PTCUT"}}}, "ppzmumu", pdhists.at("ppzLL")}},
    {"ppztata", {data_sets.at("ppztata"), *ppztataAnalysis, {{size_2, {"Particle_TauIndicies"}}}, "ppztata", pdhists.at("ppzLL")}},
    {"ppzttee_exh", {data_sets.at("ppzttee_exh"), *ppzllAnalysis, {{size_2, {"Particle_ElectronIndicies_PTCUT"}}, {size_0_ui, {"Particle_MuonIndicies_PTCUT"}}, "ppzttee_exh", pdhists.at("ppzLL")}}},
    {"ppzttmumu_exh", {data_sets.at("ppzttmumu_exh"), *ppzllAnalysis, {{size_2, {"Particle_MuonIndicies_PTCUT"}}, {size_0_ui, {"Particle_ElectronIndicies_PTCUT"}}, "ppzttmumu_exh", pdhists.at("ppzLL")}}},
    {"ppztttata_exh", {data_sets.at("ppztttata_exh"), *ppztataAnalysis, {{size_2, {"Particle_TauIndicies"}}}, "ppztttata_exh", pdhists.at("ppzLL")}}

  };

};

#endif