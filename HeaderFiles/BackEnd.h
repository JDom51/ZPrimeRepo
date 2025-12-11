#ifndef BACKEND
#define BACKEND
#include<sstream>
#include<filesystem>
#include<fstream>
#include<iostream>
#include<memory>
#include<chrono>
#include<future>
#include<THStack.h>
#include<TH2D.h>
#include<TLegend.h>

#include"DataStructures.h"
#include"LFuncs.h"
#include"DataBase.h"

using DataStructs::FrameAndData;
using DataStructs::SelectionCut;
using DataStructs::HistInfo;
using std::cout;


namespace BackEnd
{
  // [Previous functions remain unchanged - load_files, initial_cut, etc.]
  vector<string> load_files(string& analysis_mode, string& cluster_mode)
  {
    vector<string> files{};
    string file{};
    for(auto file :  DataBase::database[analysis_mode].ifiles)
    {
      size_t delimiter{file.find(",")};
      cout<<"/gluster/data/atlas/" + file.substr(0, delimiter) + "/DATA/" + file.substr(delimiter + 1) + ".root\n";
      files.push_back(DataBase::pre_path.at(cluster_mode) + "/gluster/data/atlas/" + file.substr(0, delimiter) + "/DATA/" + file.substr(delimiter + 1) + ".root");
    }
    return files;
  }
  
  void initial_cut(FrameAndData& fd, vector<DataStructs::SelectionCut>& selcuts){
    for(auto sel_cut : selcuts){
      fd.Filter(sel_cut.selection_cut, sel_cut.variables, sel_cut.variables[0]);
    }
  }

  // Structure to hold different histogram types
  struct HistogramCollection {
    vector<ROOT::RDF::RResultPtr<TH1D>> hist1d;
    vector<ROOT::RDF::RResultPtr<TH2D>> hist2d;
    vector<vector<ROOT::RDF::RResultPtr<TH1D>>> stacked_hists; // Each element is a vector of histograms for one stack
    vector<ROOT::RDF::RResultPtr<TH1D>> data_hists; // Data histograms to overlay as crosses
  };

  void save_histograms_extended(std::string outputFileName, 
                                 HistogramCollection& hists,
                                 const std::vector<HistInfo>& hist_info)
  {
    // Output ROOT file
    std::filesystem::path dir = "/gluster/data/atlas/jdombrowski/Histograms";
    std::filesystem::create_directories(dir);
    std::string file_path = (dir / (outputFileName + ".root")).string();
    TFile out_file(file_path.c_str(), "RECREATE");

    auto start = std::chrono::high_resolution_clock::now();
    cout << "Start: \n";
    
    // Materialize all 1D histograms
    for (auto& hptr : hists.hist1d) {
        hptr.GetValue();
    }
    
    // Materialize all 2D histograms
    for (auto& hptr : hists.hist2d) {
        hptr.GetValue();
    }
    
    // Materialize all stacked histograms
    for (auto& stack : hists.stacked_hists) {
        for (auto& hptr : stack) {
            hptr.GetValue();
        }
    }
    
    // Materialize all data histograms
    for (auto& hptr : hists.data_hists) {
        hptr.GetValue();
    }

    auto mid = std::chrono::high_resolution_clock::now();
    cout << "Mid: " << std::chrono::duration<double>(mid-start).count() << "\n";
    out_file.cd();

    // Write histograms based on their mode
    size_t idx_1d = 0;
    size_t idx_2d = 0;
    size_t idx_stacked = 0;
    size_t idx_data = 0;
    
    for (size_t i = 0; i < hist_info.size(); ++i) {
        const auto& h = hist_info[i];
        
        if (h.mode == 0) {
            // Normal 1D histogram
            TH1D* hist = hists.hist1d[idx_1d].GetPtr();
            hist->SetName(h.name.c_str());
            hist->GetXaxis()->SetTitle((h.x_axis + " (" + h.units + ")").c_str());
            hist->GetYaxis()->SetTitle("Events");
            hist->Write(h.name.c_str());
            idx_1d++;
        } 
        else if (h.mode == 1) {
            // Stacked histogram
            THStack* stack = new THStack(h.name.c_str(), h.name.c_str());
            TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            
            // Define colors for stacking (you can customize these)
            const int colors[] = {kRed, kBlue, kGreen+2, kOrange, kMagenta, kCyan, kYellow+2, kViolet};
            const int num_colors = 8;
            
            // Determine starting index based on data_vs_mc flag
            size_t start_idx = h.data_vs_mc ? 1 : 0;
            
            // Add MC histograms to stack
            for (size_t j = start_idx; j < hists.stacked_hists[idx_stacked].size(); ++j) {
                TH1D* hist = hists.stacked_hists[idx_stacked][j].GetPtr();
                string hist_name = h.name + "_" + h.columns[j];
                hist->SetName(hist_name.c_str());
                int color_idx = h.data_vs_mc ? (j - 1) : j; // Adjust color index if first is data
                hist->SetFillColor(colors[color_idx % num_colors]);
                hist->SetLineColor(colors[color_idx % num_colors]);
                stack->Add(hist);
                legend->AddEntry(hist, h.columns[j].c_str(), "f");
                hist->Write(); // Write individual histograms as well
            }
            
            stack->Write(h.name.c_str());
            
            // Handle data histogram if data_vs_mc is true
            TH1D* data_hist = nullptr;
            if (h.data_vs_mc) {
                data_hist = hists.data_hists[idx_data].GetPtr();
                string data_hist_name = h.name + "_data";
                data_hist->SetName(data_hist_name.c_str());
                data_hist->SetMarkerStyle(20); // Filled circle
                data_hist->SetMarkerSize(0.8);
                data_hist->SetMarkerColor(kBlack);
                data_hist->SetLineColor(kBlack);
                legend->AddEntry(data_hist, h.columns[0].c_str(), "p"); // "p" for points
                data_hist->Write();
                idx_data++;
            }
            
            legend->Write((h.name + "_legend").c_str());
            
            // Create and write a canvas with the stack and data drawn
            TCanvas* c = new TCanvas((h.name + "_canvas").c_str(), h.name.c_str(), 800, 600);
            
            // Draw stack first
            stack->Draw("HIST");
            stack->GetXaxis()->SetTitle((h.x_axis + " (" + h.units + ")").c_str());
            stack->GetYaxis()->SetTitle("Events");
            
            // Overlay data as points with error bars if present
            if (h.data_vs_mc && data_hist) {
                data_hist->Draw("E SAME"); // E for error bars, SAME to overlay
            }
            
            legend->Draw();
            c->Write();
            delete c;
            
            idx_stacked++;
        }
        else if (h.mode == 2) {
            // 2D histogram
            TH2D* hist = hists.hist2d[idx_2d].GetPtr();
            hist->SetName(h.name.c_str());
            hist->GetXaxis()->SetTitle((h.x_axis + " (" + h.units + ")").c_str());
            hist->GetYaxis()->SetTitle((h.y_axis + " (" + h.units_y + ")").c_str());
            hist->GetZaxis()->SetTitle("Events");
            hist->Write(h.name.c_str());
            idx_2d++;
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    cout << "end: " << std::chrono::duration<double>(end-mid).count() << "\n";
    out_file.Close();
    
    std::cout << "Computation: " 
              << std::chrono::duration<double>(mid-start).count() << "s\n";
    std::cout << "Writing: " 
              << std::chrono::duration<double>(end-mid).count() << "s\n";
    std::cout << "Saved " << (idx_1d + idx_2d + idx_stacked) 
              << " histograms to " << file_path << std::endl;
  }

  void save_hist_info(string outputFileName, const vector<HistInfo>& histograms, string& cluster_mode){
    ROOT::EnableImplicitMT();
    std::stringstream complete_file_name{}; 
    complete_file_name << DataBase::pre_path.at(cluster_mode) 
                       << "/gluster/data/atlas/jdombrowski/HistogramInfo/" 
                       << outputFileName << ".txt";
    std::ofstream histograms_file{complete_file_name.str()};
    std::stringstream save_stringstream{};
    save_stringstream << "{";
    
    for(auto h : histograms)
    {
      save_stringstream << "{" << h.name << ", " << h.x_axis << ", ";
      
      // Save column info based on mode
      if (h.mode == 0) {
        save_stringstream << h.columns[0];
      } else if (h.mode == 1) {
        save_stringstream << "[";
        for (size_t i = 0; i < h.columns.size(); ++i) {
          save_stringstream << h.columns[i];
          if (i < h.columns.size() - 1) save_stringstream << ";";
        }
        save_stringstream << "]";
      } else if (h.mode == 2) {
        save_stringstream << h.columns[0] << ";" << h.columns[1];
      }
      
      save_stringstream << ", " << h.nbins << ", " << h.lbound << ", " 
                       << h.ubound << ", " << h.units << ", " 
                       << h.take_point_five << ", mode=" << h.mode;
      
      if (h.mode == 2) {
        save_stringstream << ", " << h.y_axis << ", " << h.nbins_y << ", " 
                         << h.lbound_y << ", " << h.ubound_y << ", " << h.units_y;
      }
      
      save_stringstream << "}, ";
    }
    
    string save_string{save_stringstream.str()};
    save_string[save_string.size()-2] = '}';
    histograms_file << save_string;
    histograms_file.close();
  }

  void plot(const string& outputFileName, FrameAndData& fd, 
            const vector<HistInfo>& his, const string& weighting) {
    if (his.empty()) return;
    
    HistogramCollection hists;
    
    // Book all histograms based on mode
    for (const auto& h : his) {
        double lower = h.lbound - 0.5 * h.take_point_five;
        double upper = h.ubound - 0.5 * h.take_point_five;
        
        if (h.mode == 0) {
            // Normal 1D histogram
            hists.hist1d.push_back(
                fd.Histo1D(h.name.c_str(), h.name.c_str(), h.nbins, lower, upper, h.columns[0].c_str())
            );
        } 
        else if (h.mode == 1) {
            // Stacked histogram - create one histogram per column
            if (h.columns.empty()) {
                cout << "Warning: Stacked histogram '" << h.name 
                     << "' has no columns specified. Skipping.\n";
                continue;
            }
            
            if (h.data_vs_mc) {
                // First column is data (not stacked)
                string data_hist_name = h.name + "_data_component";
                hists.data_hists.push_back(
                    fd.Histo1D(data_hist_name.c_str(), data_hist_name.c_str(),
                              h.nbins, lower, upper, h.columns[0].c_str())
                );
                
                // Remaining columns are MC (stacked)
                vector<ROOT::RDF::RResultPtr<TH1D>> stack_histograms;
                for (size_t i = 0; i < h.columns.size(); ++i) {
                    string hist_name = h.name + "_component_" + std::to_string(i);
                    stack_histograms.push_back(
                        fd.Histo1D(hist_name.c_str(), hist_name.c_str(), 
                                  h.nbins, lower, upper, h.columns[i].c_str())
                    );
                }
                hists.stacked_hists.push_back(stack_histograms);
            } else {
                // All columns are stacked
                vector<ROOT::RDF::RResultPtr<TH1D>> stack_histograms;
                for (size_t i = 0; i < h.columns.size(); ++i) {
                    string hist_name = h.name + "_component_" + std::to_string(i);
                    stack_histograms.push_back(
                        fd.Histo1D(hist_name.c_str(), hist_name.c_str(), 
                                  h.nbins, lower, upper, h.columns[i].c_str())
                    );
                }
                hists.stacked_hists.push_back(stack_histograms);
            }
        }
        else if (h.mode == 2) {
            // 2D histogram
            if (h.columns.size() < 2) {
                cout << "Warning: 2D histogram '" << h.name 
                     << "' requires 2 columns. Skipping.\n";
                continue;
            }
            
            double lower_y = h.lbound_y;
            double upper_y = h.ubound_y;
            
            hists.hist2d.push_back(
                fd.Histo2D(h.name.c_str(), h.name.c_str(),
                          h.nbins, lower, upper,
                          h.nbins_y, lower_y, upper_y,
                          h.columns[0].c_str(), h.columns[1].c_str())
            );
        }
    }
    
    // Trigger computation in parallel
    cout << "Computing " << (hists.hist1d.size() + hists.hist2d.size() + hists.stacked_hists.size()) 
         << " histogram(s)..." << endl;
    
    // Apply scaling to 1D histograms
    double scale_factor = DataBase::weighting.at(weighting);
    for (auto& hptr : hists.hist1d) {
        hptr.GetValue();
        hptr->Scale(scale_factor);
    }
    
    // Apply scaling to stacked histograms
    for (auto& stack : hists.stacked_hists) {
        for (auto& hptr : stack) {
            hptr.GetValue();
            hptr->Scale(scale_factor);
        }
    }
    
    // Apply scaling to 2D histograms
    for (auto& hptr : hists.hist2d) {
        hptr.GetValue();
        hptr->Scale(scale_factor);
    }
    
    // Apply scaling to data histograms
    for (auto& hptr : hists.data_hists) {
        hptr.GetValue();
        hptr->Scale(scale_factor);
    }
    
    // Save to file
    save_histograms_extended(weighting + "_" + outputFileName, hists, his);
  }

  // [Rest of the pre_process_columns and other functions remain unchanged]
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
    // OS tau jet indicies
    fd.Define("Jet_TauTagOSIndicies", LFuncs::get_os_tau_indicies, {"Jet_TauTagPT", "Jet_TauTagCharge"});

    fd.Define("Jet_TauTagOSMass", LFuncs::use_indicies, {"Jet_TauTagMass", "Jet_TauTagOSIndicies"});
    fd.Define("Jet_TauTagOSPT", LFuncs::use_indicies, {"Jet_TauTagPT", "Jet_TauTagOSIndicies"});
    fd.Define("Jet_TauTagOSEta", LFuncs::use_indicies, {"Jet_TauTagEta", "Jet_TauTagOSIndicies"});
    fd.Define("Jet_TauTagOSPhi", LFuncs::use_indicies, {"Jet_TauTagPhi", "Jet_TauTagOSIndicies"});
    fd.Define("Jet_TauTagOSCharge", LFuncs::use_indicies_int, {"Jet_TauTagCharge", "Jet_TauTagOSIndicies"});
    fd.Define("Jet_TauOSNum", LFuncs::get_size<unsigned int>, {"Jet_TauTagOSIndicies"});



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
    fd.Define("TruthMatchIndicies", LFuncs::get_truth_match_indicies, {"TruthMatchTau", "DeltaRJetTau"});
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
    // forces root to compute all the defines here and now.
    auto dummy = fd.node.Count();
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
  
  void cache_columns(FrameAndData& fd, vector<string> cached_columns){
    cout << "Caching processed data..." << endl;
    fd.node.Cache(cached_columns);
  }
};

#endif
