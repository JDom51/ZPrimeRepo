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
  // fd.Filter(at_least_two, {"Jet_TauTagPT"}, "at_least_two");
}

void DataBase::Analyse(DataStructs::FrameAndData& fd)
{
  // cout << fd.node.GetColumnType("MissingET.MET") << "\n";
  // auto size_of = [](ROOT::VecOps::RVec<unsigned int> tautag_indicies){return tautag_indicies.size() == 2;};
  // auto get_pt2 = [](ROOT::VecOps::RVec<Float_t> met, ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> phi_e){
    // ROOT::VecOps::RVec<Float_t> pt2{};
    // pt2.push_back(met[0] * ( sin(phi_e[0]) - cos(phi_e) * tan(phi[0]) ) / ( sin(phi[1]) - cos(phi[1]) * tan(phi_e[0]) ) );
    // return met[0] * ( sin(phi_e[0]) - cos(phi_e[0]) * tan(phi[0]) ) / ( sin(phi[1]) - cos(phi[1]) * tan(phi[0]) );
    // Float_t frac1 = met[0] * cos(phi_e[0]) / ( sin(theta[1]) * cos(phi[1]) );
    // Float_t frac2 = ( tan(phi_e[0]) - tan(phi[0]) ) / ( tan(phi[1]) - tan(phi[0]) );
    // return frac1 * frac2;
  // };
  // auto get_pt1 = [](ROOT::VecOps::RVec<Float_t> met, ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> phi_e, Float_t pt2){
    // ROOT::VecOps::RVec<Float_t> pt1{};
    // pt1.push_back((met[0] * cos( phi_e[0] ) - pt2[0] * cos( phi[1] ) ) / ( cos( phi[0] ) ) );
    // return (met[0] * cos( phi_e[0] ) - pt2 * cos( phi[1] ) ) / ( cos( phi[0] ) );
    // Float_t frac1 = met[0] * cos(phi_e[0]) / ( sin(theta[0]) * cos(phi[0]) );
    // Float_t frac2 = 1 - ( tan(phi_e[0]) - tan(phi[0]) ) / ( tan(phi[1]) - tan(phi[0]) );
    // return frac1 * frac2;
  // };
  // auto get_theta = [](ROOT::VecOps::RVec<Float_t> eta){
  //   eta[0] = 2 * atan( exp( -eta[0] ) );      
  //   eta[1] = 2 * atan( exp( -eta[1] ) );
  //   return eta;
  // };
  // auto get_met = [](ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> phi_e, Float_t pt1, Float_t pt2){
  //   return ( pt1 * cos(phi[0]) + pt2 * cos(phi[1]) ) / cos(phi_e[0]);
  // };
  // auto get_inv_mass = [](ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> eta){
  //   return sqrt(2 * pt[0] * pt[1] * ( cosh(eta[0] - eta[1]) - cos(phi[0] - phi[1]) ) );
  // };
  // auto add = [](ROOT::VecOps::RVec<Float_t> pt, Float_t pt1, Float_t pt2){
  //   pt[0] = pt[0] + pt1;
  //   pt[1] = pt[1] + pt2;
  //   return pt;
  // };
  // auto filter_open_met = [](ROOT::VecOps::RVec<Float_t> met_eta, ROOT::VecOps::RVec<Float_t> tau_eta)
  // {
  //   // Useless
  //   return abs(met_eta[0]) > abs(tau_eta[0]) && abs(met_eta[0]) > abs(tau_eta[1]);
  // };


  // auto filter_phi = [](ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> phi_e)
  // {
  //   if(phi[0] < -M_PI/2 && phi[1] > 0)
  //   {
  //     if(phi_e[0] < 0){phi_e[0] += 2*M_PI;}
      
  //     phi[0] += 2* M_PI;
  //   }
  //   if(phi[0] > 0 && phi[1] < -M_PI/2) //3.145192/2
  //   {
  //     if(phi_e[0] < 0){phi_e[0] += 2*M_PI;}
  //     phi[1] += 2* M_PI;
  //   }
  //   if(phi[0] < phi[1]){return phi[0] < phi_e[0] && phi_e[0] < phi[1];}
  //   else if(phi[0] > phi[1]){return phi[1] < phi_e[0] && phi_e[0] < phi[0];}
  //   else{return phi[0] == phi_e[0];}
  // };
  
  // auto get_x2 = [](ROOT::VecOps::RVec<Float_t> met, ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> phi_e, ROOT::VecOps::RVec<Float_t> phi){
  //   // Float_t frac1 = (met[0] * cos(phi_e[0]) )/ (pt[1] * cos(phi[1]));
  //   // Float_t frac2 = ( tan(phi_e[0]) - tan(phi[0]) ) / (tan(phi[1]) - tan(phi[0]));
  //   Float_t frac1 = ( met[0] * sin(phi_e[0]) - met[0] * cos(phi_e[0]) * tan(phi[0]) ) / pt[1];
  //   Float_t frac2 = 1 / ( sin(phi[1]) - cos(phi[1]) * tan(phi[0]) ); 
  //     return 1 / (frac1 * frac2 + 1);};
  // auto get_x1 = [](ROOT::VecOps::RVec<Float_t> met, ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> phi_e, ROOT::VecOps::RVec<Float_t> phi){
  //   Float_t frac1 = ( met[0] * cos(phi_e[0]) ) / ( pt[0] * cos(phi[0]) ); 
  //   Float_t frac2 = 1 - ( tan(phi_e[0]) - tan(phi[0]) ) / ( tan(phi[1]) - tan(phi[0]) );
  //   return 1 / (frac1 * frac2 + 1);};
  // auto x_boundaries = [](Float_t x){
  //   return 0 < x && x <= 1;};

  // auto my_way = [](ROOT::VecOps::RVec<Float_t> phi)
  // {
  //   return (phi[0] > M_PI/2 && phi[1] < -M_PI/2) || (phi[0] < -M_PI/2 && phi[1] > M_PI/2);
  // };
  // auto get_angle_between = [](ROOT::VecOps::RVec<Float_t> px, ROOT::VecOps::RVec<Float_t> py, ROOT::VecOps::RVec<Float_t> pz)
  // {
  //   Float_t numerator = px[0] * px[1] + py[0] * py[1];
  //   Float_t denominator = sqrt( pow( px[0], 2) + pow( py[0], 2)) * sqrt( pow( px[1], 2) + pow( py[1], 2));
  //   return acos( numerator / denominator );
  // };
  // auto greater_zero = [](Float_t npt1, Float_t npt2)
  // {
  //   return npt1 >= 0 && npt2 >= 0;
  // };
  // auto pid_check = [](ROOT::VecOps::RVec<int> pids)
  // {
  //   size_t count{0};
  //   for(auto pid : pids)
  //   {
  //     if(abs(pid) == 15){count++;}
  //   }
  //   return count == 2;
  // };


  
  // // auto no_b2b = [](ROOT::VecOps::RVec<Float_t> phi) {
  // // 
  // // //   return abs(cos(phi[0] - phi[1])) < 0.95;
  // // };

  // auto wraparound_fix = [](ROOT::VecOps::RVec<Float_t> phi_jet, ROOT::VecOps::RVec<Float_t> eta_jet, ROOT::VecOps::RVec<Float_t> pt_jet, ROOT::VecOps::RVec<Float_t> phi_met, ROOT::VecOps::RVec<Float_t> eta_met, ROOT::VecOps::RVec<Float_t> met)
  // {
  //   Float_t px1{ pt_jet[0] * cos(phi_jet[0]) };
  //   Float_t py1{ pt_jet[0] * sin(phi_jet[0]) };
  //   Float_t pz1{ pt_jet[0] * sinh(eta_jet[0]) };
  //   Float_t p1{ static_cast<Float_t>(sqrt (  pow(px1, 2) + pow(py1, 2) + pow(pz1, 2) ) ) };
  //   Float_t px2{ pt_jet[1] * cos( phi_jet[1] ) };
  //   Float_t py2{ pt_jet[1] * sin( phi_jet[1] )};
  //   Float_t pz2{ pt_jet[1] * sinh( eta_jet[1] ) };
  //   Float_t p2{ static_cast<Float_t>(sqrt( pow(px2, 2) + pow(py2, 2) + pow(pz2, 2) )) };
  //   Float_t ex{ met[0] * cos(phi_met[0] ) };
  //   Float_t ey{ met[0] * sin(phi_met[0] ) };
  //   Float_t ez{ met[0] * sinh(eta_met[0] )};
  //   Float_t em{ static_cast<Float_t>(sqrt( pow(ex,2) + pow(ey, 2) + pow(ez, 2) ) )};
      //  this is wrong
  //   Float_t jet_angle = acos( (px1 * px2 + py1 * py2 + pz1 * pz2) / (p1 * p2) );
  //   Float_t jet_met_1_angle = acos( (px1 * ex + py1 * ey + pz1 *ez) / (p1 * em) );
  //   Float_t jet_met_2_angle = acos( (px2 * ex + py2 * ey + pz2 * ez) / (p2 * em) );

  //   return abs(jet_met_1_angle + jet_met_2_angle - jet_angle ) / jet_angle  < 0.1;
  // };

  // fd.Filter(size_of, {"Jet_TauTagIndicies"});
  
   // fd.Filter(filter_open_met, {"MissingET.Eta", "Jet_TauTagEta"});
  // // fd.Define("x2", get_x2, {"MissingET.MET", "Jet_TauTagPT", "MissingET.Phi", "Jet_TauTagPhi"});
  // // fd.Define("x1", get_x1, {"MissingET.MET", "Jet_TauTagPT", "MissingET.Phi", "Jet_TauTagPhi"});
  // // fd.Filter(x_boundaries, {"x2"});
  // // fd.Filter(x_boundaries, {"x1"});
  // // fd.Filter(filter_phi, {"Jet_TauTagPhi", "MissingET.Phi"});
  // // fd.Filter(my_way, {"Jet_TauTagPhi"});
  // // fd.Filter(no_b2b, {"Jet_TauTagPhi"});
  // // fd.Define("NeutrinoTheta", get_theta, {"Jet_TauTagEta"});
  fd.Define("AngleBetweenMET", LFuncs::met_jet_ang, {"MissingET.Phi", "Jet_DTauTagPhi"});
  fd.Define("Neutrino_PT2", LFuncs::get_col_neutrinopt2, {"MissingET.MET", "Jet_TauTagPhi", "MissingET.Phi"});
  fd.Define("Neutrino_PT1", LFuncs::get_col_neutrinopt1, {"MissingET.MET", "Jet_TauTagPhi", "MissingET.Phi", "Neutrino_PT2"});
  // fd.Filter(greater_zero, {"Neutrino_PT1", "Neutrino_PT2"});
  //fd.Define("MET_TEST", get_met, {"Jet_TauTagPhi", "MissingET.Phi", "Neutrino_PT1", "Neutrino_PT2"});
  fd.Define("Jet_TauTagNeutrinoPT", LFuncs::add_col_pt, {"Jet_TauTagPT", "Neutrino_PT1", "Neutrino_PT2"});
  fd.Define("TauJetInvMass", LFuncs::inv_mass_ml, {"Jet_TauTagPT", "Jet_TauTagPhi", "Jet_TauTagEta"});
  fd.Define("TauJetInvMassWithNeutrino", LFuncs::inv_mass_ml, {"Jet_TauTagNeutrinoPT", "Jet_TauTagPhi", "Jet_TauTagEta"});
  fd.Define("TauJetDeltaR", LFuncs::get_DeltaR, {"Jet_TauTagPhi", "Jet_TauTagEta"});

  // fd.Filter(SC::met_angle_diff, {"AngleBetweenMET"}, "AngleMET2times10to-6");


  // // fd.Filter(pid_check, {"Particle.PID"});
}
void DataBase::TruthAnalysis(DataStructs::FrameAndData& fd){
  // truth analysis was all about matching the generator level particles to the jets and calculating the relevant quantities such
  // as the Delta R the Delta Phi and inv mass with and without neutrinos

}
void DataBase::GeneratorLevelTauAnalysis(DataStructs::FrameAndData& fd){
  // generator level tau analysis was purely on focused on the generator information of the taus and no reco level variables
  fd.Define("GenTauInvMass", LFuncs::inv_mass_ml, {"Tau_PT", "Tau_Phi", "Tau_Eta"});
  fd.node.Foreach([](ROOT::VecOps::RVec<Float_t> tau_pt, ROOT::VecOps::RVec<Float_t> tau_phi, ROOT::VecOps::RVec<Float_t> tau_eta){cout<<"TauPT: "<<tau_pt<< " TauEta: " << tau_eta << " TauPhi: " << tau_phi << " \n";}, {"Tau_PT", "Tau_Phi", "Tau_Eta"});

}
void DataBase::SelectionCutAnalysis(DataStructs::FrameAndData& fd){
  // selection cut was all about restricting the inv mass region and then figuring out what the optimal selection cut was
  // Need to restrict IMass range to 80-100
  // then generate columns with the number of events for pt of specific cut.
  // for loop for range and for loop for all events.  thats selcutalgo i think
  // then need to save the example thing and save the new csv file
  // get_integral(fd, column_name, selection_cut, sel_cut_values, resolution)
  //get invmass

  auto wraparound_fix = [](ROOT::VecOps::RVec<Float_t> phi_jet, ROOT::VecOps::RVec<Float_t> pt_jet, ROOT::VecOps::RVec<Float_t> phi_met, ROOT::VecOps::RVec<Float_t> met)
  {
    Float_t px1{ pt_jet[0] * cos(phi_jet[0]) };
    Float_t py1{ pt_jet[0] * sin(phi_jet[0]) };
    Float_t p1{ static_cast<Float_t>(sqrt (  pow(px1, 2) + pow(py1, 2) ) ) };
    Float_t px2{ pt_jet[1] * cos( phi_jet[1] ) };
    Float_t py2{ pt_jet[1] * sin( phi_jet[1] )};
    Float_t p2{ static_cast<Float_t>(sqrt( pow(px2, 2) + pow(py2, 2))) };
    Float_t ex{ met[0] * cos(phi_met[0] ) };
    Float_t ey{ met[0] * sin(phi_met[0] ) };
    Float_t em{ static_cast<Float_t>(sqrt( pow(ex,2) + pow(ey, 2) ) )};
      //  this is wrong
    Float_t jet_angle = acos( (px1 * px2 + py1 * py2) / (p1 * p2) );
    Float_t jet_met_1_angle = acos( (px1 * ex + py1 * ey) / (p1 * em) );
    Float_t jet_met_2_angle = acos( (px2 * ex + py2 * ey) / (p2 * em) );
// 
    return abs(jet_met_1_angle + jet_met_2_angle - jet_angle );
  };
  // fd.Define("NeutrinoPT2", LFuncs::get_co/ // fd.Filter(filter_open_met, {"MissingET.Eta", "Jet_TauTagEta"});
  // // fd.Define("x2", get_x2, {"MissingET.MET", "Jet_TauTagPT", "MissingET.Phi", "Jet_TauTagPhi"});
  // // fd.Define("x1", get_x1, {"MissingET.MET", "Jet_TauTagPT", "MissingET.Phi", "Jet_TauTagPhi"});
  // // fd.Filter(x_boundaries, {"x2"});
  // // fd.Filter(x_boundaries, {"x1"});
  // // fd.Filter(filter_phi, {"Jet_TauTagPhi", "MissingET.Phi"});
  // // fd.Filter(my_way, {"Jet_TauTagPhi"});
  // // fd.Filter(no_b2b, {"Jet_TauTagPhi"});
  // // fd.Define("NeutrinoTheta", get_theta, {"Jet_TauTagEta"});
  // fd.Define("AngleBetween", get_angle_between, {"DJet_TauTagPx", "DJet_TauTagPy", "DcJet_TauTagPz"});
  // fd.Define("Neutrino_PT2", get_pt2, {"MissingET.MET", "Jet_TauTagPhi", "MissingET.Phi"});
  // fd.Define("Neutrino_PT1", get_pt1, {"MissingET.MET", "Jet_TauTagPhi", "MissingET.Phi", "Neutrino_PT2"});
  // fd.Filter(wraparound_fix, {"Jet_TauTagPhi", "Jet_TauTagEta", "Jet_TauTagPT", "MissingET.Phi", "MissingET.Eta", "MissingET.MET"});
  // // fd.Filter(greater_zero, {"Neutrino_PT1", "Neutrino_PT2"});
  // //fd.Define("MET_TEST", get_met, {"Jet_TauTagPhi", "MissingET.Phi", "Neutrino_PT1", "Neutrino_PT2"});
  // fd.Define("Jet_TauTagNeutrinoPT", LFuncs::add_col_pt, {"Jet_TauTagPT", "Neutrino_PT1", "Neutrino_PT2"});
  // fd.Define("TauJetInvMass", get_inv_mass, {"Jet_TauTagPT", "Jet_TauTagPhi", "Jet_TauTagEta"});
  // fd.Define("TauJetInvMassWithNeutrino", get_inv_mass, {"Jet_TauTagNeutrinoPT", "Jet_TauTagPhi", "Jet_TauTagEta"});l_neutrinopt2, {"MissingET.MET", "Jet_DTauTagPhi", "MissingET.Phi"});
  

  fd.Define("NeutrinoPT2", LFuncs::get_col_neutrinopt2, {"MissingET.MET", "Jet_DTauTagPhi", "MissingET.Phi"});
  fd.Define("NeutrinoPT1", LFuncs::get_col_neutrinopt1, {"MissingET.MET", "Jet_DTauTagPhi", "MissingET.Phi", "NeutrinoPT2"});
  fd.Define("Jet_DTauTagNeutrinoPT", LFuncs::add_col_pt, {"Jet_DTauTagPT", "NeutrinoPT1", "NeutrinoPT2"});
  fd.Define("TruthMatchDeltaPhi", LFuncs::get_delta_phi, {"Jet_DTauTagPT", "Jet_DTauTagPhi"});
  fd.Define("AngleBetweenMET", LFuncs::met_jet_ang, {"MissingET.Phi", "Jet_DTauTagPhi"});
  // fd.Define("AngleBetweenMET", wraparound_fix, {"Jet_DTauTagPhi", "Jet_DTauTagPT", "MissingET.Phi", "MissingET.MET"});
  fd.Define("TauJetInvMass", LFuncs::inv_mass_ml, {"Jet_DTauTagPT", "Jet_DTauTagEta", "Jet_DTauTagPhi"});
  fd.Define("TauJetInvMassWithNeutrino", LFuncs::inv_mass_ml, {"Jet_DTauTagNeutrinoPT", "Jet_DTauTagPhi", "Jet_DTauTagEta"});
  fd.Define("GenTauInvMass", LFuncs::inv_mass_ml, {"Tau_PT", "Tau_Phi", "Tau_Eta"});

  // fd.Define("TauTagDeltaR", LFuncs::get_DeltaR, {"Jet_TauTagPhi", "Jet_TauTagEta"});
  // TRUTH JET DEFINITIONS.
  // fd.Define("TruthMatchDeltaPhiTRUTHJET", LFuncs::get_delta_phi, {"Jet_TruthTauMatchPT", "Jet_TruthTauMatchPhi"});
  // fd.Define("AngleBetweenMETTRUTHJET", LFuncs::met_jet_ang, {"MissingET.Phi", "Jet_TruthTauMatchPhi"});
  // fd.Define("Jet_TruthTauMatchNeutrinoPT", LFuncs::add_col_pt, {"Jet_TruthTauMatchPT", "NeutrinoPT1", "NeutrinoPT2"});
  // fd.Define("TauJetInvMassTRUTHJET", LFuncs::inv_mass_ml, {"Jet_TruthTauMatchPT", "Jet_TruthTauMatchEta", "Jet_TruthTauMatchPhi"});
  // fd.Define("TauJetInvMassWithNeutrinoTRUTHJET", LFuncs::inv_mass_ml, {"Jet_TruthTauMatchNeutrinoPT", "Jet_TruthTauMatchPhi", "Jet_TruthTauMatchEta"});


  
  fd.Filter(SC::met_angle_diff, {"AngleBetweenMET"}, "METbetweenJets");
  fd.Filter(SC::deltaRcut0p3, {"DeltaRJetTauSel"}, "DeltaRCut0p3");
  // fd.Filter(SC::pt_tau_cut, {"Jet_DTauTagPT"}, "jetptg20");
  // fd.Filter(SC::gen_inv_mass_g20, {"GenTauInvMass"}, "GenTauInvMassG20");
  // fd.Filter(SC::deltaRcut0p2, {"DeltaRJetTauSel"}, "DeltaRCut0p2");
  // fd.Filter(SC::delta_phi_1, {"TruthMatchDeltaPhi"}, "DeltaPhi1");
  // fd.Filter(SC::delta_phi_1p2, {"TruthMatchDeltaPhi"}, "DeltaPhi1p2");
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
  string outputFileName = DataBase::fnames[analysis_mode];
  // HistInfo: {name of hist, x axis title, Column name, nbins, lbound, ubound, units, offset by .5}

  // Analysis code
  cout<<"Creating RDF\n";
  FrameAndData fd{FrameAndData(mode, files)};
  cout<<"PreProcessColumns'\n";
  BackEnd::pre_process_columns(fd);
  BackEnd::initial_cut(fd, DataBase::broad_cuts[analysis_mode]);
  // Analyse
  cout<<"Analysing\n";
  DataBase::analysis_funcs[analysis_mode](fd);
  cout<<"Plotting\n";
  BackEnd::plot(outputFileName + fd.save_string, fd, DataBase::histograms[analysis_mode]);
  BackEnd::save_hist_info(outputFileName, DataBase::histograms[analysis_mode]);
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
