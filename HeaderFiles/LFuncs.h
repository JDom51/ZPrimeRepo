#ifndef LFUNCS
#define LFUNCS
#include<vector>
#include<string>
#include<cmath>
using std::string;
using std::vector;

Float_t get_index(ROOT::VecOps::RVec<Float_t> possible_pairs, vector<int> indicies)
{
  // unsure whether setting default value of -1 will create problems later on but hopefully not.
  Float_t minimum{1000};
  Float_t index{-1};
  for(int i{0}; i < possible_pairs.size(); i++)
  {
    if(possible_pairs[i] < minimum && std::find(indicies.begin(), indicies.end(), i) == indicies.end())
    {
      minimum = possible_pairs[i];
      index = i;
    }
  }
  return index;
}

namespace LFuncs
{
    ROOT::VecOps::RVec<Float_t> add_col_pt(ROOT::VecOps::RVec<Float_t> pt, Float_t pt1, Float_t pt2){
    pt[0] = pt[0] + pt1;
    pt[1] = pt[1] + pt2;
    return pt;
  };
  int create_one(ROOT::VecOps::RVec<Float_t> mas)
  {
    return 1;
  }
  int create_zero(ROOT::VecOps::RVec<Float_t> mas)
  {
    return 0;
  }
  int create_minus_one(ROOT::VecOps::RVec<Float_t> mas)
  {
    return -1;
  }

  ROOT::VecOps::RVec<Float_t> get_alpha(ROOT::VecOps::RVec<Float_t> etas){
    ROOT::VecOps::RVec<Float_t> alpha{vector<Float_t>(etas.size())};
    for(int i{0}; i < etas.size(); i++){
      alpha[i] = atan(1/sinh(etas[i]));
    }
    return alpha;
  }
  ROOT::VecOps::RVec<unsigned int> get_indicies_int(ROOT::VecOps::RVec<int> int_vec, int number)
  {
    ROOT::VecOps::RVec<unsigned int> indicies{};
    for(int i{0}; i < int_vec.size(); i++)\
    {
      if(int_vec[i]  == number){indicies.push_back(i);}
    }
    return indicies;
  }
  ROOT::VecOps::RVec<unsigned int> get_indicies_greater_than_zero_ui(ROOT::VecOps::RVec<unsigned int> boolean_vec)
  {
    ROOT::VecOps::RVec<unsigned int> indicies{};
    for(int i{0}; i < boolean_vec.size(); i++)
    {
      if(boolean_vec[i] > 0){indicies.push_back(i);}
    }
    return indicies;
  };
  ROOT::VecOps::RVec<unsigned int> get_indicies_is_one(ROOT::VecOps::RVec<unsigned int> boolean_vec)
  {
    ROOT::VecOps::RVec<unsigned int> indicies{};
    for(int i{0}; i < boolean_vec.size(); i++)
    {
      if(boolean_vec[i] == 1){indicies.push_back(i);}
    }
    return indicies;
  };
  ROOT::VecOps::RVec<Float_t> use_indicies(ROOT::VecOps::RVec<Float_t> var, ROOT::VecOps::RVec<unsigned int> indicies)
  {
    ROOT::VecOps::RVec<Float_t> new_var{};
    for(auto i : indicies)
    {
      new_var.push_back(var[i]);
    }
    return new_var;
  };
  ROOT::VecOps::RVec<Float_t> not_use_indicies(ROOT::VecOps::RVec<Float_t> var, ROOT::VecOps::RVec<unsigned int> indicies){
    ROOT::VecOps::RVec<Float_t> new_var{};
    for(int i{0}; i < var.size(); i++){
      if(std::find(indicies.begin(), indicies.end(), i) == indicies.end()){
        new_var.push_back(var[i]);
      }
    }
    return new_var;
  }
  ROOT::VecOps::RVec<int> use_indicies_int(ROOT::VecOps::RVec<int> var, ROOT::VecOps::RVec<unsigned int> indicies)
  {
    ROOT::VecOps::RVec<Float_t> new_var{};
    for(auto i : indicies)
    {
      new_var.push_back(var[i]);
    }
    return new_var;
  };
  ROOT::VecOps::RVec<Float_t> use_neither_indicies(ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<unsigned int> btag, ROOT::VecOps::RVec<unsigned int> tautag)
  {
    ROOT::VecOps::RVec<Float_t> not_bt_pt{};
    for(int i{0}; i < pt.size(); i++)
    {
      if(std::find(btag.begin(), btag.end(), i) == btag.end() && std::find(tautag.begin(), tautag.end(), i) == tautag.end())
      {
        not_bt_pt.push_back(pt[i]);
      };
    }
    return not_bt_pt;
  };
  Float_t get_angle_between_jets(Float_t phi_1, Float_t phi_2, Float_t eta_1, Float_t eta_2)
  {
    return static_cast<Float_t>(acos( ( cos(phi_1) * cos(phi_2) + sin(phi_1) * sin(phi_2) + sinh(eta_1) * sinh(eta_2) ) / ( sqrt(1 + pow(sinh(eta_1), 2) )  * sqrt(1 + pow(sinh(eta_2), 2)) ) ));
  };
  template<typename Type> int get_size(ROOT::VecOps::RVec<Type> vect)
  {
    return vect.size();
  }
  ROOT::VecOps::RVec<unsigned int> get_g10_pt_indicies(ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> eta){
    ROOT::VecOps::RVec<unsigned int> indicies{};
    for(int i{0}; i < pt.size(); i++)
    {
      if((pt[i] <= 10.0) * (0.00) +
                                           (abs(eta[i]) <= 1.5) * (pt[i] > 10.0)  * (0.95) +
                         (abs(eta[i]) > 1.5 && abs(eta[i]) <= 2.5) * (pt[i] > 10.0)  * (0.85) +
                         (abs(eta[i]) > 2.5)                                   * (0.00)){
        indicies.push_back(i);
      }
    }
    return  indicies;
  }
  int get_neither_size(ROOT::VecOps::RVec<unsigned int> reference, ROOT::VecOps::RVec<unsigned int> indicies1, ROOT::VecOps::RVec<unsigned int> indicies2)
  {
    int count{0};
    for(int i{0}; i < reference.size(); i++)
    {
      if(std::find(indicies1.begin(), indicies1.end(), i) == indicies1.end() && std::find(indicies2.begin(), indicies2.end(), i) == indicies2.end()){count++;} // need to change this 
    }
    return count;
  }
  Float_t inv_mass_ml(ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> eta, ROOT::VecOps::RVec<Float_t> phi)
  {
    // if(130 <= sqrt(2 * pt[0] * pt[1] * ( cosh(eta[0] - eta[1]) - cos(phi[0] - phi[1]) ) ) && sqrt(2 * pt[0] * pt[1] * ( cosh(eta[0] - eta[1]) - cos(phi[0] - phi[1]) ) ) <= 170){
      // cout << "IMass: " << sqrt(2 * pt[0] * pt[1] * ( cosh(eta[0] - eta[1]) - cos(phi[0] - phi[1]) ) ) << " PT: " << pt << " Eta: " << eta << " Phi: " << phi <<"\n";
    // }
    return sqrt(2 * pt[0] * pt[1] * ( cosh(eta[0] - eta[1]) - cos(phi[0] - phi[1]) ) );
  }
  Float_t inv_mass_pt(ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> eta, ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> m)
  {
    Float_t px{0};
    Float_t py{0};
    Float_t pz{0};
    Float_t tot_e{0};
    Float_t tot_px{0};
    Float_t tot_py{0};
    Float_t tot_pz{0};
    for(int i{0}; i < m.size(); i++)
    {
      px = pt[i] * cos(phi[i]);
      py = pt[i] * sin(phi[i]);
      pz = pt[i] * sinh(eta[i]);
      tot_e += sqrt(pow(m[i], 2) + pow(px, 2) + pow(py, 2) + pow(pz, 2));
      tot_px += px;
      tot_py += py;
      tot_pz += pz;
    }
  
    return sqrt(pow(tot_e,2) - pow(tot_px, 2) - pow(tot_py, 2) - pow(tot_pz, 2));
  }
  Float_t inv_mass_xyz(ROOT::VecOps::RVec<Float_t> E, ROOT::VecOps::RVec<Float_t> px, ROOT::VecOps::RVec<Float_t> py, ROOT::VecOps::RVec<Float_t> pz)
  {
    Float_t tot_E{0};
    Float_t tot_px{0};
    Float_t tot_py{0};
    Float_t tot_pz{0};
    for(int i{0}; i < E.size(); i++)
    {
      tot_E += E[i];
      tot_px += px[i];
      tot_py += py[i];
      tot_pz += pz[i];
    }
    return sqrt(pow(tot_E,2) - pow(tot_px, 2) - pow(tot_py, 2) - pow(tot_pz, 2));
  }
  ROOT::VecOps::RVec<Float_t> get_px(ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> phi)
  {
    ROOT::VecOps::RVec<Float_t> px{ROOT::VecOps::RVec<Float_t>(pt.size(), 0.0f)};
    for(int i{0}; i < pt.size(); i++)
    {
      px[i] = pt[i] * cos(phi[i]);
    }
    return px;
  }
  ROOT::VecOps::RVec<Float_t> get_py(ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> phi)
  {
    ROOT::VecOps::RVec<Float_t> py{ROOT::VecOps::RVec<Float_t>(pt.size(), 0.0f)};
    for(int i{0}; i < pt.size(); i++)
    {
      py[i] = pt[i] * sin(phi[i]);
    }
    return py;
  }
  ROOT::VecOps::RVec<Float_t> get_pz(ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<Float_t> eta)
  {
    ROOT::VecOps::RVec<Float_t> pz{ROOT::VecOps::RVec<Float_t>(pt.size(), 0.0f)};
    for(int i{0}; i < pt.size(); i++)
    {
      pz[i] = pt[i] * sinh(eta[i]);
    }
    return pz;
  }


  ROOT::VecOps::RVec<Float_t> get_difference_between_vectors(ROOT::VecOps::RVec<Float_t> v1, ROOT::VecOps::RVec<Float_t>& v2){
    for(int i{0}; i< v1.size(); i++){
      v1[i] -= v2[i];
    }
    return v1;
  }
  ROOT::VecOps::RVec<Float_t> get_normalised_by_v1_difference_between_vectors(ROOT::VecOps::RVec<Float_t> v1, ROOT::VecOps::RVec<Float_t>& v2){
    for(int i{0}; i< v1.size(); i++){
      v1[i] = (v1[i]-v2[i])/v1[i];
    }
    return v1;
  }

  Float_t get_col_neutrinopt1(ROOT::VecOps::RVec<Float_t> met, ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> phi_e, Float_t pt2){
    return (met[0] * cos( phi_e[0] ) - pt2 * cos( phi[1] ) ) / ( cos( phi[0] ) );
  }

  Float_t get_col_neutrinopt2(ROOT::VecOps::RVec<Float_t> met, ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> phi_e){
    // return met[0] * ( sin(phi_e[0]) - cos(phi_e[0]) * tan(phi[0]) ) / ( sin(phi[1]) - cos(phi[1]) * tan(phi[0]) );
    // cout << "MET: " << met << " Phi: " << phi << " phi_e: " << phi_e <<"\n";
    return met[0] * ( - cos(phi_e[0]) * tan(phi[0]) + sin(phi_e[0]) ) / (- cos(phi[1]) * tan(phi[0]) + sin(phi[1]));
  }

  ROOT::VecOps::RVec<unsigned int> get_gen_indicies(ROOT::VecOps::RVec<unsigned int> tau_tag, ROOT::VecOps::RVec<unsigned int> gen_tag)
  {
    ROOT::VecOps::RVec<unsigned int> gen_indicies{};
    for(int i{0}; i < tau_tag.size(); i++)
    {
      if(tau_tag[i] == gen_tag[i] && tau_tag[i] == 1){gen_indicies.push_back(i);}
    }
    return gen_indicies;
  }
  Float_t get_delta_phi(ROOT::VecOps::RVec<Float_t> phi)
  {
    // Float_t px1 = pt[0] * cos(phi[0]);
    // Float_t py1 = pt[0] * sin(phi[0]);
    // Float_t p1 = sqrt(px1 * px1 + py1 * py1);
    // Float_t px2 = pt[1] * cos(phi[1]);
    // Float_t py2 = pt[1] * sin(phi[1]);
    // Float_t p2 = sqrt(px2 * px2 + py2 * py2);
    // return abs( acos( (px1 * px2 + py1 * py2) / (p1 * p2) ) ) ;
    return abs( acos( ( cos(phi[0]) * cos( phi[1] )  + sin(phi[0]) * sin(phi[1]) ) ) );

  };
  Float_t get_delta_phi_special(Float_t phi1, Float_t phi2)
  {
    return abs( acos( cos(phi1) * cos( phi2 )  + sin(phi1) * sin(phi2) ) );
  }

  ROOT::VecOps::RVec<unsigned int> get_tau_indicies(ROOT::VecOps::RVec<int> pid, ROOT::VecOps::RVec<int> status, ROOT::VecOps::RVec<Float_t> particle_pt)
  {
    ROOT::VecOps::RVec<unsigned int> tau_indicies{};
    for(int i{0}; i < pid.size(); i++)
    {
      if(abs(pid[i]) == 15 && status[i] == 2){tau_indicies.push_back(i);}
    }
    return tau_indicies;
  }

  ROOT::VecOps::RVec<unsigned int> get_electron_indicies(ROOT::VecOps::RVec<int> pid, ROOT::VecOps::RVec<int> status){
    ROOT::VecOps::RVec<unsigned int> electron_indicies{};
    for(int i{0}; i < pid.size(); i++){
      if(abs(pid[i]) == 11){electron_indicies.push_back(i);}
    }
    return electron_indicies;
  }
  ROOT::VecOps::RVec<unsigned int> get_muon_indicies(ROOT::VecOps::RVec<int> pid, ROOT::VecOps::RVec<int> status){
    ROOT::VecOps::RVec<unsigned int> muon_indicies{};
    for(int i{0};i < pid.size(); i++){
      if(abs(pid[i]) == 13){muon_indicies.push_back(i);}
    }
    return muon_indicies;
  }
  
  ROOT::VecOps::RVec<unsigned int> get_tau_neutrino_indicies(ROOT::VecOps::RVec<int> pid, ROOT::VecOps::RVec<int> status, ROOT::VecOps::RVec<unsigned int> intermediate_indicies){
    ROOT::VecOps::RVec<unsigned int> tau_neutrino_indicies{};
    ROOT::VecOps::RVec<unsigned int> final_tau_neutrino_indicies{};
    for(int i{0}; i < pid.size(); i++){
      if(abs(pid[i]) == 16 && status[i] == 1){tau_neutrino_indicies.push_back(i);}
    }
    for(auto i : intermediate_indicies){
      final_tau_neutrino_indicies.push_back(tau_neutrino_indicies[i]);
    }
    return final_tau_neutrino_indicies;
  }
  ROOT::VecOps::RVec<unsigned int> get_intermediate_tau_indicies(ROOT::VecOps::RVec<unsigned int> tau_indicies, ROOT::VecOps::RVec<Float_t> particle_pt){
    ROOT::VecOps::RVec<unsigned int> new_indicies{};
    for(int i{0}; i < tau_indicies.size(); i++){
      if(particle_pt[tau_indicies[i]] > 10){
        new_indicies.push_back(i);
      }
    }
    return new_indicies;
  }
  ROOT::VecOps::RVec<unsigned int> get_cut_tau_indicies(ROOT::VecOps::RVec<unsigned int> tau_indicies, ROOT::VecOps::RVec<unsigned int> selected_indicies){
    ROOT::VecOps::RVec<Float_t> new_indicies{};
    for(int i{0}; i < selected_indicies.size(); i++){
      new_indicies.push_back(tau_indicies[selected_indicies[i]]);
    }
    return new_indicies;
  }

  ROOT::VecOps::RVec<Float_t> get_truth_met_met(ROOT::VecOps::RVec<Float_t> tau_neutrino_pt, ROOT::VecOps::RVec<Float_t> tau_neutrino_phi){
    Float_t truth_met_x{0};
    Float_t truth_met_y{0};
    for(int i{0}; i < tau_neutrino_pt.size(); i++){
      truth_met_x += tau_neutrino_pt[i] * cos(tau_neutrino_phi[i]);
      truth_met_y += tau_neutrino_pt[i] * sin(tau_neutrino_phi[i]);
    }
    return {sqrt(truth_met_x * truth_met_x + truth_met_y * truth_met_y)};
  }
  ROOT::VecOps::RVec<Float_t> get_truth_met_phi(ROOT::VecOps::RVec<Float_t> tau_neutrino_pt, ROOT::VecOps::RVec<Float_t> tau_neutrino_phi){
    Float_t truth_met_x{0};
    Float_t truth_met_y{0};
    for(int i{0}; i < 2; i++){
      truth_met_x += tau_neutrino_pt[i] * cos(tau_neutrino_phi[i]);
      truth_met_y += tau_neutrino_pt[i] * sin(tau_neutrino_phi[i]);
    }
    return {atan(truth_met_y/truth_met_x)};
  }
  ROOT::VecOps::RVec<Float_t> get_delta_r_1(ROOT::VecOps::RVec<Float_t> eta_tau, ROOT::VecOps::RVec<Float_t> eta_jet, ROOT::VecOps::RVec<Float_t> phi_tau, ROOT::VecOps::RVec<Float_t> phi_jet)
  {
    ROOT::VecOps::RVec<Float_t> delta_r{};
    int size{};
    for(int j{0}; j < eta_jet.size(); j++)
    {
      delta_r.push_back( sqrt( pow( eta_tau[0] - eta_jet[j] ,2) + pow( get_delta_phi_special(phi_tau[0], phi_jet[j]) , 2) ) );
    }
    return delta_r;
  }
  Float_t met_jet_ang(ROOT::VecOps::RVec<Float_t> met_phi, ROOT::VecOps::RVec<Float_t> tau_phi){
    Float_t taumet1 = get_delta_phi_special(met_phi[0], tau_phi[0]);
    Float_t taumet2 = get_delta_phi_special(met_phi[0], tau_phi[1]);
    Float_t tautau12 = get_delta_phi_special(tau_phi[0], tau_phi[1]);
    return abs(taumet1 + taumet2 - tautau12); 
  }
  ROOT::VecOps::RVec<Float_t> get_delta_r_2(ROOT::VecOps::RVec<Float_t> eta_tau, ROOT::VecOps::RVec<Float_t> eta_jet, ROOT::VecOps::RVec<Float_t> phi_tau, ROOT::VecOps::RVec<Float_t> phi_jet)
  {
    ROOT::VecOps::RVec<Float_t> delta_r{};
    int size{}; 
    for(int j{0}; j < eta_jet.size(); j++)
    {
      delta_r.push_back( sqrt( pow( eta_tau[1] - eta_jet[j] ,2) + pow( get_delta_phi_special(phi_tau[1], phi_jet[j]) , 2) ) );
    }
    return delta_r;
  }
  Float_t get_DeltaR(ROOT::VecOps::RVec<Float_t> phi, ROOT::VecOps::RVec<Float_t> eta){
    return sqrt(pow(phi[0] - phi[1], 2) + pow(eta[0] - eta[1], 2));
  }
  ROOT::VecOps::RVec<unsigned int> get_delta_r_indicies(ROOT::VecOps::RVec<Float_t>& delta_r_1, ROOT::VecOps::RVec<float_t>& delta_r_2)
  {
    Float_t minimum_value{10000};
    ROOT::VecOps::RVec<unsigned int> indicies{0, 0};
    if(delta_r_1.size() == 1 || delta_r_2.size() == 1){return {};}

    for(int i{0}; i < delta_r_1.size(); i++)
    {
      for(int j{0}; j < delta_r_2.size(); j++)
      {
        if(i == j){continue;}
        if(delta_r_1[i] + delta_r_2[j] < minimum_value)
        {
          minimum_value = delta_r_1[i] + delta_r_2[j];
          indicies[0] = i;
          indicies[1] = j;
        }
      }
      // this line means that there is only 1 delta r for1 and 2 so there would only be a single 
    }
    // cout<<indicies[0] << " " << indicies[1]  << " " << delta_r_1 << " " << delta_r_2 << "\n";
    if(indicies[0] == indicies[1]){return {};}
    return indicies;
  };
  ROOT::VecOps::RVec<unsigned int> get_delta_r_indicies_truthjet(ROOT::VecOps::RVec<Float_t>& delta_r_1, ROOT::VecOps::RVec<float_t>& delta_r_2, ROOT::VecOps::RVec<Float_t>& tau_jet_weight)
  {
    Float_t minimum_value{10000};
    ROOT::VecOps::RVec<unsigned int> indicies{0, 0};
    if(delta_r_1.size() == 1 || delta_r_2.size() == 1){return {};}

    for(int i{0}; i < delta_r_1.size(); i++)
    {
      for(int j{0}; j < delta_r_2.size(); j++)
      {
        if(i == j){continue;}
        if(delta_r_1[i]/tau_jet_weight[i] + delta_r_2[j]/tau_jet_weight[j] < minimum_value)
        {
          minimum_value = delta_r_1[i] + delta_r_2[j];
          indicies[0] = i;
          indicies[1] = j;
        }
      }
      // this line means that there is only 1 delta r for1 and 2 so there would only be a single 
    }
    // cout<<indicies[0] << " " << indicies[1]  << " " << delta_r_1 << " " << delta_r_2 << "\n";
    if(indicies[0] == indicies[1]){return {};}
    return indicies;
  };
  ROOT::VecOps::RVec<Float_t> use_delr_indicies(ROOT::VecOps::RVec<Float_t> delr1, ROOT::VecOps::RVec<Float_t> delr2, ROOT::VecOps::RVec<unsigned int> delr_indicies)
  {
    return {delr1[delr_indicies[0]], delr2[delr_indicies[1]]};
  }
} // LFuncs

#endif 