#ifndef SC
#define SC
namespace SC
{
  // Broad selection cut is defined as a basic cut needed to start analysis.

  bool semileptonic_cut(RV<int> tau_charge, RV<int> lep_charge){
    return (tau_charge.size() == 1 && lep_charge.size() == 1) && (tau_charge[0] + lep_charge[0] == 0);
  }

  bool opp_sign(RV<int> charge){
    return charge[0] + charge[1] == 0;
  }
  bool size_2(ROOT::VecOps::RVec<unsigned int> indicies){
    return indicies.size() == 2;
  }
  bool size_0(ROOT::VecOps::RVec<Float_t> e_pt, ROOT::VecOps::RVec<Float_t> mu_pt){
    return e_pt.size() == 0 && mu_pt.size() == 0;
  }
  bool norm_truth_reco_diff_cut(ROOT::VecOps::RVec<Float_t> diff){
    return abs(diff[0]) < 0.3;
  }
  bool inv_mass_region_z(Float_t mass){
    return 80 <= mass && mass <= 100;
  }
  bool met_angle_diff(Float_t angle_diff){
    return angle_diff < 0.025;
  }
  bool met_angle_diff_fine(Float_t angle_diff){
    return angle_diff < 2 * pow(10, -6);
  }
  bool met_angle_diff_less_fine(Float_t angle_diff){
    return angle_diff < 3 * pow(10, -6);
  }
  bool gen_inv_mass_g20(Float_t mass){
    return mass > 20;
  }
  bool met_g_30(ROOT::VecOps::RVec<Float_t> met)
  {
    return met[0] > 30;
  }
  bool greater_than_0(Float_t pt)
  {
    return pt > 0;
  }
  bool delta_phi_cut_less_2p9(Float_t delta_phi)
  {
    return delta_phi < 2.9;
  }
  
  bool delta_phi_cut_less(Float_t delta_phi){
    return delta_phi < 1.5;
  }
  bool delta_phi_1p2(Float_t delta_phi){
    return delta_phi < 1.2;
  }
  bool delta_phi_1p6(Float_t delta_phi){
    return delta_phi < 1.6;
  }
  bool delta_phi_2(Float_t delta_phi){
    return delta_phi < 2;
  }
  bool delta_phi_1(Float_t delta_phi){
    return delta_phi < 1;
  }
  bool delta_phi_cut_moreq(Float_t delta_phi){
    return delta_phi >=1.5;
  }
  bool jet_number_cut(int size){
    return size > 5;
  }
  bool  deltaRcut0p2(ROOT::VecOps::RVec<Float_t> DeltaR){
    bool result{DeltaR.size() != 0};
    for(auto DR : DeltaR){
      result *= DR  < 0.2;
    }
    return result;
  }
  bool  deltaRcut0p3(ROOT::VecOps::RVec<Float_t> DeltaR){
    bool result{DeltaR.size() != 0};
    for(auto DR : DeltaR){
      result *= DR  < 0.3;
    }
    return result;
  }
  bool truth_and_fake_cut(ROOT::VecOps::RVec<Float_t> DeltaR){
    bool truth{false};
    bool fake{false};
    for(int i{0}; i < 2; i++){
      Float_t DR = DeltaR[i];
      if(!truth && DR < 0.3){truth = true;}
      if(!fake && DR > 0.3){fake = true;}
    }
    return fake && truth;
  }
  bool all_fake(ROOT::VecOps::RVec<Float_t> DeltaR){
    bool result{DeltaR.size() != 0};
    for(auto DR: DeltaR){
      result *= DR > 0.3;
    }
    return result;
  }

  bool pt_tau_cut(ROOT::VecOps::RVec<Float_t> pt){
    return pt[0] > 20 && pt[1] > 20;
  }
  bool pt_g_10_cut(ROOT::VecOps::RVec<Float_t> pt){
    bool result{true};
    for(auto _ : pt) {result *= _ > 10;}
    return result;
  }
  bool y_lt_x(Float_t y, Float_t x){
    return y < x;
  }
  bool y_gt_x(Float_t y, Float_t x){
    return y > x;
  }

  bool y_lt_x_inout(Float_t y, Float_t x, Float_t inv_mass, bool in)
  {
    if(in){
      return y < x && 80 < inv_mass && inv_mass < 100;
    } else{
      return y < x && (80 > inv_mass || inv_mass > 100);
    }
  }
  bool y_gt_x_inout(Float_t y, Float_t x, Float_t inv_mass, bool in)
  {
    if(in){
      return y > x && (80 < inv_mass && inv_mass < 100);
    } else{
      return y > x && (80 > inv_mass || inv_mass > 100);
    }
  }
};
#endif