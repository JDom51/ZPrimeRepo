#ifndef SC
#define SC
namespace SC
{
  // Broad selection cut is defined as a basic cut needed to start analysis.
  bool size_2(ROOT::VecOps::RVec<unsigned int> indicies){
    return indicies.size() == 2;
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
    return size > 4;
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
  bool pt_tau_cut(ROOT::VecOps::RVec<Float_t> pt){
    return pt[0] > 20 && pt[1] > 20;
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