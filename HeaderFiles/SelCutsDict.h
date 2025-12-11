
#ifndef SCD
#define SCD
#include"SelectionCuts.h"
#include"DataStructures.h"
#include<map>
#include<string>
#include<functional>
#include<variant>
template<typename T>
using RV = ROOT::VecOps::RVec<T>;using std::string;
using std::map;
// using DataStructs::GeneralSelectionCut;
using namespace SC;


struct GeneralSelectionCut {
    using FuncVariant = std::variant<
        std::function<bool(RV<int>, RV<int>)>,
        std::function<bool(RV<int>)>,
        std::function<bool(RV<unsigned int>)>,
        std::function<bool(RV<Float_t>, RV<Float_t>)>,
        std::function<bool(RV<Float_t>)>,
        std::function<bool(Float_t)>,
        std::function<bool(int)>,
        std::function<bool(Float_t, Float_t)>,
        std::function<bool(Float_t, Float_t, Float_t, bool)>
    >;
    
    FuncVariant cut_function;
    
    // Constructor that wraps function in std::function
    template<typename Func>
    GeneralSelectionCut(Func func) {
        cut_function = std::function(func);  // Explicit conversion
    }
    
    GeneralSelectionCut() = default;
    
    template<typename... Args>
    bool operator()(Args&&... args) const {
        return std::visit([&](auto&& func) -> bool {
            return func(std::forward<Args>(args)...);
        }, cut_function);
    }

};

// Helper function to create cuts with explicit type deduction
template<typename R, typename... Args>
GeneralSelectionCut make_cut(R(*func)(Args...)) {
    GeneralSelectionCut cut;
    cut.cut_function = std::function<R(Args...)>(func);
    return cut;
}


std::map<std::string, GeneralSelectionCut> sel_cuts;

// Two RVec<int> arguments
sel_cuts["semileptonic"] = make_cut(semileptonic_cut);

// Single RVec<int> argument
sel_cuts["opp_sign"] = make_cut(opp_sign);

// Single RVec<unsigned int> argument
sel_cuts["size_2"] = make_cut(size_2);
sel_cuts["size_0_ui"] = make_cut(size_0_ui);

// Two RVec<Float_t> arguments
sel_cuts["size_0"] = make_cut(size_0);

// Single RVec<Float_t> argument
sel_cuts["norm_truth_reco_diff"] = make_cut(norm_truth_reco_diff_cut);
sel_cuts["met_g_30"] = make_cut(met_g_30);
sel_cuts["tau_pt_cut"] = make_cut(tau_pt_cut);
sel_cuts["lep_pt_l_25"] = make_cut(lep_pt_l_25);
sel_cuts["deltaR_0p2"] = make_cut(deltaRcut0p2);
sel_cuts["deltaR_0p3"] = make_cut(deltaRcut0p3);
sel_cuts["truth_and_fake"] = make_cut(truth_and_fake_cut);
sel_cuts["all_fake"] = make_cut(all_fake);
sel_cuts["pt_tau"] = make_cut(pt_tau_cut);
sel_cuts["pt_g_10"] = make_cut(pt_g_10_cut);

// Single Float_t argument
sel_cuts["inv_mass_z"] = make_cut(inv_mass_region_z);
sel_cuts["met_angle_diff"] = make_cut(met_angle_diff);
sel_cuts["met_angle_diff_fine"] = make_cut(met_angle_diff_fine);
sel_cuts["met_angle_diff_less_fine"] = make_cut(met_angle_diff_less_fine);
sel_cuts["gen_inv_mass_g20"] = make_cut(gen_inv_mass_g20);
sel_cuts["greater_than_0"] = make_cut(greater_than_0);
sel_cuts["delta_phi_2p9"] = make_cut(delta_phi_cut_less_2p9);
sel_cuts["delta_phi_1p5"] = make_cut(delta_phi_cut_less);
sel_cuts["delta_phi_1p2"] = make_cut(delta_phi_1p2);
sel_cuts["delta_phi_1p6"] = make_cut(delta_phi_1p6);
sel_cuts["delta_phi_2"] = make_cut(delta_phi_2);
sel_cuts["delta_phi_1"] = make_cut(delta_phi_1);
sel_cuts["delta_phi_moreq_1p5"] = make_cut(delta_phi_cut_moreq);

// Single int argument
sel_cuts["jet_num_g5"] = make_cut(jet_number_cut);
sel_cuts["jet_num_g4"] = make_cut(num_jet_4);

// Two Float_t arguments
sel_cuts["y_lt_x"] = make_cut(y_lt_x);
sel_cuts["y_gt_x"] = make_cut(y_gt_x);

// Four arguments (Float_t, Float_t, Float_t, bool)
sel_cuts["y_lt_x_inout"] = make_cut(y_lt_x_inout);
sel_cuts["y_gt_x_inout"] = make_cut(y_gt_x_inout);


#endif