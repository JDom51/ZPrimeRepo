#include<vector>
#include<string>
#include<cmath>
using std::string;
using std::vector;


namespace LFuncs
{
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


  template<typename Type> int get_size(ROOT::VecOps::RVec<Type> vect)
  {
    return vect.size();
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
  };


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
  };
} // LFuncs
