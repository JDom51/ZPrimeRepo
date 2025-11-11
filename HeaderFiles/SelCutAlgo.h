#ifndef SEL_CUT_ALGO
#define SEL_CUT_ALGO

#include"DataStructures.h"

using DataStructs::Integral;
using DataStructs::SelCutFunc;
//RDF.Take(column_name) takes out the values into an actual variable i can work with.

namespace SCAlgo
{
  // get mc data 
  // get integral of distribution with a specific selection cut
  // calculate metric at that value of selection cut.
  // repeat for all values of the selection cut.
  template<typename T> vector<vector<Integral<T>>> get_integral(FrameAndData& fd, string column_name, string inv_mass_column_name, SelCutFunc selection_cut, vector<T> sel_cut_values, size_t resolution){
    vector<Integral<T>> integrals_in{};
    vector<Integral<T>> integrals_out{};
    auto cutting_variable{fd.node.Take<T>(column_name)};
    auto inv_mass{fd.node.Take<Float_t>(inv_mass_column_name)};
    // cout<<"Data size: "<<cutting_variable->size()<<"\n";
    // cout<<"IMassSize: "<<inv_mass->size()<<"\n";
    T dif{(sel_cut_values[1] - sel_cut_values[0]) / resolution};
    for(auto sel_cut_value{sel_cut_values[0]}; sel_cut_value < sel_cut_values[1]; sel_cut_value += dif){
      cout << sel_cut_value << "\n";
      Integral<T> integral_in{sel_cut_value};
      Integral<T> integral_out{sel_cut_value};
      for(int i{0}; i < cutting_variable->size(); i++){
      //   for(auto value : cutting_variable->at(i)){
          T value = cutting_variable->at(i);  
          // cout<<i<<" "<<value<<"\n";
          integral_in.add(selection_cut(value, sel_cut_value, inv_mass->at(i), true));
          integral_out.add(selection_cut(value, sel_cut_value, inv_mass->at(i), false));
      //   }
      }
      integrals_in.push_back(integral_in);
      integrals_out.push_back(integral_out);
    }
  // return {};
  return {integrals_in, integrals_out};
  }

};
#endif
