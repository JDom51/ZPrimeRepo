#ifndef DATA_STRUCTURES
#define DATA_STRUCTURES

#include<string>
#include<vector>
#include<map>

using std::string;
using std::vector;
using std::map;


namespace DataStructs
{
  struct FrameAndData
  {
    ROOT::RDataFrame frame;
    ROOT::RDF::RNode node;
    string save_string;
    
    FrameAndData(string mode, vector<string> files) : frame{mode, files}, node{frame}, save_string{} {}
    ROOT::RDF::RResultPtr<TH1D> Histo1D(string hist_id, string hist_name, int nbins, double lbound, double ubound, string plotting_column)
    {
      return node.Histo1D({hist_id.c_str(), hist_name.c_str(), nbins, lbound, ubound}, plotting_column);
    }
    template<typename T>
    void Define(string column_name, T lambda_func, vector<string> var_names){
      node = node.Define(column_name, lambda_func, var_names);
    }
    template<typename Type> void Filter(Type func, vector<string> column_names, string cut_name)
    {
      node = node.Filter(func, column_names);
      save_string += " " + cut_name;
      cout<<save_string<<"\n";
    }
  };
  struct HistInfo
  {
    string name;
    string x_axis;
    vector<string> columns;
    int nbins;
    double lbound;
    double ubound;
    string units;
    bool take_point_five;
    unsigned int mode; // 0 : histo1d, 1 : stacked histo1d, 2: scatter plot 
    HistInfo(string nam, string x_ax, string m_column, int bins, double lowbound, double upbound, string unit) : 
      name{nam}, x_axis{x_ax}, columns{{m_column}}, nbins{bins}, lbound{lowbound}, ubound{upbound}, units{unit}, take_point_five{false}, mode{0} {}
    HistInfo(string nam, string x_ax, string m_column, int bins, double lowbound, double upbound, string unit, bool m_take_point_five) : 
      name{nam}, x_axis{x_ax},  columns{{m_column}}, nbins{bins}, lbound{lowbound}, ubound{upbound}, units{unit}, take_point_five{m_take_point_five}, mode{0} {}
    HistInfo(string nam, string x_ax, vector<string> m_columns, int bins, double lowbound, double upbound, string unit, unsigned int m_mode) : 
      name{nam}, x_axis{x_ax}, columns{m_columns}, nbins{bins}, lbound{lowbound}, ubound{upbound}, units{unit}, take_point_five{false}, mode{m_mode} {}
    void update_column(vector<string> override_column){
      columns=override_column;
    }
  };
  struct HistInfoPair{
    HistInfo info;
    vector<string> override_column;
    HistInfoPair(HistInfo m_info, vector<string> m_override_column={}) : info{m_info}, override_column{m_override_column}{
      if(m_override_column.size() != 0){
        info.update_column(override_column);
      }
    };
  };
  typedef bool (*BroadCut)(ROOT::VecOps::RVec<unsigned int> indicies);
  typedef bool (*SelCutFunc)(Float_t y, Float_t x, Float_t m, bool in);

  struct SelectionCut{
    BroadCut selection_cut;
    vector<string> variables;
    SelectionCut(BroadCut m_selection_cut, vector<string> m_variables) : selection_cut{m_selection_cut}, variables{m_variables} {};
  };
  
  // template<typename T> using GenSelectionCut = bool (*)(vector<T> columns);

  template<typename T> struct Integral{
    T value;
    size_t integral;
    Integral(T sel_cut_value) : value{sel_cut_value}, integral{0} {};
    void add(bool result){integral += result;}
  };

  typedef void (*AnalysisFunc)(FrameAndData& fd);
  struct DataBase{
    vector<string> ifiles;
    AnalysisFunc a_func;
    vector<SelectionCut> sel_cuts;
    string fname;
    vector<HistInfo> hists;
  };
};
#endif