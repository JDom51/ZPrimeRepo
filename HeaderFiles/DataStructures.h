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
    FrameAndData(ROOT::RDataFrame&& rdf) : frame{rdf}, node{frame} {}
    ROOT::RDF::RResultPtr<TH1D> Histo1D(string hist_id, string hist_name, int nbins, double lbound, double ubound, string plotting_column)
    {
      return node.Histo1D({hist_id.c_str(), hist_name.c_str(), nbins, lbound, ubound}, plotting_column);
    }
    ROOT::RDF::RResultPtr<TH2D> Histo2D(string hist_id, string hist_name, 
                                         int nbins_x, double lbound_x, double ubound_x,
                                         int nbins_y, double lbound_y, double ubound_y,
                                         string x_column, string y_column)
    {
      return node.Histo2D({hist_id.c_str(), hist_name.c_str(), 
                           nbins_x, lbound_x, ubound_x,
                           nbins_y, lbound_y, ubound_y}, 
                          x_column, y_column);
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
    unsigned int mode; // 0: normal 1D hist, 1: stacked 1D hist, 2: 2D histogram
    bool data_vs_mc; // For mode 1: if true, first column is data (crosses), rest are MC (stacked)
    
    // For 2D histograms
    string y_axis;
    int nbins_y;
    double lbound_y;
    double ubound_y;
    string units_y;
    
    // Constructor for normal 1D histogram (mode == 0)
    HistInfo(string nam, string x_ax, string m_column, int bins, double lowbound, double upbound, string unit) : 
      name{nam}, x_axis{x_ax}, columns{{m_column}}, nbins{bins}, lbound{lowbound}, ubound{upbound}, 
      units{unit}, take_point_five{false}, mode{0}, data_vs_mc{false},
      y_axis{""}, nbins_y{0}, lbound_y{0}, ubound_y{0}, units_y{""} {}
    
    // Constructor for normal 1D histogram with offset option (mode == 0)
    HistInfo(string nam, string x_ax, string m_column, int bins, double lowbound, double upbound, 
             string unit, bool m_take_point_five) : 
      name{nam}, x_axis{x_ax}, columns{{m_column}}, nbins{bins}, lbound{lowbound}, ubound{upbound}, 
      units{unit}, take_point_five{m_take_point_five}, mode{0}, data_vs_mc{false},
      y_axis{""}, nbins_y{0}, lbound_y{0}, ubound_y{0}, units_y{""} {}
    
    // Constructor for stacked histogram (mode == 1)
    HistInfo(string nam, string x_ax, vector<string> m_columns, int bins, double lowbound, 
             double upbound, string unit, unsigned int m_mode) : 
      name{nam}, x_axis{x_ax}, columns{m_columns}, nbins{bins}, lbound{lowbound}, ubound{upbound}, 
      units{unit}, take_point_five{false}, mode{m_mode}, data_vs_mc{false},
      y_axis{""}, nbins_y{0}, lbound_y{0}, ubound_y{0}, units_y{""} {}
    
    // Constructor for stacked histogram with data vs MC option (mode == 1)
    HistInfo(string nam, string x_ax, vector<string> m_columns, int bins, double lowbound, 
             double upbound, string unit, unsigned int m_mode, bool m_data_vs_mc) : 
      name{nam}, x_axis{x_ax}, columns{m_columns}, nbins{bins}, lbound{lowbound}, ubound{upbound}, 
      units{unit}, take_point_five{false}, mode{m_mode}, data_vs_mc{m_data_vs_mc},
      y_axis{""}, nbins_y{0}, lbound_y{0}, ubound_y{0}, units_y{""} {}
    
    // Constructor for 2D histogram (mode == 2)
    HistInfo(string nam, string x_ax, string y_ax, string x_column, string y_column,
             int bins_x, double lowbound_x, double upbound_x, string unit_x,
             int bins_y, double lowbound_y, double upbound_y, string unit_y) :
      name{nam}, x_axis{x_ax}, columns{{x_column, y_column}}, 
      nbins{bins_x}, lbound{lowbound_x}, ubound{upbound_x}, units{unit_x},
      take_point_five{false}, mode{2}, data_vs_mc{false},
      y_axis{y_ax}, nbins_y{bins_y}, lbound_y{lowbound_y}, ubound_y{upbound_y}, units_y{unit_y} {}
    
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
