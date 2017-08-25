#include <Rcpp.h>

#include <neurostr/io/parser_dispatcher.h>
#include <neurostr/core/log.h>
#include <neurostr/core/neuron.h>
#include <neurostr/core/node.h>

#include <neurostr/selector/selector.h>
#include <neurostr/selector/universal_selector.h>
#include <neurostr/selector/neurite_selector.h>
#include <neurostr/selector/node_selector.h>

#include <neurostr/measure/universal_measure.h>
#include <neurostr/measure/branch_measure.h>
#include <neurostr/measure/neurite_measure.h>
#include <neurostr/measure/node_measure.h>
#include <neurostr/measure/aggregate.h>
#include <neurostr/measure/measure_operations.h>

namespace ns = neurostr::selector;
namespace nm = neurostr::measure;

using namespace Rcpp;

// [[Rcpp::plugins(cpp14)]]
template <typename T>
using const_selector_reference = std::reference_wrapper<const T>;
using const_node_reference  = const_selector_reference<neurostr::Node>;

//Get terminal nodes when the neurites had been cut with box cutter
std::vector<const_node_reference> get_terminal_cut_nodes(const neurostr::Neurite& n)
{
  //Get those nodes that hasnt been cut
  auto non_cut = ns::diff_selector_factory(ns::neurite_node_selector, ns::compose_selector(ns::property_exists_factory<neurostr::Node>("cut"), ns::neurite_node_selector));
  std::vector<const_node_reference> selection = non_cut(n);

  std::vector<const_node_reference> new_terminals;

  bool is_terminal;
  std::vector<const_node_reference> descendants;
  int num_new_terminals = 0;
  int i = 0;
  for(auto it = selection.begin(); it != selection.end(); ++it, ++i)
  {
    descendants = ns::node_descendants(*it);
    is_terminal = true; //If it has not descendants is terminal
    for(auto it2 = descendants.begin(); it2 != descendants.end(); ++it2)
    {
      is_terminal = is_terminal && it2->get().properties.exists("cut"); //If all its descendants were cut it is terminal
    }
    if(is_terminal)
    {
      new_terminals.emplace_back(*it);
      if(descendants.size() > 0)
      {
        num_new_terminals++;
      }
    }
  }

  return(new_terminals);
}


//Get the coordinates of the nodes
NumericMatrix get_node_coords(std::vector<const_node_reference> terminal_nodes)
{

  auto n = std::distance(terminal_nodes.begin(), terminal_nodes.end());

  // Create nx3 matrix
  NumericMatrix m (n,3);

  // Fill it
  int i = 0;
  for (auto it = terminal_nodes.begin(); it != terminal_nodes.end(); ++it, ++i) {
    m(i, 0) = it->get().x();
    m(i, 1) = it->get().y();
    m(i, 2) = it->get().z();
  }

  return(m);
}


// [[Rcpp::export]]
NumericMatrix c_get_cut_nodes(std::string json_info) {
  std::ifstream aux;
  rapidjson::Document doc;
  doc.Parse(json_info.c_str());
  neurostr::io::JSONParser* p = new neurostr::io::JSONParser(aux);

  std::vector<const_node_reference> nodes_without_descendant;//Nodes that not have descendant, i.e., are terminal. It is used when box cutter has been applied in the neurite.
  //std::vector<const_node_reference> terminal_nodes;//Terminal nodes of neurite before cutting

  auto r = p->read_string("reconstruction",doc);

  delete p;
  bool first = true;

  for(auto n_it = r->begin(); n_it != r->end(); ++n_it){
    neurostr::Neuron& n = *n_it;
    n.center();
    for(auto it = n.begin_neurite(); it != n.end_neurite(); ++it){
      if(first)
      {
        //terminal_nodes=get_terminal_nodes(*it);
        nodes_without_descendant=get_terminal_cut_nodes(*it);
        first=false;
      }else{
        //auto neurite_terminal_nodes=get_terminal_nodes(*it);
        auto neurite_non_cut_nodes=get_terminal_cut_nodes(*it);
        //terminal_nodes.insert(terminal_nodes.end(),neurite_terminal_nodes.begin(),neurite_terminal_nodes.end());
        nodes_without_descendant.insert(nodes_without_descendant.end(),neurite_non_cut_nodes.begin(),neurite_non_cut_nodes.end());

      }
    }
  } // End neuron for

  Rcout<<nodes_without_descendant.size();
  NumericMatrix terminal_nodes_coords=get_node_coords(nodes_without_descendant);

  return(terminal_nodes_coords);
}
