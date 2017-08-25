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

//Get a std::vector with the terminal nodes of the neurite
std::vector<const_node_reference> get_terminal_nodes(const neurostr::Neurite& n)
{
  std::vector<const_node_reference> selection=ns::neurite_terminal_selector(n);

  return(selection);
}

//Get the coordinates of the nodes
NumericMatrix get_node_info(std::vector<const_node_reference> terminal_nodes)
{

  auto n = std::distance(terminal_nodes.begin(), terminal_nodes.end());

  // Create nx3 matrix
  NumericMatrix m (n,4);

  // Fill it
  int i = 0;
  for (auto it = terminal_nodes.begin(); it != terminal_nodes.end(); ++it, ++i) {
    m(i, 0) = it->get().x();
    m(i, 1) = it->get().y();
    m(i, 2) = it->get().z();
    m(i, 3) = it->get().id();
  }

  return(m);
}


// [[Rcpp::export]]
NumericMatrix c_get_terminal_nodes(std::string json_info) {
  std::ifstream aux;
  rapidjson::Document doc;
  doc.Parse(json_info.c_str());
  neurostr::io::JSONParser* p = new neurostr::io::JSONParser(aux);

  //Save all the terminal nodes
  std::vector<const_node_reference> terminal_nodes;

  auto r = p->read_string("reconstruction",doc);

  delete p;
  bool first = true;

  for(auto n_it = r->begin(); n_it != r->end(); ++n_it){
    neurostr::Neuron& n = *n_it;
    n.center();
    /*neurostr::Node::id_type my_id=55;
    Rcout<< n.find(my_id).node()->id()<<std::endl;*/

    for(auto it = n.begin_neurite(); it != n.end_neurite(); ++it){
      if(first)
      {
        terminal_nodes=get_terminal_nodes(*it);
        first=false;
      }else{
        auto neurite_terminal_nodes=get_terminal_nodes(*it);
        terminal_nodes.insert(terminal_nodes.end(),neurite_terminal_nodes.begin(),neurite_terminal_nodes.end());
      }
    }
  } // End neuron for

  NumericMatrix terminal_nodes_coords=get_node_info(terminal_nodes);

  return(terminal_nodes_coords);
}
