#include <Rcpp.h>

#include <neurostr/io/parser_dispatcher.h>
#include <neurostr/core/log.h>
#include <neurostr/core/neuron.h>

#include <neurostr/measure/universal_measure.h>
#include <neurostr/measure/branch_measure.h>
#include <neurostr/measure/neurite_measure.h>
#include <neurostr/measure/node_measure.h>
#include <neurostr/measure/aggregate.h>
#include <neurostr/measure/measure_operations.h>

#include <neurostr/selector/neurite_selector.h>

#include "include/escape_string.hpp"

#include <sstream>

namespace ns = neurostr::selector;
namespace nm = neurostr::measure;

using namespace Rcpp;

// [[Rcpp::plugins(cpp14)]]
std::string escape_string(const std::string& s){
  return "\""+s+"\"";
}
std::string escape_string(const char *c){
  return escape_string(std::string(c));
}

std::map<std::string, std::vector<float>> get_neurite_measures(const neurostr::Neurite& n,  const std::vector<std::string>& markers){

  std::map<std::string, std::vector<float>> m; // measures

  // Aux vector for single values
  std::vector<float> aux;

  // Number of bifurcations
  float nbifs = ns::neurite_bifurcation_selector(n).size();
  aux.push_back(nbifs);
  m.emplace( "N_bifurcations", aux );

  // Number of branches
  aux.clear();
  aux.push_back(n.size());
  m.emplace( "N_branches", aux );

  // Number of nodes
  float nnodes = ns::neurite_node_selector(n).size();
  aux.clear();
  aux.push_back(nnodes);
  m.emplace( "N_nodes", aux );

  // Node (compartment) length
  m.emplace( "node_length" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_length_to_parent))(n));

  //Node local orientation, azimuth and elevation

  std::vector<std::pair<float, float> > angles = nm::selectorMeasureCompose(ns::neurite_node_selector, nm::measureEach(nm::node_local_orientation))(n);

  std::vector <float> azimuth (angles.size());
  std::vector <float> elevation (angles.size());

  for(int i=0;i<angles.size();i++)
  {
    azimuth[i]=angles.at(i).first;
    elevation[i]=angles.at(i).second;
  }

  m.emplace( "azimuth_angle", azimuth );
  m.emplace( "elevation_angle", elevation );

  auto node_order = nm::selectorMeasureCompose(ns::neurite_node_selector,
                                                 nm::measureEach(nm::node_order))(n);
  aux.clear();
  aux.insert(aux.end(),node_order.begin(),node_order.end());

  m.emplace( "node_order", aux);

  auto node_num_descendant = nm::selectorMeasureCompose(ns::neurite_node_selector,
                                               nm::measureEach(nm::desc_count))(n);
  aux.clear();
  aux.insert(aux.end(),node_num_descendant.begin(),node_num_descendant.end());

  m.emplace( "node_num_descendant", aux);

  // Neurite box volume
  aux.clear();
  aux.push_back(nm::selectorMeasureCompose(ns::neurite_node_selector,
                                           nm::box_volume)(n));
  m.emplace( "box_volume", aux);

  // Node euclidean distance to root
  m.emplace( "node_root_dist", nm::selectorMeasureCompose(ns::neurite_node_selector,
                                            nm::measureEach(nm::node_distance_to_root))(n));

  // Node path distance to root
  m.emplace( "node_root_path", nm::selectorMeasureCompose(ns::neurite_node_selector,
                                            nm::measureEach(nm::node_path_to_root))(n));

  // Branch length
  m.emplace( "branch_length", nm::selectorMeasureCompose(ns::neurite_branch_selector,
                                            nm::measureEach(nm::branch_length))(n));

  // Branch volume
  m.emplace( "branch_volume", nm::selectorMeasureCompose(ns::neurite_branch_selector,
                                            nm::measureEach(nm::selectorMeasureCompose(ns::branch_node_selector,
                                                                                       nm::measureEachAggregate( nm::node_volume,
                                                                                                                 nm::aggregate::sum_aggr_factory<float,float>(0.0)))))(n));

  // Branch surface
  m.emplace( "branch_surface", nm::selectorMeasureCompose(ns::neurite_branch_selector,
                                            nm::measureEach(nm::selectorMeasureCompose(ns::branch_node_selector,
                                                                                       nm::measureEachAggregate( nm::node_compartment_surface,
                                                                                                                 nm::aggregate::sum_aggr_factory<float,float>(0.0)))))(n));

  // Terminal branch length
  m.emplace( "terminal_branch_length", nm::selectorMeasureCompose(ns::neurite_terminal_branch_selector,
                                            nm::measureEach(nm::branch_length)) (n));

  // Terminal branch order

  // Auxiliar - branch order outputs integers
  auto orders = nm::selectorMeasureCompose(ns::neurite_terminal_branch_selector,
                                           nm::measureEach(nm::branch_order))(n);
  aux.clear();
  aux.insert(aux.end(),orders.begin(),orders.end());

  m.emplace( "terminal_branch_order", aux);

  // Terminal nodes distance to root
  m.emplace( "terminal_nodes_root_dist", nm::selectorMeasureCompose(ns::neurite_terminal_selector,
                                            nm::measureEach(nm::node_distance_to_root)) (n));

  // Terminal nodes path to root
  m.emplace( "terminal_nodes_root_path", nm::selectorMeasureCompose(ns::neurite_terminal_selector,
                                            nm::measureEach(nm::node_path_to_root)) (n));

  // Branch tortuosity
  m.emplace( "branch_tortuosity", nm::selectorMeasureCompose(ns::neurite_branch_selector,
                                            nm::measureEach(nm::tortuosity))(n));

  // Hillman taper rate
  m.emplace( "hill_taper_rate", nm::selectorMeasureCompose(ns::neurite_branch_selector,
                                            nm::measureEach(nm::taper_rate_hillman))(n));
  // Burker taper rate
  m.emplace( "burker_taper_rate", nm::selectorMeasureCompose(ns::neurite_branch_selector,
                                            nm::measureEach(nm::taper_rate_burker))(n));

  // Branch fractal dimension
  m.emplace( "branch_fractal_dimension", nm::selectorMeasureCompose( ns::neurite_branch_selector,
                                             nm::measureEach( nm::branch_fractal_dim))(n));

  /** Bifurcation measures **/
  if(n.size() > 1){

    // Local bifurcation angle
    m.emplace( "local_bifurcation_angle", nm::selectorMeasureCompose(ns::neurite_non_terminal_branches,
                                              nm::measureEach(nm::local_bifurcation_angle)) (n));

    // Local tilt angle
    m.emplace( "local_tilt_angle", nm::selectorMeasureCompose(ns::neurite_non_terminal_branches,
                                              nm::measureEach(nm::local_tilt_angle)) (n));

    // Local torque angle
    m.emplace( "local_torque_angle", nm::selectorMeasureCompose(ns::neurite_non_terminal_branches,
                                              nm::measureEach(nm::local_torque_angle)) (n));

    // Remote bifurcation angle
    m.emplace( "remote_bifurcation_angle", nm::selectorMeasureCompose(ns::neurite_non_terminal_branches,
                                              nm::measureEach(nm::remote_bifurcation_angle)) (n));

    // Remote tilt angle
    m.emplace( "remote_tilt_angle", nm::selectorMeasureCompose(ns::neurite_non_terminal_branches,
                                              nm::measureEach(nm::remote_tilt_angle)) (n));

    // Remote torque angle
    m.emplace( "remote_torque_angle", nm::selectorMeasureCompose(ns::neurite_non_terminal_branches,
                                              nm::measureEach(nm::remote_torque_angle)) (n));

    // Child diameter ratio
    m.emplace( "child_diam_ratio", nm::selectorMeasureCompose(ns::neurite_non_terminal_branches,
                                              nm::measureEach(nm::child_diam_ratio)) (n));
    // Partition asymmetry
    m.emplace( "partition_asymmetry", nm::selectorMeasureCompose(ns::neurite_non_terminal_branches,
                                              nm::measureEach( nm::partition_asymmetry )) (n));
  } else {
    NSTR_LOG_(info, std::string("Neurite ") + std::to_string(n.id()) + " doesn't have any bifuractions. Bifurcations measures are skipped.." );
  }

  // MARKERS

  for(auto mark_name = markers.begin(); mark_name != markers.end(); ++mark_name){

    aux.clear();
    aux.push_back(static_cast<float>(ns::neurite_marker_selector(*mark_name)(n).size()));
    m.emplace(std::string("marker_count_")+*mark_name, aux);
  }

  return m;

}

void print_neurite_id(const neurostr::Neurite& n, std::ostream& os){
  os << escape_string("neuron") << " : " << escape_string(n.neuron().id()) << ", ";
  os << escape_string("neurite") << " : " << n.id() << ", ";
  os << escape_string("neurite_type") << " : ";
  if(n.type() == neurostr::NeuriteType::kAxon){
    os << escape_string("Axon");
  } else if(n.type() == neurostr::NeuriteType::kApical){
    os << escape_string("Apical");
  } else if(n.type() == neurostr::NeuriteType::kDendrite){
    os << escape_string("Dendrite");
  } else {
    os << escape_string("Unknown");
  }
}

// Note: This should be done with rapidjson
void print_vector_measures(std::map<std::string, std::vector<float>>& m ,
                           std::ostream& os ){
  bool first = true;
  // Measures json element
  os << escape_string("measures") << " : { ";

  // Print each measure
  for(auto it = m.begin(); it!=m.end();++it ){

    // Values vector reference
    std::vector<float>& v = it->second;

    // Remove nans (FIX! nans shouldnt appear here)
    for(auto val = v.begin(); val != v.end(); ++val){
      if( std::isnan(*val) ){
        val = v.erase(val)-1;
        NSTR_LOG_(info, std::string("Nan value removed in measure ") + it->first );
      }
    }

    // If values vector is not empty
    if(v.size() > 0){
      if(first){
        first = false;
      } else {
        os << ", ";
      }

      // Print key
      os << escape_string(it->first) << " : " ;

      // Print value vector
      os << "[ " ;
      os << v.front();

      // Print rest
      for( auto val = std::next(v.begin(),1); val != v.end(); ++val){
        os << ", " << *val;
      }

      // Close array
      os << "]";

    } // End if vector is empty
  } // End for loop
  os << " }"; // Close measures
}

void print_neurite_measures(const neurostr::Neurite& n, const std::vector<std::string>& markers, std::ostream& os){
  os << "{" ;
  // Print neurite ID
  print_neurite_id(n,os);
  os << ", ";

  // Get measures
  auto m = get_neurite_measures(n,markers);

  // Print them
  print_vector_measures(m,os);

  // End obj
  os << "}";
}

// [[Rcpp::export]]
std::string c_neurite_feature_extractor(std::string json_info, bool omitapical, bool omitaxon, bool omitdend, bool correct) {

  std::vector<std::string> markers;
  std::ostringstream oss;

  std::ifstream aux;
  rapidjson::Document doc;
  doc.Parse(json_info.c_str());

  neurostr::io::JSONParser* p = new neurostr::io::JSONParser(aux);
  auto r = p->read_string("reconstruction",doc);

  delete p;
  bool first = true;
  oss << "[" << std::endl;
  // For each neuron
  for(auto n_it = r->begin(); n_it != r->end(); ++n_it){
    neurostr::Neuron& n = *n_it;

    if(omitapical) n.erase_apical();
    if(omitaxon) n.erase_axon();
    if(omitdend) n.erase_dendrites();
    if(correct) n.correct();

    for(auto it = n.begin_neurite(); it != n.end_neurite(); ++it){
      if(!first){
        oss << " , ";
      }
      first = false;
      print_neurite_measures(*it, markers, oss);
    }
  } // End neuron for

  oss << "]" << std::endl;

  std::string features = oss.str();


  return features;
}


