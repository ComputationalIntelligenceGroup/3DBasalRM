#include <Rcpp.h>

#include <neurostr/io/parser_dispatcher.h>
#include <neurostr/core/log.h>
#include <neurostr/core/neuron.h>

#include <neurostr/core/geometry.h>
#include <neurostr/measure/universal_measure.h>
#include <neurostr/measure/branch_measure.h>
#include <neurostr/measure/neurite_measure.h>
#include <neurostr/measure/node_measure.h>
#include <neurostr/measure/aggregate.h>
#include <neurostr/measure/measure_operations.h>

#include <neurostr/selector/neurite_selector.h>

#include "include/escape_string.hpp"
#include "include/remove_Z_jumps.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <sstream>

namespace ns = neurostr::selector;
namespace nm = neurostr::measure;

using namespace Rcpp;
namespace bg =        boost::geometry;
using point_type =    bg::model::point<float, 3, bg::cs::cartesian>;

Eigen::MatrixXf c_prcomp(std::vector<point_type> v)
{
  auto n_elem=std::distance(v.begin(),v.end());
  Eigen::MatrixXf m(n_elem,3);

  int i=0;
  for(auto it = v.begin(); it != v.end(); ++it,i++){
    m(i, 0) = it->get<0>();
    m(i, 1) = it->get<1>();
    m(i, 2) = it->get<2>();
  }

  Eigen::MatrixXf mean=m.colwise().mean();
  Eigen::MatrixXf p_centered=m-mean;

  Eigen::JacobiSVD<Eigen::MatrixXf> svd=p_centered.jacobiSvd(Eigen::ComputeFullU);
  Eigen::MatrixXf U=svd.matrixU();


  Eigen::MatrixXf S=U.col(2);
Rcout<<"Valores: "<<S(0)<<" "<<S(1)<<" "<<S(2)<<std::endl;
  return(S);
}

std::vector<point_type> soma_node_positions(const neurostr::Neuron &n) {
  std::vector<point_type> v;
  auto n_elem=std::distance(n.begin_soma(),n.end_soma());
  v.reserve(n_elem);
  for(auto it = n.begin_soma(); it != n.end_soma(); ++it)
    v.push_back(it->position());
  return v;
}


// [[Rcpp::plugins(cpp14)]]
std::string escape_string(const std::string& s){
  return "\""+s+"\"";
}
std::string escape_string(const char *c){
  return escape_string(std::string(c));
}

std::map<std::string, std::vector<float>> get_neurite_measures_prueba(const neurostr::Neurite& n,  const std::vector<std::string>& markers){

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
  //Rcout<<"X: "<<n.neuron().up().get<0>()<<"Y: "<<n.neuron().up().get<1>()<<"Z: "<<n.neuron().up().get<2>()<<std::endl;
  auto list_nodes=ns::neurite_node_selector(n);
  auto soma_nodes=soma_node_positions(n.neuron());
  Rcout<<c_prcomp(soma_nodes);
  if(list_nodes.size()>115)
  {
  auto my_node=list_nodes.at(115);
  auto node_parent=my_node.get().parent();
  auto basis=node_parent.local_basis(node_parent.parent(),n.neuron().up());
  auto pos = node_parent.vectorTo(my_node);
  Rcout<<"Parent_parent_vector: "<<std::endl;
  Rcout<<"X="<<node_parent.position().get<0>()<<" Y="<<node_parent.position().get<1>()<<" Z="<<node_parent.position().get<2>()<<std::endl;
  Rcout<<"X="<<node_parent.parent().position().get<0>()<<" Y="<<node_parent.parent().position().get<1>()<<" Z="<<node_parent.parent().position().get<2>()<<std::endl;

  Rcout<<"Basis"<<std::endl;
  Rcout<<"Vector 1: X="<<basis.at(0).get<0>()<<" Y="<<basis.at(0).get<1>()<<" Y="<<basis.at(0).get<2>()<<std::endl;
  Rcout<<"Vector 2: X="<<basis.at(1).get<0>()<<" Y="<<basis.at(1).get<1>()<<" Y="<<basis.at(1).get<2>()<<std::endl;
  Rcout<<"Vector 3: X="<<basis.at(2).get<0>()<<" Y="<<basis.at(2).get<1>()<<" Y="<<basis.at(2).get<2>()<<std::endl;

  Rcout<<"Parent Vector"<<std::endl;
  Rcout<<"X="<<pos.get<0>()<<" Y="<<pos.get<1>()<<" Y="<<pos.get<2>()<<std::endl;

  std::pair<float,float> angulos=neurostr::geometry::local_orientation(pos, basis);
  Rcout<<"Azimuth: "<<angulos.first<<"  Elevacion: "<<angulos.second<<std::endl;
  }
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

// [[Rcpp::export]]
std::vector<std::map<std::string, std::vector<float>>> c_neurite_feature_extractor_prueba(std::string json_info, bool omitapical, bool omitaxon, bool omitdend, bool correct, bool remove_zjumps) {

  std::vector<std::map<std::string, std::vector<float>>> m;
  std::vector<std::string> markers;

  std::ifstream aux;
  rapidjson::Document doc;
  doc.Parse(json_info.c_str());

  neurostr::io::JSONParser* p = new neurostr::io::JSONParser(aux);
  auto r = p->read_string("reconstruction",doc);

  delete p;
  bool first = true;
  // For each neuron
  for(auto n_it = r->begin(); n_it != r->end(); ++n_it){
    neurostr::Neuron& n = *n_it;

    if(omitapical) n.erase_apical();
    if(omitaxon) n.erase_axon();
    if(omitdend) n.erase_dendrites();
    if(correct) n.correct();
    if(remove_zjumps) remove_Z_jumps_neuron(n);

    for(auto it = n.begin_neurite(); it != n.end_neurite(); ++it){
      if(!first){
      }
      first = false;
      auto f=get_neurite_measures_prueba(*it,markers);
      m.emplace_back(f);
    }
  } // End neuron for

  return m;
}


