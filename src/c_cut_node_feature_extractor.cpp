/* This code computes features for each node. It is used to characterize the dendrite at node level.
 * We combine this measures to obtain a model that represents the dendrite basal arborization.
 */
#include <Rcpp.h>

#include <neurostr/io/parser_dispatcher.h>
#include <neurostr/core/log.h>
#include <neurostr/core/neuron.h>

#include <neurostr/core/geometry.h>

#include "include/escape_string.hpp"
#include "include/node_measures_extension.hpp"
#include "include/remove_Z_jumps.hpp"

#include <sstream>

namespace ns = neurostr::selector;
namespace nm = neurostr::measure;

using namespace Rcpp;
namespace bg =        boost::geometry;
using point_type =    bg::model::point<float, 3, bg::cs::cartesian>;

// [[Rcpp::plugins(cpp14)]]
std::string escape_string(const std::string& s){
  return "\""+s+"\"";
}
std::string escape_string(const char *c){
  return escape_string(std::string(c));
}

template <typename T>
using const_selector_reference = std::reference_wrapper<const T>;
using const_node_reference  = const_selector_reference<neurostr::Node>;


std::map<std::string, std::vector<float>> get_node_measures_by_id(neurostr::Neuron& n, NumericVector id_nodes){
  std::map<std::string, std::vector<float>> m; // measures

  // Aux vector for single values
  std::vector<float> aux;

  //Search nodes by id
  std::vector<const_node_reference> neuron_nodes;
  for(int i=0;i<id_nodes.length();++i)
  {
    neurostr::Node::id_type id_node=id_nodes.at(i);
    auto it_node=n.find(id_node);
    if(it_node.begin()!=it_node.end())
    {
      neuron_nodes.emplace_back(*it_node.node());
    }
  }

  //Definition of variables to compute
  std::vector<float> compartment_length;
  std::vector<float> node_to_brach_root_dist;
  std::vector<float> length_to_branch_root;
  std::vector<float> fract_dim_to_branch_root;
  std::vector<float> tortuosity;

  std::pair<float, float> angles;
  std::vector <float> azimuth_angle;
  std::vector <float> elevation_angle;

  std::vector<float> node_order;
  std::vector<float> node_num_descendant;
  std::vector<float> distance_to_root;
  std::vector<float> path_to_root;

  std::vector<float> subtree_box_volume;
  std::vector<float> subtree_length;
  std::vector<float> subtree_no_bif;
  std::vector<float> subtree_terminals;
  std::vector<float> subtree_width;
  std::vector<float> subtree_height;
  std::vector<float> subtree_depth;
  std::vector<float> subtree_max_order;
  std::vector<float> subtree_min_order;
  std::vector<float> subtree_max_length;
  std::vector<float> subtree_min_length;

  //For each node compute its measures
  for(auto it = neuron_nodes.begin(); it != neuron_nodes.end(); ++it)
  {
    compartment_length.emplace_back(nm::node_length_to_parent(*it));
    node_to_brach_root_dist.emplace_back(nm::node_distance_to_branch_root(*it));
    length_to_branch_root.emplace_back(nm::node_length_to_branch_root(*it));
    fract_dim_to_branch_root.emplace_back(nm::node_fract_dim_to_branch_root(*it));
    tortuosity.emplace_back(nm::node_branch_tortuosity(*it));

    angles=nm::node_local_orientation(*it);
    azimuth_angle.emplace_back(angles.first);
    elevation_angle.emplace_back(angles.second);

    node_order.emplace_back(nm::node_order(*it));
    node_num_descendant.emplace_back(nm::desc_count(*it));
    distance_to_root.emplace_back(nm::node_distance_to_root(*it));
    path_to_root.emplace_back(nm::node_path_to_root(*it));

    subtree_box_volume.emplace_back(nm::node_subtree_box_volume(*it));
    subtree_length.emplace_back(nm::node_subtree_length(*it));
    subtree_no_bif.emplace_back(nm::node_subtree_bifurcations(*it));
    subtree_terminals.emplace_back(nm::node_subtree_terminals(*it));
    subtree_width.emplace_back(nm::node_subtree_width(*it));
    subtree_height.emplace_back(nm::node_subtree_height(*it));
    subtree_depth.emplace_back(nm::node_subtree_depth(*it));
    subtree_max_order.emplace_back(nm::node_subtree_max_order(*it));
    subtree_min_order.emplace_back(nm::node_subtree_min_order(*it));
    subtree_max_length.emplace_back(nm::node_subtree_max_length(*it));
    subtree_min_length.emplace_back(nm::node_subtree_min_length(*it));
  }

  m.emplace("compartment_length",compartment_length);
  m.emplace("node_to_brach_root_dist",node_to_brach_root_dist);
  m.emplace("length_to_branch_root",length_to_branch_root);
  m.emplace("fract_dim_to_branch_root",fract_dim_to_branch_root);
  m.emplace("tortuosity",tortuosity);
  m.emplace("azimuth_angle",azimuth_angle);
  m.emplace("elevation_angle",elevation_angle);
  m.emplace("node_order",node_order);
  m.emplace("node_num_descendant",node_num_descendant);
  m.emplace("distance_to_root",distance_to_root);
  m.emplace("path_to_root",path_to_root);
  m.emplace("subtree_box_volume",subtree_box_volume);
  m.emplace("subtree_length",subtree_length);
  m.emplace("subtree_no_bif",subtree_no_bif);
  m.emplace("subtree_terminals",subtree_terminals);
  m.emplace("subtree_width",subtree_width);
  m.emplace("subtree_height",subtree_height);
  m.emplace("subtree_height",subtree_height);
  m.emplace("subtree_depth",subtree_depth);
  m.emplace("subtree_max_order",subtree_max_order);
  m.emplace("subtree_min_order",subtree_min_order);
  m.emplace("subtree_max_length",subtree_max_length);
  m.emplace("subtree_min_length",subtree_min_length);

  return m;

}

// [[Rcpp::export]]
std::vector<std::map<std::string, std::vector<float>>> c_cut_node_feature_extractor(std::string json_info, NumericVector id_nodes, bool omitapical, bool omitaxon, bool omitdend, bool correct, bool remove_zjumps) {

  std::vector<std::map<std::string, std::vector<float>>> m;
  std::vector<std::string> markers;

  //Read reconstruction from JSON
  std::ifstream aux;
  rapidjson::Document doc;
  doc.Parse(json_info.c_str());
  neurostr::io::JSONParser* p = new neurostr::io::JSONParser(aux);
  auto r = p->read_string("reconstruction",doc);
  delete p;

  // Compute the measures of the nodes which its id has been introduced by the user
  for(auto n_it = r->begin(); n_it != r->end(); ++n_it){
    neurostr::Neuron& n = *n_it;
    if(omitapical) n.erase_apical();
    if(omitaxon) n.erase_axon();
    if(omitdend) n.erase_dendrites();
    if(correct) n.correct();
    if(remove_zjumps) remove_Z_jumps_neuron(n);
    auto f = get_node_measures_by_id(n,id_nodes);
    m.emplace_back(f);
  } // End neuron for

  return m;
}


