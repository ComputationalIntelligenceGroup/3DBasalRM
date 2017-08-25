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



std::map<std::string, std::vector<float>> get_node_measures(const neurostr::Neurite& n,  const std::vector<std::string>& markers){
  std::map<std::string, std::vector<float>> m; // measures

  // Aux vector for single values
  std::vector<float> aux;

  // Node (compartment) length
  m.emplace( "compartment_length" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_length_to_parent))(n));

  // Node distance in a straight line to the root of the branch
  m.emplace( "node_to_brach_root_dist" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_distance_to_branch_root))(n));

  // Node length to branch root
  m.emplace( "length_to_branch_root" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_length_to_branch_root))(n));
  // Fractal distance to branch root
  m.emplace( "fract_dim_to_branch_root" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_fract_dim_to_branch_root))(n));

  //Tortuosity from the node to the root of the branch
  m.emplace( "tortuosity" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_branch_tortuosity))(n));

  //Partition asymmetry from the node to the root of the branch
  //m.emplace( "partition_asymmetry" , nm::selectorMeasureCompose(ns::neurite_node_selector,
  //                                           nm::measureEach(nm::node_branch_partition_asymmetry))(n));

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

  //Node centrifugal order
  auto node_order = nm::selectorMeasureCompose(ns::neurite_node_selector,
                                               nm::measureEach(nm::node_order))(n);
  aux.clear();
  aux.insert(aux.end(),node_order.begin(),node_order.end());

  m.emplace( "node_order", aux);

  //Number of descendants of the node
  auto node_num_descendant = nm::selectorMeasureCompose(ns::neurite_node_selector,
                                                        nm::measureEach(nm::desc_count))(n);
  aux.clear();
  aux.insert(aux.end(),node_num_descendant.begin(),node_num_descendant.end());

  m.emplace( "node_num_descendant", aux);

  // Node euclidean distance to root
  m.emplace( "distance_to_root", nm::selectorMeasureCompose(ns::neurite_node_selector,
                                            nm::measureEach(nm::node_distance_to_root))(n));

  // Node path distance to root
  m.emplace( "path_to_root", nm::selectorMeasureCompose(ns::neurite_node_selector,
                                            nm::measureEach(nm::node_path_to_root))(n));

  //Subtree measures
  m.emplace( "subtree_box_volume" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_box_volume))(n));

  m.emplace( "subtree_length" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_length))(n));

  m.emplace( "subtree_no_bif" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_bifurcations))(n));

  m.emplace( "subtree_terminals" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_terminals))(n));

  m.emplace( "subtree_width" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_width))(n));

  m.emplace( "subtree_height" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_height))(n));

  m.emplace( "subtree_depth" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_depth))(n));

  //m.emplace( "node_max_distance_bw_nodes" , nm::selectorMeasureCompose(ns::neurite_node_selector,
  //                                           nm::measureEach(nm::node_max_distance_bw_nodes))(n));

  m.emplace( "subtree_max_order" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_max_order))(n));

  m.emplace( "subtree_min_order" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_min_order))(n));

  m.emplace( "subtree_max_length" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_max_length))(n));

  m.emplace( "subtree_min_length" , nm::selectorMeasureCompose(ns::neurite_node_selector,
                                             nm::measureEach(nm::node_subtree_min_length))(n));

  /*Compute descendants azimuth and elevation angles and length*/
  auto list_nodes=ns::neurite_node_selector(n);//Get all nodes
  std::vector <float> desc_azimuth (list_nodes.size(),0);
  std::vector <float> desc_elevation (list_nodes.size(),0);
  std::vector <float> desc_length (list_nodes.size(),0);
  std::vector <float> desc_node_num_desc (list_nodes.size(),0);

  std::vector <float> desc_azimuth2 (list_nodes.size(),0);
  std::vector <float> desc_elevation2 (list_nodes.size(),0);
  std::vector <float> desc_length2 (list_nodes.size(),0);
  std::vector <float> desc_node_num_desc2 (list_nodes.size(),0);
  std::vector <float> desc_longer (list_nodes.size(),1);

  int num_desc=0;
  for(auto i=0;i<list_nodes.size();++i)
  {
    auto descendants=ns::node_descendants(list_nodes.at(i));
    num_desc=descendants.size();
    if(num_desc>0)
    {
//Ordenar desc de acuerdo al angulo azimuth
      for(int j=0;j<num_desc;++j)
      {
        auto desc_angles=nm::node_local_orientation(descendants.at(j));
        auto desc_l=nm::node_length_to_parent(descendants.at(j));
        auto desc_num_desc=nm::desc_count(descendants.at(j));
        if(j==0)
        {
          desc_azimuth[i]=desc_angles.first;
          desc_elevation[i]=desc_angles.second;
          desc_length[i]=desc_l;
          desc_node_num_desc[i]=desc_num_desc;
        }else{
          if(desc_azimuth[i]>desc_angles.first)
          {
            desc_azimuth2[i]=desc_angles.first;
            desc_elevation2[i]=desc_angles.second;
            desc_length2[i]=desc_l;
            desc_node_num_desc2[i]=desc_num_desc;
            if(desc_length2[i]>desc_length[i])
            {
              desc_longer[i]=2;
            }
          }else{
            desc_azimuth2[i]=desc_azimuth[i];
            desc_elevation2[i]=desc_elevation[i];
            desc_length2[i]=desc_length[i];
            desc_node_num_desc2[i]=desc_node_num_desc[i];

            desc_azimuth[i]=desc_angles.first;
            desc_elevation[i]=desc_angles.second;
            desc_length[i]=desc_l;
            desc_node_num_desc[i]=desc_num_desc;
          }
        }
      }
      if(desc_length2[i]>desc_length[i])
      {
        desc_longer[i]=2;
      }
    }
  }

  m.emplace("desc_azimuth_angle", desc_azimuth );
  m.emplace("desc_elevation_angle", desc_elevation);
  m.emplace("desc_length", desc_length);
  m.emplace("desc_azimuth_angle2", desc_azimuth2);
  m.emplace("desc_elevation_angle2", desc_elevation2);
  m.emplace("desc_length2", desc_length2);
  m.emplace("desc_longer", desc_longer);

  // MARKERS
  for(auto mark_name = markers.begin(); mark_name != markers.end(); ++mark_name){

    aux.clear();
    aux.push_back(static_cast<float>(ns::neurite_marker_selector(*mark_name)(n).size()));
    m.emplace(std::string("marker_count_")+*mark_name, aux);
  }

  return m;

}

// [[Rcpp::export]]
std::vector<std::map<std::string, std::vector<float>>> c_node_feature_extractor(std::string json_info, bool omitapical, bool omitaxon, bool omitdend, bool correct, bool remove_zjumps) {

  std::vector<std::map<std::string, std::vector<float>>> m;
  std::vector<std::string> markers;

  std::ifstream aux;
  rapidjson::Document doc;
  doc.Parse(json_info.c_str());
  neurostr::io::JSONParser* p = new neurostr::io::JSONParser(aux);
  auto r = p->read_string("reconstruction",doc);

  delete p;


  // For each neuron
  for(auto n_it = r->begin(); n_it != r->end(); ++n_it){
    neurostr::Neuron& n = *n_it;
    if(omitapical) n.erase_apical();
    if(omitaxon) n.erase_axon();
    if(omitdend) n.erase_dendrites();
    if(correct) n.correct();
    if(remove_zjumps) remove_Z_jumps_neuron(n);


    for(auto it = n.begin_neurite(); it != n.end_neurite(); ++it){
     auto f=get_node_measures(*it,markers);
      m.emplace_back(f);
    }

  } // End neuron for

  return m;
}


