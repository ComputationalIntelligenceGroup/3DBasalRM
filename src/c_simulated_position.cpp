/* This code computes features for each node. It is used to characterize the dendrite at node level.
 * We combine this measures to obtain a model that represents the dendrite basal arborization.
 */
#include <Rcpp.h>

#include <neurostr/io/parser_dispatcher.h>
#include <neurostr/core/log.h>
#include <neurostr/core/neuron.h>
#include <neurostr/io/JSONWriter.h>
#include <neurostr/selector/neuron_selector.h>

#include <neurostr/core/geometry.h>

#include "include/escape_string.hpp"
#include "include/node_measures_extension.hpp"
#include "include/remove_Z_jumps.hpp"

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/error/en.h>

#include <sstream>

namespace ns = neurostr::selector;
namespace nm = neurostr::measure;
using namespace Rcpp;
namespace bg =        boost::geometry;
using point_type =    bg::model::point<float, 3, bg::cs::cartesian>;
template <typename T>
using const_selector_reference = std::reference_wrapper<const T>;
using const_node_reference  = const_selector_reference<neurostr::Node>;

//Cast spherical coordinates to cartesian
point_type spherical2cartesian(float azimuth, float elevation, float radius, const std::array<point_type,3>& basis,int id_parent, point_type position)
{
  point_type cartesian_positions;
  point_type cartesian_coords;
  cartesian_positions.set<0>(std::cos(azimuth)*std::sin(elevation));
  cartesian_positions.set<1>(std::sin(azimuth)*std::sin(elevation));
  cartesian_positions.set<2>(std::cos(elevation));

  point_type base0;
  point_type base1;
  point_type base2;

  //Transpose basis matrix
  base0.set<0>(basis[0].get<0>());base0.set<1>(basis[1].get<0>());base0.set<2>(basis[2].get<0>());
  base1.set<0>(basis[0].get<1>());base1.set<1>(basis[1].get<1>());base1.set<2>(basis[2].get<1>());
  base2.set<0>(basis[0].get<2>());base2.set<1>(basis[1].get<2>());base2.set<2>(basis[2].get<2>());

  //Compute the placement of the points
  cartesian_coords.set<0>(radius*bg::dot_product(base0, cartesian_positions)+position.get<0>());
  cartesian_coords.set<1>(radius*bg::dot_product(base1, cartesian_positions)+position.get<1>());
  cartesian_coords.set<2>(radius*bg::dot_product(base2, cartesian_positions)+position.get<2>());

  return(cartesian_coords);
}

//Given a neuron and the spherical positions of new nodes add new nodes to the neuron
NumericVector get_node_coordinates(neurostr::Neuron& n, NumericMatrix new_nodes, NumericVector is_bifurcation,int init_id){
  NumericVector id_sim_nodes(new_nodes.rows());

  //For each new node
  for(int i=0;i<new_nodes.rows();++i)
  {
    //Get the id of the parent of the new node. Then get the parent using its id
    neurostr::Node::id_type id_node=new_nodes(i,new_nodes.cols()-1);
    auto it_parent_node=n.find(id_node);

    //If the reference of the parent is not empty, i.e., the new node has a parent
    if(it_parent_node.begin()!=it_parent_node.end())
    {
      //Place the new nodes in the neuron space
      auto grandparent=ns::node_parent(*it_parent_node.node());
      auto basis = it_parent_node->local_basis(grandparent, it_parent_node->branch().neurite().neuron().up());
      auto cartesian_coords=spherical2cartesian(new_nodes(i,0),new_nodes(i,1),new_nodes(i,2), basis, new_nodes(i,new_nodes.cols()-1), it_parent_node->position());

      //If we do not the next id to assign search the max id in the neuron
      if(init_id==0)
      {
        //Get all the nodes in the neuron and search for the max id
        std::vector<const_node_reference> neuron_nodes=ns::neuron_node_selector(n);
        neurostr::Node::id_type max_id_node=0;
        for(auto it=neuron_nodes.begin();it!=neuron_nodes.end();++it)
        {
          if(max_id_node<it->get().id())
          {
            max_id_node=it->get().id();
          }
        }
        init_id=max_id_node;
      }

      //Generate the new node
      id_sim_nodes(i)=init_id+i+1;
      neurostr::Node new_node_obj(id_sim_nodes(i),cartesian_coords,grandparent.radius());

      if(is_bifurcation(i)==0)//If it is not a bifurcation then insert the new node at the end of the branch
      {
        it_parent_node->branch().push_back(new_node_obj);
      }else{//If the node is a bifurcation create a new branch and add the node at the beginning of it
        neurostr::Neurite::branch_iterator pos=it_parent_node->branch().neurite().find(it_parent_node->branch());
        auto newpos=it_parent_node->branch().neurite().append_branch(pos,neurostr::Branch());

        newpos->order(pos->order()+1);
        newpos->root(*it_parent_node.node());
        newpos->push_back(new_node_obj);
      }
    }
  }
  return (id_sim_nodes);
}

//Add the simulated nodes to the neuron and return the ids of the new nodes
// [[Rcpp::export]]
List c_simulated_position(std::string json_info, NumericMatrix new_nodes, NumericVector is_bifurcation, int init_id, bool omitapical, bool omitaxon, bool omitdend, bool correct, bool remove_zjumps) {
  NumericVector id_sim_nodes;

  std::ifstream aux;
  std::ostringstream oss;

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
    id_sim_nodes=get_node_coordinates(n,new_nodes, is_bifurcation,init_id);
  } // End neuron for

  neurostr::io::JSONWriter writer(oss);
  writer.write(*r);

  std::string simulated_neuron = oss.str();

  return (List::create(simulated_neuron,id_sim_nodes));
}


