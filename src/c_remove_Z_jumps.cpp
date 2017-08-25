#include "include/remove_Z_jumps.hpp"
#include <Rcpp.h>
#include <typeinfo>
#include <unistd.h>
namespace ns = neurostr::selector;

namespace bg =        boost::geometry;
using point_type =    bg::model::point<float, 3, bg::cs::cartesian>;

//This function removes the jumps over the Z axis of a neurite
void remove_Z_jumps_neurite(neurostr::Neurite& n)
{
  auto branches = ns::neurite_branch_selector(n);
  std::vector<neurostr::Node> tmp;
  point_type direction;

  for(auto it = n.begin_branch(); it != n.end_branch(); ++it)
  {
    if(it->size() > 1)
    {
      tmp.insert(tmp.end(), it->begin(), it->end());
      for(auto node = it->begin(); node != it->end();)
      {
        ns::node_parent(*node);
        if(node->valid_parent())
        {
          direction = node->parent().vectorTo(*node);
          if((direction.get<0>() == 0) && (direction.get<1>() == 0) && (direction.get<2>() != 0))
          {
            node = it->erase(node);
          }else{
            ++node;
          }//end IF
        }else{
          ++node;
        }//end IF
      }//end FOR
    }//end IF
  }//end FOR
}

//Remove the jumps over Z axis of a neuron
void remove_Z_jumps_neuron(neurostr::Neuron& n)
{
  for(auto it = n.begin_neurite(); it != n.end_neurite(); ++it){
    remove_Z_jumps_neurite(*it);
  }
}



