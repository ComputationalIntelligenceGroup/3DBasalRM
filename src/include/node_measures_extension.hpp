#ifndef NEUROSTR_MEASURE_NODE_EXTENSION_MEASURE_H_
#define NEUROSTR_MEASURE_NODE_EXTENSION_MEASURE_H_

/* This code extends the measures integrated in the raw version of neurostr.
 * Most of them compute a measure for the nodes.
 */

#include <neurostr/core/node.h>
#include <neurostr/core/branch.h>
#include <neurostr/core/neuron.h>

#include <neurostr/measure/universal_measure.h>
#include <neurostr/measure/branch_measure.h>
#include <neurostr/measure/neurite_measure.h>
#include <neurostr/measure/node_measure.h>
#include <neurostr/measure/aggregate.h>
#include <neurostr/measure/measure_operations.h>

#include <neurostr/selector/selector.h>
#include <neurostr/selector/node_selector.h>
#include <neurostr/selector/neurite_selector.h>
#include <neurostr/selector/neuron_selector.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>

namespace neurostr {
namespace measure {

using node_iterator = typename std::vector<selector::node_reference>::iterator;
using const_node_iterator = typename std::vector<selector::const_node_reference>::iterator;
/** Get node properties */

/* Parent-related measures */

// Get distance in straight line from the node to the root of the branch
static const auto node_distance_to_branch_root = [](const Node& n) -> float {
  if(n.valid_branch())
  {
    if(n.branch().has_root())
    {
      auto branch_root = n.branch().root();
      return (n.distance(branch_root));
    }else{
      return (0);
    }
  }
};

// Get length along the branch from the node to the root of the branch
static const auto node_length_to_branch_root = [](const Node& n) -> float {
  if(n.valid_branch())
  {
    float len = n.branch().has_root() ? n.branch().first().distance(n.branch().root()) : 0.;
    if(*n.branch().begin() != n)//If the node is the first elementof the branch finish
    {
      for(auto it = (n.branch().begin() + 1); *it != n; ++it)//Compute distance between each pair of consecutive nodes
      {
        len += it->distance(*(it-1));
      }

      len+=n.parent().distance(n);//Compute distance from the node to its parent
    }
      return(len);
  }else{
    return(std::numeric_limits<float>::quiet_NaN());
  }
};

static const auto node_fract_dim_to_branch_root = [](const Node& n) -> float {
  if(n.valid_branch())
  {
    // Get all the nodes in the branch and insert the root
    auto nodes=neurostr::selector::branch_node_selector(n.branch());

    // If the branch has root include it for the computation of fractal distance and increase the position
    // of the index of the input node in the vector
    int idx_node = 0;
    if(n.branch().has_root())
    {
      nodes.emplace(nodes.begin(), n.branch().root());
      ++idx_node;
    }

    //Get the index of the input node in the vector
    for(auto it = n.branch().begin(); *it != n; ++it) ++idx_node;

    //If the input node is not the last of the branch, remove all the nodes after it
    if((idx_node+1) < nodes.size())
    {
      nodes.erase(nodes.begin() + (idx_node + 1),nodes.end());
    }
    return(node_set_fractal_dim(nodes.begin(), nodes.end()));
  }else{
    return(std::numeric_limits<float>::quiet_NaN());
  }
};

// Get tortuosity from the root to the node
static const auto node_branch_tortuosity = [](const Node& n) -> float {
  if(n.valid_branch())
  {
    auto tmp_branch = std::make_unique<neurostr::Branch>();
    for(auto it = n.branch().begin(); *it != n; ++it)
    {
      tmp_branch->insert(tmp_branch->begin(), *it);
    }
    tmp_branch->insert(tmp_branch->begin(), n);

    return(tortuosity(*tmp_branch));
  }else{
    return(std::numeric_limits<float>::quiet_NaN());
  }
};

// Get partition asymmetry
static const auto node_branch_partition_asymmetry = [](const Node& n) -> float {
  if(n.valid_branch())
  {
    return(partition_asymmetry(n.branch()));
  }else{
    return(std::numeric_limits<float>::quiet_NaN());
  }
};

/*  NODE MEASUREMENTS
    BASED ON NEURITE COMPUTATIONS */

//Given a node return all the nodes of the neurite
static const std::vector<selector::const_node_reference> node_subtree(const neurostr::Node& n){
  std::vector<selector::const_node_reference> nodes;
  if(n.valid_branch())
  {
    if(n.branch().valid_neurite())
    {
      auto order=node_order(n);

      if(order>0)
      {
        //Get all the branches with equal or lower value than the input node
        std::vector<selector::const_branch_reference> branches;
        std::vector<selector::const_branch_reference> tmp_branch;
        for(int i = 0; i <= order; ++i)
        {
          tmp_branch = selector::compose_selector(selector::branch_order_filter_factory(i), selector::neurite_branch_selector)(n.branch().neurite());
          branches.insert(branches.end(), tmp_branch.begin(), tmp_branch.end());
        }

        std::vector<selector::const_node_reference> tmp_nodes;
        for(auto it = branches.begin(); it != branches.end(); ++it)
        {
          tmp_nodes = selector::branch_node_selector(it->get());
          nodes.insert(nodes.end(), tmp_nodes.begin(), tmp_nodes.end());
        }
      }
    }
  }
  return(nodes);
};

static const auto node_subtree_box_volume = [](const neurostr::Node& n) -> float {
  std::vector<selector::const_node_reference> nodes = node_subtree(n);
  if(!nodes.empty())
  {
    return(box_volume(nodes.begin(), nodes.end()));
  }else{
    return(0.0);
  }
};

static const auto node_subtree_length = [](const neurostr::Node& n) -> float {
  float subtree_length = 0.0;
  std::vector<selector::const_node_reference> nodes = node_subtree(n);

  if(!nodes.empty())
  {
    for(auto it = nodes.begin(); it != nodes.end(); ++it)
    {
      if(it->get().valid_parent())
      {
        subtree_length += it->get().distance(it->get().parent());
      }
    }
  }

  return(subtree_length);
};

static const auto node_subtree_bifurcations = [](const neurostr::Node& n) -> float {
  float num_descendants = 0.0;
  std::vector<selector::const_node_reference> nodes = node_subtree(n);

  if(!nodes.empty())
  {
    auto order = node_order(n);
    for(auto it = nodes.begin(); it != nodes.end(); ++it)
    {
      if((order >= node_order(*it)) && (selector::node_descendants(*it).size()==2)) ++num_descendants;
    }
  }

  return(num_descendants);
};

static const auto node_subtree_terminals = [](const neurostr::Node& n) -> float {
  float num_terminals = 0.0;
  std::vector<selector::const_node_reference> nodes=node_subtree(n);//Get all the nodes of same or lower order

  if(!nodes.empty())//If it is not the root of the neurite
  {
    auto order = node_order(n);

    //For each node get its descendants and check if they are in the same order
    for(auto it = nodes.begin(); it != nodes.end(); ++it)
    {
      if((order >= node_order(*it)) && selector::node_descendants(*it).empty())
      {
        ++num_terminals;
      }else if((order == node_order(*it)) && (selector::node_self(*it) == it->get().branch().last()))
      {
        ++num_terminals;
      }
    }
  }
  return(num_terminals);
};

static const float compute_dimension(const std::vector<selector::const_node_reference> nodes, const int position){
  float dim_length = 0.0;

  if(!nodes.empty())//If it is not the root of the neurite
  {
    // Create nx3 matrix
    Eigen::Matrix<float, Eigen::Dynamic, 3> m;
    auto n = std::distance(nodes.begin(), nodes.end());
    m.resize(n, 3);

    int i = 0;
    for (auto it = nodes.begin(); it != nodes.end(); ++it, ++i) {
      m(i, 0) = it->get().x();
      m(i, 1) = it->get().y();
      m(i, 2) = it->get().z();
    }

    // Then perform jacobisvd
    Eigen::JacobiSVD<decltype(m)> svd(m, Eigen::ComputeFullV);

    Eigen::Matrix<float, Eigen::Dynamic, 3> aux = m * svd.matrixV();

    // Eigen allocates matrices in column-major

    auto data = aux.data();
    std::vector<float> dimensions;
    for (int i = 0; i < 3; ++i) {
      float minv = std::numeric_limits<float>::max();
      float maxv = std::numeric_limits<float>::min();
      int tmp = n * i;
      for (int j = 0; j < n; ++j) {
        if (data[tmp + j] > maxv) maxv = data[tmp + j];
        if (data[tmp + j] < minv) minv = data[tmp + j];
      }
      dimensions.emplace_back(maxv - minv);
    }
    std::sort(dimensions.begin(), dimensions.end());
    dim_length = dimensions.at(position);
  }
  return dim_length;
};

static const auto node_subtree_depth = [](const neurostr::Node& n) -> float {
  std::vector<selector::const_node_reference> nodes = node_subtree(n);
  return(compute_dimension(nodes, 0));
};

static const auto node_subtree_height= [](const neurostr::Node& n) -> float {
  std::vector<selector::const_node_reference> nodes = node_subtree(n);
  return(compute_dimension(nodes, 1));
};

static const auto node_subtree_width = [](const neurostr::Node& n) -> float {
  std::vector<selector::const_node_reference> nodes = node_subtree(n);
  return(compute_dimension(nodes, 2));
};

static const auto node_max_distance_bw_nodes = [](const neurostr::Node& n) -> float {
  std::vector<selector::const_node_reference> nodes = node_subtree(n);
  float max_distance = 0;
  if(!nodes.empty())//If it is not the root of the neurite
  {
    float tmp_distance = 0;
    for(auto it = nodes.begin(); it != (nodes.end() - 1); ++it)
    {
      for(auto it2 = it + 1; it2 != nodes.end(); ++it2)
      {
        tmp_distance = it->get().distance(*it2);
        if(max_distance < tmp_distance)
        {
          max_distance = tmp_distance;
        }
      }
    }
  }

  return(max_distance);
};


static const std::vector<selector::const_node_reference> get_terminal_nodes(const std::vector<selector::const_node_reference> nodes) {
  std::vector<selector::const_node_reference> terminal_nodes;
  if(!nodes.empty())//If it is not the root of the neurite
  {
    //For each node get its descendants and check if they are in the same o
    for(auto it = nodes.begin(); it != nodes.end(); ++it)
    {
      if(selector::node_descendants(*it).empty())
      {
        terminal_nodes.insert(terminal_nodes.end(), *it);
      }
    }
  }
  return(terminal_nodes);
};

static const auto node_subtree_max_order = [](const neurostr::Node& n) -> float {
  std::vector<selector::const_node_reference> nodes = node_subtree(n);
  nodes = get_terminal_nodes(nodes);
  float max_order = 0;

  if(!nodes.empty())
  {

    for(auto it = nodes.begin(); it != nodes.end(); ++it)
    {
      if((max_order < node_order(*it)) && (selector::node_descendants(*it).empty())) max_order = node_order(*it);
    }
  }

  return(max_order);
};

static const auto node_subtree_min_order = [](const neurostr::Node& n) -> float {
  std::vector<selector::const_node_reference> nodes = node_subtree(n);
  nodes = get_terminal_nodes(nodes);
  float min_order = 0;
  auto order = node_order(n);
  bool first = true;

  if(!nodes.empty())
  {
    for(auto it = nodes.begin(); it != nodes.end(); ++it)
    {
      //If the node does not have descendants
      if(selector::node_descendants(*it).empty())
      {
        if(first)
        {
          min_order = node_order(*it);
          first = false;
        }else if(min_order > node_order(*it))
        {
          min_order = node_order(*it);
        }
      }
    }
  }

  return(min_order);
};

static const auto node_subtree_max_length = [](const neurostr::Node& n) -> float {
  std::vector<selector::const_node_reference> nodes = node_subtree(n);
  nodes = get_terminal_nodes(nodes);
  float max_length = 0;

  if(!nodes.empty())
  {
    float path_to_root;
    for(auto it = nodes.begin(); it != nodes.end(); ++it)
    {
      path_to_root = node_path_to_root(*it);
      if(max_length > path_to_root) max_length = path_to_root;
    }
  }

  return(max_length);
};

static const auto node_subtree_min_length = [](const neurostr::Node& n) -> float {
  std::vector<selector::const_node_reference> nodes = node_subtree(n);
  nodes = get_terminal_nodes(nodes);
  float min_length = 0;

  if(!nodes.empty())
  {
    float path_to_root;
    for(auto it = nodes.begin(); it != nodes.end(); ++it)
    {
      path_to_root = node_path_to_root(*it);
      if((min_length < path_to_root) && (min_length != 0)) min_length = path_to_root;
    }
  }

  return(min_length);
};

} // Measure ns
} // Neurostr ns
#endif
